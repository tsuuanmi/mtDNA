#!/usr/bin/env python3
"""Web UI cho mtDNA batch pipeline.

Thay cho việc sửa tay `.env` + `scripts/batch_pipeline.sh` trước mỗi lần chạy:
form trên trình duyệt -> dựng sẵn dòng lệnh -> gọi `scripts/batch_pipeline.sh`
với các cờ tương ứng, rồi stream log về trang.

Chỉ dùng thư viện chuẩn của Python, cố ý không dùng FastAPI: venv của repo là
venv production đang chạy thật, không cài thêm gói vào đó.

    python3 web/server.py [--port 8765] [--host 0.0.0.0]
"""

from __future__ import annotations

import argparse
import json
import os
import queue
import re
import shlex
import signal
import subprocess
import threading
import time
import uuid
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse, parse_qs

REPO_ROOT = Path(__file__).resolve().parent.parent
WEB_DIR = REPO_ROOT / "web"
RUNS_DIR = WEB_DIR / "runs"
PRESETS_FILE = WEB_DIR / "presets.json"
BATCH_SCRIPT = "scripts/batch_pipeline.sh"

# Hai định dạng batch đang dùng: MS_DDMMYY_NNN (mới) và YYYYMMDD_mtDNA_NN (cũ).
BATCH_RE = re.compile(r"^(MS_\d{6}_\d{3}|\d{8}_mtDNA_\d+)$")
SHEET_ID_RE = re.compile(r"^[A-Za-z0-9_-]{20,}$")
VALID_STEPS = {"1", "2", "3", "4", "5", "6"}
# Beyond this the per-batch worker share drops to nothing useful and the run
# spends its time in contention rather than work.
MAX_JOBS = 32
VALID_PIPELINES = {"tracy", "blastn"}

MAX_LINES = 5000  # số dòng log giữ trong bộ nhớ cho mỗi lần chạy


class ConfigError(ValueError):
    """Người dùng nhập sai — trả về 400 kèm thông báo tiếng Việt."""


# =============================================================================
# .env
# =============================================================================

def load_env() -> dict[str, str]:
    """Đọc `.env` ở gốc repo. Chỉ để hiển thị giá trị mặc định trên form."""
    env: dict[str, str] = {}
    path = REPO_ROOT / ".env"
    if not path.exists():
        return env
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip()
        if value[:1] in {'"', "'"} and value[-1:] == value[:1] and len(value) > 1:
            value = value[1:-1]
        else:
            # Bỏ chú thích cuối dòng, ví dụ: GOOGLE_SHEET_ID=abc #Main api
            value = value.split(" #", 1)[0].split("\t#", 1)[0].strip()
        env[key] = value
    return env


def load_presets() -> dict:
    return json.loads(PRESETS_FILE.read_text(encoding="utf-8"))


# Mỗi nhiệm vụ một Google Sheet riêng. Giữ ID trong .env chứ không trong
# presets.json: presets.json được git theo dõi, .env thì không.
SHEET_KEYS = (
    ("main", "GOOGLE_SHEET_ID", "Sheet chính (batch thường quy)"),
    ("training", "GOOGLE_SHEET_ID_TRAINING", "Sheet training"),
    ("rerun", "GOOGLE_SHEET_ID_RERUN", "Sheet rerun"),
)


def load_sheets(env: dict[str, str]) -> dict[str, dict[str, str]]:
    """Các sheet đã khai trong .env, theo khoá ngắn."""
    sheets = {}
    for key, var, label in SHEET_KEYS:
        value = env.get(var, "").strip()
        if value:
            sheets[key] = {"key": key, "label": label, "id": value, "env_var": var}
    return sheets


SHEET_PLACEHOLDER = "***"


def redact_cmd(cmd: list[str]) -> list[str]:
    """Thay giá trị sau --sheet-id bằng dấu che.

    ID sheet không được xuất hiện ở bất cứ đâu người dùng nhìn thấy: dòng lệnh
    xem trước, dòng đầu của log, hay file web/runs/*.log. Nó nằm trong .env và
    chỉ đi thẳng tới tiến trình con.
    """
    out = list(cmd)
    for i, part in enumerate(out):
        if part == "--sheet-id" and i + 1 < len(out):
            out[i + 1] = SHEET_PLACEHOLDER
    return out


def display_cmd(cmd: list[str]) -> str:
    parts = [shlex.quote(c) for c in cmd]
    for i, part in enumerate(parts):
        if part == "--sheet-id" and i + 1 < len(parts):
            parts[i + 1] = SHEET_PLACEHOLDER
    return " ".join(parts)


def redact_cfg(cfg: dict) -> dict:
    out = dict(cfg)
    if out.get("sheet_id"):
        out["sheet_id"] = SHEET_PLACEHOLDER
    return out


# =============================================================================
# Kiểm tra dữ liệu vào và dựng dòng lệnh
# =============================================================================

def _clean_dir(value, field: str, must_exist: bool) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if not text.startswith("/"):
        raise ConfigError(f"{field} phải là đường dẫn tuyệt đối, nhận được: {text!r}")
    if "\x00" in text or "\n" in text:
        raise ConfigError(f"{field} chứa ký tự không hợp lệ")
    if must_exist and not Path(text).is_dir():
        raise ConfigError(f"{field} không tồn tại: {text}")
    return text


def parse_request(payload: dict) -> dict:
    """Chuẩn hoá + kiểm tra payload từ form. Ném ConfigError nếu sai."""
    steps = [str(s).strip() for s in payload.get("steps") or []]
    steps = [s for s in steps if s]
    if not steps:
        raise ConfigError("Chưa chọn bước nào để chạy.")
    bad = [s for s in steps if s not in VALID_STEPS]
    if bad:
        raise ConfigError(f"Bước không hợp lệ: {', '.join(bad)} (chỉ nhận 1-6)")
    steps = sorted(set(steps))

    pipeline = str(payload.get("pipeline") or "tracy").strip()
    if pipeline not in VALID_PIPELINES:
        raise ConfigError(f"Engine phải là tracy hoặc blastn, nhận được {pipeline!r}")

    raw_batches = payload.get("batches")
    if isinstance(raw_batches, str):
        raw_batches = re.split(r"[,\s]+", raw_batches)
    batches: list[str] = []
    for item in raw_batches or []:
        bid = str(item).strip()
        if not bid:
            continue
        if not BATCH_RE.match(bid):
            raise ConfigError(
                f"Mã batch không hợp lệ: {bid!r}. "
                "Định dạng đúng: MS_110826_015 hoặc 20241223_mtDNA_12"
            )
        if bid not in batches:
            batches.append(bid)
    if not batches:
        raise ConfigError("Chưa nhập mã batch nào.")

    # 0 means "as many as there are batches", which is the sensible default:
    # nothing is gained by holding batches back, and the core cap already keeps
    # the total work bounded however many run at once.
    jobs_raw = payload.get("jobs", 0)
    try:
        jobs = int(jobs_raw)
    except (TypeError, ValueError):
        raise ConfigError(f"Số batch chạy cùng lúc phải là số nguyên, nhận được {jobs_raw!r}") from None
    if not 0 <= jobs <= MAX_JOBS:
        raise ConfigError(f"Số batch chạy cùng lúc phải trong khoảng 0–{MAX_JOBS}, nhận được {jobs}")
    if jobs == 0:
        jobs = min(len(batches), MAX_JOBS)
    else:
        # More jobs than batches is harmless but pointless; trim so the command
        # shown in the preview matches what will actually happen.
        jobs = min(jobs, len(batches))

    # Free text, no accounts. It exists so a person looking at the run list can
    # tell whose run they are about to stop, not to authenticate anyone.
    owner = re.sub(r"[^\w .,()@-]", "", str(payload.get("owner") or "").strip(), flags=re.UNICODE)[:40]

    rerun = bool(payload.get("rerun"))
    sequencher_dir = _clean_dir(payload.get("sequencher_dir"), "Thư mục Sequencher", must_exist=True)
    results_dir = _clean_dir(payload.get("results_dir"), "Thư mục kết quả", must_exist=False)

    env = load_env()
    sheets = load_sheets(env)

    # `sheet` là khoá ngắn (main/training/rerun); `sheet_id` là ID gõ tay.
    sheet_key = str(payload.get("sheet") or "").strip()
    sheet_id = str(payload.get("sheet_id") or "").strip()
    if sheet_key:
        if sheet_key not in sheets:
            known = ", ".join(sheets) or "(chưa khai sheet nào trong .env)"
            raise ConfigError(f"Không có sheet tên {sheet_key!r}. Đang có: {known}")
        sheet_id = sheets[sheet_key]["id"]
    if sheet_id and not SHEET_ID_RE.match(sheet_id):
        raise ConfigError("Google Sheet ID không hợp lệ.")

    upload = payload.get("upload")
    if upload is not None:
        upload = bool(upload)

    # Bật upload thì phải chỉ rõ sheet. Không rơi về GOOGLE_SHEET_ID trong .env:
    # đó là sheet chính, và một lượt FASTA hay rerun lỡ bật upload sẽ ghi đè
    # worksheet thật của batch thường quy.
    if upload and not sheet_id:
        raise ConfigError(
            "Bật upload thì phải chọn sheet ghi vào. "
            "Chọn ở ô 'Ghi vào sheet nào', hoặc tắt upload."
        )

    return {
        "steps": steps,
        "pipeline": pipeline,
        "batches": batches,
        "rerun": rerun,
        "jobs": jobs,
        "owner": owner,
        "sequencher_dir": sequencher_dir,
        "results_dir": results_dir,
        "sheet_id": sheet_id,
        "sheet_key": sheet_key,
        "upload": upload,
    }


def build_command(cfg: dict) -> list[str]:
    """Dựng argv cho batch_pipeline.sh. Luôn là list, không bao giờ qua shell."""
    cmd = ["bash", BATCH_SCRIPT]
    if cfg["rerun"]:
        cmd.append("--rerun")
    cmd += ["-s", ",".join(cfg["steps"])]
    cmd += ["-p", cfg["pipeline"]]
    cmd += ["-b", ",".join(cfg["batches"])]
    if cfg["jobs"] > 1:
        cmd += ["-j", str(cfg["jobs"])]
    if cfg["sequencher_dir"]:
        cmd += ["--sequencher-dir", cfg["sequencher_dir"]]
    if cfg["results_dir"]:
        cmd += ["--results-dir", cfg["results_dir"]]
    if cfg["sheet_id"]:
        cmd += ["--sheet-id", cfg["sheet_id"]]
    if cfg["upload"] is True:
        cmd.append("--upload")
    elif cfg["upload"] is False:
        cmd.append("--no-upload")
    return cmd


# =============================================================================
# Soát dữ liệu đầu vào trên đĩa
# =============================================================================

def _match_batch_files(root: Path, batch_id: str, suffix: str) -> list[str]:
    """Tìm file thuộc đúng batch này.

    Neo vào trọn mã batch và bắt buộc có dấu phân cách (_ . hoặc khoảng trắng)
    ngay sau — nếu chỉ so khớp chuỗi con thì batch 20241222_mtDNA_12 sẽ kéo
    nhầm cả mtDNA_120, mtDNA_121.
    """
    if not root.is_dir():
        return []
    hits: list[str] = []
    patterns = (f"{batch_id}{suffix}", f"{batch_id}_*{suffix}",
                f"{batch_id}.*{suffix}", f"{batch_id} *{suffix}")
    for pattern in patterns:
        for path in root.rglob(pattern):
            if path.is_file() and "@eaDir" not in path.parts:
                hits.append(str(path))
    return sorted(set(hits))


def inspect_inputs(cfg: dict, env: dict[str, str]) -> list[dict]:
    """Với mỗi batch, cho biết sẽ dùng ZIP nào và có bao nhiêu file AB1."""
    lab_dir = Path(env.get("LAB_DATA_DIR", ""))
    seq_dir = Path(cfg["sequencher_dir"] or env.get("SEQUENCHER_DIR", ""))
    needs_ab1 = "2" in cfg["steps"]
    needs_zip = "4" in cfg["steps"] or cfg["rerun"]

    report = []
    for batch in cfg["batches"]:
        # "blockers" khác "warnings": blocker là thứ chắc chắn làm lượt chạy
        # thất bại, không phải điều đáng lưu ý.
        item: dict = {"batch": batch, "zips": [], "ab1_count": 0,
                      "warnings": [], "blockers": []}

        zips = _match_batch_files(seq_dir, batch, ".zip")
        item["zips"] = zips
        if needs_zip:
            if not zips:
                item["warnings"].append(f"Không thấy ZIP Sequencher nào trong {seq_dir}")
            elif len(zips) > 1:
                # Rerun chạy riêng từng thư mục nguồn, nên nhiều ZIP ở nhiều thư
                # mục là bình thường. Nhiều ZIP trong CÙNG một thư mục mới đáng lo:
                # chúng bị giải nén chung, và nếu trùng mẫu thì QC loại sạch.
                folders = Counter(str(Path(z).parent) for z in zips)
                if not cfg["rerun"]:
                    item["warnings"].append(
                        f"Có {len(zips)} ZIP khớp — tất cả sẽ được giải nén chung")
                elif max(folders.values()) > 1:
                    item["warnings"].append(
                        f"Có thư mục chứa nhiều hơn một ZIP của batch này. Chúng được "
                        f"giải nén chung — nếu hai ZIP có cùng mẫu thì QC loại hết và "
                        f"batch hỏng. Kiểm tra trước khi chạy.")
                else:
                    item["warnings"].append(
                        f"Có {len(zips)} ZIP ở {len(folders)} thư mục — mỗi thư mục chạy "
                        f"riêng và ra một tab riêng")

        if needs_ab1 and lab_dir.is_dir():
            # pipeline.sh khớp khi BẤT KỲ thư mục cha nào mang đúng mã batch;
            # thư mục ở đây lồng 2 tầng nên phải đếm theo tập file duy nhất,
            # nếu cộng dồn từng thư mục sẽ đếm trùng.
            ab1_dirs = {d for pattern in (f"{batch}*", f"*/{batch}*")
                        for d in lab_dir.glob(pattern) if d.is_dir()}
            files: set[str] = set()
            for d in ab1_dirs:
                files.update(str(f) for f in d.rglob("*.ab1"))
            count = len(files)
            item["ab1_count"] = count
            if count == 0:
                item["warnings"].append(f"Bước 2 (prepare) được chọn nhưng không thấy file .ab1 nào trong {lab_dir}")

        metadata = Path(env.get("DATA_DIR", "")) / "metadata" / f"{batch}.xlsx"
        item["metadata"] = str(metadata) if metadata.is_file() else ""
        if cfg["rerun"] and not item["metadata"]:
            # Lời cũ ("sẽ được copy ở bước prepare") sai: danh sách mẫu rerun được
            # lập ở giai đoạn 1, trước khi bước 2 chạy. Và từ bản mới của
            # batch_pipeline.sh, thiếu metadata làm cả lượt thoát ngay từ đầu.
            item["blockers"].append(
                "Thiếu file metadata Excel. Danh sách mẫu rerun được lập trước khi "
                "bước prepare chạy nên bước 2 không cứu được, và cả lượt sẽ dừng "
                "ngay từ đầu — không batch nào chạy.")

        report.append(item)
    return report


# =============================================================================
# Quản lý tiến trình
# =============================================================================

# =============================================================================
# Tóm tắt tiến độ
# =============================================================================

STEP_LABELS = {
    "1": "Kiểm tra Positive Control",
    "2": "Chuẩn bị dữ liệu AB1",
    "3": "Chạy engine phân tích",
    "4": "So sánh với Sequencher",
    "5": "Xuất FASTA",
    "6": "Sinh báo cáo PDF",
}
# Cùng số 2 nhưng khác việc hẳn: rerun không gom AB1 từ đầu mà thay read lab
# vừa giải trình tự lại vào đúng vùng. Gọi chung một tên là hiểu nhầm.
RERUN_STEP_LABELS = {"2": "Thay read lab đã chạy lại"}

RE_STEP_RUN = re.compile(r"\[Step ([1-6])/6")
RE_STEP_SKIP = re.compile(r"Skipping Step ([1-6])")
RE_BATCH_START = re.compile(r"\[([A-Za-z0-9_]+)\] Running pipeline")
# batch_pipeline.sh prefixes every line of a batch's output with its id, which
# is the only way to attribute a step once batches run concurrently and their
# streams interleave.
RE_BATCH_PREFIX = re.compile(r"^\[([A-Za-z0-9_]+)\] ")
# log_message viết "<timestamp> | INFO | … - [batch] …", không phải đầu dòng,
# nên tên batch phải moi từ giữa câu.
RE_NO_SAMPLES = re.compile(r"\[([A-Za-z0-9_]+)\] No rerun samples found")
RE_PREP_ABORT = re.compile(r"sample list for batch ([A-Za-z0-9_]+)")
# apply_reads tóm tắt mỗi batch một dòng; cộng dồn để hiện tổng cả lượt.
RE_READS_APPLIED = re.compile(
    r"([A-Za-z0-9_]+): thay (\d+) read, thêm (\d+) read, xoá \d+ read gốc, trên (\d+) mẫu")
RE_READS_NO_DIR = re.compile(r"No rerun AB1 directory .* for batch ([A-Za-z0-9_]+)")
RE_READS_FAILED = re.compile(r"Could not apply rerun reads for batch ([A-Za-z0-9_]+)")
RE_LEVEL = re.compile(r"\|\s*(WARNING|ERROR|CRITICAL)\s*\|")

MAX_PROBLEMS = 40


class Progress:
    """Dựng bảng tiến độ từ log, thay cho việc đổ nguyên log ra màn hình.

    Các mốc được nhận ra từ chính dòng log mà pipeline vẫn in; không phải thêm
    thiết bị đo nào vào script.
    """

    def __init__(self, cfg: dict, upload: bool):
        self.items: list[dict] = []
        self.problems: list[str] = []
        self._current_batch = cfg["batches"][0] if cfg["batches"] else ""
        self._tabs = 0

        labels = {**STEP_LABELS, **(RERUN_STEP_LABELS if cfg["rerun"] else {})}
        if cfg["rerun"]:
            self._add("prep", "Chọn mẫu rerun từ ZIP", "")
        for batch in cfg["batches"]:
            for step in cfg["steps"]:
                self._add(f"{batch}:{step}", labels[step], batch)
        self._add("merge", "Gộp kết quả so sánh", "")
        if upload:
            self._add("upload", "Đưa lên Google Sheet", "")

    def _add(self, key: str, label: str, batch: str) -> None:
        self.items.append({"key": key, "label": label, "batch": batch,
                           "status": "pending", "detail": ""})

    def _find(self, key: str) -> dict | None:
        return next((i for i in self.items if i["key"] == key), None)

    def _set(self, key: str, status: str, detail: str = "") -> bool:
        item = self._find(key)
        if item is None:
            return False
        # So cả detail chứ không chỉ status: bước upload ở lại "running" suốt
        # nhiều tab, chỉ có tên tab đổi.
        if item["status"] == status and (not detail or item["detail"] == detail):
            return False
        item["status"] = status
        if detail:
            item["detail"] = detail
        return True

    def _finish_running(self, batch: str | None = None) -> bool:
        """Mark running steps done. Limited to one batch when given."""
        changed = False
        for item in self.items:
            if item["status"] == "running" and (batch is None or item["batch"] == batch):
                item["status"] = "done"
                changed = True
        return changed

    def feed(self, line: str) -> bool:
        """Cập nhật theo một dòng log. Trả về True nếu bảng có thay đổi."""
        changed = False

        # Prefer the line's own prefix over the last batch seen: with --jobs > 1
        # the lines arrive interleaved, and "last one wins" attributes steps to
        # whichever batch happened to log most recently.
        prefix = RE_BATCH_PREFIX.match(line)
        batch = prefix.group(1) if prefix else self._current_batch

        if RE_LEVEL.search(line) or "Traceback" in line:
            if len(self.problems) < MAX_PROBLEMS:
                self.problems.append(line)
                changed = True

        match = RE_BATCH_START.search(line)
        if match:
            self._current_batch = match.group(1)
            changed |= self._set("prep", "done")
            return True

        match = RE_STEP_SKIP.search(line)
        if match:
            return changed | self._set(f"{batch}:{match.group(1)}", "skipped")

        match = RE_STEP_RUN.search(line)
        if match:
            # Only close out steps of the same batch. Finishing every running
            # step would mark another batch's in-flight work as done.
            changed |= self._finish_running(batch)
            return changed | self._set(f"{batch}:{match.group(1)}", "running")

        # Thay read xong ở giai đoạn 1, mỗi batch một dòng tóm tắt.
        match = RE_READS_APPLIED.search(line)
        if match:
            changed |= self._set("prep", "done")
            return changed | self._set(
                f"{match.group(1)}:2", "done",
                f"thay {match.group(2)}, thêm {match.group(3)} read / {match.group(4)} mẫu")
        match = RE_READS_NO_DIR.search(line)
        if match:
            return changed | self._set(f"{match.group(1)}:2", "skipped", "không thấy thư mục AB1 lab")
        match = RE_READS_FAILED.search(line)
        if match:
            return changed | self._set(f"{match.group(1)}:2", "failed", "vẫn chạy, nhưng bằng read gốc")
        if "Rerun units of work" in line:
            return changed | self._set("prep", "done")
        if "Cannot prepare rerun sample list" in line:
            match = RE_PREP_ABORT.search(line)
            who = match.group(1) if match else ""
            return changed | self._set(
                "prep", "failed",
                f"thiếu metadata: {who}" if who else "thiếu metadata")
        if "Running in RERUN mode" in line:
            return changed | self._set("prep", "running")
        if "Rerun ZIP sample summary" in line:
            # ví dụ: "... 12 unique sample(s) in ZIP, 2 found in raw TXT"
            tail = line.split("summary for batch", 1)[-1].strip()
            return changed | self._set("prep", "running", tail[:120])
        if "No rerun samples found" in line:
            # Chỉ bỏ qua đúng batch đó. Trước đây đánh dấu skipped cho MỌI batch
            # còn chờ: một batch không có mẫu nào là cả bảng tiến độ trắng xoá.
            match = RE_NO_SAMPLES.search(line)
            target = match.group(1) if match else batch
            changed |= self._set("prep", "running")
            for item in self.items:
                if item["status"] == "pending" and item["batch"] == target:
                    item["status"] = "skipped"
                    changed = True
            return changed
        if "Merging comparison results for batches" in line:
            changed |= self._finish_running()
            return changed | self._set("merge", "running")
        if "Merging comparison results for folder" in line:
            # Rerun ghi mỗi thư mục nguồn một tab, nên gộp và upload xen kẽ
            # nhau. Đây thuộc về bước upload, không phải quay lại bước gộp.
            folder = line.split("for folder", 1)[-1].split(":")[0].strip()
            changed |= self._set("merge", "done")
            return changed | self._set("upload", "running", f"đang gộp {folder}"[:80])
        if "Google Sheets upload disabled" in line:
            return changed | self._set("upload", "skipped", "đã tắt cho lần chạy này")
        if "credentials file not found" in line:
            return changed | self._set("upload", "skipped", "không có file credentials")
        if line.startswith("Uploading ") or " - Uploading " in line:
            changed |= self._set("merge", "done")
            worksheet = line.split("Uploading", 1)[-1].strip().split(" ")[0]
            self._tabs += 1
            return changed | self._set(
                "upload", "running", f"{worksheet} (tab {self._tabs})"[:80])
        if "Batch pipeline completed successfully" in line:
            return changed | self._finish_running()

        return changed

    def finalize(self, status: str) -> None:
        for item in self.items:
            if item["status"] == "running":
                item["status"] = "done" if status == "success" else "failed"
            elif item["status"] == "pending":
                item["status"] = "done" if status == "success" else "skipped"

    def snapshot(self) -> dict:
        return {"items": [dict(i) for i in self.items], "problems": list(self.problems)}


class Run:
    def __init__(self, run_id: str, cmd: list[str], cfg: dict):
        self.id = run_id
        self.cmd = cmd
        self.cfg = cfg
        self.owner = cfg.get("owner", "")
        self.cores = 0
        self.started_at = time.time()
        self.finished_at: float | None = None
        self.status = "running"          # running | success | failed | stopped
        self.returncode: int | None = None
        self.lines: list[str] = []
        self.subscribers: set[queue.Queue] = set()
        self.proc: subprocess.Popen | None = None
        self.log_path = RUNS_DIR / f"{run_id}.log"
        self._lock = threading.Lock()
        upload = cfg["upload"]
        if upload is None:
            upload = not cfg["rerun"]
        self.progress = Progress(cfg, bool(upload))

    def summary(self) -> dict:
        return {
            "id": self.id,
            "cmd": redact_cmd(self.cmd),
            "cmd_display": display_cmd(self.cmd),
            "cfg": redact_cfg(self.cfg),
            "status": self.status,
            "returncode": self.returncode,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "owner": self.owner,
            "cores": self.cores,
            "batches": list(self.cfg["batches"]),
            "started_display": datetime.fromtimestamp(self.started_at).strftime("%H:%M:%S %d/%m"),
            "line_count": len(self.lines),
            "progress": self.progress.snapshot(),
        }

    def emit(self, line: str) -> None:
        with self._lock:
            self.lines.append(line)
            progress_changed = self.progress.feed(line)
            if len(self.lines) > MAX_LINES:
                del self.lines[: len(self.lines) - MAX_LINES]
            payload = {"line": line}
            if progress_changed:
                payload["progress"] = self.progress.snapshot()
            dead = []
            for sub in self.subscribers:
                try:
                    sub.put_nowait(payload)
                except queue.Full:
                    dead.append(sub)
            for sub in dead:
                self.subscribers.discard(sub)

    def snapshot(self) -> list[str]:
        with self._lock:
            return list(self.lines)

    def attach(self, start: int = 0) -> tuple[list[str], dict, queue.Queue]:
        """Lấy log đã có và đăng ký nhận log mới trong cùng một lần giữ khoá.

        Nếu tách làm hai bước, dòng phát ra đúng khe giữa chúng sẽ bị mất
        (subscribe sau) hoặc bị lặp (subscribe trước).
        """
        q: queue.Queue = queue.Queue(maxsize=MAX_LINES)
        with self._lock:
            backlog = self.lines[start:] if start < len(self.lines) else []
            snapshot = self.progress.snapshot()
            self.subscribers.add(q)
        return backlog, snapshot, q

    def unsubscribe(self, q: queue.Queue) -> None:
        with self._lock:
            self.subscribers.discard(q)


class CoreBudget:
    """Tracks how many cores the running jobs have been given.

    Every run gets the same share, and the total is allowed past the core
    count. Rationing was tried and traded away: a run's worker counts are
    fixed at startup and cannot be shrunk for a newcomer, so sharing out a
    budget meant whoever started later ran slower for the whole of the earlier
    run. Oversubscribing costs some scheduler churn instead, which is the
    cheaper of the two.

    The numbers are still tracked so the UI can show when the machine is
    overcommitted.
    """

    def __init__(self, total: int, per_run: int):
        self.total = max(1, total)
        self.per_run = max(1, per_run)
        self.reserved: dict[str, int] = {}

    def used(self) -> int:
        return sum(self.reserved.values())

    def free(self) -> int:
        return max(0, self.total - self.used())

    def reserve(self, run_id: str) -> int:
        self.reserved[run_id] = self.per_run
        return self.per_run

    def release(self, run_id: str) -> None:
        self.reserved.pop(run_id, None)


class RunManager:
    """Holds every run, several of them live at once.

    Concurrency is limited per batch rather than globally. Two batches never
    touch the same paths -- everything a batch writes carries its id -- but the
    same batch run twice would have each wiping the other's tool output
    mid-flight, so that is refused.
    """

    HISTORY_LIMIT = 30

    def __init__(self):
        self.lock = threading.Lock()
        self.runs: list[Run] = []          # newest first, active and finished
        self.budget = CoreBudget(
            total=int(os.environ.get("MTDNA_TOTAL_CORES") or os.cpu_count() or 8),
            per_run=int(os.environ.get("MTDNA_MAX_CORES") or 128),
        )

    def active(self) -> list[Run]:
        return [r for r in self.runs if r.status == "running"]

    def is_busy(self) -> bool:
        return bool(self.active())

    def get(self, run_id: str) -> Run | None:
        return next((r for r in self.runs if r.id == run_id), None)

    def _busy_batches(self) -> dict[str, Run]:
        return {b: r for r in self.active() for b in r.cfg["batches"]}

    def start(self, cmd: list[str], cfg: dict) -> Run:
        with self.lock:
            busy = self._busy_batches()
            clashing = [b for b in cfg["batches"] if b in busy]
            if clashing:
                other = busy[clashing[0]]
                owner = f" ({other.owner})" if other.owner else ""
                raise ConfigError(
                    f"Batch {', '.join(clashing)} đang được chạy ở một lượt khác{owner}. "
                    "Đợi lượt đó xong, hoặc bỏ batch đó ra."
                )

            run_id = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:4]
            cores = self.budget.reserve(run_id)
            cmd = [*cmd, "--max-cores", str(cores)]
            run = Run(run_id, cmd, cfg)
            run.cores = cores
            RUNS_DIR.mkdir(parents=True, exist_ok=True)

            env = os.environ.copy()
            venv_bin = REPO_ROOT / ".venv" / "bin"
            if venv_bin.is_dir():
                env["PATH"] = f"{venv_bin}:{env.get('PATH', '')}"
            env["PYTHONUNBUFFERED"] = "1"

            run.proc = subprocess.Popen(
                cmd,
                cwd=str(REPO_ROOT),
                env=env,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                start_new_session=True,   # để dừng được cả cây tiến trình con
            )
            self.runs.insert(0, run)
            # Trim only finished runs; a live one must stay reachable however
            # many have piled up behind it.
            finished = [r for r in self.runs if r.status != "running"]
            for stale in finished[self.HISTORY_LIMIT:]:
                self.runs.remove(stale)

        threading.Thread(target=self._pump, args=(run,), daemon=True).start()
        return run

    def _pump(self, run: Run) -> None:
        header = f"$ {display_cmd(run.cmd)}"
        run.emit(header)
        with run.log_path.open("w", encoding="utf-8") as log:
            log.write(header + "\n")
            assert run.proc is not None and run.proc.stdout is not None
            for line in run.proc.stdout:
                line = line.rstrip("\n")
                run.emit(line)
                log.write(line + "\n")
                log.flush()
        run.returncode = run.proc.wait()
        MANAGER.budget.release(run.id)
        run.finished_at = time.time()
        if run.status != "stopped":
            run.status = "success" if run.returncode == 0 else "failed"
        run.progress.finalize(run.status)
        elapsed = int(run.finished_at - run.started_at)
        run.emit(f"__END__ {run.status} rc={run.returncode} ({elapsed}s)")

    def stop(self, run_id: str | None = None) -> bool:
        """Stop one run. Without an id, stops the newest running one."""
        if run_id:
            run = self.get(run_id)
        else:
            run = next(iter(self.active()), None)
        if run is None or run.status != "running" or run.proc is None:
            return False
        run.status = "stopped"
        try:
            os.killpg(os.getpgid(run.proc.pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            run.proc.terminate()
        return True


MANAGER = RunManager()


# =============================================================================
# HTTP
# =============================================================================

class Handler(BaseHTTPRequestHandler):
    server_version = "mtdna-ui"
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):  # bớt ồn trên terminal
        if "/api/stream" not in (self.path or ""):
            print(f"[{self.log_date_time_string()}] {fmt % args}", flush=True)

    # -- helpers ------------------------------------------------------------
    def _send_json(self, payload: dict, status: int = 200) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, path: Path, content_type: str) -> None:
        if not path.is_file():
            self._send_json({"error": "Không tìm thấy file"}, 404)
            return
        body = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _read_json(self) -> dict:
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            return {}
        if length > 1_000_000:
            raise ConfigError("Yêu cầu quá lớn")
        raw = self.rfile.read(length)
        try:
            data = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ConfigError(f"JSON không hợp lệ: {exc}") from exc
        if not isinstance(data, dict):
            raise ConfigError("Payload phải là object JSON")
        return data

    # -- routes -------------------------------------------------------------
    def do_GET(self) -> None:
        route = urlparse(self.path).path
        try:
            if route in ("/", "/index.html"):
                self._send_file(WEB_DIR / "index.html", "text/html; charset=utf-8")
            elif route == "/api/config":
                self._send_json(self._config_payload())
            elif route == "/api/status":
                self._send_json(self._status_payload())
            elif route == "/api/stream":
                self._stream()
            elif route == "/api/log":
                params = parse_qs(urlparse(self.path).query)
                run = MANAGER.get((params.get("id") or [""])[0])
                if run is None:
                    self._send_json({"error": "Không tìm thấy lần chạy này"}, 404)
                else:
                    self._send_json({"lines": run.snapshot(), "run": run.summary()})
            else:
                self._send_json({"error": "Không có route này"}, 404)
        except ConfigError as exc:
            self._send_json({"error": str(exc)}, 400)
        except BrokenPipeError:
            pass
        except Exception as exc:  # noqa: BLE001 - đừng để server chết vì 1 request
            self._send_json({"error": f"Lỗi máy chủ: {exc}"}, 500)

    def do_POST(self) -> None:
        route = urlparse(self.path).path
        try:
            if route == "/api/preview":
                cfg = parse_request(self._read_json())
                env = load_env()
                self._send_json({
                    "cmd": redact_cmd(build_command(cfg)),
                    "cmd_display": display_cmd(build_command(cfg)),
                    "effective": self._effective(cfg, env),
                    "inputs": inspect_inputs(cfg, env),
                })
            elif route == "/api/run":
                cfg = parse_request(self._read_json())
                run = MANAGER.start(build_command(cfg), cfg)
                self._send_json({"run": run.summary()})
            elif route == "/api/stop":
                body = self._read_json()
                stopped = MANAGER.stop(str(body.get("id") or "") or None)
                self._send_json({"stopped": stopped})
            else:
                self._send_json({"error": "Không có route này"}, 404)
        except ConfigError as exc:
            self._send_json({"error": str(exc)}, 400)
        except BrokenPipeError:
            pass
        except Exception as exc:  # noqa: BLE001
            self._send_json({"error": f"Lỗi máy chủ: {exc}"}, 500)

    # -- payloads -----------------------------------------------------------
    def _config_payload(self) -> dict:
        env = load_env()
        presets = load_presets()
        rerun_root = Path("/mnt/nas/bca/mtDNA/Sequencher_temp/rerun")
        rerun_dirs = []
        if rerun_root.is_dir():
            rerun_dirs = sorted(
                (str(p) for p in rerun_root.iterdir() if p.is_dir() and p.name != "@eaDir"),
                reverse=True,
            )
        sheets = load_sheets(env)
        return {
            "steps": presets["steps"],
            "presets": presets["presets"],
            # Chỉ gửi tên. ID không rời khỏi server, kể cả dạng đã che.
            "sheets": [{"key": v["key"], "label": v["label"]} for v in sheets.values()],
            "env": {
                "sequencher_dir": env.get("SEQUENCHER_DIR", ""),
                "results_dir": env.get("RESULTS_DIR", ""),
                "lab_data_dir": env.get("LAB_DATA_DIR", ""),
                "has_credentials": (REPO_ROOT / env.get("CREDENTIALS_FILE", "x")).is_file(),
            },
            "rerun_dirs": rerun_dirs,
        }

    def _effective(self, cfg: dict, env: dict[str, str]) -> dict:
        upload = cfg["upload"]
        if upload is None:
            upload = not cfg["rerun"]

        sheets = load_sheets(env)
        sheet_id = cfg["sheet_id"] or env.get("GOOGLE_SHEET_ID", "")
        # Nói rõ tên sheet chứ không chỉ "form hay .env": ghi nhầm batch rerun
        # đè lên sheet chính là lỗi đã từng xảy ra.
        by_id = {v["id"]: v["label"] for v in sheets.values()}
        sheet_label = by_id.get(sheet_id, "sheet lạ (ID gõ tay)" if sheet_id else "(chưa có)")

        return {
            "sequencher_dir": cfg["sequencher_dir"] or env.get("SEQUENCHER_DIR", ""),
            "results_dir": cfg["results_dir"] or env.get("RESULTS_DIR", ""),
            "upload": upload,
            "sheet_label": sheet_label,
            "jobs": cfg["jobs"],
        }

    def _status_payload(self) -> dict:
        active = MANAGER.active()
        return {
            "busy": bool(active),
            "active": [r.summary() for r in active],
            "recent": [r.summary() for r in MANAGER.runs if r.status != "running"][:10],
            "busy_batches": sorted({b for r in active for b in r.cfg["batches"]}),
            "cores": {
                "total": MANAGER.budget.total,
                "per_run": MANAGER.budget.per_run,
                "used": MANAGER.budget.used(),
                "free": MANAGER.budget.free(),
                "overcommitted": MANAGER.budget.used() > MANAGER.budget.total,
            },
        }

    def _stream(self) -> None:
        params = parse_qs(urlparse(self.path).query)
        run = MANAGER.get((params.get("id") or [""])[0]) or next(iter(MANAGER.active()), None)
        if run is None:
            self._send_json({"error": "Chưa có lần chạy nào"}, 404)
            return

        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream; charset=utf-8")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "close")
        self.end_headers()

        try:
            start = max(0, int((params.get("from") or ["0"])[0]))
        except ValueError:
            start = 0

        backlog, snapshot, sub = run.attach(start)
        try:
            # Tiến độ đi trước để trang vẽ được ngay, kể cả khi nối lại giữa chừng.
            self._sse({"progress": snapshot})
            for line in backlog:
                self._sse({"line": line})
            if run.status != "running" and not any(l.startswith("__END__") for l in backlog):
                self._sse({"line": f"__END__ {run.status} rc={run.returncode}",
                           "progress": run.progress.snapshot()})
                return
            while True:
                try:
                    payload = sub.get(timeout=15)
                except queue.Empty:
                    self.wfile.write(b": ping\n\n")
                    self.wfile.flush()
                    if run.status != "running":
                        return
                    continue
                self._sse(payload)
                if payload.get("line", "").startswith("__END__"):
                    return
        except (BrokenPipeError, ConnectionResetError):
            pass
        finally:
            run.unsubscribe(sub)

    def _sse(self, payload: dict) -> None:
        body = json.dumps(payload, ensure_ascii=False)
        self.wfile.write(f"data: {body}\n\n".encode("utf-8"))
        self.wfile.flush()


def main() -> None:
    presets = load_presets()
    parser = argparse.ArgumentParser(description="Web UI cho mtDNA batch pipeline")
    parser.add_argument("--host", default="0.0.0.0", help="mặc định 0.0.0.0 (mở trong LAN)")
    parser.add_argument("--port", type=int, default=presets.get("port", 8765))
    args = parser.parse_args()

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    server.daemon_threads = True
    print(f"mtDNA pipeline UI  ->  http://{args.host}:{args.port}")
    print(f"Repo: {REPO_ROOT}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nĐang tắt…")
        MANAGER.stop()
        server.shutdown()


if __name__ == "__main__":
    main()
