#!/usr/bin/env python3
"""Tự kiểm web UI.

Chạy server riêng trên cổng trống, gọi lần lượt các endpoint, kiểm tra cả
đường đúng lẫn đường sai. Mọi lần chạy thật đều ghi vào thư mục tạm và bị
ép `--no-upload`; test sẽ tự hỏng nếu có bất kỳ lệnh nào định upload lên
Google Sheet.

    python3 web/selftest.py                     # chỉ kiểm API, ~5 giây
    python3 web/selftest.py --full              # thêm một lượt chạy pipeline thật
    python3 web/selftest.py --tmp-dir /duong/dan
"""

from __future__ import annotations

import argparse
import json
import os
import re
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

PASS, FAIL = [], []


def check(name: str, ok: bool, detail: str = "") -> bool:
    (PASS if ok else FAIL).append(name)
    mark = "  ok  " if ok else " FAIL "
    print(f"[{mark}] {name}" + (f"  -- {detail}" if detail and not ok else ""))
    return ok


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def call(base: str, path: str, payload=None) -> tuple[int, dict]:
    url = base + path
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json"} if data else {},
        method="POST" if data is not None else "GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        body = exc.read().decode()
        try:
            return exc.code, json.loads(body)
        except json.JSONDecodeError:
            return exc.code, {"error": body[:200]}


# =============================================================================

def test_progress_model() -> None:
    """Bảng tiến độ dựng từ log — kiểm trực tiếp, không cần chạy pipeline."""
    sys.path.insert(0, str(REPO_ROOT / "web"))
    import server  # noqa: PLC0415

    cfg = {"batches": ["MS_270726_003"], "steps": ["2", "3", "4", "5"], "rerun": False,
           "upload": True, "pipeline": "tracy", "sequencher_dir": "", "results_dir": "",
           "sheet_id": "x", "sheet_key": "main"}
    log = [
        "| INFO | - [MS_270726_003] Running pipeline (normal mode)",
        "| INFO | - Skipping Step 1: Positive Control validation",
        "| INFO | - [Step 2/6 - prepare] Organizing AB1 files",
        "| INFO | - [Step 3/6 - automate] Running Tracy-only pipeline",
        "| WARNING | - Sequencher ZIP TXT has sample not found in raw TXT",
        "| INFO | - [Step 4/6 - manual] Extracting manual data",
        "| INFO | - Skipping Step 5: Generate FASTA files",
        "| INFO | - Merging comparison results for batches: MS_270726_003",
        "| INFO | - Uploading Batch_MS_270726_003 to sheet 1poUNuBm…",
        "| INFO | - Batch pipeline completed successfully!",
    ]

    ok = server.Progress(cfg, upload=True)
    for line in log:
        ok.feed(line)
    ok.finalize("success")
    by_label = {i["label"]: i for i in ok.snapshot()["items"]}
    check("tiến độ: bước đã chạy đánh dấu xong",
          by_label["Chạy engine phân tích"]["status"] == "done")
    check("tiến độ: bước bị bỏ đánh dấu skipped",
          by_label["Xuất FASTA"]["status"] == "skipped")
    check("tiến độ: có mục upload khi bật upload", "Đưa lên Google Sheet" in by_label)
    check("tiến độ: upload ghi rõ worksheet",
          by_label["Đưa lên Google Sheet"]["detail"].startswith("Batch_"))
    check("tiến độ: cảnh báo được gom lại", len(ok.snapshot()["problems"]) == 1)

    bad = server.Progress(cfg, upload=True)
    for line in log[:4]:
        bad.feed(line)
    bad.feed("| ERROR | - Traceback (most recent call last):")
    bad.finalize("failed")
    items = {i["label"]: i["status"] for i in bad.snapshot()["items"]}
    check("tiến độ: hỏng thì bước đang chạy đánh dấu failed",
          items["Chạy engine phân tích"] == "failed", str(items))
    check("tiến độ: hỏng thì các bước sau không bị coi là xong",
          all(v != "done" for k, v in items.items() if k != "Chuẩn bị dữ liệu AB1"),
          str(items))
    check("tiến độ: lỗi được gom lại", len(bad.snapshot()["problems"]) == 1)

    no_up = server.Progress(dict(cfg, upload=False), upload=False)
    check("tiến độ: không upload thì không có mục upload",
          all(i["label"] != "Đưa lên Google Sheet" for i in no_up.snapshot()["items"]))


COMPARISON_HEADERS = [
    "Sample ID", "Batch", "Analyzed Range (Pipeline)", "Variants (Pipeline)",
    "Variants Unique (Pipeline)", "Analyzed Range (Sequencher)", "Variants (Sequencher)",
    "Variants Unique (Sequencher)", "Variants Count (Sequencher)", "Concordant", "Flag",
    "Flagged (Pipeline)", "Flagged (Sequencher)",
]


def test_progress_rerun_tabs() -> None:
    """Rerun ghi nhiều tab: bảng tiến độ phải đi theo tab mới nhất."""
    sys.path.insert(0, str(REPO_ROOT / "web"))
    import server  # noqa: PLC0415

    cfg = {"batches": ["A_1", "B_2"], "steps": ["3", "4"], "rerun": True,
           "upload": True, "pipeline": "blastn", "sequencher_dir": "", "results_dir": "",
           "sheet_id": "x", "sheet_key": "rerun"}
    prog = server.Progress(cfg, upload=True)
    for line in [
        "| INFO | - Merging comparison results for folder Rerun_260618_T: A_1",
        "| INFO | - Uploading Rerun_260618_T to the configured sheet",
        "| INFO | - Merging comparison results for folder Rerun_260801_RT_013: B_2",
        "| INFO | - Uploading Rerun_260801_RT_013 to the configured sheet",
    ]:
        prog.feed(line)
    by_label = {i["label"]: i for i in prog.snapshot()["items"]}
    upload = by_label["Đưa lên Google Sheet"]

    # Lấy đúng tên tab: log thật là "Uploading X to the configured sheet",
    # cắt theo " to sheet" sẽ nuốt cả phần đuôi vào tên tab.
    check("rerun: tiến độ bám tab mới nhất",
          upload["detail"].startswith("Rerun_260801_RT_013"), upload["detail"])
    check("rerun: đếm đủ số tab đã ghi", "(tab 2)" in upload["detail"], upload["detail"])
    check("rerun: gộp theo thư mục không kéo ngược về bước gộp",
          by_label["Gộp kết quả so sánh"]["status"] == "done")


def test_input_blockers() -> None:
    """Thiếu metadata khi rerun là thứ giết cả lượt — phải báo đỏ, không lẫn cảnh báo."""
    sys.path.insert(0, str(REPO_ROOT / "web"))
    import server  # noqa: PLC0415

    with tempfile.TemporaryDirectory() as tmpdir:
        data = Path(tmpdir) / "data"
        (data / "metadata").mkdir(parents=True)
        (data / "metadata" / "CO_META_001.xlsx").write_bytes(b"x")
        env = {"DATA_DIR": str(data), "LAB_DATA_DIR": str(tmpdir), "SEQUENCHER_DIR": tmpdir}
        base = {"steps": ["3", "4"], "rerun": True, "sequencher_dir": tmpdir}

        thieu = server.inspect_inputs({**base, "batches": ["KHONG_META_001"]}, env)[0]
        co = server.inspect_inputs({**base, "batches": ["CO_META_001"]}, env)[0]

    check("thiếu metadata khi rerun bị xếp là blocker", len(thieu["blockers"]) == 1,
          str(thieu["blockers"]))
    check("blocker nói rõ cả lượt sẽ dừng",
          "không batch nào chạy" in " ".join(thieu["blockers"]), str(thieu["blockers"]))
    check("blocker KHÔNG nằm lẫn trong warnings",
          not any("metadata" in w for w in thieu["warnings"]), str(thieu["warnings"]))
    check("có metadata thì không bị chặn", not co["blockers"], str(co["blockers"]))

    # Lượt chạy chết ở giai đoạn 1 thì bảng tiến độ phải nêu tên batch.
    cfg = {"batches": ["A_1"], "steps": ["3", "4"], "rerun": True, "upload": False,
           "pipeline": "blastn", "sequencher_dir": "", "results_dir": "",
           "sheet_id": "", "sheet_key": "rerun"}
    prog = server.Progress(cfg, upload=False)
    prog.feed("| ERROR | - ERROR: Cannot prepare rerun sample list for batch A_1")
    item = next(i for i in prog.snapshot()["items"] if i["key"] == "prep")
    check("tiến độ: báo hỏng ở bước chọn mẫu", item["status"] == "failed", item["status"])
    check("tiến độ: nêu tên batch thiếu metadata", "A_1" in item["detail"], item["detail"])


def test_progress_rerun_reads() -> None:
    """Bước 2 của rerun là thay read lab, không phải gom AB1 như chạy thường."""
    sys.path.insert(0, str(REPO_ROOT / "web"))
    import server  # noqa: PLC0415

    base = {"steps": ["2", "3", "4"], "upload": False, "sequencher_dir": "",
            "results_dir": "", "sheet_id": "", "sheet_key": "rerun"}
    rerun_cfg = {**base, "batches": ["B_1", "B_2"], "rerun": True, "pipeline": "blastn"}
    thuong_cfg = {**base, "batches": ["B_1"], "rerun": False, "pipeline": "tracy"}

    def buoc2(prog, batch="B_1"):
        return next(i for i in prog.snapshot()["items"] if i["key"] == f"{batch}:2")

    check("chạy thường: bước 2 vẫn là gom AB1",
          buoc2(server.Progress(thuong_cfg, upload=False))["label"] == "Chuẩn bị dữ liệu AB1")
    check("rerun: bước 2 đổi tên thành thay read",
          buoc2(server.Progress(rerun_cfg, upload=False))["label"] == "Thay read lab đã chạy lại")

    prog = server.Progress(rerun_cfg, upload=False)
    for line in [
        "| INFO | - B_1: thay 4 read, thêm 0 read, xoá 4 read gốc, trên 1 mẫu",
        "| INFO | - B_2: thay 9 read, thêm 2 read, xoá 9 read gốc, trên 5 mẫu",
        "| WARNING | - 3 mẫu không có AB1 chạy lại, sẽ dùng read gốc: X, Y, Z",
    ]:
        prog.feed(line)
    check("rerun: bước 2 của từng batch báo đúng số của batch đó",
          buoc2(prog, "B_1")["detail"] == "thay 4, thêm 0 read / 1 mẫu", buoc2(prog, "B_1")["detail"])
    check("rerun: batch thứ hai không bị batch đầu ghi đè",
          buoc2(prog, "B_2")["detail"] == "thay 9, thêm 2 read / 5 mẫu", buoc2(prog, "B_2")["detail"])
    check("rerun: mẫu thiếu AB1 vào mục cảnh báo",
          any("không có AB1 chạy lại" in w for w in prog.snapshot()["problems"]))

    thieu = server.Progress(rerun_cfg, upload=False)
    thieu.feed("| INFO | - WARNING: No rerun AB1 directory at /x for batch B_1; using original reads")
    check("rerun: không có thư mục AB1 thì bước 2 báo bỏ qua",
          buoc2(thieu)["status"] == "skipped", buoc2(thieu)["status"])

    hong = server.Progress(rerun_cfg, upload=False)
    hong.feed("| INFO | - WARNING: Could not apply rerun reads for batch B_1; using original reads")
    check("rerun: thay read hỏng thì bước 2 báo hỏng",
          buoc2(hong)["status"] == "failed", buoc2(hong)["status"])


def test_theme() -> None:
    """Nền sáng/tối: mọi màu phải là biến, và hai bảng màu phải đủ như nhau."""
    html = (REPO_ROOT / "web" / "index.html").read_text()
    css = html.split("<style>", 1)[1].split("</style>", 1)[0]

    def block(marker: str) -> str:
        start = css.index("{", css.index(marker))
        depth = 0
        for k in range(start, len(css)):
            if css[k] == "{":
                depth += 1
            elif css[k] == "}":
                depth -= 1
                if depth == 0:
                    return css[start:k]
        return ""

    def defs(b: str) -> set[str]:
        return {m.group(1) for m in re.finditer(r"(--[a-z0-9-]+)\s*:", b)}

    root = defs(block(":root{")) - {"--mono"}
    media = defs(block("@media (prefers-color-scheme: light)"))
    light = defs(block(':root[data-theme="light"]'))
    used = {m.group(1) for m in re.finditer(r"var\((--[a-z0-9-]+)\)", css)}

    check("nền: có bảng màu theo hệ điều hành", bool(media))
    check("nền: có bảng màu khi bấm nút", bool(light))
    check("nền: biến nào dùng cũng được khai ở :root",
          not (used - root - {"--mono"}), str(sorted(used - root - {"--mono"})[:4]))
    # Thiếu một biến ở bảng sáng thì nó rơi về màu tối -- chữ tối trên nền tối.
    check("nền: bảng sáng (hệ điều hành) khai đủ màu",
          not (root - media), str(sorted(root - media)[:4]))
    check("nền: bảng sáng (nút bấm) khai đủ màu",
          not (root - light), str(sorted(root - light)[:4]))
    check("nền: hai bảng sáng không lệch nhau",
          not (media ^ light), str(sorted(media ^ light)[:4]))

    # Màu viết thẳng trong rule sẽ không đổi theo nền.
    hard = [ln.strip() for ln in css.splitlines()
            if re.search(r"^\s*[.#a-z@][^{]*\{[^}]*:#[0-9a-fA-F]{3,6}", ln)]
    check("nền: không còn màu viết cứng trong rule", not hard, str(hard[:2]))

    check("nút đổi nền có trên trang", 'id="btntheme"' in html)
    check("không lộ đường dẫn máy chủ trên trang",
          "repo_root" not in html and 'id="repo"' not in html)


def test_sheet_formatting() -> None:
    """Định dạng worksheet — dựng request rồi soi, không gọi tới Google."""
    sys.path.insert(0, str(REPO_ROOT))
    try:
        from src.modules.sheets.uploader import (  # noqa: PLC0415
            UNIFORM_COLUMN_WIDTH_PX, build_format_requests,
        )
    except ImportError as exc:
        check("nạp được module uploader", False, str(exc))
        return

    sheet_id = 424242
    n_cols = len(COMPARISON_HEADERS)
    reqs = build_format_requests(sheet_id, n_cols, n_rows=91)
    kinds = [k for r in reqs for k in r]

    check("đóng băng hàng tiêu đề", "updateSheetProperties" in kinds)
    check("bật bộ lọc", "setBasicFilter" in kinds)
    check("trả chiều cao dòng lại cho Google tự co", "autoResizeDimensions" in kinds)

    # Mọi cột một bề rộng, đặt bằng một request duy nhất trải hết các cột.
    widths = [r for r in reqs
              if r.get("updateDimensionProperties", {}).get("range", {})
              .get("dimension") == "COLUMNS"]
    check("chỉ một request đặt bề rộng cho tất cả cột", len(widths) == 1, str(len(widths)))
    if widths:
        rng = widths[0]["updateDimensionProperties"]["range"]
        check("request bề rộng phủ đúng mọi cột",
              rng["startIndex"] == 0 and rng["endIndex"] == n_cols, str(rng))
        check("dùng bề rộng mặc định khi không truyền",
              widths[0]["updateDimensionProperties"]["properties"]["pixelSize"]
              == UNIFORM_COLUMN_WIDTH_PX)

    custom = build_format_requests(sheet_id, n_cols, 91, column_width=160)
    custom_px = next(r["updateDimensionProperties"]["properties"]["pixelSize"] for r in custom
                     if r.get("updateDimensionProperties", {}).get("range", {})
                     .get("dimension") == "COLUMNS")
    check("--column-width đổi được bề rộng", custom_px == 160, str(custom_px))

    def fmt_of(start_row: int) -> dict:
        for r in reqs:
            if r.get("repeatCell", {}).get("range", {}).get("startRowIndex") == start_row:
                return r["repeatCell"]["cell"]["userEnteredFormat"]
        return {}

    check("ô dữ liệu WRAP để chữ chảy dọc", fmt_of(1).get("wrapStrategy") == "WRAP")
    check("ô dữ liệu canh trên, không canh giữa dọc",
          fmt_of(1).get("verticalAlignment") == "TOP")
    check("tiêu đề WRAP và canh giữa",
          fmt_of(0).get("wrapStrategy") == "WRAP"
          and fmt_of(0).get("horizontalAlignment") == "CENTER")
    check("tiêu đề in đậm có nền",
          fmt_of(0).get("textFormat", {}).get("bold") is True
          and "backgroundColor" in fmt_of(0))

    problems = []

    def walk(node) -> None:
        if isinstance(node, dict):
            if "sheetId" in node:
                if node["sheetId"] != sheet_id:
                    problems.append(f"sheetId sai: {node['sheetId']}")
                if node.get("endColumnIndex", 0) > n_cols:
                    problems.append(f"vượt số cột: {node}")
                if ("startRowIndex" in node and "endRowIndex" in node
                        and node["endRowIndex"] <= node["startRowIndex"]):
                    problems.append(f"range rỗng: {node}")
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    walk(reqs)
    check("mọi range hợp lệ và trỏ đúng worksheet", not problems, str(problems[:2]))

    # Bảng rỗng vẫn phải dựng được request hợp lệ, không sinh range âm.
    empty = build_format_requests(sheet_id, n_cols, n_rows=0)
    auto = next(r for r in empty if "autoResizeDimensions" in r)
    check("bảng rỗng không sinh range rỗng",
          auto["autoResizeDimensions"]["dimensions"]["endIndex"] >= 1)

    # fields phải kể tên đúng những khoá được đặt, nếu không Sheets bỏ qua im lặng.
    mismatched = []
    for r in reqs:
        if "repeatCell" not in r:
            continue
        fields = r["repeatCell"]["fields"]
        for key in r["repeatCell"]["cell"]["userEnteredFormat"]:
            if key not in fields:
                mismatched.append(key)
    check("fields kể đủ khoá được đặt", not mismatched, str(mismatched))


def test_core_budget() -> None:
    """Cấp nhân cho từng lượt — kiểm trực tiếp, không cần chạy pipeline."""
    sys.path.insert(0, str(REPO_ROOT / "web"))
    import server  # noqa: PLC0415

    b = server.CoreBudget(total=256, per_run=128)
    shares = [b.reserve(f"run{i}") for i in range(4)]
    check(f"mọi lượt được cấp bằng nhau: {shares}", len(set(shares)) == 1, str(shares))
    check("mỗi lượt đúng 128 nhân", shares[0] == 128, str(shares[0]))

    # Vượt số nhân là chấp nhận được, nhưng phải nói ra để còn biết.
    check("báo được khi cấp quá số nhân thật", b.used() > b.total,
          f"{b.used()}/{b.total}")
    check("phần trống không tụt xuống âm", b.free() == 0, str(b.free()))

    for i in range(4):
        b.release(f"run{i}")
    check("trả lại đủ nhân khi các lượt kết thúc", b.used() == 0, str(b.used()))
    check("lượt sau vẫn được cấp đầy", b.reserve("moi") == 128)


def test_config(base: str) -> dict:
    status, cfg = call(base, "/api/config")
    check("GET /api/config trả về 200", status == 200, str(status))
    check("có đủ 6 bước", len(cfg.get("steps", [])) == 6)
    check("có ít nhất 1 quy trình sẵn", len(cfg.get("presets", [])) >= 1)
    check("danh sách sheet chỉ có tên, không có ID",
          all(set(x) == {"key", "label"} for x in cfg.get("sheets", [])),
          str(cfg.get("sheets")))
    return cfg


def test_index(base: str) -> None:
    try:
        with urllib.request.urlopen(base + "/", timeout=10) as resp:
            html = resp.read().decode()
        check("GET / trả về trang HTML", resp.status == 200 and "<title>" in html)
        check("trang có bảng tiến độ", 'id="plist"' in html)
        check("trang có log chi tiết (thu gọn)", "<details" in html and 'id="log"' in html)
        check("form và màn tiến độ tách riêng",
              'id="view-form"' in html and 'id="view-run"' in html)
    except Exception as exc:  # noqa: BLE001
        check("GET / trả về trang HTML", False, str(exc))


def test_rejects(base: str) -> None:
    """Dữ liệu vào sai phải bị chặn ở server, không phải chỉ ở giao diện."""
    cases = [
        ("bỏ trống bước", {"steps": [], "batches": "MS_110826_001"}),
        ("bước ngoài 1-6", {"steps": ["7"], "batches": "MS_110826_001"}),
        ("mã batch sai định dạng", {"steps": ["4"], "batches": "batch_bat_ky"}),
        ("bỏ trống batch", {"steps": ["4"], "batches": ""}),
        ("engine lạ", {"steps": ["4"], "pipeline": "bwa", "batches": "MS_110826_001"}),
        ("đường dẫn tương đối", {"steps": ["4"], "batches": "MS_110826_001",
                                 "sequencher_dir": "khong/tuyet/doi"}),
        ("thư mục không tồn tại", {"steps": ["4"], "batches": "MS_110826_001",
                                   "sequencher_dir": "/khong/he/ton/tai"}),
        ("chèn lệnh shell", {"steps": ["4"], "batches": "MS_110826_001; rm -rf /"}),
        ("chèn lệnh qua backtick", {"steps": ["4"], "batches": "MS_110826_001`id`"}),
        ("xuống dòng trong đường dẫn", {"steps": ["4"], "batches": "MS_110826_001",
                                        "sequencher_dir": "/tmp\nrm -rf /"}),
        ("sheet id sai định dạng", {"steps": ["4"], "batches": "MS_110826_001", "sheet_id": "x"}),
    ]
    for name, payload in cases:
        status, body = call(base, "/api/preview", payload)
        check(f"chặn: {name}", status == 400 and "error" in body, f"status={status}")


def test_presets(base: str, cfg: dict) -> None:
    """Mỗi quy trình sẵn phải dựng ra dòng lệnh chạy được."""
    for p in cfg["presets"]:
        payload = {
            "steps": p["steps"], "pipeline": p["pipeline"], "batches": "MS_270726_003",
            "rerun": p["rerun"], "sequencher_dir": p["sequencher_dir"],
            "results_dir": p["results_dir"], "sheet": p.get("sheet", ""),
            "upload": p["upload"] if isinstance(p["upload"], bool) else None,
        }
        status, body = call(base, "/api/preview", payload)
        if not check(f"quy trình '{p['label']}' dựng được lệnh", status == 200,
                     body.get("error", "")):
            continue
        cmd = body["cmd"]
        check(f"  '{p['id']}': lệnh là mảng tham số", isinstance(cmd, list) and cmd[0] == "bash")
        check(f"  '{p['id']}': đúng các bước", cmd[cmd.index("-s") + 1] == ",".join(sorted(p["steps"])))
        check(f"  '{p['id']}': cờ rerun khớp", ("--rerun" in cmd) == p["rerun"])
        if p["upload"] is False:
            check(f"  '{p['id']}': có --no-upload", "--no-upload" in cmd)
        if p["upload"] is True:
            check(f"  '{p['id']}': có --upload", "--upload" in cmd)
            check(f"  '{p['id']}': có --sheet-id đi kèm", "--sheet-id" in cmd)


def test_upload_guard(base: str) -> None:
    """Rerun không được tự ý upload — đây là lỗi từng làm hỏng sheet production."""
    status, body = call(base, "/api/preview", {
        "steps": ["3", "4"], "pipeline": "blastn", "batches": "MS_010426_001",
        "rerun": True, "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_temp/rerun",
        "upload": None,
    })
    if check("rerun: xem trước chạy được", status == 200, body.get("error", "")):
        check("rerun mặc định KHÔNG upload", body["effective"]["upload"] is False)

    status, body = call(base, "/api/preview", {
        "steps": ["2", "3", "4"], "pipeline": "tracy", "batches": "MS_270726_003",
        "rerun": False, "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_temp",
        "upload": None,
    })
    if check("thường quy: xem trước chạy được", status == 200, body.get("error", "")):
        check("thường quy mặc định CÓ upload", body["effective"]["upload"] is True)


def test_no_sheet_id_leak(base: str, cfg: dict) -> None:
    """Không ID sheet nào được xuất hiện trong bất cứ phản hồi nào.

    Kiểm bằng cách lấy ID thật từ .env rồi tìm nguyên văn trong phần thân trả
    về, thay vì tin vào việc chỗ nào đó đã che — chỉ cần bỏ sót một chỗ là lộ.
    """
    real_ids = [v for v in (load_env_value(var) for _, var, _ in
                            (("main", "GOOGLE_SHEET_ID", ""),
                             ("training", "GOOGLE_SHEET_ID_TRAINING", ""),
                             ("rerun", "GOOGLE_SHEET_ID_RERUN", ""))) if v]
    if not check("đọc được ID sheet từ .env để dò", bool(real_ids)):
        return

    probes: list[tuple[str, str]] = []
    for endpoint in ("/api/config", "/api/status"):
        _, body = call(base, endpoint)
        probes.append((endpoint, json.dumps(body, ensure_ascii=False)))
    try:
        with urllib.request.urlopen(base + "/", timeout=10) as resp:
            probes.append(("/ (trang HTML)", resp.read().decode()))
    except Exception as exc:  # noqa: BLE001
        check("đọc được trang HTML để dò", False, str(exc))

    for preset in cfg["presets"]:
        _, body = call(base, "/api/preview", {
            "steps": preset["steps"], "pipeline": preset["pipeline"],
            "batches": "MS_270726_003", "rerun": preset["rerun"],
            "sequencher_dir": preset["sequencher_dir"],
            "results_dir": preset["results_dir"], "sheet": preset["sheet"],
            "upload": preset["upload"],
        })
        probes.append((f"/api/preview ({preset['id']})", json.dumps(body, ensure_ascii=False)))

    for where, text in probes:
        hit = next((i for i in real_ids if i in text), None)
        check(f"không lộ ID sheet ở {where}", hit is None,
              f"tìm thấy {hit[:6] if hit else ''}…")

    # File log trên đĩa cũng không được giữ ID: nó tồn lâu hơn phiên làm việc.
    logs = sorted((REPO_ROOT / "web" / "runs").glob("*.log"))
    if logs:
        leaked = [
            log.name for log in logs
            if any(i in log.read_text(errors="replace") for i in real_ids)
        ]
        check(f"không lộ ID sheet trong {len(logs)} file web/runs/*.log",
              not leaked, str(leaked[:3]))

    # Dòng lệnh vẫn phải có --sheet-id, chỉ là giá trị đã bị che.
    _, body = call(base, "/api/preview", {
        "steps": ["4"], "batches": "MS_270726_003", "sheet": "main", "upload": True,
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_temp",
    })
    cmd = body.get("cmd", [])
    if check("xem trước vẫn dựng được lệnh có --sheet-id", "--sheet-id" in cmd, str(body.get("error"))):
        value = cmd[cmd.index("--sheet-id") + 1]
        check("giá trị --sheet-id đã bị che", value == "***", value)


def test_sheets(base: str, cfg: dict) -> None:
    """Mỗi nhiệm vụ phải ghi vào đúng sheet của nó."""
    keys = {x["key"] for x in cfg.get("sheets", [])}
    check("có khai sheet trong .env", bool(keys), str(keys))

    seen_ids = set()
    for sheet in cfg.get("sheets", []):
        status, body = call(base, "/api/preview", {
            "steps": ["4"], "batches": "MS_270726_003", "sheet": sheet["key"],
            "upload": True, "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_temp",
        })
        if not check(f"sheet '{sheet['key']}': chọn được", status == 200, body.get("error", "")):
            continue
        eff = body["effective"]
        check(f"sheet '{sheet['key']}': xem trước gọi đúng tên sheet",
              eff["sheet_label"] == sheet["label"],
              f"{eff['sheet_label']} != {sheet['label']}")
        seen_ids.add(eff["sheet_label"])
    check("các sheet phân biệt nhau", len(seen_ids) == len(keys), str(seen_ids))

    status, body = call(base, "/api/preview", {
        "steps": ["4"], "batches": "MS_270726_003", "sheet": "khong_ton_tai",
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_temp",
    })
    check("chặn: tên sheet không có thật", status == 400)

    # Cái bẫy chính: bật upload mà không chọn sheet thì trước đây rơi về sheet
    # chính trong .env — một lượt FASTA hay rerun sẽ ghi đè batch thường quy.
    status, body = call(base, "/api/preview", {
        "steps": ["5"], "batches": "MS_270726_003", "sheet": "", "upload": True,
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher",
    })
    check("chặn: bật upload mà không chọn sheet", status == 400, f"status={status}")


def test_preset_sheets(base: str, cfg: dict) -> None:
    """Preset rerun/training không được trỏ vào sheet chính."""
    by_key = {x["key"]: x for x in cfg.get("sheets", [])}
    main = by_key.get("main", {}).get("label")
    for preset in cfg["presets"]:
        if preset["id"] in ("rerun", "training"):
            check(f"preset '{preset['id']}' dùng sheet riêng, không phải sheet chính",
                  preset.get("sheet") not in ("", "main", None),
                  f"sheet={preset.get('sheet')!r}")
            check(f"preset '{preset['id']}' có upload (API riêng đang dùng thật)",
                  preset["upload"] is True, f"upload={preset['upload']!r}")
        if preset["id"] == "fasta":
            check("preset 'fasta' không gắn sheet nào", not preset.get("sheet"))
    check("preset thường quy dùng sheet chính",
          next(p["sheet"] for p in cfg["presets"] if p["id"] == "automate") == "main")
    check("biết ID sheet chính để đối chiếu", bool(main), str(main))


def test_no_crosstalk(base: str, cfg: dict) -> None:
    """Quy trình phụ upload thì phải upload vào sheet của nó, không phải sheet chính."""
    main_label = next((x["label"] for x in cfg["sheets"] if x["key"] == "main"), None)
    if not check("biết sheet chính để đối chiếu", bool(main_label)):
        return
    for preset in cfg["presets"]:
        if preset["id"] in ("automate", "compare") or preset["upload"] is not True:
            continue
        status, body = call(base, "/api/preview", {
            "steps": preset["steps"], "pipeline": preset["pipeline"],
            "batches": "MS_270726_003", "rerun": preset["rerun"],
            "sequencher_dir": preset["sequencher_dir"],
            "results_dir": preset["results_dir"], "sheet": preset["sheet"],
            "upload": True,
        })
        if not check(f"'{preset['id']}': xem trước chạy được", status == 200,
                     body.get("error", "")):
            continue
        eff = body["effective"]
        check(f"'{preset['id']}': upload thật sự bật", eff["upload"] is True)
        check(f"'{preset['id']}': KHÔNG ghi vào sheet chính",
              eff["sheet_label"] != main_label,
              f"đang trỏ vào {eff['sheet_label']}")


def test_inputs(base: str) -> None:
    """Số liệu xem trước phải khớp với `find` mà pipeline thật dùng."""
    status, body = call(base, "/api/preview", {
        "steps": ["2", "4"], "pipeline": "tracy", "batches": "MS_270726_003",
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_training", "upload": False,
    })
    if not check("soát dữ liệu vào chạy được", status == 200, body.get("error", "")):
        return
    item = body["inputs"][0]
    check("khớp đúng 1 zip cho batch", len(item["zips"]) == 1, str(item["zips"]))

    lab = Path("/mnt/nas/bca/mtDNA/Lab/Raw_data")
    if lab.is_dir():
        out = subprocess.run(
            ["find", str(lab), "-type", "f", "-name", "*.ab1",
             "-path", "*/MS_270726_003[_. ]*/*"],
            capture_output=True, text=True,
        )
        real = len([x for x in out.stdout.splitlines() if x])
        check(f"số file ab1 khớp find thật ({item['ab1_count']} vs {real})",
              item["ab1_count"] == real)


def wait_for_run(base: str, run_id: str, timeout: int = 900) -> dict | None:
    """Chờ một lượt cụ thể kết thúc, trả về bản tóm tắt của chính nó."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        _, st = call(base, "/api/status")
        if not any(r["id"] == run_id for r in st.get("active", [])):
            break
        time.sleep(3)
    _, st = call(base, "/api/status")
    for bucket in ("active", "recent"):
        for run in st.get(bucket, []):
            if run["id"] == run_id:
                return run
    return None


def test_concurrency(base: str, tmp: Path) -> None:
    """Nhiều lượt chạy song song được, trừ khi trùng batch."""
    # Cố ý chọn bước 5 (fasta): chỉ sinh file cục bộ. Bước 6 có gọi API portal
    # nên không dùng cho test, kể cả khi bị dừng ngay sau đó.
    def payload(batches: str, owner: str, out: str) -> dict:
        return {
            "steps": ["5"], "pipeline": "tracy", "batches": batches, "owner": owner,
            "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_training",
            "results_dir": str(tmp / out), "sheet": "", "upload": False,
        }

    status, first = call(base, "/api/run", payload("MS_270726_003", "Nguoi A", "conc-a"))
    if not check("lượt 1 được nhận", status == 200, first.get("error", "")):
        return
    check("lượt 1 ghi nhận người chạy", first["run"]["owner"] == "Nguoi A", str(first["run"]["owner"]))

    # Batch khác -> chạy song song được.
    status, second = call(base, "/api/run", payload("MS_270726_009", "Nguoi B", "conc-b"))
    check("lượt 2 batch KHÁC chạy song song được", status == 200, second.get("error", ""))

    # Cùng batch -> từ chối, vì hai lượt sẽ xoá output của nhau giữa chừng.
    status, third = call(base, "/api/run", payload("MS_270726_003", "Nguoi C", "conc-c"))
    check("lượt 3 batch TRÙNG bị từ chối", status == 400, f"status={status}")
    if status == 400:
        check("lời từ chối nêu tên batch và người đang chạy",
              "MS_270726_003" in third.get("error", "") and "Nguoi A" in third.get("error", ""),
              third.get("error", ""))

    _, st = call(base, "/api/status")
    check("trạng thái liệt kê đủ các lượt đang chạy", len(st.get("active", [])) >= 2,
          str(len(st.get("active", []))))
    check("trạng thái nêu batch nào đang bận",
          "MS_270726_003" in st.get("busy_batches", []), str(st.get("busy_batches")))

    # Dừng từng lượt theo id, không phải "lượt hiện tại".
    for run in list(st.get("active", [])):
        call(base, "/api/stop", {"id": run["id"]})

    deadline = time.time() + 60
    while time.time() < deadline:
        _, st = call(base, "/api/status")
        if not st.get("active"):
            break
        time.sleep(0.5)
    _, st = call(base, "/api/status")
    check("dừng được từng lượt theo id", not st.get("active"),
          str([r["owner"] for r in st.get("active", [])]))


def test_stream_resume(base: str, run_id: str, total: int) -> None:
    """?from=N chỉ gửi phần còn thiếu, không phát lại từ đầu."""
    import urllib.request as ur
    keep = max(0, total - 3)
    req = ur.Request(f"{base}/api/stream?id={run_id}&from={keep}")
    got: list[str] = []
    saw_progress = False
    try:
        with ur.urlopen(req, timeout=20) as resp:
            for raw in resp:
                text = raw.decode().strip()
                if not text.startswith("data: "):
                    continue
                # Gói SSE không phải lúc nào cũng có "line": gói đầu là tiến độ.
                payload = json.loads(text[6:])
                if "progress" in payload:
                    saw_progress = True
                line = payload.get("line")
                if line is None:
                    continue
                got.append(line)
                if line.startswith("__END__") or len(got) > total:
                    break
    except Exception as exc:  # noqa: BLE001
        check("stream nối lại được", False, str(exc))
        return
    check("stream gửi tiến độ ngay khi nối vào", saw_progress)
    check(f"stream ?from={keep} chỉ gửi phần thiếu ({len(got)} dòng, không phải {total})",
          0 < len(got) <= total - keep + 1, f"nhận {len(got)}")


def _uploader_stub(tmp: Path, name: str) -> tuple[Path, Path]:
    """A `python` that swallows uploader.py calls and forwards everything else.

    Placed ahead on PATH, it records what the upload *would* have sent and
    returns 0, so the dispatch can be checked without a packet reaching Google.
    """
    stub_dir = tmp / f"stub-{name}"
    stub_dir.mkdir(parents=True, exist_ok=True)
    record = tmp / f"uploader-calls-{name}.txt"
    record.write_text("")

    real_python = REPO_ROOT / ".venv" / "bin" / "python"
    if not real_python.exists():
        real_python = Path(sys.executable)

    stub = stub_dir / "python"
    stub.write_text(
        "#!/bin/bash\n"
        'if [[ "${1:-}" == *sheets/uploader.py ]]; then\n'
        f'    printf "%s\\n" "$*" >> "{record}"\n'
        "    exit 0\n"
        "fi\n"
        f'exec "{real_python}" "$@"\n'
    )
    stub.chmod(0o755)
    return stub_dir, record


def expected_rerun_folders(sequencher_dir: str, batch: str) -> set[str]:
    """Thư mục cha trực tiếp của mỗi ZIP thuộc batch — mỗi cái là một worksheet.

    Bám đúng luật tìm ZIP của batch_pipeline.sh: khớp trọn batch id rồi phải có
    dấu ngăn, để 20241222_mtDNA_12 không vơ luôn ZIP của mtDNA_120.
    """
    root = Path(sequencher_dir.rstrip("/"))
    folders = set()
    for zip_path in root.rglob("*.zip"):
        name = zip_path.name
        if name == f"{batch}.zip" or (
            name.startswith(batch) and len(name) > len(batch) and name[len(batch)] in "_. "
        ):
            folders.add(zip_path.parent.name)
    return folders


def check_upload_dispatch(
    tmp: Path, cfg: dict, *, label: str, sheet_key: str, env_var: str,
    batch: str, results: Path, sequencher_dir: str, rerun: bool,
) -> None:
    """Upload có được gọi không, và gọi vào sheet nào — không gửi gì lên Google."""
    if not results.is_dir():
        print(f"       (bỏ qua {label}: chưa có kết quả để upload)")
        return
    if sheet_key not in {x["key"] for x in cfg["sheets"]}:
        print(f"       (bỏ qua {label}: .env chưa khai sheet {sheet_key})")
        return

    sheet_id = load_env_value(env_var)
    if not check(f"{label}: đọc được ID sheet từ .env", bool(sheet_id)):
        return

    stub_dir, record = _uploader_stub(tmp, sheet_key)
    env = os.environ.copy()
    env["PATH"] = f"{stub_dir}:{env.get('PATH', '')}"

    cmd = ["bash", "scripts/batch_pipeline.sh"]
    if rerun:
        cmd.append("--rerun")
    cmd += ["-s", "5", "-p", "blastn" if rerun else "tracy", "-b", batch,
            "--sequencher-dir", sequencher_dir, "--results-dir", str(results),
            "--sheet-id", sheet_id, "--upload"]

    proc = subprocess.run(
        cmd, cwd=str(REPO_ROOT), env=env, capture_output=True, text=True, timeout=900
    )
    if not check(f"{label}: batch_pipeline chạy xong", proc.returncode == 0,
                 (proc.stderr or proc.stdout)[-400:]):
        return

    calls = [c for c in record.read_text().splitlines() if c.strip()]
    check(f"{label}: uploader.py được gọi ({len(calls)} lần)", len(calls) >= 1)
    def uploaded_to(name: str) -> bool:
        return any(f"-w {name} -r" in c or c.rstrip().endswith(f"-w {name}") for c in calls)

    if rerun:
        # Rerun không dùng All_Batches: mỗi thư mục nguồn là một đợt rerun
        # riêng, gộp chung lại thì cùng một mẫu ra hai dòng không phân biệt được.
        check(f"{label}: KHÔNG có worksheet All_Batches",
              not uploaded_to("All_Batches"), str(calls)[:220])
    else:
        check(f"{label}: có worksheet All_Batches", uploaded_to("All_Batches"), str(calls)[:220])
    if rerun:
        # Rerun gom vài mẫu từ nhiều batch, nên không tạo tab cho từng batch.
        # Bảng gộp được tách theo thư mục cha TRỰC TIẾP của ZIP: một rerun
        # thường trỏ vào cả cây .../rerun, mỗi đợt rerun một thư mục con, nên
        # lấy tên gốc cây sẽ dồn tất cả vào một tab.
        folders = expected_rerun_folders(sequencher_dir, batch)
        for folder in sorted(folders):
            check(f"{label}: có worksheet mang tên thư mục '{folder}'",
                  uploaded_to(folder), str(calls)[:220])
        check(f"{label}: KHÔNG tạo worksheet riêng cho từng batch",
              not any(f"Batch_{batch}" in c for c in calls), str(calls)[:200])
        check(f"{label}: số lần upload = số thư mục nguồn",
              len(calls) == len(folders), f"{len(calls)} vs {len(folders)}")
    else:
        check(f"{label}: có worksheet của batch",
              any(f"Batch_{batch}" in c for c in calls), str(calls)[:200])
    check(f"{label}: mọi lời gọi trỏ vào sheet {sheet_key}",
          all(sheet_id in c for c in calls),
          "\n".join(c.replace(sheet_id, f"<{sheet_key.upper()}>") for c in calls))

    if sheet_key != "main":
        main_id = load_env_value("GOOGLE_SHEET_ID")
        check(f"{label}: KHÔNG lời gọi nào chạm sheet chính",
              bool(main_id) and not any(main_id in c for c in calls))


def test_fasta_export(base: str, tmp: Path) -> None:
    """Xuất FASTA (bước 5) chạy thật, ghi vào thư mục tạm, không upload."""
    batch = "MS_270726_003"
    results = tmp / "results"
    if not (results / "manual_pipeline" / "regenerate" / batch).is_dir():
        print("       (bỏ qua FASTA: cần test_real_run chạy trước)")
        return

    payload = {
        "steps": ["5"], "pipeline": "tracy", "batches": batch, "rerun": False,
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_training",
        "results_dir": str(results), "sheet": "", "upload": False,
    }
    status, body = call(base, "/api/preview", payload)
    if not check("FASTA: xem trước chạy được", status == 200, body.get("error", "")):
        return
    if body["effective"]["upload"] is not False or "--no-upload" not in body["cmd"]:
        check("FASTA: đã chặn upload", False, "lệnh vẫn có thể upload — dừng test")
        return
    check("FASTA: đã chặn upload", True)
    check("FASTA: ghi vào thư mục tạm", str(results) in body["cmd"])

    print(f"       → đang xuất FASTA thật vào {results / 'fasta' / batch}")
    status, body = call(base, "/api/run", payload)
    if not check("FASTA: được nhận", status == 200, body.get("error", "")):
        return
    run_id = body["run"]["id"]

    cur = wait_for_run(base, run_id)
    if not check("FASTA: kết thúc thành công", cur and cur["status"] == "success",
                 f"status={cur and cur['status']} rc={cur and cur['returncode']}"):
        _, log = call(base, f"/api/log?id={run_id}")
        for line in log.get("lines", [])[-15:]:
            print("        |", line[:170])
        return

    fasta_dir = results / "fasta" / batch
    files = sorted(fasta_dir.glob("*.fasta")) if fasta_dir.is_dir() else []
    if not check(f"FASTA: sinh ra file .fasta ({len(files)} file)", bool(files), str(fasta_dir)):
        return

    # File rỗng hoặc chỉ có header thì coi như hỏng, dù bước chạy xong.
    text = files[0].read_text()
    check("FASTA: file có header vùng HV", text.lstrip().startswith(">"), text[:40])
    bases = "".join(l for l in text.splitlines() if not l.startswith(">"))
    check(f"FASTA: có chuỗi thật ({len(bases)} base)", len(bases) > 100, str(len(bases)))
    check("FASTA: chỉ gồm ký tự nucleotide",
          set(bases.upper()) <= set("ACGTNRYSWKMBDHV-"), str(sorted(set(bases.upper()))[:8]))


def load_env_value(key: str) -> str:
    """Đọc một biến trong .env, không qua server."""
    path = REPO_ROOT / ".env"
    if not path.exists():
        return ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if line.startswith(f"{key}=") and not line.startswith("#"):
            value = line.split("=", 1)[1].strip().strip('"').strip("'")
            return value.split(" #", 1)[0].strip()
    return ""


def test_real_run(base: str, tmp: Path) -> str | None:
    """Chạy thật, ghi hết vào thư mục tạm, tuyệt đối không upload."""
    results = tmp / "results"
    payload = {
        "steps": ["3", "4"], "pipeline": "tracy", "batches": "MS_270726_003",
        "rerun": False,
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_training",
        "results_dir": str(results),
        "upload": False,
    }
    status, body = call(base, "/api/preview", payload)
    if not check("chạy thật: xem trước chạy được", status == 200, body.get("error", "")):
        return None

    # Chốt chặn: không cho test chạm vào Google Sheet trong bất kỳ hoàn cảnh nào.
    if body["effective"]["upload"] is not False or "--no-upload" not in body["cmd"]:
        check("chạy thật: đã chặn upload", False, "lệnh vẫn có thể upload — dừng test")
        return None
    check("chạy thật: đã chặn upload", True)
    check("chạy thật: ghi vào thư mục tạm", str(results) in body["cmd"])

    print(f"       → đang chạy pipeline thật, kết quả vào {results}")
    started = time.time()
    status, body = call(base, "/api/run", payload)
    if not check("chạy thật: được nhận", status == 200, body.get("error", "")):
        return None
    run_id = body["run"]["id"]

    cur = wait_for_run(base, run_id)
    elapsed = int(time.time() - started)
    ok = check(f"chạy thật: kết thúc thành công ({elapsed}s)",
               cur and cur["status"] == "success",
               f"status={cur and cur['status']} rc={cur and cur['returncode']}")
    if not ok:
        _, log = call(base, f"/api/log?id={run_id}")
        for line in log.get("lines", [])[-15:]:
            print("        |", line[:170])
        return run_id

    merged = results / "modules" / "comparison" / "all_batches_comparison.xlsx"
    check("chạy thật: có file so sánh tổng hợp", merged.is_file(), str(merged))
    per_batch = results / "modules" / "comparison" / "MS_270726_003" / "MS_270726_003.xlsx"
    check("chạy thật: có file so sánh của batch", per_batch.is_file(), str(per_batch))
    return run_id


def test_rerun_scope(base: str, tmp: Path) -> None:
    """Rerun chỉ được phân tích đúng số mẫu trong danh sách rerun.

    Hồi quy cho lỗi: `-l` trước đây chỉ được đọc ở bước 2, mà SOP rerun chạy
    `-s 3,4` nên bước 2 bị bỏ qua — cả batch bị phân tích lại thay vì vài mẫu.
    """
    batch = "MS_010426_001"
    results = tmp / "rerun"
    payload = {
        "steps": ["3", "4"], "pipeline": "blastn", "batches": batch, "rerun": True,
        "sequencher_dir": "/mnt/nas/bca/mtDNA/Sequencher_temp/rerun",
        "results_dir": str(results), "upload": False,
    }
    status, body = call(base, "/api/preview", payload)
    if not check("rerun: xem trước chạy được", status == 200, body.get("error", "")):
        return
    if not body["inputs"][0]["zips"]:
        print("       (bỏ qua: không còn ZIP rerun cho batch mẫu)")
        return
    if body["effective"]["upload"] is not False or "--no-upload" not in body["cmd"]:
        check("rerun: đã chặn upload", False, "lệnh vẫn có thể upload — dừng test")
        return
    check("rerun: đã chặn upload", True)

    print(f"       → đang chạy rerun thật, kết quả vào {results}")
    status, body = call(base, "/api/run", payload)
    if not check("rerun: được nhận", status == 200, body.get("error", "")):
        return
    run_id = body["run"]["id"]

    cur = wait_for_run(base, run_id)
    if not check("rerun: kết thúc thành công",
                 cur and cur["status"] == "success",
                 f"status={cur and cur['status']} rc={cur and cur['returncode']}"):
        _, log = call(base, f"/api/log?id={run_id}")
        for line in log.get("lines", [])[-15:]:
            print("        |", line[:170])
        return

    # Đơn vị việc của rerun là cặp (batch, thư mục nguồn): cùng một batch nằm ở
    # hai thư mục là hai lần phân tích Sequencher khác nhau, phải chạy tách ra
    # thì bản này mới không đè bản kia.
    pairs_file = results / "rerun_lists" / "pairs.tsv"
    if not check("rerun: có bảng cặp (batch, thư mục)", pairs_file.is_file(), str(pairs_file)):
        return
    pairs = [ln.split("\t") for ln in pairs_file.read_text().splitlines() if ln.strip()]
    if not check("rerun: sinh ra ít nhất một cặp việc", len(pairs) >= 1, str(pairs)):
        return

    for pair_batch, tab, _folder in pairs:
        rerun_list = results / "rerun_lists" / tab / f"{pair_batch}_rerun.txt"
        expected = len([x for x in rerun_list.read_text().splitlines() if x.strip()]) \
            if rerun_list.is_file() else -1
        comparison = (results / "by_folder" / tab / "modules" / "comparison"
                      / pair_batch / f"{pair_batch}.tsv")
        if not check(f"rerun: có file so sánh cho '{tab}'", comparison.is_file(), str(comparison)):
            return
        rows = len([x for x in comparison.read_text().splitlines() if x.strip()]) - 1
        check(f"rerun: '{tab}' chỉ so sánh đúng {expected} mẫu trong danh sách rerun (thấy {rows})",
              expected > 0 and rows == expected)

    # Rerun không gộp All_Batches: các thư mục là những đợt rerun riêng biệt.
    merged_dir = results / "modules" / "comparison"
    check("rerun: KHÔNG tạo All_Batches",
          not (merged_dir / "all_batches_comparison.xlsx").exists())
    tabs = sorted({tab for _b, tab, _f in pairs})
    for tab in tabs:
        check(f"rerun: có bảng gộp riêng cho '{tab}'",
              (merged_dir / f"comparison_{tab}.xlsx").is_file(),
              str(merged_dir / f"comparison_{tab}.xlsx"))


# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(description="Tự kiểm web UI của pipeline")
    parser.add_argument("--full", action="store_true", help="chạy thêm một lượt pipeline thật")
    parser.add_argument("--tmp-dir", default=None, help="thư mục tạm để chứa kết quả test")
    args = parser.parse_args()

    # Ghi ra file thì Python gom khối; ép line-buffer để theo dõi được tiến độ.
    sys.stdout.reconfigure(line_buffering=True)

    tmp = Path(args.tmp_dir) if args.tmp_dir else Path(tempfile.mkdtemp(prefix="mtdna-selftest-"))
    tmp.mkdir(parents=True, exist_ok=True)

    port = free_port()
    base = f"http://127.0.0.1:{port}"
    print(f"Server test : {base}")
    print(f"Thư mục tạm : {tmp}\n")

    proc = subprocess.Popen(
        [sys.executable, "web/server.py", "--host", "127.0.0.1", "--port", str(port)],
        cwd=str(REPO_ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
    )
    try:
        for _ in range(60):
            try:
                urllib.request.urlopen(base + "/api/status", timeout=1)
                break
            except Exception:  # noqa: BLE001
                if proc.poll() is not None:
                    print("Server chết ngay khi khởi động:\n" + (proc.stderr.read() or ""))
                    return 1
                time.sleep(0.2)
        else:
            print("Server không lên sau 12 giây")
            return 1

        test_progress_model()
        test_progress_rerun_tabs()
        test_theme()
        test_progress_rerun_reads()
        test_input_blockers()
        test_sheet_formatting()
        test_core_budget()
        cfg = test_config(base)
        test_index(base)
        test_rejects(base)
        test_presets(base, cfg)
        test_upload_guard(base)
        test_no_sheet_id_leak(base, cfg)
        test_sheets(base, cfg)
        test_preset_sheets(base, cfg)
        test_no_crosstalk(base, cfg)
        test_inputs(base)
        test_concurrency(base, tmp)

        if args.full:
            test_rerun_scope(base, tmp)
            check_upload_dispatch(
                tmp, cfg, label="upload rerun", sheet_key="rerun",
                env_var="GOOGLE_SHEET_ID_RERUN", batch="MS_010426_001",
                results=tmp / "rerun",
                sequencher_dir="/mnt/nas/bca/mtDNA/Sequencher_temp/rerun", rerun=True,
            )
            run_id = test_real_run(base, tmp)
            test_fasta_export(base, tmp)
            check_upload_dispatch(
                tmp, cfg, label="upload training", sheet_key="training",
                env_var="GOOGLE_SHEET_ID_TRAINING", batch="MS_270726_003",
                results=tmp / "results",
                sequencher_dir="/mnt/nas/bca/mtDNA/Sequencher_training", rerun=False,
            )
            if run_id:
                _, log = call(base, f"/api/log?id={run_id}")
                test_stream_resume(base, run_id, len(log.get("lines", [])))
        else:
            print("\n(bỏ qua lượt chạy pipeline thật — thêm --full nếu muốn)")
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()

    print(f"\n{'=' * 60}")
    print(f"Đạt: {len(PASS)}   Hỏng: {len(FAIL)}")
    if FAIL:
        for name in FAIL:
            print("  hỏng:", name)
    print(f"Kết quả test nằm ở: {tmp}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
