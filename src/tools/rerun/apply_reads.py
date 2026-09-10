#!/usr/bin/env python3
"""Thay read gốc bằng read đã giải trình tự lại, theo từng vùng.

Lab chạy lại một mẫu khi read đầu tiên quá xấu, và thường chỉ chạy lại đúng
vùng hỏng chứ không chạy lại cả mẫu -- gần 80% mẫu rerun chỉ có 1-2 read. Nên
đây không phải "thay cả mẫu": khoá thay thế là (mã mẫu, vùng, chiều), vùng nào
có read mới thì vùng đó bị thay, các vùng còn lại giữ nguyên read gốc.

Phải xoá read gốc chứ không chỉ copy read mới vào: tên file mã hoá cả giếng,
ngày và số thứ tự nên hai lần chạy không bao giờ trùng tên, mà tracy thì gom
mọi .ab1 có chứa mã mẫu trong tên. Để cả hai lại là đưa đúng read xấu -- cái
vừa phải chạy lại vì nó xấu -- vào phân tích.
"""
from __future__ import annotations

import json
import re
import shutil
import sys
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from loguru import logger

MANIFEST_NAME = ".rerun_applied.json"

# Chiều F/R đứng trước hoặc sau vị trí mồi: HV1F, HV1-15930F, HV1R-16410.
# Vị trí mồi không tham gia khớp -- HV1-15930F vẫn là "HV1 chiều F".
_REGION = re.compile(r"^(HV\d)(?:[-_][^FR]*)?([FR])(?:[-_].*)?$", re.IGNORECASE)


def batch_token(batch_id: str) -> str:
    """Phần định danh dùng để nhận thư mục bên Lab.

    Thư mục lab mang NGÀY CHẠY LẠI chứ không phải ngày batch
    (``20260714_mtDNA_427_2 mẫu`` là của batch ``20260610_mtDNA_427``), nên
    không neo được vào trọn mã batch. Bỏ trường ngày, giữ phần còn lại.
    """
    head, _, rest = batch_id.partition("_")
    return rest if (len(head) == 8 and head.isdigit() and rest) else batch_id


def folder_matches(name: str, token: str) -> bool:
    """Tên thư mục có mang đúng token này không -- không phải chuỗi con.

    Chặn hai đầu: mtDNA_32 không được trúng mtDNA_320 hay mtDNA_328, mà
    mtDNA_32_2 mẫu thì vẫn trúng.
    """
    return re.search(rf"(?<![0-9A-Za-z]){re.escape(token)}(?![0-9])", name) is not None


@dataclass(frozen=True)
class Read:
    """Một file AB1 đã tách được tên."""

    path: Path
    sample: str
    region: str  # HV1 | HV2 | HV3
    direction: str  # F | R
    date: str  # YYYYMMDD

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.sample, self.region, self.direction)


def parse_ab1_name(path: Path) -> Read | None:
    """``<giếng>_<ngày>_<mã mẫu>_<vùng>_<số>.ab1`` -> Read, hoặc None.

    Tách từ phải sang vì mã mẫu có thể chứa dấu gạch dưới (``257715_2``,
    ``LN_26_AB6721``), nên đếm trường từ trái sẽ cắt nhầm. Vùng chiếm một hoặc
    hai trường (``HV1R`` so với ``HV2_29F``).
    """
    parts = path.stem.split("_")
    if len(parts) < 4 or len(parts[1]) != 8 or not parts[1].isdigit() or not parts[-1].isdigit():
        return None
    for span in (1, 2):
        if len(parts) < 3 + span:
            continue
        matched = _REGION.match("_".join(parts[-1 - span : -1]))
        if not matched:
            continue
        sample = "_".join(parts[2 : -1 - span])
        if not sample:
            return None
        return Read(
            path=path,
            sample=sample,
            region=matched.group(1).upper(),
            direction=matched.group(2).upper(),
            date=parts[1],
        )
    return None


def collect_rerun_reads(
    rerun_root: Path, wanted: set[str], token: str | None = None
) -> tuple[dict[tuple[str, str, str], Read], list[Path]]:
    """Read mới nhất cho từng (mẫu, vùng, chiều), cùng danh sách file không tách được.

    Mới nhất tính riêng cho từng vùng, không phải cho cả mẫu: một mẫu có thể
    được chạy lại HV2R hôm 15/7 rồi HV3F hôm 31/7, và cả hai đều phải giữ.
    """
    best: dict[tuple[str, str, str], Read] = {}
    unparsed: list[Path] = []
    # Chỉ nhận read nằm trong thư mục mang đúng định danh batch. Mã mẫu một
    # mình chưa đủ: 6730 mã nằm ở nhiều batch, nên không có ràng buộc này thì
    # read chạy lại của batch khác vẫn lọt vào.
    roots = [rerun_root]
    if token:
        roots = [d for d in rerun_root.iterdir()
                 if d.is_dir() and folder_matches(d.name, token)]
    for root in roots:
        for path in root.rglob("*.ab1"):
            read = parse_ab1_name(path)
            if read is None:
                # Chỉ báo file thuộc mẫu đang quan tâm; cả cây có hàng nghìn
                # file của batch khác, kêu hết thì không ai đọc.
                if any(f"_{sample}_" in path.name for sample in wanted):
                    unparsed.append(path)
                continue
            if read.sample not in wanted:
                continue
            current = best.get(read.key)
            if current is None or read.date > current.date:
                best[read.key] = read
    return best, unparsed


def index_originals(raw_dir: Path, wanted: set[str]) -> tuple[dict[tuple[str, str, str], list[Path]], list[Path]]:
    """Read gốc trong thư mục làm việc, gom theo cùng khoá."""
    found: dict[tuple[str, str, str], list[Path]] = defaultdict(list)
    unparsed: list[Path] = []
    if not raw_dir.is_dir():
        return found, unparsed
    for path in sorted(raw_dir.glob("*.ab1")):
        read = parse_ab1_name(path)
        if read is None:
            if any(f"_{sample}_" in path.name for sample in wanted):
                unparsed.append(path)
            continue
        if read.sample in wanted:
            found[read.key].append(path)
    return found, unparsed


def apply(
    batch_id: str,
    samples: set[str],
    rerun_root: Path,
    raw_dir: Path,
    *,
    dry_run: bool = False,
) -> dict:
    token = batch_token(batch_id)
    replacements, bad_new = collect_rerun_reads(rerun_root, samples, token)
    originals, bad_old = index_originals(raw_dir, samples)

    # Một mẫu có tên file không đọc được thì không thể biết read nào thuộc vùng
    # nào, mà đè nhầm vùng còn tệ hơn không đè. Bỏ qua nguyên mẫu đó và nói rõ.
    skipped: set[str] = set()
    for path in bad_old + bad_new:
        for sample in samples:
            if f"_{sample}_" in path.name:
                skipped.add(sample)
                logger.error(f"Không đọc được tên file, bỏ qua mẫu {sample}: {path.name}")

    replaced: list[dict] = []
    added: list[dict] = []
    removed: list[str] = []

    for key in sorted(replacements):
        sample, region, direction = key
        if sample in skipped:
            continue
        new = replacements[key]
        old_paths = originals.get(key, [])
        record = {
            "sample": sample,
            "region": f"{region}{direction}",
            "new_file": new.path.name,
            "new_date": new.date,
            "replaced_files": [p.name for p in old_paths],
        }
        if not dry_run:
            for old in old_paths:
                old.unlink()
            raw_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(new.path, raw_dir / new.path.name)
        removed.extend(p.name for p in old_paths)
        (replaced if old_paths else added).append(record)

    touched = sorted({r["sample"] for r in replaced + added})
    manifest = {
        "batch_id": batch_id,
        "applied_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "rerun_root": str(rerun_root),
        "samples": touched,
        "replaced": replaced,
        "added": added,
        "skipped_samples": sorted(skipped),
    }

    if not dry_run and (replaced or added):
        # Bước 2 thường sẽ rm -rf thư mục này rồi nạp lại từ Lab/Raw_data, mà ổ
        # đó không có read rerun -- tức là mang đúng read xấu quay về. File kê
        # khai để bước 2 phát hiện và cảnh báo.
        (raw_dir / MANIFEST_NAME).write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    logger.info(
        f"{batch_id}: thay {len(replaced)} read, thêm {len(added)} read, "
        f"xoá {len(removed)} read gốc, trên {len(touched)} mẫu"
        + (f", bỏ qua {len(skipped)} mẫu" if skipped else "")
        + (" (thử, chưa ghi gì)" if dry_run else "")
    )
    return manifest


def main() -> int:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--batch-id", required=True)
    parser.add_argument("--samples", type=Path, required=True,
                        help="file danh sách mã mẫu cần chạy lại, mỗi dòng một mã")
    parser.add_argument("--rerun-root", type=Path, required=True,
                        help="thư mục AB1 lab chạy lại, ví dụ Lab/Raw_data/rerun")
    parser.add_argument("--raw-dir", type=Path, required=True,
                        help="thư mục làm việc của batch, ví dụ science/data/raw/<batch>")
    parser.add_argument("--dry-run", action="store_true", help="chỉ báo, không đụng file")
    parser.add_argument("--report", type=Path, default=None, help="ghi kê khai ra đây (JSON)")
    args = parser.parse_args()

    samples = {
        line.strip()
        for line in args.samples.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    }
    if not samples:
        logger.warning(f"Danh sách mẫu rỗng: {args.samples}")
        return 0
    if not args.rerun_root.is_dir():
        logger.error(f"Không thấy thư mục AB1 rerun: {args.rerun_root}")
        return 1

    manifest = apply(args.batch_id, samples, args.rerun_root, args.raw_dir, dry_run=args.dry_run)

    missing = sorted(samples - set(manifest["samples"]) - set(manifest["skipped_samples"]))
    if missing:
        # Không phải lỗi: mẫu có thể nằm trong ZIP Sequencher mà lab chưa chạy
        # lại read. Khi đó phân tích dùng read gốc, và cần nói rõ ra.
        logger.warning(
            f"{len(missing)} mẫu không có AB1 chạy lại, sẽ dùng read gốc: {', '.join(missing[:10])}"
            + (" …" if len(missing) > 10 else "")
        )
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
