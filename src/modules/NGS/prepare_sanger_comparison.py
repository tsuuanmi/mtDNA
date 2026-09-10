"""Filter and normalize Sanger samples requested by the FIS manifest."""

import csv
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any

from loguru import logger

from src.core.models import Tool, Variant
from src.modules.NGS.fis_sanger_common import (
    build_comparison_sample,
    filter_variants_to_intervals,
    load_json_object,
    normalize_intervals,
    normalize_variant,
    write_sample_batch,
)

REQUIRED_MANIFEST_COLUMNS = {"Comparison_Key", "Sanger_Sample", "Match_Status"}


def arg_parser() -> Namespace:
    """Parse Sanger preparation arguments."""
    parser = ArgumentParser(description="Prepare selected Sanger samples for shared comparison")
    parser.add_argument("-s", "--sanger-json", type=Path, required=True, help="Full Sanger batch JSON")
    parser.add_argument("-m", "--manifest-file", type=Path, required=True, help="FIS-to-Sanger manifest TSV")
    parser.add_argument("-o", "--output-file", type=Path, required=True, help="Prepared Sanger batch JSON")
    parser.add_argument("--batch-id", help="Batch ID stored in prepared Samples")
    return parser.parse_args()


def read_matched_manifest(manifest_path: Path) -> list[dict[str, str]]:
    """Read matched manifest rows and validate their comparison keys."""
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames or [])
        if missing:
            msg = f"Manifest {manifest_path} is missing columns: {', '.join(sorted(missing))}"
            raise ValueError(msg)
        rows = [dict(row) for row in reader if row.get("Match_Status") == "Matched"]

    seen: set[str] = set()
    for row in rows:
        comparison_key = row.get("Comparison_Key", "")
        if not comparison_key:
            msg = f"Manifest {manifest_path} contains an empty comparison key"
            raise ValueError(msg)
        if comparison_key in seen:
            msg = f"Duplicate comparison key in manifest: {comparison_key}"
            raise ValueError(msg)
        seen.add(comparison_key)
    return rows


def _flatten_grouped_variants(raw_record: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten Sanger snp/insertion/deletion groups."""
    grouped = raw_record.get("variants")
    if isinstance(grouped, list):
        return [entry for entry in grouped if isinstance(entry, dict)]
    if not isinstance(grouped, dict):
        return []
    return [
        entry
        for group_name in ("snps", "insertions", "deletions")
        for entry in grouped.get(group_name, [])
        if isinstance(entry, dict)
    ]


def prepare_sanger_comparison(
    sanger_json: Path,
    manifest_file: Path,
    output_file: Path,
    batch_id: str | None = None,
) -> None:
    """Prepare only the Sanger records referenced by matched manifest rows."""
    sanger_data = load_json_object(sanger_json)
    manifest_rows = read_matched_manifest(manifest_file)
    samples = []

    for manifest_row in manifest_rows:
        comparison_key = manifest_row["Comparison_Key"]
        sanger_sample = manifest_row["Sanger_Sample"]
        raw_record = sanger_data.get(sanger_sample)
        if not isinstance(raw_record, dict):
            msg = f"Matched Sanger sample {sanger_sample!r} is missing from {sanger_json}"
            raise KeyError(msg)

        intervals = normalize_intervals(raw_record.get("intervals"))
        variants: list[Variant] = []
        for raw_variant in _flatten_grouped_variants(raw_record):
            variant = normalize_variant(raw_variant, sanger_sample)
            if variant is not None:
                variants.append(variant)
        variants = filter_variants_to_intervals(variants, intervals)

        sanger_batch = str(raw_record.get("batch") or manifest_row.get("Sanger_Batch") or "")
        samples.append(
            build_comparison_sample(
                sample_id=comparison_key,
                variants=variants,
                source_tool=Tool.SANGER,
                intervals=intervals,
                batch_id=batch_id or sanger_batch or None,
                information={
                    "sanger_sample_id": sanger_sample,
                    "sanger_batch": sanger_batch or None,
                    "fis_sample_id": manifest_row.get("FIS_Sample") or comparison_key,
                },
            )
        )

    write_sample_batch(samples, output_file)
    logger.info("Prepared selected Sanger samples: {}", len(samples))


def main() -> None:
    """Run selected-Sanger preparation."""
    args = arg_parser()
    prepare_sanger_comparison(
        sanger_json=args.sanger_json,
        manifest_file=args.manifest_file,
        output_file=args.output_file,
        batch_id=args.batch_id,
    )


if __name__ == "__main__":
    main()
