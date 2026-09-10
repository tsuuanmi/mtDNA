"""Map canonical FIS Samples to Sanger comparison keys."""

import csv
from argparse import ArgumentParser, Namespace
from collections import Counter
from pathlib import Path
from typing import Any, NamedTuple

import pandas as pd
from loguru import logger

from src.core.models import Tool
from src.core.sample import load_sample_batch
from src.modules.NGS.fis_sanger_common import (
    build_comparison_sample,
    filter_variants_to_intervals,
    load_json_object,
    normalize_intervals,
    write_sample_batch,
)
from src.tools.fis.utils import clean_sample_id

MAPPING_MT_COLUMN = "mtDNA"
MAPPING_STR_COLUMN = "STR"
CORRECTION_COLUMNS = ["CE_Sample_ID", "Sanger_Sample", "Sanger_Batch", "ID_Change_Note"]
MANIFEST_FIELDS = [
    "Comparison_Key",
    "FIS_Sample",
    "FIS_sample_id_base",
    "Sanger_Sample",
    "Sanger_Batch",
    "Mapping_Source",
    "Match_Status",
]


class SangerMapping(NamedTuple):
    """Resolved Sanger target and mapping method."""

    sample_id: str | None
    source: str


def arg_parser() -> Namespace:
    """Parse FIS preparation arguments."""
    parser = ArgumentParser(description="Map canonical FIS samples for shared comparison")
    parser.add_argument("-f", "--fis-json", type=Path, required=True, help="Canonical FIS batch JSON")
    parser.add_argument("-s", "--sanger-json", type=Path, required=True, help="Full Sanger batch JSON")
    parser.add_argument("-t", "--mapping-file", type=Path, required=True, help="MT-to-STR/HID mapping TSV")
    parser.add_argument(
        "-u", "--correction-file", type=Path, required=True, help="Authoritative CE-to-Sanger correction TSV"
    )
    parser.add_argument("-o", "--output-file", type=Path, required=True, help="Prepared FIS batch JSON")
    parser.add_argument("-m", "--manifest-file", type=Path, required=True, help="FIS-to-Sanger manifest TSV")
    parser.add_argument("--batch-id", help="Batch ID stored in prepared Samples")
    return parser.parse_args()


def parse_sample_mapping(mapping_path: Path) -> dict[str, str]:
    """Parse MT and STR/HID aliases into canonical Sanger IDs."""
    frame = pd.read_csv(mapping_path, sep="\t", dtype=str)
    missing = {MAPPING_MT_COLUMN, MAPPING_STR_COLUMN} - set(frame.columns)
    if missing:
        msg = f"Mapping file {mapping_path} is missing required columns: {', '.join(sorted(missing))}"
        raise ValueError(msg)
    mapping: dict[str, str] = {}
    for mt_value, str_value in frame[[MAPPING_MT_COLUMN, MAPPING_STR_COLUMN]].itertuples(index=False):
        mt_id, str_id = clean_sample_id(mt_value), clean_sample_id(str_value)
        if not mt_id or not str_id:
            continue
        for alias in {mt_id, str_id}:
            if alias in mapping and mapping[alias] != mt_id:
                msg = f"Mapping alias {alias!r} points to multiple IDs"
                raise ValueError(msg)
            mapping[alias] = mt_id
    return mapping


def load_sanger_corrections(
    correction_path: Path, sanger_data: dict[str, Any], fis_sample_ids: set[str]
) -> dict[str, str]:
    """Load and validate documented corrections relevant to this FIS batch."""
    frame = pd.read_csv(correction_path, sep="\t", dtype=str, keep_default_na=False)
    missing = set(CORRECTION_COLUMNS) - set(frame.columns)
    if missing:
        msg = f"Correction file {correction_path} is missing required columns: {', '.join(sorted(missing))}"
        raise ValueError(msg)
    frame = frame.loc[frame["ID_Change_Note"] != "", CORRECTION_COLUMNS]
    if frame.empty:
        msg = f"Correction file has no documented changes: {correction_path}"
        raise ValueError(msg)
    if (frame[CORRECTION_COLUMNS[:-1]] == "").any().any():
        msg = f"Correction file has empty documented values: {correction_path}"
        raise ValueError(msg)
    if frame["CE_Sample_ID"].duplicated().any():
        msg = f"Correction file has duplicate CE samples: {correction_path}"
        raise ValueError(msg)
    corrections: dict[str, str] = {}
    relevant = frame.loc[frame["CE_Sample_ID"].isin(fis_sample_ids)]
    for fis_sample, sanger_sample, declared_batch, _ in relevant.itertuples(index=False):
        target = sanger_data.get(sanger_sample)
        if not isinstance(target, dict):
            msg = f"Correction target {sanger_sample!r} is missing from Sanger JSON"
            raise KeyError(msg)
        actual_batch = str(target.get("batch") or "")
        if declared_batch != actual_batch and not declared_batch.startswith(f"{actual_batch}_"):
            msg = f"Correction target {sanger_sample!r} has incompatible batch {declared_batch!r}"
            raise ValueError(msg)
        corrections[fis_sample] = sanger_sample
    return corrections


def resolve_sanger_mapping(
    sample_id: str, base_id: str, sanger_ids: set[str], mapping: dict[str, str], corrections: dict[str, str]
) -> SangerMapping:
    """Resolve documented correction, exact ID, then longest alias."""
    if sample_id in corrections:
        return SangerMapping(corrections[sample_id], "correction")
    clean_id = base_id.replace("_L02", "").replace("_L01", "").strip("_")
    if clean_id in sanger_ids:
        return SangerMapping(clean_id, "exact")
    aliases = [alias for alias in mapping if clean_id == alias or clean_id.startswith(f"{alias}_")]
    return SangerMapping(mapping[max(aliases, key=len)], "alias") if aliases else SangerMapping(None, "unresolved")


def write_manifest(rows: list[dict[str, str]], output_path: Path) -> None:
    """Write the comparison manifest."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare_fis_comparison(
    fis_json: Path,
    sanger_json: Path,
    mapping_file: Path,
    correction_file: Path,
    output_file: Path,
    manifest_file: Path,
    batch_id: str | None = None,
) -> None:
    """Map canonical FIS Samples and filter them to Sanger intervals."""
    fis_samples = load_sample_batch(fis_json)
    sanger_data = load_json_object(sanger_json)
    mapping = parse_sample_mapping(mapping_file)
    fis_ids = {sample.sample_id for sample in fis_samples}
    corrections = load_sanger_corrections(correction_file, sanger_data, fis_ids)
    prepared = []
    rows: list[dict[str, str]] = []
    for sample in fis_samples:
        metadata = (sample.information or {}).get("fis_metadata") or {}
        base_id = str(metadata.get("FIS_sample_id_base") or sample.sample_id)
        resolved = resolve_sanger_mapping(sample.sample_id, base_id, set(sanger_data), mapping, corrections)
        target = sanger_data.get(resolved.sample_id) if resolved.sample_id else None
        found = isinstance(target, dict)
        intervals = normalize_intervals(target.get("intervals")) if found else None
        variants = filter_variants_to_intervals(sample.variants, intervals) if found else list(sample.variants)
        prepared.append(
            build_comparison_sample(
                sample_id=sample.sample_id,
                variants=variants,
                source_tool=Tool.FIS,
                intervals=intervals,
                batch_id=batch_id,
                information={
                    "fis_metadata": metadata,
                    "sanger_sample_id": resolved.sample_id,
                    "mapping_source": resolved.source,
                    "sanger_record_found": found,
                },
            )
        )
        rows.append(
            {
                "Comparison_Key": sample.sample_id,
                "FIS_Sample": sample.sample_id,
                "FIS_sample_id_base": base_id,
                "Sanger_Sample": resolved.sample_id or "",
                "Sanger_Batch": str(target.get("batch", "")) if found else "",
                "Mapping_Source": resolved.source,
                "Match_Status": "Matched" if found else "Not found",
            }
        )
    write_sample_batch(prepared, output_file)
    write_manifest(rows, manifest_file)
    sources = Counter(row["Mapping_Source"] for row in rows)
    logger.success(
        "Prepared FIS samples: total={}, matched={}, corrections={}",
        len(prepared),
        sum(row["Match_Status"] == "Matched" for row in rows),
        sources["correction"],
    )


def main() -> None:
    args = arg_parser()
    prepare_fis_comparison(
        args.fis_json,
        args.sanger_json,
        args.mapping_file,
        args.correction_file,
        args.output_file,
        args.manifest_file,
        args.batch_id,
    )


if __name__ == "__main__":
    main()
