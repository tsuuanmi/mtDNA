"""Add FIS-specific metadata to standard shared comparison results."""

from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from src.core.models import Sample, Variant
from src.core.sample import load_sample_batch
from src.core.variants import iupac_bases_compatible, iupac_to_bases, pos_base
from src.modules.NGS.fis_sanger_common import load_json_object, normalize_intervals
from src.tools.fis.utils import select_profile

FIS_NOMENCLATURE_COLUMN = "FIS_Nomenclature correction genotype"
FIS_REVIEW_SPECIAL_INDEL_POSITION = 16193
TRANSFORMED_COLUMN_MAP = {
    "Nomenclature correction genotype": FIS_NOMENCLATURE_COLUMN,
    "Nomenclature QC": "FIS_Nomenclature QC",
    "Coverage pass rate(%)": "FIS_Coverage pass rate(%)",
    "Failed position": "FIS_Failed position",
    "# N in consensus seq": "FIS_# N in consensus seq",
}
FINAL_COLUMNS = [
    "FIS_Sample",
    "Sanger_Sample",
    "Sanger_Batch",
    "FIS_Variants",
    *TRANSFORMED_COLUMN_MAP,
    "Analyzed_Intervals",
    "Sanger_Variants",
    "Concordant",
    "FIS_Correct",
    "FIS_Unique",
    "Sanger_Unique",
    "FIS_QC",
    "FIS_Flag",
    "FIS_Flag_Level",
    "FIS_Flag_Info",
    "FIS_Flag_Reasons",
    "Sanger_Flag",
    "Sanger_Flag_Reasons",
]


def arg_parser() -> Namespace:
    """Parse final enrichment arguments."""
    parser = ArgumentParser(description="Enrich shared FIS-Sanger comparison results")
    parser.add_argument("-c", "--comparison-json", type=Path, required=True, help="Shared comparison JSON")
    parser.add_argument("-f", "--fis-json", type=Path, required=True, help="FIS analyzer JSON")
    parser.add_argument(
        "-t",
        "--fis-transformed-tsv",
        type=Path,
        required=True,
        help="Transformed FIS TSV containing nomenclature and coverage columns",
    )
    parser.add_argument("-s", "--sanger-json", type=Path, required=True, help="Prepared Sanger JSON")
    parser.add_argument("-m", "--manifest-file", type=Path, required=True, help="FIS-to-Sanger manifest TSV")
    parser.add_argument("-o", "--output-file", type=Path, required=True, help="Final enriched TSV")
    return parser.parse_args()


def _clean_display(value: object, default: str = "") -> str:
    """Convert missing metadata values to a stable display string."""
    if value is None or str(value).strip().lower() in {"", "nan", "none", "<na>"}:
        return default
    return str(value)


def _fis_metadata_lookup(samples: list[Sample]) -> dict[str, dict[str, Any]]:
    """Index canonical FIS Sample metadata by public sample ID."""
    lookup: dict[str, dict[str, Any]] = {}
    for sample in samples:
        metadata = (sample.information or {}).get("fis_metadata")
        if not isinstance(metadata, dict):
            msg = f"Canonical FIS sample {sample.sample_id!r} has no fis_metadata"
            raise TypeError(msg)
        if sample.sample_id in lookup:
            msg = f"Duplicate FIS metadata key: {sample.sample_id}"
            raise ValueError(msg)
        lookup[sample.sample_id] = metadata
    return lookup


def _transformed_fis_lookup(transformed_tsv: Path) -> dict[str, dict[str, str]]:
    """Load required per-sample nomenclature and coverage fields from transformed FIS TSV."""
    transformed = pd.read_csv(transformed_tsv, sep="\t", dtype=str, keep_default_na=False)
    required_columns = {"FIS_Sample", *TRANSFORMED_COLUMN_MAP.values()}
    missing_columns = required_columns - set(transformed.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        msg = f"Transformed FIS TSV {transformed_tsv} is missing columns: {missing}"
        raise ValueError(msg)
    duplicated = transformed.loc[transformed["FIS_Sample"].duplicated(), "FIS_Sample"].tolist()
    if duplicated:
        msg = f"Transformed FIS TSV contains duplicate samples: {', '.join(duplicated[:5])}"
        raise ValueError(msg)
    return {
        str(row["FIS_Sample"]): {column: str(row[source]) for column, source in TRANSFORMED_COLUMN_MAP.items()}
        for row in transformed.to_dict(orient="records")
    }


def _format_intervals(raw_intervals: object) -> str:
    """Format prepared Sanger intervals like the legacy enriched report."""
    intervals = normalize_intervals(raw_intervals)
    if not intervals:
        return "N/A"
    ranges = [f"{start}-{end}" for region_ranges in intervals.values() for start, end in region_ranges]
    return ", ".join(ranges) if ranges else "None"


def _profile(metadata: dict[str, Any]) -> str:
    """Return the complete corrected FIS profile for display."""
    profile, _ = select_profile(metadata.get(FIS_NOMENCLATURE_COLUMN), _clean_display(metadata.get("FIS_Variants")))
    return _clean_display(profile, "None")


def _split_flag_reasons(value: object) -> list[str]:
    """Split a standard comparison flag cell and remove empty sentinels."""
    text = _clean_display(value)
    if not text or text in {"No", "None"}:
        return []
    return [reason.strip() for reason in text.split("; ") if reason.strip()]


def _sanger_flag_data(comparison_row: dict[str, Any]) -> tuple[str, str]:
    """Combine standard Sanger sample and variant flag columns."""
    reasons = []
    seen: set[str] = set()
    for column in ("Sample Flags (sanger)", "Variant Flags (sanger)"):
        for reason in _split_flag_reasons(comparison_row.get(column)):
            if reason not in seen:
                reasons.append(reason)
                seen.add(reason)
    return ("Y" if reasons else "N", ", ".join(reasons))


def _variant_display_key(variant: Variant) -> str:
    """Return the simplified display token used by core comparison rows."""
    allele = "DEL" if variant.seq == "-" else variant.seq
    return f"{variant.pos}{allele}"


def _is_16193_indel(variant: Variant) -> bool:
    """Return whether a variant is the FIS-review special 16193 indel."""
    return pos_base(variant.pos) == FIS_REVIEW_SPECIAL_INDEL_POSITION and (
        variant.seq == "-" or isinstance(variant.pos, str)
    )


def _is_reference_covering_iupac(variant: Variant) -> bool:
    """Return whether an ambiguous FIS call includes its rCRS reference base."""
    return len(iupac_to_bases(variant.seq)) > 1 and iupac_bases_compatible(variant.ref, variant.seq)


def _unique_tokens(comparison_row: dict[str, Any], column: str) -> list[str]:
    """Split a standard core unique-variant display field."""
    value = _clean_display(comparison_row.get(column))
    return [] if value in {"", "None", "Error"} else value.split()


def _fis_correct(comparison_row: dict[str, Any] | None, variants: list[Variant]) -> str:
    """Refine core concordance for FIS-review special and IUPAC calls."""
    concordant = _clean_display((comparison_row or {}).get("Concordant"), "N/A")
    if concordant != "N" or comparison_row is None:
        return concordant
    if _unique_tokens(comparison_row, "Variants Unique (sanger)"):
        return "N"
    by_display = {_variant_display_key(variant): variant for variant in variants}
    fis_unique = _unique_tokens(comparison_row, "Variants Unique (fis)")
    if fis_unique and all(
        (variant := by_display.get(token)) is not None
        and (_is_16193_indel(variant) or _is_reference_covering_iupac(variant))
        for token in fis_unique
    ):
        return "Y"
    return "N"


def _build_final_row(
    manifest_row: dict[str, str],
    metadata: dict[str, Any],
    fis_variants: list[Variant],
    transformed_row: dict[str, str],
    comparison_row: dict[str, Any] | None,
    sanger_record: dict[str, Any] | None,
) -> dict[str, Any]:
    """Join one standard comparison row to FIS and Sanger metadata."""
    matched = comparison_row is not None and manifest_row["Match_Status"] == "Matched"
    sanger_flag, sanger_reasons = _sanger_flag_data(comparison_row or {})
    resolved_sanger_sample = manifest_row["Sanger_Sample"]
    sanger_sample_display = (
        resolved_sanger_sample
        if matched
        else f"{resolved_sanger_sample} (Not found)"
        if resolved_sanger_sample
        else "Not found"
    )

    return {
        "FIS_Sample": manifest_row["FIS_Sample"],
        "Sanger_Sample": sanger_sample_display,
        "Sanger_Batch": manifest_row["Sanger_Batch"] if matched else "N/A",
        "Analyzed_Intervals": _format_intervals((sanger_record or {}).get("intervals")) if matched else "N/A",
        "FIS_Variants": _profile(metadata),
        **transformed_row,
        "Sanger_Variants": _clean_display((comparison_row or {}).get("Variants (sanger)"), "None"),
        "Concordant": _clean_display((comparison_row or {}).get("Concordant"), "N/A"),
        "FIS_Correct": _fis_correct(comparison_row, fis_variants),
        "FIS_Unique": _clean_display((comparison_row or {}).get("Variants Unique (fis)"), "None"),
        "Sanger_Unique": _clean_display((comparison_row or {}).get("Variants Unique (sanger)"), "None"),
        "FIS_QC": _clean_display(metadata.get("FIS_QC")),
        "FIS_Flag": _clean_display(metadata.get("FIS_Flag"), "N"),
        "FIS_Flag_Level": metadata.get("FIS_Flag_Level", 0),
        "FIS_Flag_Info": _clean_display(metadata.get("FIS_Flag_Info")),
        "FIS_Flag_Reasons": _clean_display(metadata.get("FIS_Flag_Reasons")),
        "Sanger_Flag": sanger_flag if matched else "N",
        "Sanger_Flag_Reasons": sanger_reasons if matched else "",
    }


def enrich_fis_sanger_comparison(
    comparison_json: Path,
    fis_json: Path,
    fis_transformed_tsv: Path,
    sanger_json: Path,
    manifest_file: Path,
    output_file: Path,
) -> None:
    """Create the final FIS-enriched TSV and Excel comparison reports."""
    comparison_data = load_json_object(comparison_json)
    sanger_data = load_json_object(sanger_json)
    manifest = pd.read_csv(manifest_file, sep="\t", dtype=str, keep_default_na=False)
    fis_samples = load_sample_batch(fis_json)
    metadata_lookup = _fis_metadata_lookup(fis_samples)
    variants_lookup = {sample.sample_id: sample.variants for sample in fis_samples}
    transformed_lookup = _transformed_fis_lookup(fis_transformed_tsv)

    summary_rows_raw = comparison_data.get("summary_rows", [])
    if not isinstance(summary_rows_raw, list):
        msg = f"Comparison JSON has invalid summary_rows: {comparison_json}"
        raise TypeError(msg)
    comparison_lookup = {
        str(row["Sample ID"]): row for row in summary_rows_raw if isinstance(row, dict) and row.get("Sample ID")
    }

    final_rows = []
    for manifest_row_raw in manifest.to_dict(orient="records"):
        manifest_row = {column: str(value) for column, value in manifest_row_raw.items()}
        comparison_key = manifest_row["Comparison_Key"]
        fis_sample = manifest_row["FIS_Sample"]
        metadata = metadata_lookup.get(fis_sample, {})
        transformed_row = transformed_lookup.get(fis_sample)
        if transformed_row is None:
            msg = f"FIS sample {fis_sample!r} is missing from {fis_transformed_tsv}"
            raise KeyError(msg)
        comparison_row = comparison_lookup.get(comparison_key)
        sanger_record_raw = sanger_data.get(comparison_key)
        sanger_record = sanger_record_raw if isinstance(sanger_record_raw, dict) else None
        final_rows.append(
            _build_final_row(
                manifest_row,
                metadata,
                variants_lookup.get(fis_sample, []),
                transformed_row,
                comparison_row,
                sanger_record,
            )
        )

    result = pd.DataFrame(final_rows, columns=FINAL_COLUMNS)
    result["_sort_key"] = result["Sanger_Sample"].apply(
        lambda sample: "zzz_not_found" if "Not found" in str(sample) else str(sample)
    )
    result = result.sort_values("_sort_key").drop(columns="_sort_key").reset_index(drop=True)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    excel_file = output_file.with_suffix(".xlsx")
    result.to_csv(output_file, sep="\t", index=False)
    result.to_excel(excel_file, index=False)
    logger.success("Saved enriched comparison: {} and {} ({} rows)", output_file, excel_file, len(result))


def main() -> None:
    """Run final FIS comparison enrichment."""
    args = arg_parser()
    enrich_fis_sanger_comparison(
        comparison_json=args.comparison_json,
        fis_json=args.fis_json,
        fis_transformed_tsv=args.fis_transformed_tsv,
        sanger_json=args.sanger_json,
        manifest_file=args.manifest_file,
        output_file=args.output_file,
    )


if __name__ == "__main__":
    main()
