"""Generate called-heteroplasmy reports from NGS comparison artifacts."""

from __future__ import annotations

import json
import re
from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from src.core.variants import parse_variant_position_allele, pos_sort_key
from src.modules.NGS.split_fis_sanger_comparison_sets import load_sets, resolve_set_id

RAW_FIS_COLUMNS = [
    "FIS_Raw_Position",
    "FIS_Raw_Marker",
    "FIS_Raw_Genotype",
    "FIS_Raw_AlleleFrequency",
    "FIS_Raw_TotalDepth",
    "FIS_Raw_Ref(Depth):Alt(Depth)",
    "FIS_Raw_QC_Info",
    "FIS_Raw_Genotype(raw)",
]
OUTPUT_COLUMNS = [
    "FIS_Sample",
    "Sanger_Sample",
    "Position",
    "Heteroplasmy_Source",
    "NGS_Genotype",
    "Sanger_Genotype",
    "True_Heteroplasmy (Y/N/U)",
    "Top_Base",
    "Top_Depth",
    "Second_Base",
    "Second_Depth",
    "Total_Depth",
    "Top_Frequency",
    "Second_Frequency",
    "Group",
    *RAW_FIS_COLUMNS,
]
FIS_COLUMNS = {
    "sample_id_base",
    "Position",
    "Marker",
    "Genotype",
    "AlleleFrequency",
    "TotalDepth",
    "Ref(Depth):Alt(Depth)",
    "QC_Info",
    "Genotype(raw)",
}
MAPPING_COLUMNS = {"FIS_Sample", "Sanger_Sample"}
CORRECTION_COLUMNS = {"CE_Sample_ID", "Sanger_Sample", "ID_Change_Note"}
COMPARISON_COLUMNS = {"Sample ID", "Position", "fis Variant", "Variant Flags (fis)"}
HETEROPLASMY_RE = re.compile(r"(?:^|;\s*)Heteroplasmy at (?P<position>\d+(?:\.\d+)?)(?:;|$)")
DEPTH_RE = re.compile(r"^(?P<first>[^()]+)\((?P<first_depth>\d+)\):(?P<second>[^()]+)\((?P<second_depth>\d+)\)$")


@dataclass(frozen=True)
class DepthCall:
    """The two reported FIS alleles and their depths."""

    first_allele: str
    first_depth: int
    second_allele: str
    second_depth: int
    total_depth: int

    def report_fields(self) -> dict[str, str | int | float]:
        """Return depth fields ordered by descending allele depth."""
        alleles = sorted(
            ((self.first_allele, self.first_depth), (self.second_allele, self.second_depth)),
            key=lambda item: item[1],
            reverse=True,
        )
        (top_base, top_depth), (second_base, second_depth) = alleles
        return {
            "Top_Base": top_base,
            "Top_Depth": top_depth,
            "Second_Base": second_base,
            "Second_Depth": second_depth,
            "Total_Depth": self.total_depth,
            "Top_Frequency": round(top_depth / self.total_depth, 4),
            "Second_Frequency": round(second_depth / self.total_depth, 4),
        }


@dataclass(frozen=True)
class FisRecord:
    """Raw FIS source values and optionally parseable biallelic depths."""

    raw_fields: dict[str, str]
    depth_call: DepthCall | None


def _require_columns(frame: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        msg = f"{path} is missing required columns: {', '.join(missing)}"
        raise ValueError(msg)


def _parse_depth_call(row: pd.Series) -> DepthCall | None:
    depth_text = str(row["Ref(Depth):Alt(Depth)"]).strip()
    match = DEPTH_RE.fullmatch(depth_text)
    if match is None:
        return None
    total_depth = int(str(row["TotalDepth"]))
    if total_depth <= 0:
        return None
    return DepthCall(
        first_allele=match["first"].strip(),
        first_depth=int(match["first_depth"]),
        second_allele=match["second"].strip(),
        second_depth=int(match["second_depth"]),
        total_depth=total_depth,
    )


def _genotype_positions(genotype: str) -> list[str]:
    """Return every position represented in a comma-separated FIS genotype."""
    positions: list[str] = []
    for call in genotype.split(","):
        position, allele = parse_variant_position_allele(call.strip())
        if not position or not allele:
            msg = f"Malformed FIS genotype: {genotype!r}"
            raise ValueError(msg)
        positions.append(position)
    return positions


def _raw_fields(row: pd.Series) -> dict[str, str]:
    return {
        "FIS_Raw_Position": str(row["Position"]),
        "FIS_Raw_Marker": str(row["Marker"]),
        "FIS_Raw_Genotype": str(row["Genotype"]),
        "FIS_Raw_AlleleFrequency": str(row["AlleleFrequency"]),
        "FIS_Raw_TotalDepth": str(row["TotalDepth"]),
        "FIS_Raw_Ref(Depth):Alt(Depth)": str(row["Ref(Depth):Alt(Depth)"]),
        "FIS_Raw_QC_Info": str(row["QC_Info"]),
        "FIS_Raw_Genotype(raw)": str(row["Genotype(raw)"]),
    }


def load_fis_records(fis_file: Path) -> dict[tuple[str, str], FisRecord]:
    """Load one raw FIS record per sample/genotype position."""
    fis = pd.read_csv(fis_file, sep="\t", dtype=str, keep_default_na=False)
    _require_columns(fis, FIS_COLUMNS, fis_file)
    records: dict[tuple[str, str], FisRecord] = {}
    for _, row in fis.iterrows():
        record = FisRecord(raw_fields=_raw_fields(row), depth_call=_parse_depth_call(row))
        for position in dict.fromkeys(_genotype_positions(str(row["Genotype"]))):
            key = (str(row["sample_id_base"]), position)
            if key in records:
                msg = f"Duplicate FIS records for {key[0]} at {key[1]}"
                raise ValueError(msg)
            records[key] = record
    return records


def load_sanger_samples(mapping_file: Path) -> dict[str, str]:
    """Load resolved Sanger sample IDs keyed by FIS sample ID."""
    mapping = pd.read_csv(mapping_file, sep="\t", dtype=str, keep_default_na=False)
    _require_columns(mapping, MAPPING_COLUMNS, mapping_file)
    duplicate_samples = mapping.loc[mapping["FIS_Sample"].duplicated(), "FIS_Sample"].tolist()
    if duplicate_samples:
        msg = f"Duplicate FIS samples in {mapping_file}: {', '.join(duplicate_samples[:5])}"
        raise ValueError(msg)
    return dict(zip(mapping["FIS_Sample"], mapping["Sanger_Sample"], strict=True))


def load_sanger_corrections(correction_file: Path) -> dict[str, str]:
    """Load documented FIS-to-Sanger ID corrections."""
    corrections = pd.read_csv(correction_file, sep="\t", dtype=str, keep_default_na=False)
    _require_columns(corrections, CORRECTION_COLUMNS, correction_file)
    corrections = corrections.loc[corrections["ID_Change_Note"] != ""]
    duplicate_samples = corrections.loc[corrections["CE_Sample_ID"].duplicated(), "CE_Sample_ID"].tolist()
    if duplicate_samples:
        msg = f"Duplicate corrected FIS samples in {correction_file}: {', '.join(duplicate_samples[:5])}"
        raise ValueError(msg)
    return dict(zip(corrections["CE_Sample_ID"], corrections["Sanger_Sample"], strict=True))


def _heteroplasmy_position(flags: str) -> str | None:
    match = HETEROPLASMY_RE.search(flags)
    return match["position"] if match else None


def _variant_is_heteroplasmy(flags: object, position: str) -> bool:
    return isinstance(flags, list) and any(_heteroplasmy_position(str(flag)) == position for flag in flags)


def _validated_sanger_sample_id(fis_sample: str, record: dict[str, Any], corrections: dict[str, str]) -> str:
    """Return the prepared Sanger ID after validating any correction."""
    information = record.get("information")
    sanger_sample = str(information.get("sanger_sample_id") or "") if isinstance(information, dict) else ""
    corrected_sample = corrections.get(fis_sample)
    if corrected_sample and sanger_sample != corrected_sample:
        msg = f"Prepared Sanger record {fis_sample!r} does not use corrected ID {corrected_sample!r}"
        raise ValueError(msg)
    return sanger_sample


def _add_sanger_heteroplasmies(
    fis_sample: str, record: dict[str, Any], heteroplasmies: dict[tuple[str, str], str]
) -> None:
    """Add the record's flagged Sanger heteroplasmies to the shared lookup."""
    variants = record.get("variants")
    variant_flags = record.get("variant_flags")
    if not isinstance(variants, list) or not isinstance(variant_flags, dict):
        return
    for variant in variants:
        if not isinstance(variant, dict):
            continue
        position = str(variant.get("pos", ""))
        ref = str(variant.get("ref", ""))
        allele = str(variant.get("seq", ""))
        flags = variant_flags.get(f"{position}|{ref}|{allele}")
        if not position or not allele or not _variant_is_heteroplasmy(flags, position):
            continue
        key = (fis_sample, position)
        if key in heteroplasmies and heteroplasmies[key] != allele:
            msg = f"Multiple Sanger heteroplasmy calls for {key[0]} at {key[1]}"
            raise ValueError(msg)
        heteroplasmies[key] = allele


def load_prepared_sanger_data(
    sanger_file: Path, corrections: dict[str, str]
) -> tuple[dict[str, str], dict[tuple[str, str], str]]:
    """Load corrected Sanger IDs and heteroplasmy calls from prepared Sanger JSON."""
    data: Any = json.loads(sanger_file.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        msg = f"Expected a JSON object in {sanger_file}"
        raise TypeError(msg)

    sample_ids: dict[str, str] = {}
    heteroplasmies: dict[tuple[str, str], str] = {}
    for raw_fis_sample, raw_record in data.items():
        if not isinstance(raw_record, dict):
            continue
        fis_sample = str(raw_fis_sample)
        sanger_sample = _validated_sanger_sample_id(fis_sample, raw_record, corrections)
        if sanger_sample:
            sample_ids[fis_sample] = sanger_sample
        _add_sanger_heteroplasmies(fis_sample, raw_record, heteroplasmies)
    return sample_ids, heteroplasmies


def load_sample_list(sample_list_file: Path) -> set[str]:
    """Load unique sample IDs from a one-ID-per-line text file."""
    sample_ids = [line.strip() for line in sample_list_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    duplicates = sorted(sample_id for sample_id, count in Counter(sample_ids).items() if count > 1)
    if duplicates:
        msg = f"Duplicate sample IDs in {sample_list_file}: {', '.join(duplicates[:5])}"
        raise ValueError(msg)
    return set(sample_ids)


def _genotype(position: str, allele: str) -> str:
    return "" if allele in {"", "None"} else f"{position}{allele}"


def _group_for_sample(sample_id: str, groups: dict[str, str]) -> str:
    matched_id = resolve_set_id(sample_id, set(groups))
    return groups[matched_id] if matched_id else ""


def _empty_raw_fields() -> dict[str, str]:
    return dict.fromkeys(RAW_FIS_COLUMNS, "")


def _write_report(report: pd.DataFrame, output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(output_file, sep="\t", index=False)


def _report_row(
    comparison_row: pd.Series,
    fis_records: dict[tuple[str, str], FisRecord],
    sanger_samples: dict[str, str],
    sanger_heteroplasmies: dict[tuple[str, str], str],
    groups: dict[str, str],
) -> dict[str, str | int | float] | None:
    """Build one report row from FIS and corrected Sanger calls."""
    sample_id = str(comparison_row["Sample ID"])
    position = str(comparison_row["Position"])
    fis_position = _heteroplasmy_position(str(comparison_row["Variant Flags (fis)"]))
    if fis_position is not None and fis_position != position:
        msg = f"FIS heteroplasmy position does not match comparison position for {sample_id!r}"
        raise ValueError(msg)
    sanger_allele = sanger_heteroplasmies.get((sample_id, position))
    if fis_position is None and sanger_allele is None:
        return None

    source = "Both" if fis_position and sanger_allele else "NGS_Only" if fis_position else "Sanger_Only"
    row: dict[str, str | int | float] = {
        "FIS_Sample": sample_id,
        "Sanger_Sample": sanger_samples.get(sample_id, ""),
        "Position": position,
        "Heteroplasmy_Source": source,
        "NGS_Genotype": _genotype(position, str(comparison_row["fis Variant"])),
        "Sanger_Genotype": _genotype(position, sanger_allele or ""),
        "True_Heteroplasmy (Y/N/U)": "U",
        "Top_Base": "",
        "Top_Depth": "",
        "Second_Base": "",
        "Second_Depth": "",
        "Total_Depth": "",
        "Top_Frequency": "",
        "Second_Frequency": "",
        "Group": _group_for_sample(sample_id, groups),
        **_empty_raw_fields(),
    }
    if fis_position:
        fis_record = fis_records.get((sample_id, position))
        if fis_record is None:
            logger.warning("No matching FIS raw record for {} at {}", sample_id, position)
        else:
            row.update(fis_record.raw_fields)
            if fis_record.depth_call is not None:
                row.update(fis_record.depth_call.report_fields())
    return row


def generate_heteroplasmy_report(
    fis_file: Path,
    comparison_file: Path,
    mapping_file: Path,
    sanger_file: Path,
    correction_file: Path,
    sets_file: Path,
    output_file: Path,
    sample_list_file: Path | None = None,
    filtered_output_file: Path | None = None,
) -> pd.DataFrame:
    """Create a full report and, optionally, a report filtered to listed samples."""
    if (sample_list_file is None) != (filtered_output_file is None):
        msg = "sample_list_file and filtered_output_file must be provided together"
        raise ValueError(msg)

    fis_records = load_fis_records(fis_file)
    corrections = load_sanger_corrections(correction_file)
    sanger_samples = load_sanger_samples(mapping_file)
    prepared_sanger_ids, sanger_heteroplasmies = load_prepared_sanger_data(sanger_file, corrections)
    sanger_samples.update(prepared_sanger_ids)
    groups = load_sets(sets_file)
    comparison = pd.read_csv(comparison_file, sep="\t", dtype=str, keep_default_na=False)
    _require_columns(comparison, COMPARISON_COLUMNS, comparison_file)

    rows: list[dict[str, str | int | float]] = []
    for _, comparison_row in comparison.iterrows():
        row = _report_row(comparison_row, fis_records, sanger_samples, sanger_heteroplasmies, groups)
        if row is not None:
            rows.append(row)

    report = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not report.empty:
        report = report.sort_values(
            by=["FIS_Sample", "Position"],
            key=lambda column: column.map(pos_sort_key) if column.name == "Position" else column,
            kind="stable",
        )
    _write_report(report, output_file)
    logger.success("Saved {} called heteroplasmy rows to {}", len(report), output_file)

    if sample_list_file is not None and filtered_output_file is not None:
        sample_ids = load_sample_list(sample_list_file)
        filtered_report = report.loc[report["FIS_Sample"].isin(sorted(sample_ids))]
        _write_report(filtered_report, filtered_output_file)
        logger.success("Saved {} filtered heteroplasmy rows to {}", len(filtered_report), filtered_output_file)
    return report


def arg_parser() -> ArgumentParser:
    """Build the command-line parser."""
    parser = ArgumentParser(description="Generate called-heteroplasmy TSV reports from NGS comparison artifacts")
    parser.add_argument("-f", "--fis-file", type=Path, required=True, help="Raw FIS.tsv artifact")
    parser.add_argument("-c", "--comparison-file", type=Path, required=True, help="Comparison variant-level TSV")
    parser.add_argument("-m", "--mapping-file", type=Path, required=True, help="Comparison sample mapping TSV")
    parser.add_argument("-r", "--sanger-file", type=Path, required=True, help="Prepared corrected Sanger JSON")
    parser.add_argument(
        "-u", "--correction-file", type=Path, required=True, help="Authoritative FIS-to-Sanger correction TSV"
    )
    parser.add_argument("-s", "--sets-file", type=Path, required=True, help="Headerless Sample_ID/batch/set TSV")
    parser.add_argument("-o", "--output-file", type=Path, required=True, help="Full output heteroplasmy TSV")
    parser.add_argument("--sample-list", type=Path, help="Optional one-sample-ID-per-line filter TXT")
    parser.add_argument("--filtered-output-file", type=Path, help="Output TSV for --sample-list samples")
    return parser


def main() -> None:
    """Run the report generator."""
    args = arg_parser().parse_args()
    generate_heteroplasmy_report(
        args.fis_file,
        args.comparison_file,
        args.mapping_file,
        args.sanger_file,
        args.correction_file,
        args.sets_file,
        args.output_file,
        args.sample_list,
        args.filtered_output_file,
    )


if __name__ == "__main__":
    main()
