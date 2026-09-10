"""FIS Excel workbooks to canonical :class:`~src.core.models.Sample` objects."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO

from src.config import get_settings
from src.core.flagging import SampleFlagger, deduplicate_sample_flags, flag_variants
from src.core.models import Sample, Tool, Variant, validate_position
from src.core.variants import parse_variant_position_allele, pos_base
from src.tools.fis.utils import (
    clean_sample_id,
    clean_value,
    format_genotypes,
    json_safe_value,
    raw_variant_row,
    select_profile,
)

R1_SHEET = "T1.BasicStatistics"
R2_SUMMARY_SHEET = "T1.MT.Summary.QC"
T2_SHEET = "T2.Genotype"
T3_SHEET = "T3.Seq"
SAMPLE_COLUMN = "SampleCode"
QC_COLUMN = "SampleQC"
NOMENCLATURE_COLUMN = "Nomenclature correction genotype"
COVERAGE_COLUMN = "Coverage pass rate(%)"
R2_RESULT_COLUMNS = (
    NOMENCLATURE_COLUMN,
    "Nomenclature QC",
    COVERAGE_COLUMN,
    "Failed position",
    "# N in consensus seq",
    "Info",
)
R2_OUTPUT_COLUMNS = {column: f"FIS_{column}" for column in R2_RESULT_COLUMNS}
MAX_COVERAGE_PERCENT = 100


@dataclass
class FisRun:
    """Parsed FIS workbooks reused by canonical and NGS artifact writers."""

    t2: pd.DataFrame
    t3: pd.DataFrame | None
    qc_fail_samples: list[str]
    r2_results: dict[str, dict[str, object | None]]
    qc_flags: dict[str, str]


def load_reference_bases(reference_path: Path) -> dict[int, str]:
    """Load rCRS bases keyed by one-based position."""
    record = next(SeqIO.parse(reference_path, "fasta"), None)
    if record is None:
        msg = f"Reference FASTA is empty: {reference_path}"
        raise ValueError(msg)
    return {index + 1: str(base).upper() for index, base in enumerate(record.seq)}


def _build_qc_info(row: pd.Series) -> str:  # noqa: C901
    """Build legacy FIS quality information from R2 summary fields."""
    coverage = float(str(clean_value(row.get(COVERAGE_COLUMN)) or MAX_COVERAGE_PERCENT))
    info = str(clean_value(row.get("Info")) or "/")
    failed = str(clean_value(row.get("Failed position")) or "/")
    rare = str(clean_value(row.get("Rare mutations")) or "/")
    n_count = int(float(str(clean_value(row.get("# N in consensus seq")) or 0)))
    sample_qc = str(clean_value(row.get("Sample QC")) or "Pass")
    parts: list[str] = []
    if coverage < MAX_COVERAGE_PERCENT:
        parts.append("Low depth")
        if failed != "/":
            parts.append(f"Failed position ({failed})")
        if info != "/":
            parts.extend(tag.strip().capitalize() for tag in info.split(";") if tag.strip().lower() != "low depth")
        elif n_count > 0:
            parts.append(f"Possible wrong consensus ({int(n_count)}N)")
    elif info != "/":
        parts.extend(tag.strip().capitalize() for tag in info.split(";") if tag.strip())
        if failed != "/":
            parts.append(f"Failed position ({failed})")
    else:
        if failed != "/":
            parts.append(f"Failed position ({failed})")
        if rare != "/":
            parts.append(f"Rare mutations ({rare})")
        if n_count > 0:
            parts.append(f"Possible wrong consensus ({int(n_count)}N)")
    if sample_qc not in {"Pass", "Check"}:
        parts.append(f"QC: {sample_qc}")
    return "; ".join(parts)


def load_run(r1_path: Path, r2_path: Path) -> FisRun:
    """Read the FIS source workbooks once."""
    r1 = pd.read_excel(r1_path, sheet_name=R1_SHEET)
    r1[SAMPLE_COLUMN] = r1[SAMPLE_COLUMN].map(clean_sample_id)
    qc_fail = r1.loc[r1[QC_COLUMN] == "QCFail", SAMPLE_COLUMN].dropna().tolist()
    summary = pd.read_excel(r2_path, sheet_name=R2_SUMMARY_SHEET)
    summary[SAMPLE_COLUMN] = summary[SAMPLE_COLUMN].map(clean_sample_id)
    r2_results: dict[str, dict[str, object | None]] = {}
    qc_flags: dict[str, str] = {}
    for _, row in summary.iterrows():
        sample = clean_sample_id(row.get(SAMPLE_COLUMN))
        if sample is None:
            continue
        if sample in r2_results:
            msg = f"Duplicate FIS R2 summary sample: {sample!r}"
            raise ValueError(msg)
        values = {column: clean_value(row.get(column)) for column in R2_RESULT_COLUMNS}
        corrected = values[NOMENCLATURE_COLUMN]
        if corrected is not None:
            values[NOMENCLATURE_COLUMN] = format_genotypes([str(corrected)])
        r2_results[sample] = values
        qc_flags[sample] = _build_qc_info(row)
    sheet_names = set(pd.ExcelFile(r2_path).sheet_names)
    t2 = pd.read_excel(r2_path, sheet_name=T2_SHEET)
    t2[SAMPLE_COLUMN] = t2[SAMPLE_COLUMN].map(clean_sample_id)
    t3 = pd.read_excel(r2_path, sheet_name=T3_SHEET) if T3_SHEET in sheet_names else None
    if t3 is not None:
        t3[SAMPLE_COLUMN] = t3[SAMPLE_COLUMN].map(clean_sample_id)
    return FisRun(
        t2=t2,
        t3=t3,
        qc_fail_samples=qc_fail,
        r2_results=r2_results,
        qc_flags=qc_flags,
    )


def profile_to_variants(
    profile: str, reference_bases: dict[int, str], *, normalize_alleles: bool = True
) -> list[Variant]:
    """Convert a selected FIS profile to reference-backed variants.

    Canonical variants use uppercase alleles. Flagging of a raw fallback profile
    may retain its original case to report low-confidence lowercase calls.
    """
    variants: list[Variant] = []
    seen: set[tuple[int | str, str, str]] = set()
    for genotype in profile.split():
        position_text, allele_text = parse_variant_position_allele(genotype)
        allele = allele_text.upper() if normalize_alleles else allele_text
        if not position_text or not allele:
            msg = f"Malformed FIS genotype: {genotype!r}"
            raise ValueError(msg)
        if allele.upper().endswith("N"):
            continue
        try:
            position = validate_position(position_text)
        except (TypeError, ValueError) as error:
            msg = f"Invalid FIS genotype position: {genotype!r}"
            raise ValueError(msg) from error
        base = pos_base(position)
        if allele.upper() == "DEL":
            ref, seq = reference_bases.get(base), "-"
        elif isinstance(position, str) and "." in position:
            ref, seq = "-", allele
        else:
            ref, seq = reference_bases.get(base), allele
        if ref is None:
            msg = f"FIS genotype position is outside the rCRS reference: {genotype!r}"
            raise ValueError(msg)
        key = (position, ref, seq)
        if key not in seen:
            variants.append(Variant(pos=position, ref=ref, seq=seq))
            seen.add(key)
    return variants


def _raw_profile(sample_rows: pd.DataFrame) -> str:
    """Format the raw T2 genotype calls for one sample."""
    values = sample_rows["Genotype"].dropna().astype(str).tolist()
    return format_genotypes(values)


def _raw_qc_flags(sample_rows: pd.DataFrame) -> str:
    """Format low-depth T2 calls retained in the NGS artifact."""
    flags: list[str] = []
    for _, row in sample_rows.iterrows():
        qc_info = row.get("QC_Info")
        if qc_info is None or "Low Depth" not in str(qc_info):
            continue
        position = str(row.get("Position", ""))
        base = position.split(":", maxsplit=1)[-1]
        depth = clean_value(row.get("Ref(Depth):Alt(Depth)"))
        flags.append(f"{base}_{depth}" if depth not in {None, "-"} else f"{base}_LD")
    return " ".join(flags)


def _full_intervals() -> dict[str, list[list[int]]]:
    """Return configured full HV intervals for FIS calls."""
    return {region: [list(span)] for region, span in get_settings().regions.REGIONS.items()}


def build_samples(run: FisRun, reference_path: Path, batch_id: str | None = None) -> dict[str, Sample]:
    """Build one canonical FIS Sample per non-QC-fail T2 sample."""
    reference_bases = load_reference_bases(reference_path)
    samples: dict[str, Sample] = {}
    for sample_id, sample_rows in run.t2.groupby(SAMPLE_COLUMN, sort=False):
        sample = str(sample_id)
        if sample in run.qc_fail_samples:
            continue
        sorted_rows = sample_rows.sort_values("Position")
        raw_profile = _raw_profile(sorted_rows)
        r2_result = run.r2_results.get(sample, {})
        selected_profile, source = select_profile(r2_result.get(NOMENCLATURE_COLUMN), raw_profile)
        variants = profile_to_variants(selected_profile, reference_bases)
        intervals = _full_intervals()
        variant_flags = flag_variants(variants)
        flagging_variants = profile_to_variants(
            selected_profile, reference_bases, normalize_alleles=source != "raw_fallback"
        )
        flagger = SampleFlagger(flagging_variants, intervals=intervals)
        _, reasons, _ = flagger.analyze()
        reasons = deduplicate_sample_flags(reasons, variant_flags)
        qc_info = run.qc_flags.get(sample, "")
        flag_status = "Y" if reasons or qc_info else "N"
        flag_level = flagger.classify(flag_status, qc_info, reasons, selected_profile)
        metadata: dict[str, object | None] = {
            "FIS_Sample": sample,
            "FIS_sample_id_base": sample,
            **{R2_OUTPUT_COLUMNS[column]: json_safe_value(r2_result.get(column)) for column in R2_RESULT_COLUMNS},
            "FIS_Variants": raw_profile,
            "FIS_Flag_Profile": selected_profile,
            "FIS_Flag_Profile_Source": source,
            "FIS_QC": _raw_qc_flags(sorted_rows),
            "FIS_Flag": flag_status,
            "FIS_Flag_Level": flag_level,
            "FIS_Flag_Info": qc_info,
            "FIS_Flag_Reasons": ", ".join(reasons),
        }
        samples[sample] = Sample(
            sample_id=sample,
            variants=variants,
            source_tool=Tool.FIS,
            intervals=intervals,
            sample_flags=[qc_info] if qc_info else [],
            variant_flags=variant_flags,
            batch_id=batch_id,
            information={"fis_metadata": metadata},
        )
    return samples


def raw_records(run: FisRun, samples: dict[str, Sample]) -> dict[str, dict[str, Any]]:
    """Return the historical raw FIS JSON shape from the shared parsed run."""
    records: dict[str, dict[str, Any]] = {}
    for sample_id, sample in samples.items():
        rows = run.t2.loc[run.t2[SAMPLE_COLUMN].astype(str) == sample_id]
        metadata = (sample.information or {}).get("fis_metadata", {})
        records[sample_id] = {
            "sample_code": sample_id,
            "variants": [raw_variant_row(row) for _, row in rows.iterrows()],
            "fis_metadata": metadata,
        }
    return records
