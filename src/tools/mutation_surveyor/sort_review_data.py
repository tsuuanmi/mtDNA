"""Post-step AB1 sorting by no-triggering-flag/triggering-flag classification.

Reads the combined final profiles to classify each LID by region review status,
then copies AB1 trace files into categorized folders for manual review workflow.

Two output folders:
    no_triggering_flag/AB1/  — AB1 traces from autopassing regions
        HV1F/         — HV1F forward traces from HV1-autopassing samples
        HV1R/         — HV1R reverse traces from HV1-autopassing samples
        HV2F/         — HV2F forward traces from HV2+HV3-autopassing samples
        HV3R/         — HV3R reverse traces from HV2+HV3-autopassing samples
    triggering_flag/AB1/     — AB1 traces from triggering-flag regions
        HV1F/         — HV1F traces from samples where HV1 has triggering flags
        HV1R/         — HV1R traces from samples where HV1 has triggering flags
        HV2F/         — HV2F traces from samples where HV2+HV3 has triggering flags
        HV3R/         — HV3R traces from samples where HV2+HV3 has triggering flags

    no_triggering_flag/AB1/ contains traces for ALL autopassing regions:
    - Samples passing both HV1 and HV2+HV3 have all 4 primer traces.
    - Samples passing only HV1 have HV1F+HV1R traces.
    - Samples passing only HV2+HV3 have HV2F+HV3R traces.

Classification uses HV1 Review Class and HV2-3 Review Class columns:
    "Auto..." (no "Review") → no triggering flag for that region
    "Review"                → has triggering flag for that region

Output includes an Excel file (no_triggering_flag_classification.xlsx) with:
    - Classification sheet: per-LID pass/fail status for HV1, HV2+HV3, and both regions
    - Summary sheet: aggregate counts of LIDs in each classification category
"""

import shutil
from enum import Enum
from pathlib import Path

import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field

from src.tools.mutation_surveyor.utils import is_control_lid, safe_str

# ---------------------------------------------------------------------------
# Primer type constants
# ---------------------------------------------------------------------------

PRIMER_TYPES = ("HV1F", "HV1R", "HV2F", "HV3R")

# Mapping from region group to individual primer types
REGION_TO_PRIMERS: dict[str, list[str]] = {
    "HV1": ["HV1F", "HV1R"],
    "HV2-3": ["HV2F", "HV3R"],
}


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


class ReviewCategory(str, Enum):  # noqa: UP042
    """Classification category for a sample based on region review status."""

    NO_TRIGGERING_FLAG = "no_triggering_flag"
    TRIGGERING_FLAG = "triggering_flag"


class LIDClassification(BaseModel):
    """Per-LID classification result."""

    lid: str = Field(..., description="Lab ID")
    hv1_no_triggering: bool = Field(..., description="Whether HV1 region has no triggering flags")
    hv23_no_triggering: bool = Field(..., description="Whether HV2+HV3 region has no triggering flags")
    category: ReviewCategory = Field(..., description="Overall classification category")


def _is_no_triggering(review_class: str, region: str = "") -> bool:
    """Check if a review class string indicates no triggering flag.

    No triggering flag: string contains "Auto" and does NOT contain "Review".
    Matches the logic in flagging.py build_region_status_flags().

    Args:
        review_class: Review class string (e.g. "Auto - no review-triggering flag", "Review").
        region: Optional region label for logging.

    Returns:
        True if no triggering flag, False if review needed.
    """
    s = safe_str(review_class)
    if not s:
        logger.warning("Empty review class for {} — treating as review", region or "unknown")
        return False
    return "Auto" in s and "Review" not in s


def classify_lids(df: pd.DataFrame) -> dict[str, LIDClassification]:
    """Classify each LID in the combined profiles DataFrame.

    Reads HV1 Review Class and HV2-3 Review Class columns to determine
    no-triggering-flag/triggering-flag status per region, then assigns a ReviewCategory.

    Args:
        df: Combined profiles LID DataFrame (from final_profiles xlsx).

    Returns:
        Dict mapping LID → LIDClassification.
    """
    classifications: dict[str, LIDClassification] = {}

    for _, row in df.iterrows():
        lid = safe_str(row.get("LID", ""))
        if not lid or is_control_lid(lid):
            continue

        hv1_review = safe_str(row.get("HV1 Review Class", ""))
        hv23_review = safe_str(row.get("HV2-3 Review Class", ""))

        hv1_no_triggering = _is_no_triggering(hv1_review, "HV1")
        hv23_no_triggering = _is_no_triggering(hv23_review, "HV2-3")

        if hv1_no_triggering and hv23_no_triggering:
            category = ReviewCategory.NO_TRIGGERING_FLAG
        else:
            category = ReviewCategory.TRIGGERING_FLAG

        classifications[lid] = LIDClassification(
            lid=lid,
            hv1_no_triggering=hv1_no_triggering,
            hv23_no_triggering=hv23_no_triggering,
            category=category,
        )

    logger.info(
        "Classified {} LIDs: {} no_triggering_flag, {} triggering_flag",
        len(classifications),
        sum(1 for c in classifications.values() if c.category == ReviewCategory.NO_TRIGGERING_FLAG),
        sum(1 for c in classifications.values() if c.category == ReviewCategory.TRIGGERING_FLAG),
    )

    return classifications


# ---------------------------------------------------------------------------
# Trace filename discovery
# ---------------------------------------------------------------------------


def _classify_primer(primer: str) -> str | None:
    """Classify a primer string into one of the known primer types.

    Returns the canonical primer key (HV1F, HV1R, HV2F, HV3R) or None if unknown.
    """
    if primer in PRIMER_TYPES:
        return primer
    if primer.startswith(("HV1", "HV2", "HV3")):
        # Could be HV1F, HV1R, HV2F, HV3R or other HV-prefixed strings
        return primer if primer in PRIMER_TYPES else None
    return None


def _collect_trace_filenames(traces_df: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    """Build per-LID trace filename lookup from Traces xlsx.

    Groups Sample Name by LID and individual primer type
    (HV1F, HV1R, HV2F, HV3R).

    Args:
        traces_df: Traces DataFrame (from HV1 or HV2-3 merged xlsx Traces sheet).

    Returns:
        Dict mapping LID → {"HV1F": [filenames...], "HV1R": [...],
                             "HV2F": [...], "HV3R": [...]}.
    """
    result: dict[str, dict[str, list[str]]] = {}

    for _, row in traces_df.iterrows():
        lid = safe_str(row.get("LID", ""))
        if not lid or is_control_lid(lid):
            continue

        sample_name = safe_str(row.get("Sample Name", ""))
        if not sample_name:
            continue

        primer = safe_str(row.get("Primer", ""))
        primer_key = _classify_primer(primer)
        if primer_key is None:
            logger.warning("Unknown primer '{}' for LID {}, skipping", primer, lid)
            continue

        if lid not in result:
            result[lid] = {p: [] for p in PRIMER_TYPES}
        result[lid][primer_key].append(sample_name)

    # Deduplicate filenames per LID per primer
    for lid_data in result.values():
        for primer in PRIMER_TYPES:
            seen: set[str] = set()
            unique: list[str] = []
            for fn in lid_data[primer]:
                if fn not in seen:
                    seen.add(fn)
                    unique.append(fn)
            lid_data[primer] = unique

    return result


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------


def _write_classification_excel(
    classifications: dict[str, LIDClassification],
    output_path: Path,
) -> None:
    """Write classification results to an Excel file with two sheets.

    Sheets:
        Classification: per-LID pass status for HV1, HV2+HV3, and both regions.
        Summary: aggregate counts for each classification category.

    Args:
        classifications: Dict mapping LID → LIDClassification.
        output_path: Path to write the xlsx file.
    """
    # Build Classification sheet data
    rows = []
    for lid, cls in sorted(classifications.items()):
        rows.append({
            "LID": lid,
            "pass HV1": cls.hv1_no_triggering,
            "pass HV2+HV3": cls.hv23_no_triggering,
            "pass both": cls.hv1_no_triggering and cls.hv23_no_triggering,
        })
    classification_df = pd.DataFrame(rows, columns=["LID", "pass HV1", "pass HV2+HV3", "pass both"])

    # Build Summary sheet data
    total = len(classifications)
    pass_hv1 = sum(1 for c in classifications.values() if c.hv1_no_triggering)
    pass_hv23 = sum(1 for c in classifications.values() if c.hv23_no_triggering)
    pass_both = sum(1 for c in classifications.values() if c.hv1_no_triggering and c.hv23_no_triggering)
    no_triggering_flag = sum(1 for c in classifications.values() if c.category == ReviewCategory.NO_TRIGGERING_FLAG)
    triggering_flag = sum(1 for c in classifications.values() if c.category == ReviewCategory.TRIGGERING_FLAG)

    summary_rows = [
        {"Metric": "Total LIDs", "Count": total},
        {"Metric": "Pass HV1", "Count": pass_hv1},
        {"Metric": "Pass HV2+HV3", "Count": pass_hv23},
        {"Metric": "Pass both", "Count": pass_both},
        {"Metric": "No triggering flag (pass both)", "Count": no_triggering_flag},
        {"Metric": "Triggering flag (fail one or both)", "Count": triggering_flag},
    ]
    summary_df = pd.DataFrame(summary_rows, columns=["Metric", "Count"])

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        classification_df.to_excel(writer, sheet_name="Classification", index=False)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)

    logger.info("Wrote classification Excel to {} ({} LIDs)", output_path, total)


class SortResult(BaseModel):
    """Result of AB1 sort operation."""

    no_triggering_flag_lids: list[str] = Field(
        default_factory=list,
        description="LIDs with no triggering flags in both regions",
    )
    triggering_flag_lids: list[str] = Field(
        default_factory=list,
        description="LIDs with triggering flags in one or both regions",
    )
    files_copied: int = Field(0, description="Number of AB1 files successfully copied")
    files_missing: int = Field(0, description="Number of AB1 files not found in raw dir")
    output_dir: str = Field("", description="Output directory path")


# ---------------------------------------------------------------------------
# Helper: copy traces for a primer type
# ---------------------------------------------------------------------------


def _copy_primer_traces(
    lid: str,
    primer: str,
    trace_filenames: dict[str, dict[str, list[str]]],
    raw_dir: Path,
    dest_dir: Path,
) -> int:
    """Copy AB1 traces for a specific LID and primer type to a destination directory.

    Args:
        lid: Lab ID.
        primer: Primer type (HV1F, HV1R, HV2F, or HV3R).
        trace_filenames: Per-LID trace filename lookup.
        raw_dir: Raw AB1 source directory.
        dest_dir: Destination directory for the traces.

    Returns:
        Number of files successfully copied.
    """
    copied = 0
    for fn in trace_filenames.get(lid, {}).get(primer, []):
        src = raw_dir / fn
        dst = dest_dir / fn
        if src.exists():
            shutil.copy2(src, dst)
            copied += 1
        else:
            logger.warning("{} trace not found in raw dir: {}", primer, src)
    return copied


# ---------------------------------------------------------------------------
# Main sort function
# ---------------------------------------------------------------------------


def sort_review_data(  # noqa: C901, PLR0912, PLR0915
    batch_id: str,
    results_dir: str | Path,
    raw_dir: str | Path,
    excels_dir: str | Path,
    output_dir: str | Path | None = None,
) -> SortResult:
    """Sort AB1 files by no-triggering-flag/triggering-flag classification.

    Reads combined final profiles to classify LIDs, then copies AB1 trace files
    into categorized folders organized by primer type:
        no_triggering_flag/AB1/HV1F/   — HV1F traces from no-triggering-flag samples
        no_triggering_flag/AB1/HV1R/   — HV1R traces from no-triggering-flag samples
        no_triggering_flag/AB1/HV2F/   — HV2F traces from no-triggering-flag samples
        no_triggering_flag/AB1/HV3R/   — HV3R traces from no-triggering-flag samples
        triggering_flag/AB1/HV1F/      — HV1F traces where HV1 has triggering flags
        triggering_flag/AB1/HV1R/      — HV1R traces where HV1 has triggering flags
        triggering_flag/AB1/HV2F/      — HV2F traces where HV2+HV3 has triggering flags
        triggering_flag/AB1/HV3R/      — HV3R traces where HV2+HV3 has triggering flags

    Samples failing both regions have traces in both regional subfolders.

    Args:
        batch_id: Batch identifier (e.g. MS_300326_004).
        raw_dir: Raw AB1 directory (e.g. DATA_DIR/raw/MS_300326_004).
        excels_dir: Directory containing xlsx files.
        output_dir: Output directory. Defaults to results_dir/data/.

    Returns:
        SortResult with classification counts and copy statistics.
    """
    results_dir = Path(results_dir)
    raw_dir = Path(raw_dir)
    if output_dir is None:
        output_dir = results_dir
    output_dir = Path(output_dir)

    excels_dir = Path(excels_dir)

    # --- Read combined final profiles for classification ---
    lid_xlsx = excels_dir / f"{batch_id}_final_profiles.xlsx"
    lid_df = pd.read_excel(lid_xlsx, sheet_name="LID")

    lid_df.columns = [str(c).strip() for c in lid_df.columns]
    classifications = classify_lids(lid_df)

    # --- Read Traces xlsx for per-LID filenames ---
    # Must read BOTH HV23 and HV1 merged files to get all trace filenames.
    # HV23_merged has HV2F/HV3R traces; HV1_merged has HV1F/HV1R traces.
    traces_dfs: list[pd.DataFrame] = []
    for xlsx_name in [f"{batch_id}_HV23_merged.xlsx", f"{batch_id}_HV1_merged.xlsx"]:
        xlsx_path = excels_dir / xlsx_name
        if xlsx_path.exists():
            df = pd.read_excel(xlsx_path, sheet_name="Traces")
            df.columns = [str(c).strip() for c in df.columns]
            traces_dfs.append(df)

    if not traces_dfs:
        msg = f"No Traces files found in {excels_dir}"
        raise FileNotFoundError(msg)

    traces_df = pd.concat(traces_dfs, ignore_index=True)
    trace_filenames = _collect_trace_filenames(traces_df)

    # --- Create output directories ---
    primer_dirs: dict[str, dict[str, Path]] = {}
    for category in ("no_triggering_flag", "triggering_flag"):
        primer_dirs[category] = {}
        for primer in PRIMER_TYPES:
            d = output_dir / category / "AB1" / primer
            d.mkdir(parents=True, exist_ok=True)
            primer_dirs[category][primer] = d

    no_triggering_flag_dir = output_dir / "no_triggering_flag"

    # --- Classify and copy ---
    no_triggering_flag_lids: list[str] = []
    triggering_flag_lids: list[str] = []
    files_copied = 0
    files_missing = 0

    for lid, cls in classifications.items():
        if cls.category == ReviewCategory.NO_TRIGGERING_FLAG:
            no_triggering_flag_lids.append(lid)
        elif cls.category == ReviewCategory.TRIGGERING_FLAG:
            triggering_flag_lids.append(lid)

        # Copy AB1 traces for autopassing regions to no_triggering_flag/
        if cls.hv1_no_triggering:
            for primer in REGION_TO_PRIMERS["HV1"]:
                n = _copy_primer_traces(
                    lid,
                    primer,
                    trace_filenames,
                    raw_dir,
                    primer_dirs["no_triggering_flag"][primer],
                )
                files_copied += n
                files_missing += len(trace_filenames.get(lid, {}).get(primer, [])) - n
        else:
            # HV1 has triggering flags — copy to triggering_flag/
            for primer in REGION_TO_PRIMERS["HV1"]:
                n = _copy_primer_traces(
                    lid,
                    primer,
                    trace_filenames,
                    raw_dir,
                    primer_dirs["triggering_flag"][primer],
                )
                files_copied += n
                files_missing += len(trace_filenames.get(lid, {}).get(primer, [])) - n

        if cls.hv23_no_triggering:
            for primer in REGION_TO_PRIMERS["HV2-3"]:
                n = _copy_primer_traces(
                    lid,
                    primer,
                    trace_filenames,
                    raw_dir,
                    primer_dirs["no_triggering_flag"][primer],
                )
                files_copied += n
                files_missing += len(trace_filenames.get(lid, {}).get(primer, [])) - n
        else:
            # HV2+HV3 has triggering flags — copy to triggering_flag/
            for primer in REGION_TO_PRIMERS["HV2-3"]:
                n = _copy_primer_traces(
                    lid,
                    primer,
                    trace_filenames,
                    raw_dir,
                    primer_dirs["triggering_flag"][primer],
                )
                files_copied += n
                files_missing += len(trace_filenames.get(lid, {}).get(primer, [])) - n

    # --- Write classification Excel ---
    classification_xlsx = no_triggering_flag_dir / "no_triggering_flag_classification.xlsx"
    _write_classification_excel(classifications, classification_xlsx)

    logger.success(
        "AB1 sort complete: {} files copied, {} missing | no_triggering_flag={}, triggering_flag={}",
        files_copied,
        files_missing,
        len(no_triggering_flag_lids),
        len(triggering_flag_lids),
    )

    return SortResult(
        no_triggering_flag_lids=sorted(no_triggering_flag_lids),
        triggering_flag_lids=sorted(triggering_flag_lids),
        files_copied=files_copied,
        files_missing=files_missing,
        output_dir=str(output_dir),
    )
