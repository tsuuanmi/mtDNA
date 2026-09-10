"""Sequencher quality control — filename and range validation plus the QC gate.

This module is the single source of truth for Sequencher filename parsing and
analysis-range validation, and owns the pre-ETL pass/fail QC gate used by the
pipeline. Region and interval helpers (``get_intervals_from_manual``,
``is_full_region``, region-group detection) remain in ``utils.py``.

QC contract:
  - ``run_qc(filename, file_path=None)`` returns a ``QCResult`` (``passed``, ``reason``).
  - Hard fail (no JSON output): no extractable LID, missing Sequencher variant
    table header when ``file_path`` is provided, **or**
    ``validate_analysis_ranges`` returns ``should_display=False`` (invalid ranges:
    out-of-bounds, wrong direction, unsorted, cross-region, prefix mismatch).
  - Soft (do NOT gate JSON): ``SAMPLE_FLAG_RANGE_MISSING`` (incomplete coverage) stays
    as a flag attached to the Sample by ``etl.process()`` for passing samples.
  - Hard fail: ``format_error`` (non-standard separator) — a malformed
    separator can silently drop a region's range, so it gates JSON.
"""

import re
from pathlib import Path
from typing import Any, NamedTuple

from loguru import logger

from src.config import get_settings

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SAMPLE_FLAG_RANGE_POSSIBLY_WRONG = "Range - Possibly wrong"
SAMPLE_FLAG_RANGE_MISSING = "Range - Missing"
SEQUENCHER_TABLE_HEADER_ERROR = "Invalid Sequencher TXT format: missing variant table header"
DUPLICATE_LID_ERROR = "Duplicate LID in Sequencher input"

# Regex pattern for valid space-separated interval pairs in filenames
# Overlapping scan for \d+-\d+ pairs used to robustly separate the LID
# from the ranges even when the LID-to-range separator is non-standard.
# The (?<![A-Za-z0-9]) lookbehind ensures the range start is preceded by a
# separator (not a letter or digit), preventing trailing LID digits from
# being mistaken for a range start (e.g. LN_26_AA0016-73 → "0016" is part
# of the LID, not a range pair).
LID_RANGE_PAIR_RE = re.compile(r"(?=(?<![A-Za-z0-9])(\d+)-(\d+))")
# Format: LID-a1-b1 a2-b2 ... an-bn (spaces only between pairs)
VALID_INTERVAL_PATTERN = re.compile(r"^\d+-\d+( \d+-\d+)*$")
# Required Sequencher variant-table header. Accepts spaces/tabs between columns.
SEQUENCHER_TABLE_HEADER_PATTERN = re.compile(r"\bPos\b\s+\bSeq\b\s+\bCon\b\s+Required\s+Edit\b")

# Regex pattern for region prefixes in Sequencher filenames
# Matches "HV1_" or "HV2-3_" at the start of a filename (before LID)
REGION_PREFIX_PATTERN = re.compile(r"^(HV1|HV2-3)_")

# ---------------------------------------------------------------------------
# FilenameParseResult — canonical return type for filename parsing
# ---------------------------------------------------------------------------


class FilenameParseResult(NamedTuple):
    """Result of parsing a manual/Sequencher result filename.

    Attributes:
        sample_id: The LID extracted from the filename (region prefix stripped).
        ranges: Either ``"FULL REGION"`` for full-coverage files, or a list
            of ``[start, end]`` pairs for partial-coverage files. An empty
            list ``[]`` indicates a format error where ranges could not be
            parsed.
        format_error: ``None`` for valid format, ``SAMPLE_FLAG_RANGE_POSSIBLY_WRONG``
            for filenames using non-space separators between interval pairs.
        region_prefix: The region prefix detected in the filename (e.g. ``"HV1"``
            or ``"HV2-3"``), or ``None`` if no prefix was found.
    """

    sample_id: str
    ranges: str | list[list[int]]
    format_error: str | None
    region_prefix: str | None = None

    @property
    def is_full_region(self) -> bool:
        """True when the file covers all 3 HV regions."""
        return self.ranges == "FULL REGION"

    @property
    def has_format_error(self) -> bool:
        """True when the filename uses a non-standard separator."""
        return self.format_error is not None


# ---------------------------------------------------------------------------
# Region prefix helpers
# ---------------------------------------------------------------------------


def extract_region_prefix(filename: str) -> tuple[str | None, str]:
    """Extract the region prefix (HV1_ or HV2-3_) from a filename.

    Sequencher filenames may include a region prefix before the LID:
        HV1_LID.TXT, HV2-3_LID.TXT, HV1_LID-ranges.TXT, HV2-3_LID-ranges.TXT

    Args:
        filename: The filename (with or without extension).

    Returns:
        Tuple of (region_prefix, remainder) where region_prefix is ``"HV1"``
        or ``"HV2-3"`` if a prefix was found, and remainder is the filename
        portion after the prefix. If no prefix is found, returns
        ``(None, original_filename)``.
    """
    match = REGION_PREFIX_PATTERN.match(filename)
    if match:
        prefix = match.group(1)
        remainder = filename[match.end() :]
        return prefix, remainder
    return None, filename


# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------


def _split_lid_and_ranges(stripped: str) -> tuple[str, str, str]:
    r"""Split a prefix-stripped filename into (sample_id, range_part, separator).

    The sample ID is the leading LID; ``range_part`` holds the interval
    pairs. The LID-to-range separator is normally ``-``, but some
    filenames use a space, dot, or underscore instead (a non-standard
    separator).

    The standard format ``LID-a1-b1 a2-b2 ... an-bn`` (hyphen between LID
    and first pair, spaces between subsequent pairs) is tried first:
    split on the first ``-`` and validate the remainder against
    ``VALID_INTERVAL_PATTERN``. This unambiguously separates the LID from
    the ranges because LIDs never contain ``-`` after prefix stripping.

    When the standard format does not match (non-standard separator or
    non-space pair separators), fall back to a robust overlapping scan
    that locates the first ``\d+-\d+`` pair whose start < end and treats
    everything before it as the LID.

    Args:
        stripped: The filename body with region prefix and extension
            removed.

    Returns:
        Tuple of (sample_id, range_part, separator). When no genuine
        range pair is found, falls back to a ``-`` split.
    """
    # --- Standard format: LID-a1-b1 a2-b2 ... an-bn ---
    # Split on the first hyphen; LIDs never contain "-" after prefix
    # stripping, so the first "-" separates the LID from the ranges.
    parts = stripped.split("-", 1)
    if len(parts) > 1:
        candidate_lid, candidate_ranges = parts[0], parts[1]
        if VALID_INTERVAL_PATTERN.match(candidate_ranges):
            return candidate_lid, candidate_ranges, "-"

    # --- Fallback: robust scan for non-standard separators ---
    for m in LID_RANGE_PAIR_RE.finditer(stripped):
        start, end = int(m.group(1)), int(m.group(2))
        if start < end:
            pos = m.start()
            # Strip the trailing separator char(s) between LID and pair.
            sample_id = re.sub(r"[^A-Za-z0-9]+$", "", stripped[:pos])
            separator = stripped[len(sample_id) : pos]
            return sample_id, stripped[pos:], separator

    # No genuine range pair found; fall back to a single hyphen split.
    if len(parts) <= 1:
        return parts[0], "", ""
    return parts[0], stripped[len(parts[0]) + 1 :], "-"


def parse_filename_ranges(filename: str) -> FilenameParseResult:
    """Parse a manual/Sequencher result filename to extract sample ID, ranges, and format validation.

    This is the canonical function that unifies:
    - ``parse_manual_filename`` (from manual_pipeline.py)
    - ``validate_filename_format`` (from src.tools.sequencher.utils)
    - ``extract_ranges_from_filename`` (from src.tools.sequencher.utils)

    Accepted standard formats:
    1. LID.TXT → FilenameParseResult(sample_id="LID", ranges="FULL REGION", format_error=None)
    2. LID-a1-b1 a2-b2 ... an-bn.TXT → FilenameParseResult(sample_id="LID", ranges=[[a1,b1],...], format_error=None)

    Any other separator (underscore, comma, etc.) between interval pairs,
    or any separator other than "-" between the LID and the first interval
    pair, → format_error = SAMPLE_FLAG_RANGE_POSSIBLY_WRONG

    Args:
        filename: The filename to parse.

    Returns:
        FilenameParseResult with sample_id, ranges, and format_error.

    Examples:
        >>> parse_filename_ranges("2513827-73-302 316-340 438-576 16024-16365.TXT")
        FilenameParseResult(
            sample_id='2513827',
            ranges=[[73, 302], [316, 340], [438, 576], [16024, 16365]],
            format_error=None,
        )

        >>> parse_filename_ranges("2510768-73-340_438-573_16024-16365.TXT")
        FilenameParseResult(
            sample_id='2510768',
            ranges=[[73, 340], [438, 573], [16024, 16365]],
            format_error='Range - Possibly wrong',
        )

        >>> parse_filename_ranges("LN_26_AA4522 73-340 438-573 16024-16365.TXT")
        FilenameParseResult(
            sample_id='LN_26_AA4522',
            ranges=[[73, 340], [438, 573], [16024, 16365]],
            format_error='Range - Possibly wrong',
        )

        >>> parse_filename_ranges("2513827.TXT")
        FilenameParseResult(sample_id='2513827', ranges='FULL REGION', format_error=None)
    """
    base_name = filename.rsplit(".", 1)[0]  # Remove extension

    # Step 1: Detect and strip region prefix (HV1_ or HV2-3_)
    region_prefix, stripped = extract_region_prefix(base_name)

    # Step 2: Split the remainder into sample ID (LID) and range part.
    # The LID-to-range separator is normally "-", but some filenames use
    # a space, dot, or underscore instead. _split_lid_and_ranges locates
    # the first genuine range pair (start < end) so the LID is extracted
    # correctly even when the separator is non-standard.
    sample_id, range_part, lid_separator = _split_lid_and_ranges(stripped)

    # No range part -> LID only, full region
    if not range_part:
        return FilenameParseResult(
            sample_id=sample_id,
            ranges="FULL REGION",
            format_error=None,
            region_prefix=region_prefix,
        )

    # Check if range part uses only spaces between pairs (valid format).
    # The LID-to-range separator must be "-"; any other separator (space,
    # dot, underscore, ...) is a non-standard separator -> flag it.
    if VALID_INTERVAL_PATTERN.match(range_part):
        pairs = range_part.split(" ")
        ranges: list[list[int]] = []
        for pair in pairs:
            start_str, end_str = pair.split("-")
            ranges.append([int(start_str), int(end_str)])

        format_error = None if lid_separator == "-" else SAMPLE_FLAG_RANGE_POSSIBLY_WRONG
        return FilenameParseResult(
            sample_id=sample_id,
            ranges=ranges,
            format_error=format_error,
            region_prefix=region_prefix,
        )

    # Invalid separator detected -- try to extract ranges anyway for partial data
    range_matches = re.findall(r"(\d+)-(\d+)", range_part)
    if range_matches:
        ranges = [[int(s), int(e)] for s, e in range_matches]
        return FilenameParseResult(
            sample_id=sample_id,
            ranges=ranges,
            format_error=SAMPLE_FLAG_RANGE_POSSIBLY_WRONG,
            region_prefix=region_prefix,
        )

    # Could not extract any ranges
    return FilenameParseResult(
        sample_id=sample_id,
        ranges=[],
        format_error=SAMPLE_FLAG_RANGE_POSSIBLY_WRONG,
        region_prefix=region_prefix,
    )


# ---------------------------------------------------------------------------
# Range validation (rules 1-6)
# ---------------------------------------------------------------------------


def _expected_regions_for_prefix(region_prefix: str) -> set[str]:
    """Return the set of HV region names expected for a given region prefix.

    Args:
        region_prefix: ``"HV1"`` or ``"HV2-3"``.

    Returns:
        Set of region names: ``{"HV1"}`` or ``{"HV2", "HV3"}``.
    """
    if region_prefix == "HV1":
        return {"HV1"}
    if region_prefix == "HV2-3":
        return {"HV2", "HV3"}
    # Unknown prefix — allow all regions (no restriction)
    return {"HV1", "HV2", "HV3"}


def _check_range_direction(ranges: list[list[int]]) -> tuple[bool, str | None]:
    """Rule 4: Range direction — for any range a-b, a < b."""
    for start, end in ranges:
        if start >= end:
            return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG
    return True, None


def _check_position_bounds(ranges: list[list[int]], valid_bounds: list[tuple[int, int]]) -> tuple[bool, str | None]:
    """Rule 2: Position bounds — all positions within valid bounds."""
    for start, end in ranges:
        if not any(lo <= start <= hi or lo <= end <= hi for lo, hi in valid_bounds):
            return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG
    return True, None


def _check_sorted_order(ranges: list[list[int]]) -> tuple[bool, str | None]:
    """Rule 1: Sorted order — a1 < b1 < a2 < b2 < ... < an < bn."""
    for i in range(len(ranges) - 1):
        if ranges[i][1] >= ranges[i + 1][0]:
            return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG
    return True, None


def _check_region_consistency(ranges: list[list[int]], regions: dict) -> tuple[bool, str | None, set[str]]:
    """Rule 3: Region consistency — each (a, b) pair belongs to the same HV region.

    Returns (is_valid, error_flag, found_regions).
    """

    def _find_region(pos: int) -> str | None:
        for region_name, (start, end) in regions.items():
            if start <= pos <= end:
                return region_name
        return None

    found_regions: set[str] = set()
    for start, end in ranges:
        start_region = _find_region(start)
        end_region = _find_region(end)
        if start_region is None or start_region != end_region:
            return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG, found_regions
        found_regions.add(start_region)
    return True, None, found_regions


def _check_prefix_consistency(
    region_prefix: str,
    found_regions: set[str],
    filename: str | None = None,
) -> tuple[bool, str | None]:
    """Rule 6: Region-prefix consistency — ranges must match the declared region prefix."""
    expected_regions = _expected_regions_for_prefix(region_prefix)
    unexpected = found_regions - expected_regions
    if unexpected:
        label = f"File {filename}: " if filename else ""
        logger.warning(
            "{}region prefix '{}' expects ranges in {}, but found ranges in {}",
            label,
            region_prefix,
            expected_regions,
            found_regions,
        )
        return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG
    return True, None


def _check_completeness(
    region_prefix: str | None,
    found_regions: set[str],
    regions: dict,
    filename: str | None = None,
) -> tuple[bool, str | None] | None:
    """Rule 5: Completeness — at least one pair in each required HV region.

    Returns None if complete, otherwise (should_display, flag_string).
    """
    label = f"File {filename}: " if filename else ""
    if region_prefix is not None:
        expected = _expected_regions_for_prefix(region_prefix)
        missing = expected - found_regions
        if missing:
            logger.info("{}missing regions for prefix '{}': {}", label, region_prefix, missing)
            return True, SAMPLE_FLAG_RANGE_MISSING
    elif len(found_regions) < len(regions):
        missing = set(regions.keys()) - found_regions
        logger.info("{}missing regions {}", label, missing)
        return True, SAMPLE_FLAG_RANGE_MISSING
    return None


def validate_analysis_ranges(
    ranges: str | list[list[int]],
    region_prefix: str | None = None,
    filename: str | None = None,
) -> tuple[bool, str | None]:
    """Validate analysis ranges from a Sequencher TXT filename.

    Checks 6 rules:
    1. Sorted order: a1 < b1 < a2 < b2 < ... < an < bn
    2. Position bounds: all positions within {73-340, 438-576, 16024-16365}
    3. Region consistency: each (a,b) pair belongs to the same HV region
    4. Range direction: for any range a-b, a < b
    5. Completeness: at least one pair in each of HV1, HV2, HV3
    6. Prefix consistency: when a region prefix is present, ranges lie in its regions

    Args:
        ranges: List of [start, end] pairs, ``"FULL REGION"`` string, or empty list.
        region_prefix: Optional region prefix ("HV1" or "HV2-3") for prefix validation.
        filename: Optional source filename, included in completeness/prefix-consistency
            log messages so it is clear which sample is affected.

    Returns:
        (should_display_variants, flag_string)
        - (True, None) for valid full-region or complete ranges
        - (True, SAMPLE_FLAG_RANGE_MISSING) when valid but missing HV region(s)
        - (False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG) when rules 1-4 are violated
    """
    regions = get_settings().regions.REGIONS

    # "FULL REGION" or empty means default full coverage — always valid
    if ranges == "FULL REGION" or not ranges:
        return True, None

    # Type narrowing: beyond this point, ranges is list[list[int]]
    if not isinstance(ranges, list):
        return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG

    return _validate_ranges_list(ranges, region_prefix, regions, filename)


def _validate_ranges_list(
    ranges: list[list[int]],
    region_prefix: str | None,
    regions: dict,
    filename: str | None = None,
) -> tuple[bool, str | None]:
    """Validate a list of range pairs against all rules (rules 1-6)."""
    valid_bounds = [
        (regions["HV2"][0], regions["HV2"][1]),  # 73-340
        (regions["HV3"][0], regions["HV3"][1]),  # 438-576
        (regions["HV1"][0], regions["HV1"][1]),  # 16024-16365
    ]

    for check in (
        _check_range_direction(ranges),
        _check_position_bounds(ranges, valid_bounds),
        _check_sorted_order(ranges),
    ):
        ok, err = check
        if not ok:
            return False, err

    ok, err, found_regions = _check_region_consistency(ranges, regions)
    if not ok:
        return False, err

    # Rule 6: Region-prefix consistency
    if region_prefix is not None:
        _, prefix_result = _check_prefix_consistency(region_prefix, found_regions, filename)
        if prefix_result == SAMPLE_FLAG_RANGE_POSSIBLY_WRONG:
            return False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG

    # Rule 5: Completeness
    return _check_completeness(region_prefix, found_regions, regions, filename) or (True, None)


# ---------------------------------------------------------------------------
# Batch-level QC helpers
# ---------------------------------------------------------------------------


def reject_duplicate_lids(qc_results: list[dict[str, Any]]) -> None:
    """Mark every file for a duplicated LID as a hard QC failure.

    Duplicate detection is batch-level because a single filename cannot know
    whether another TXT in the same input directory resolves to the same LID.
    The function mutates the QC report rows in place so downstream task
    selection and the written ``qc_report.json`` stay consistent.

    Args:
        qc_results: QC report rows containing at least ``sample_id``,
            ``filename``, ``status``, and ``reason`` keys.
    """
    files_by_lid: dict[str, list[str]] = {}
    for result in qc_results:
        sample_id = str(result.get("sample_id") or "")
        if not sample_id:
            continue
        filename = str(result.get("filename") or "")
        files_by_lid.setdefault(sample_id, []).append(filename)

    for sample_id, filenames in files_by_lid.items():
        if len(filenames) <= 1:
            continue

        joined_filenames = ", ".join(sorted(filenames))
        reason = f"{DUPLICATE_LID_ERROR}: {sample_id} appears in multiple TXT files ({joined_filenames})"
        for result in qc_results:
            if result.get("sample_id") == sample_id:
                result["status"] = "fail"
                result["reason"] = reason


# ---------------------------------------------------------------------------
# TXT format validation
# ---------------------------------------------------------------------------


def has_sequencher_variant_table_header(file_path: str | Path) -> bool:
    """Return True when a TXT file contains the required variant table header.

    Sequencher variant TXT output includes a table with these columns, in order::

        Pos       Seq       Con       Required Edit

    Files without this header are not valid variant TXT exports and should not
    enter ETL, because the parser would otherwise treat prose lines as data and
    emit an empty or misleading variant set.
    """
    try:
        with Path(file_path).open(encoding="utf-8") as f:
            return any(SEQUENCHER_TABLE_HEADER_PATTERN.search(line) for line in f)
    except UnicodeDecodeError:
        with Path(file_path).open(encoding="latin-1") as f:
            return any(SEQUENCHER_TABLE_HEADER_PATTERN.search(line) for line in f)
    except OSError as exc:
        logger.warning("Unable to read Sequencher TXT for QC ({}): {}", file_path, exc)
        return False


# ---------------------------------------------------------------------------
# QC gate
# ---------------------------------------------------------------------------


class QCResult(NamedTuple):
    """Result of the pre-ETL Sequencher QC gate.

    Attributes:
        sample_id: The LID extracted from the filename (empty string when no
            LID could be extracted).
        passed: True when the sample passes QC and should produce JSON.
        reason: ``None`` on pass; a human-readable failure reason on fail.
    """

    sample_id: str
    passed: bool
    reason: str | None


def run_qc(filename: str, file_path: str | Path | None = None) -> QCResult:
    """Run the pre-ETL quality-control gate on a Sequencher TXT file.

    Hard fail (no JSON output) when:
      - no LID can be extracted from the filename, **or**
      - ``file_path`` is provided and the TXT is missing the required
        Sequencher variant table header, **or**
      - ``validate_analysis_ranges`` returns ``should_display=False`` (invalid
        ranges: out-of-bounds, wrong direction, unsorted, cross-region, prefix
        mismatch), **or**
      - ``format_error`` (non-standard separator) — a malformed separator can
        silently drop a region's range, so it is unsafe to emit.

    Soft (do NOT gate JSON — stays as a flag on the Sample for passing samples):
      - ``SAMPLE_FLAG_RANGE_MISSING`` (incomplete coverage)

    Args:
        filename: The TXT filename to validate.
        file_path: Optional path to the TXT content. When provided, QC also
            requires the Sequencher variant table header columns: ``Pos``,
            ``Seq``, ``Con``, and ``Required Edit``.

    Returns:
        A ``QCResult`` with the extracted ``sample_id``, ``passed`` flag, and a
        ``reason`` (``None`` on pass).
    """
    parse_result = parse_filename_ranges(filename)

    if not parse_result.sample_id:
        return QCResult(
            sample_id="",
            passed=False,
            reason="No sample ID (LID) could be extracted from filename",
        )

    if file_path is not None and not has_sequencher_variant_table_header(file_path):
        return QCResult(
            sample_id=parse_result.sample_id,
            passed=False,
            reason=SEQUENCHER_TABLE_HEADER_ERROR,
        )

    # Full-region files are valid after TXT content/header validation.
    if parse_result.is_full_region:
        return QCResult(sample_id=parse_result.sample_id, passed=True, reason=None)

    # A non-standard separator is a HARD QC failure. A malformed separator
    # can silently drop a region's range (e.g. "16024.16365" instead of
    # "16024-16365"), producing incomplete coverage that is not safe to emit.
    if parse_result.has_format_error:
        return QCResult(
            sample_id=parse_result.sample_id,
            passed=False,
            reason=f"Filename format error (non-standard separator) {parse_result.format_error}",
        )

    should_display, flag = validate_analysis_ranges(
        parse_result.ranges,
        region_prefix=parse_result.region_prefix,
        filename=filename,
    )
    if not should_display:
        return QCResult(
            sample_id=parse_result.sample_id,
            passed=False,
            reason=flag or SAMPLE_FLAG_RANGE_POSSIBLY_WRONG,
        )

    # Passing sample (ranges valid, no format error). Soft flags like
    # Range - Missing are attached by etl.process() for passing samples.
    return QCResult(sample_id=parse_result.sample_id, passed=True, reason=None)
