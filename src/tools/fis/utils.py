"""FIS workbook and profile parsing helpers."""

from typing import Any

import pandas as pd

from src.core.variants import is_position_in_intervals, parse_variant_position_allele

MISSING_VALUE_TEXT = {"", "nan", "none", "<na>"}
MISSING_PROFILE_VALUES = MISSING_VALUE_TEXT | {"-", "/"}


def clean_value(value: object) -> object | None:
    """Return None for spreadsheet missing values."""
    if value is None or str(value).strip().lower() in MISSING_VALUE_TEXT:
        return None
    return value


def clean_sample_id(value: object | None) -> str | None:
    """Normalize an FIS sample identifier."""
    cleaned = clean_value(value)
    if cleaned is None:
        return None
    text = str(cleaned).strip()
    if "." in text and text.replace(".", "").isdigit():
        return str(int(float(text)))
    return text


def format_genotypes(values: list[str]) -> str:
    """Deduplicate and sort comma-separated FIS genotype calls."""
    genotypes = [genotype.strip() for value in values for genotype in value.split(",") if genotype.strip()]
    unique = list(dict.fromkeys(genotypes))
    unique.sort(key=genotype_position)
    return " ".join(unique)


def genotype_position(genotype: str) -> int:
    """Return the integer base coordinate for an FIS genotype token."""
    position, _ = parse_variant_position_allele(genotype)
    try:
        return int(position.split(".", maxsplit=1)[0])
    except ValueError:
        return 0


def parse_genotypes(profile: object | None) -> list[str]:
    """Split a space-separated FIS profile into genotype tokens."""
    cleaned = clean_value(profile)
    if cleaned is None:
        return []
    return [token.strip() for token in str(cleaned).split() if token.strip()]


def select_profile(corrected_profile: object | None, raw_profile: str) -> tuple[str, str]:
    """Return corrected profile when available, otherwise raw profile."""
    corrected = clean_value(corrected_profile)
    if corrected is None or str(corrected).strip().lower() in MISSING_PROFILE_VALUES:
        return raw_profile, "raw_fallback"
    return format_genotypes([str(corrected)]), "corrected"


def filter_genotypes_by_intervals(genotypes: list[str], intervals: list[tuple[int, int]]) -> list[str]:
    """Keep genotype tokens located in the supplied intervals."""
    return [
        genotype
        for genotype in genotypes
        if (position := parse_variant_position_allele(genotype)[0]) and is_position_in_intervals(position, intervals)
    ]


def json_safe_value(value: object) -> object | None:
    """Convert common pandas values to JSON-compatible values."""
    cleaned = clean_value(value)
    if cleaned is None:
        return None
    if isinstance(cleaned, pd.Timestamp):
        return cleaned.isoformat()
    item = getattr(cleaned, "item", None)
    return item() if callable(item) else cleaned


def row_position(value: object) -> int | None:
    """Parse a T2 position value such as ``MT:73``."""
    if value is None or str(value).strip().lower() in MISSING_PROFILE_VALUES:
        return None
    text = str(value)
    return int(text.split(":", maxsplit=1)[1] if ":" in text else text)


def _integer_or_none(value: object) -> int | None:
    """Convert one scalar spreadsheet integer or return None."""
    cleaned = clean_value(value)
    return int(float(str(cleaned))) if cleaned is not None else None


def raw_variant_row(row: pd.Series) -> dict[str, Any]:
    """Retain the raw T2 columns used by the historical NGS artifact."""
    return {
        "pos": row_position(row.get("Position")),
        "Marker": row.get("Marker"),
        "Genotype": clean_value(row.get("Genotype")),
        "AlleleFrequency": clean_value(row.get("AlleleFrequency")),
        "TotalDepth": _integer_or_none(row.get("TotalDepth")),
        "Ref(Depth):Alt(Depth)": clean_value(row.get("Ref(Depth):Alt(Depth)")),
        "QC_Info": clean_value(row.get("QC_Info")),
        "Genotype(raw)": clean_value(row.get("Genotype(raw)")),
    }
