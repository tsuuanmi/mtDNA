"""
Shared constants and helpers for the Mutation Surveyor module.

Contains IUPAC maps, region constants, and utility functions
used by flagging, ETL, pipeline, and sort_review_data modules.

Position type is defined in src.core.models — import Position from there.
"""

import re

import pandas as pd

from src.core.models import Position
from src.core.variants import normalize_position

# ---------------------------------------------------------------------------
# IUPAC ambiguity codes (reverse map for concordance)
# ---------------------------------------------------------------------------
UNORDERED_MAP: dict[str, str] = {
    "".join(sorted("AT")): "W",
    "".join(sorted("AG")): "R",
    "".join(sorted("AC")): "M",
    "".join(sorted("GT")): "K",
    "".join(sorted("CT")): "Y",
    "".join(sorted("CG")): "S",
    "".join(sorted("AGT")): "D",
    "".join(sorted("ACT")): "H",
    "".join(sorted("ACG")): "V",
    "".join(sorted("CGT")): "B",
}

# ---------------------------------------------------------------------------
# String / display helpers
# ---------------------------------------------------------------------------


def safe_str(x: object) -> str:
    """Convert to stripped string; NaN/None → empty string."""
    if isinstance(x, float | int):
        return "" if pd.isna(x) else str(x).strip()
    if x is None:
        return ""
    return str(x).strip()


# ---------------------------------------------------------------------------
# Variant token parsing
# ---------------------------------------------------------------------------


def parse_variant_tokens(text: object) -> list[tuple[Position, str, str]]:
    """Parse a space-separated variant string into (position, allele, token) triples.

    Positions follow the canonical Position type:
    - int for base positions (e.g. 309, 16569)
    - str for insertion positions (e.g. "309.1", "217.10")

    Args:
        text: Space-separated variant tokens (e.g. "152C 309.1C 16189DEL")

    Returns:
        List of (position, allele, original_token) triples.
    """
    if not isinstance(text, str):
        if isinstance(text, float) and pd.isna(text):
            return []
        return []
    if not text:
        return []

    s = safe_str(text)
    if not s:
        return []

    tokens: list[tuple[Position, str, str]] = []
    for part in s.split():
        # Match: position + allele (e.g. "152C", "309.1C", "16189DEL", "16189het_DEL")
        m = re.match(r"^(\d+(?:\.\d+)?)([A-Za-z]+(?:_[Hh][Ee][Tt])?(?:_[Dd][Ee][Ll])?)$", part)
        if m:
            pos = normalize_position(m.group(1))
            allele = m.group(2).replace("_HET", "_het")
            tokens.append((pos, allele, part))
            continue

        # Fallback: position only (e.g. bare number or partial match)
        m2 = re.match(r"^(\d+(?:\.\d+)?)", part)
        if m2:
            pos = normalize_position(m2.group(1))
            allele = part[m2.end() :] if m2.end() < len(part) else ""
            allele = allele.replace("_HET", "_het")
            tokens.append((pos, allele, part))

    return tokens


# ---------------------------------------------------------------------------
# Control LID detection (shared across inventory_qc, control_qc, merge_profiles)
# ---------------------------------------------------------------------------

_CONTROL_PC_RE = re.compile(r"^PC\d+_.+$", re.IGNORECASE)
_CONTROL_NTC_RE = re.compile(r"^NTC\d+$", re.IGNORECASE)

EXPECTED_PRIMERS: list[str] = ["HV1F", "HV1R", "HV2F", "HV3R"]


def is_control_lid(lid: object) -> bool:
    """Return True for control LIDs (PC or NTC) expected in MS reports.

    Matches:
      - NTC1, NTC2, ... (exact NTC + digits)
      - PC1_mtDNA_247547, PC2_mtDNA_247547, ... (PC + digits + underscore + suffix)

    Does NOT match strings that merely start with "PC" or "NTC" (e.g. "PCABC").
    """
    s = safe_str(lid)
    return _CONTROL_NTC_RE.fullmatch(s) is not None or _CONTROL_PC_RE.fullmatch(s) is not None
