"""Region-grouped views of canonical final batch JSON.

The regional table is a display projection of the final regenerate JSON. Every
sample field and every variant dict is preserved verbatim; only the
``variants`` container changes from type groups to region groups.
"""

import json
from pathlib import Path
from typing import Any

from src.core.variants import get_hv_region_for_position, pos_sort_key

REGION_ORDER = ("HV1", "HV2", "HV3")
TYPE_ORDER = ("snps", "insertions", "deletions")


def _flatten_type_grouped_variants(variants: object) -> list[dict[str, Any]]:
    """Flatten canonical type-grouped variants without changing row values."""
    if not isinstance(variants, dict):
        return []

    rows: list[dict[str, Any]] = []
    for variant_type in TYPE_ORDER:
        group = variants.get(variant_type, [])
        if isinstance(group, list):
            rows.extend(dict(row) for row in group if isinstance(row, dict))
    return rows


def build_region_grouped_table(final_batch: dict[str, Any]) -> dict[str, Any]:
    """Group final canonical variants by HV region without changing their values."""
    table: dict[str, Any] = {}
    for sample_id, sample in final_batch.items():
        if not isinstance(sample, dict):
            table[sample_id] = sample
            continue

        groups: dict[str, list[dict[str, Any]]] = {region: [] for region in REGION_ORDER}
        for variant in _flatten_type_grouped_variants(sample.get("variants")):
            pos = variant.get("pos")
            region = get_hv_region_for_position(pos) if isinstance(pos, (int, str)) else None
            if region not in groups:
                msg = f"Final variant at position {pos!r} is outside HV1/HV2/HV3"
                raise ValueError(msg)
            groups[region].append(variant)
        for variants in groups.values():
            variants.sort(key=lambda variant: pos_sort_key(variant["pos"]))

        table[sample_id] = {**sample, "variants": groups}
    return table


def write_region_grouped_table(statistic_fullbatch_path: Path, output_path: Path) -> None:
    """Read final ``statistic_fullbatch.json`` and write its region-grouped view."""
    final_batch = json.loads(statistic_fullbatch_path.read_text(encoding="utf-8"))
    table = build_region_grouped_table(final_batch)
    output_path.write_text(json.dumps(table, indent=2), encoding="utf-8")
