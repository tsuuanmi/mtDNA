"""Deterministic, full-precision population reports without person metadata."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING

from loguru import logger

from src.modules.statistics.population.models import (
    DataIssue,
    HaplotypePopulationStat,
    SampleHaplotype,
    VariantPopulationStat,
)

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pydantic import BaseModel

    from src.modules.statistics.population.models import PopulationResult


def _write_tsv(path: Path, rows: Iterable[BaseModel], model: type[BaseModel]) -> None:
    """Keep model field order; encode nested structures as deterministic JSON."""
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(model.model_fields), delimiter="\t")
        writer.writeheader()
        for row in rows:
            values = row.model_dump(mode="json")
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in values.items()
                }
            )


def _executive_summary(result: PopulationResult) -> str:
    summary = result.summary
    variants = summary.variants["all"]
    haplotypes = summary.haplotypes

    def number(value: float | None, *, percent: bool = False) -> str:
        if value is None:
            return "unavailable"
        if percent:
            return f"{100 * value:.{summary.options.presentation_precision}f}%"
        return f"{value:.{summary.options.presentation_precision}f}"

    lines = [
        "Vietnamese mtDNA Dataset Summary",
        f"Analysis timestamp: {summary.analysis_timestamp}",
        "Scope: current cohort only; frequencies are not estimates for the general population.",
        f"Regions analyzed: {', '.join(summary.regions_analyzed)}",
        f"Region definitions: {json.dumps(summary.region_definitions, sort_keys=True)}",
        f"Samples analyzed: {summary.n_samples_analyzed}",
        f"Complete / incomplete combined profiles: "
        f"{summary.n_samples_complete_haplotype} / {summary.n_samples_incomplete_haplotype}",
        f"Samples with data issues: {len({row.sample_id for row in result.data_issues})} "
        f"(see mtDNA_population_data_issues.tsv)",
        "",
        "Variant diversity (carrier frequency among callable samples):",
        f"Distinct variants: {variants.n_distinct_variants}",
        f"Singleton / doubleton variants: {variants.n_singleton_variants} / {variants.n_doubleton_variants}",
        f"Variants with at most two carriers: {variants.n_variants_carrier_count_le_2}",
        f"Variants with frequency <0.1% / <0.5% / <1%: "
        f"{variants.n_variants_af_lt_0_001} / {variants.n_variants_af_lt_0_005} / {variants.n_variants_af_lt_0_01}",
        "",
        "Combined HV1+HV2+HV3 haplotypes:",
        f"Distinct / singleton / recurrent haplotypes: {haplotypes.n_distinct_haplotypes} / "
        f"{haplotypes.n_singleton_haplotypes} / {haplotypes.n_recurrent_haplotypes}",
        f"Individuals with a unique haplotype: {haplotypes.n_individuals_with_unique_haplotype} "
        f"({number(haplotypes.fraction_individuals_with_unique_haplotype, percent=True)})",
        f"Most common haplotype: {haplotypes.most_common_haplotype_count} individuals "
        f"({number(haplotypes.most_common_haplotype_frequency, percent=True)})",
        f"Empirical random match probability: {number(haplotypes.empirical_random_match_probability)}",
        f"Discrimination power: {number(haplotypes.discrimination_power)}",
        f"Corrected haplotype diversity: {number(haplotypes.haplotype_diversity)}",
        f"Squared haplotype frequency sum (with replacement): {number(haplotypes.haplotype_frequency_squared_sum)}",
        "",
        "Reference database explicit negatives (observed=false only):",
    ]
    databases = set(summary.reference_not_observed_counts)
    for variant in result.variants:
        databases.update(variant.reference_annotations)
    if not databases:
        lines.append("  unavailable: no reference database annotations supplied")
    for database in sorted(databases):
        negative_count = sum(
            annotation.observed is False
            for variant in result.variants
            if (annotation := variant.reference_annotations.get(database)) is not None
        )
        lines.append(f"  {database}: {negative_count} explicitly not observed variants")
    lines.extend(
        [
            "Missing or unknown annotations do not establish absence or novelty.",
            "Regional metrics, additional thresholds, methods and provenance are in the JSON/TSV reports.",
            "",
            "mtDNA is maternally inherited: maternal relatives can share the same haplotype.",
            "These cohort statistics characterize maternal lineages, not unique individual identification.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_reports(result: PopulationResult, output_dir: str | Path) -> list[Path]:
    """Write seven reports into a new/empty directory; summary hashes other outputs."""
    output_dir = Path(output_dir)
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        msg = f"Population report output directory must be new or empty: {output_dir}"
        raise ValueError(msg)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "mtDNA_population_summary.json"
    tsv_reports = [
        ("mtDNA_variant_population_frequencies.tsv", result.variants, VariantPopulationStat),
        ("mtDNA_haplotype_frequencies.tsv", result.haplotypes, HaplotypePopulationStat),
        ("mtDNA_region_haplotype_frequencies.tsv", result.region_haplotypes, HaplotypePopulationStat),
        ("mtDNA_sample_haplotypes.tsv", result.sample_haplotypes, SampleHaplotype),
        ("mtDNA_population_data_issues.tsv", result.data_issues, DataIssue),
    ]
    paths: list[Path] = []
    for filename, rows, model in tsv_reports:
        path = output_dir / filename
        _write_tsv(path, rows, model)
        logger.info("Population reports: wrote {}", path)
        paths.append(path)
    executive_path = output_dir / "mtDNA_population_executive_summary.txt"
    with executive_path.open("x", encoding="utf-8") as handle:
        handle.write(_executive_summary(result))
    logger.info("Population reports: wrote {}", executive_path)
    paths.append(executive_path)
    output_hashes: dict[str, str] = {}
    for path in paths:
        with path.open("rb") as handle:
            output_hashes[path.name] = hashlib.file_digest(handle, "sha256").hexdigest()
    result.summary.output_sha256 = output_hashes
    with summary_path.open("x", encoding="utf-8") as handle:
        json.dump(
            result.summary.model_dump(mode="json"),
            handle,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        handle.write("\n")
    logger.info("Population reports: wrote {}", summary_path)
    return [summary_path, *paths]
