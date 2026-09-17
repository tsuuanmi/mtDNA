"""Typed population-analysis contracts; upstream Sample remains authoritative."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from src.core.models import Sample  # noqa: TC001 - Pydantic resolves this type at runtime

REGIONS = ("HV1", "HV2", "HV3")
COMBINED = "HV1+HV2+HV3"


class AnalysisOptions(BaseModel):
    confidence_level: float | None = Field(default=0.95, gt=0, lt=1)
    presentation_precision: int = Field(default=6, ge=0, le=15)
    additional_frequency_thresholds: list[float] = Field(default_factory=list)


class PopulationCohort(BaseModel):
    samples: list[Sample]
    input_sha256: str | None = None
    eligibility_basis: str = "All records in the supplied upstream-approved cohort export"


class ReferenceVariantAnnotation(BaseModel):
    observed: bool | None = None
    reference_frequency: float | None = Field(default=None, ge=0, le=1)
    reference_population: str | None = None
    database_version: str | None = None
    provenance: str | None = None


class FrequencyInterval(BaseModel):
    confidence_level: float | None = None
    frequency_ci_lower: float | None = None
    frequency_ci_upper: float | None = None
    frequency_ci_method: Literal["wilson"] | None = None


class VariantPopulationStat(FrequencyInterval):
    variant: str
    region: str
    position: int | str
    reference: str
    alternate: str
    carrier_count: int
    callable_sample_count: int
    population_frequency: float
    is_singleton: bool
    is_doubleton: bool
    carrier_count_le_2: bool
    frequency_lt_0_001: bool
    frequency_lt_0_005: bool
    frequency_lt_0_01: bool
    additional_frequency_flags: dict[str, bool] = Field(default_factory=dict)
    reference_annotations: dict[str, ReferenceVariantAnnotation] = Field(default_factory=dict)


class SampleHaplotype(BaseModel):
    sample_id: str
    hv1_haplotype: str | None
    hv2_haplotype: str | None
    hv3_haplotype: str | None
    combined_haplotype: str | None
    haplotype_id: str | None
    complete_hv_profile: bool
    region_complete: dict[str, bool]
    incomplete_reasons: list[str]


class DataIssue(BaseModel):
    """One fixable upstream representation problem, named for the affected sample."""

    sample_id: str
    issue: str
    region: str
    details: str


class HaplotypePopulationStat(FrequencyInterval):
    scope: str
    haplotype_id: str
    canonical_haplotype: str
    count: int
    frequency: float
    is_singleton: bool


class HaplotypeDiversity(BaseModel):
    n_complete_profiles: int
    n_distinct_haplotypes: int
    n_singleton_haplotypes: int
    n_recurrent_haplotypes: int
    n_individuals_with_unique_haplotype: int
    fraction_individuals_with_unique_haplotype: float | None
    most_common_haplotype_count: int
    most_common_haplotype_frequency: float | None
    empirical_random_match_probability: float | None
    discrimination_power: float | None
    haplotype_frequency_squared_sum: float | None
    haplotype_diversity: float | None


class VariantDiversity(BaseModel):
    n_distinct_variants: int
    n_singleton_variants: int
    n_doubleton_variants: int
    n_variants_carrier_count_le_2: int
    n_variants_af_lt_0_001: int
    n_variants_af_lt_0_005: int
    n_variants_af_lt_0_01: int
    additional_frequency_counts: dict[str, int] = Field(default_factory=dict)


class PopulationSummary(BaseModel):
    analysis_timestamp: str
    schema_version: str = "1.0"
    software_version: str
    input_sha256: str | None
    reference_sha256: str
    region_configuration_sha256: str
    region_definitions: dict[str, list[int]]
    regions_analyzed: list[str]
    eligibility_basis: str
    n_samples_analyzed: int
    n_samples_complete_haplotype: int
    n_samples_incomplete_haplotype: int
    data_issue_counts: dict[str, int] = Field(default_factory=dict)
    variants: dict[str, VariantDiversity]
    haplotypes: HaplotypeDiversity
    region_haplotypes: dict[str, HaplotypeDiversity]
    reference_not_observed_counts: dict[str, int]
    options: AnalysisOptions
    methods: dict[str, str]
    output_sha256: dict[str, str] = Field(default_factory=dict)


class PopulationResult(BaseModel):
    summary: PopulationSummary
    variants: list[VariantPopulationStat]
    haplotypes: list[HaplotypePopulationStat]
    region_haplotypes: list[HaplotypePopulationStat]
    sample_haplotypes: list[SampleHaplotype]
    data_issues: list[DataIssue] = Field(default_factory=list)
