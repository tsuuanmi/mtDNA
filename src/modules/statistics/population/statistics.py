"""Pure cohort frequency and haplotype diversity statistics."""

from math import sqrt
from statistics import NormalDist

from src.modules.statistics.population.models import FrequencyInterval, HaplotypeDiversity


def _validate_count(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        msg = "Counts must be integers"
        raise TypeError(msg)
    if value < 0:
        msg = "Counts must be nonnegative"
        raise ValueError(msg)


def frequency_interval(count: int, total: int, confidence_level: float | None) -> FrequencyInterval:
    """Return a two-sided Wilson score interval, or null bounds when unavailable."""
    _validate_count(count)
    _validate_count(total)
    if count > total:
        msg = "Count must not exceed total"
        raise ValueError(msg)
    if confidence_level is not None and (isinstance(confidence_level, bool) or not 0 < confidence_level < 1):
        msg = "Confidence level must be strictly between zero and one"
        raise ValueError(msg)
    if confidence_level is None or total == 0:
        return FrequencyInterval(confidence_level=confidence_level)
    # The lower tail avoids rounding (1 + confidence_level) / 2 to one.
    z_squared = NormalDist().inv_cdf((1 - confidence_level) / 2) ** 2
    frequency = count / total
    denominator = 1 + z_squared / total
    center = (frequency + z_squared / (2 * total)) / denominator
    radius = sqrt(frequency * (1 - frequency) / total + z_squared / (4 * total**2)) * sqrt(z_squared)
    radius /= denominator
    return FrequencyInterval(
        confidence_level=confidence_level,
        frequency_ci_lower=0.0 if count == 0 else max(0.0, center - radius),
        frequency_ci_upper=1.0 if count == total else min(1.0, center + radius),
        frequency_ci_method="wilson",
    )


def diversity(counts: list[int]) -> HaplotypeDiversity:
    """Summarize positive haplotype counts; zero-count categories are unobserved.

    Pairwise matching samples without replacement is distinct from the squared
    frequency sum. Corrected diversity and discrimination power equal one minus
    that pairwise probability, and require at least two complete profiles.
    """
    for count in counts:
        _validate_count(count)
    observed = [count for count in counts if count > 0]
    total = sum(observed)
    singletons = observed.count(1)
    largest = max(observed, default=0)
    squared_sum = sum(count * count for count in observed) / (total * total) if total else None
    pairwise = sum(count * (count - 1) for count in observed) / (total * (total - 1)) if total > 1 else None
    corrected = 1 - pairwise if pairwise is not None else None
    return HaplotypeDiversity(
        n_complete_profiles=total,
        n_distinct_haplotypes=len(observed),
        n_singleton_haplotypes=singletons,
        n_recurrent_haplotypes=len(observed) - singletons,
        n_individuals_with_unique_haplotype=singletons,
        fraction_individuals_with_unique_haplotype=singletons / total if total else None,
        most_common_haplotype_count=largest,
        most_common_haplotype_frequency=largest / total if total else None,
        empirical_random_match_probability=pairwise,
        discrimination_power=corrected,
        haplotype_frequency_squared_sum=squared_sum,
        haplotype_diversity=corrected,
    )
