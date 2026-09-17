"""Pure population statistics edge cases and Wilson score checks."""

import math

import pytest

from src.modules.statistics.population.statistics import diversity, frequency_interval


@pytest.mark.parametrize(
    ("counts", "distinct", "singletons", "recurrent", "largest", "squared", "pairwise"),
    [
        ([1, 1, 1, 1], 4, 4, 0, 1, 0.25, 0.0),
        ([4], 1, 0, 1, 4, 1.0, 1.0),
        ([3, 2, 1], 3, 1, 2, 3, 14 / 36, 8 / 30),
    ],
)
def test_diversity_populations(counts, distinct, singletons, recurrent, largest, squared, pairwise):
    result = diversity(counts)
    total = sum(counts)
    assert result.n_complete_profiles == total
    assert result.n_distinct_haplotypes == distinct
    assert result.n_singleton_haplotypes == singletons
    assert result.n_recurrent_haplotypes == recurrent
    assert result.n_individuals_with_unique_haplotype == singletons
    assert result.fraction_individuals_with_unique_haplotype == pytest.approx(singletons / total)
    assert result.most_common_haplotype_count == largest
    assert result.most_common_haplotype_frequency == pytest.approx(largest / total)
    assert result.haplotype_frequency_squared_sum == pytest.approx(squared)
    assert result.empirical_random_match_probability == pytest.approx(pairwise)
    assert result.discrimination_power == pytest.approx(1 - pairwise)
    assert result.haplotype_diversity == pytest.approx(total / (total - 1) * (1 - squared))


@pytest.mark.parametrize("counts", [[], [0, 0]])
def test_empty_diversity(counts):
    result = diversity(counts)
    assert result.n_complete_profiles == 0
    assert result.n_distinct_haplotypes == 0
    assert result.n_singleton_haplotypes == 0
    assert result.n_recurrent_haplotypes == 0
    assert result.n_individuals_with_unique_haplotype == 0
    assert result.most_common_haplotype_count == 0
    assert result.fraction_individuals_with_unique_haplotype is None
    assert result.most_common_haplotype_frequency is None
    assert result.haplotype_frequency_squared_sum is None
    assert result.empirical_random_match_probability is None
    assert result.discrimination_power is None
    assert result.haplotype_diversity is None


def test_one_profile_diversity():
    result = diversity([1])
    assert result.n_complete_profiles == 1
    assert result.n_distinct_haplotypes == 1
    assert result.n_singleton_haplotypes == 1
    assert result.n_recurrent_haplotypes == 0
    assert result.n_individuals_with_unique_haplotype == 1
    assert result.fraction_individuals_with_unique_haplotype == 1
    assert result.most_common_haplotype_count == 1
    assert result.most_common_haplotype_frequency == 1
    assert result.haplotype_frequency_squared_sum == 1
    assert result.empirical_random_match_probability is None
    assert result.discrimination_power is None
    assert result.haplotype_diversity is None


def test_wilson_known_interval():
    result = frequency_interval(5, 10, 0.95)
    assert result.confidence_level == 0.95
    assert result.frequency_ci_method == "wilson"
    assert result.frequency_ci_lower == pytest.approx(0.236593090512564)
    assert result.frequency_ci_upper == pytest.approx(0.763406909487436)


@pytest.mark.parametrize(("count", "lower", "upper"), [(0, 0, 0.277532799862889), (10, 0.722467200137111, 1)])
def test_wilson_boundary_counts(count, lower, upper):
    result = frequency_interval(count, 10, 0.95)
    assert result.frequency_ci_lower == pytest.approx(lower)
    assert result.frequency_ci_upper == pytest.approx(upper)


@pytest.mark.parametrize(("count", "total", "level"), [(0, 0, 0.95), (1, 2, None), (0, 0, None)])
def test_unavailable_interval(count, total, level):
    result = frequency_interval(count, total, level)
    assert result.confidence_level == level
    assert result.frequency_ci_lower is None
    assert result.frequency_ci_upper is None
    assert result.frequency_ci_method is None


@pytest.mark.parametrize("level", [0, 1, -0.1, 1.1, math.nan, math.inf, True])
def test_invalid_confidence_level(level):
    with pytest.raises(ValueError, match="Confidence level"):
        frequency_interval(1, 2, level)


@pytest.mark.parametrize(("count", "total"), [(-1, 2), (1, -2), (3, 2)])
def test_invalid_count_range(count, total):
    with pytest.raises(ValueError, match="Count"):
        frequency_interval(count, total, None)


@pytest.mark.parametrize("count", [True, 1.5, "1"])
def test_noninteger_counts(count):
    with pytest.raises(TypeError, match="integers"):
        frequency_interval(count, 2, 0.95)
    with pytest.raises(TypeError, match="integers"):
        diversity([count])


def test_negative_diversity_count():
    with pytest.raises(ValueError, match="nonnegative"):
        diversity([2, -1])


def test_extreme_valid_confidence_level():
    result = frequency_interval(1, 2, math.nextafter(1.0, 0.0))
    assert result.frequency_ci_lower is not None
    assert result.frequency_ci_upper is not None
    assert 0 <= result.frequency_ci_lower < result.frequency_ci_upper <= 1
