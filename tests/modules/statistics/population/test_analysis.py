"""Synthetic approved-cohort regression tests; never read production sample data."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import UTC, datetime
from typing import TYPE_CHECKING

import pytest
from Bio import SeqIO

from src.config import get_settings
from src.core.models import Sample, Tool, Variant
from src.core.sample import generate_sequence
from src.modules.statistics.population.analysis import analyze_population
from src.modules.statistics.population.cohort import load_cohort
from src.modules.statistics.population.models import (
    AnalysisOptions,
    PopulationCohort,
    PopulationResult,
    PopulationSummary,
    ReferenceVariantAnnotation,
)
from src.modules.statistics.population.reporting import write_reports

if TYPE_CHECKING:
    from pathlib import Path

TIMESTAMP = datetime(2025, 1, 2, 3, 4, 5, tzinfo=UTC)


@pytest.fixture(scope="module")
def reference() -> str:
    return str(SeqIO.read(get_settings().directories.ref / "rCRS.fasta", "fasta").seq).upper()


def _variant(reference: str, position: int, alternate: str | None = None) -> Variant:
    ref = reference[position - 1]
    return Variant(pos=position, ref=ref, seq=alternate or next(base for base in "ACGT" if base != ref))


def _sample(sample_id: str, variants: list[Variant] | None = None) -> Sample:
    sample = Sample(
        sample_id=sample_id,
        source_tool=Tool.UNKNOWN,
        variants=variants or [],
        intervals={region: [list(bounds)] for region, bounds in get_settings().regions.REGIONS.items()},
    )
    generated = generate_sequence(sample)
    sample.hv1 = "".join(generated.hv_seqs["HV1"]["seq"])
    sample.hv2 = "".join(generated.hv_seqs["HV2"]["seq"])
    sample.hv3 = "".join(generated.hv_seqs["HV3"]["seq"])
    return sample


def _analyze(samples: list[Sample], options: AnalysisOptions | None = None) -> PopulationResult:
    return analyze_population(PopulationCohort(samples=samples), analysis_timestamp=TIMESTAMP, options=options)


def _copies(sample: Sample, count: int) -> list[Sample]:
    return [sample.model_copy(update={"sample_id": f"{sample.sample_id}-{index:04d}"}) for index in range(count)]


def _export(path: Path, samples: list[Sample]) -> None:
    path.write_text(
        json.dumps({sample.sample_id: sample.model_dump(mode="json") for sample in samples}), encoding="utf-8"
    )


@pytest.mark.parametrize("counts", [(1, 1, 1, 1), (4,), (3, 2, 1)])
def test_complete_haplotype_count_integration(reference: str, counts: tuple[int, ...]) -> None:
    samples = [
        sample
        for index, count in enumerate(counts)
        for sample in _copies(_sample(f"group-{index}", [_variant(reference, 100 + index)]), count)
    ]
    result = _analyze(samples)
    metrics = result.summary.haplotypes
    n = sum(counts)
    rmp = sum(count * (count - 1) for count in counts) / (n * (n - 1))
    assert sorted(row.count for row in result.haplotypes) == sorted(counts)
    assert result.summary.n_samples_complete_haplotype == n
    assert metrics.n_distinct_haplotypes == len(counts)
    assert metrics.n_singleton_haplotypes == counts.count(1)
    assert metrics.n_recurrent_haplotypes == sum(count > 1 for count in counts)
    assert metrics.n_individuals_with_unique_haplotype == counts.count(1)
    assert metrics.empirical_random_match_probability == pytest.approx(rmp)
    assert metrics.discrimination_power == pytest.approx(1 - rmp)
    assert metrics.haplotype_diversity == pytest.approx(1 - rmp)
    assert sum(row.frequency for row in result.haplotypes) == pytest.approx(1)


def test_callable_denominator_is_eight_not_ten(reference: str) -> None:
    samples = _copies(_sample("carrier", [_variant(reference, 100)]), 2)
    samples += _copies(_sample("reference"), 6)
    samples += _copies(_sample("unresolved", [_variant(reference, 100, "N")]), 2)
    result = _analyze(samples)
    (row,) = result.variants
    assert result.summary.n_samples_analyzed == 10
    assert (row.carrier_count, row.callable_sample_count, row.population_frequency) == (2, 8, 0.25)
    assert row.is_doubleton
    assert row.carrier_count_le_2
    assert not row.is_singleton
    assert result.summary.n_samples_incomplete_haplotype == 2


def test_numeric_variant_sort_and_identical_duplicates(reference: str) -> None:
    variants = [
        _variant(reference, 16100),
        Variant(pos="309.10", ref="-", seq="T"),
        _variant(reference, 450),
        Variant(pos="309.9", ref="-", seq="C"),
        _variant(reference, 100),
    ]
    first = _sample("z", [*variants, variants[-1], variants[1]])
    second = _sample("a", list(reversed(variants)))
    result = _analyze([first, second])
    assert [row.position for row in result.variants] == [100, "309.9", "309.10", 450, 16100]
    assert all(row.carrier_count == 2 for row in result.variants)
    assert len(result.haplotypes) == 1
    assert result.haplotypes[0].count == 2
    assert [row.sample_id for row in result.sample_haplotypes] == ["a", "z"]
    assert result == _analyze([second, first])


def test_incomplete_hv3_excluded_only_from_combined_and_hv3(reference: str) -> None:
    complete = _sample("complete", [_variant(reference, 100)])
    incomplete = complete.model_copy(update={"sample_id": "incomplete"}, deep=True)
    assert incomplete.intervals is not None
    incomplete.intervals["HV3"] = []
    result = _analyze([complete, incomplete])
    assert result.summary.n_samples_complete_haplotype == 1
    assert result.summary.n_samples_incomplete_haplotype == 1
    assert result.summary.region_haplotypes["HV1"].n_complete_profiles == 2
    assert result.summary.region_haplotypes["HV2"].n_complete_profiles == 2
    assert result.summary.region_haplotypes["HV3"].n_complete_profiles == 1
    mapping = result.sample_haplotypes[1]
    assert mapping.hv1_haplotype is not None
    assert mapping.hv3_haplotype is None
    assert mapping.combined_haplotype is None
    assert mapping.haplotype_id is None
    assert "HV3:missing_intervals" in mapping.incomplete_reasons
    assert result.variants[0].callable_sample_count == 2


def test_hv1_identity_does_not_override_hv3_difference(reference: str) -> None:
    result = _analyze([_sample("ref"), _sample("changed", [_variant(reference, 450)])])
    assert result.summary.region_haplotypes["HV1"].n_distinct_haplotypes == 1
    assert result.summary.region_haplotypes["HV3"].n_distinct_haplotypes == 2
    assert result.summary.haplotypes.n_distinct_haplotypes == 2
    assert len({row.haplotype_id for row in result.sample_haplotypes}) == 2


@pytest.mark.parametrize("count", [0, 1])
def test_empty_and_single_profile_have_undefined_pairwise_metrics(count: int) -> None:
    result = _analyze(_copies(_sample("reference"), count))
    assert result.summary.n_samples_analyzed == count
    assert result.summary.haplotypes.n_complete_profiles == count
    assert result.summary.haplotypes.empirical_random_match_probability is None
    assert result.summary.haplotypes.discrimination_power is None
    assert result.summary.haplotypes.haplotype_diversity is None
    assert len(result.haplotypes) == count
    assert result.variants == []


def test_missing_intervals_are_not_inferred_from_resolved_sequences() -> None:
    sample = _sample("no-intervals")
    sample.intervals = None
    result = _analyze([sample])
    assert result.summary.n_samples_complete_haplotype == 0
    assert result.haplotypes == []
    assert result.sample_haplotypes[0].incomplete_reasons == [
        "HV1:missing_intervals",
        "HV2:missing_intervals",
        "HV3:missing_intervals",
    ]


def test_misreferenced_variant_never_counts_and_blocks_only_its_position(reference: str) -> None:
    wrong_ref = next(base for base in "ACGT" if base != reference[99])
    corrupt = _sample("corrupt", [Variant(pos=100, ref=wrong_ref, seq="C")])
    result = _analyze([corrupt, _sample("reference")])
    assert result.variants == []
    assert result.summary.n_samples_complete_haplotype == 1
    mapping = result.sample_haplotypes[0]
    assert not mapping.complete_hv_profile
    assert "HV2:unverifiable_variant_representation" in mapping.incomplete_reasons


def test_conflicting_and_reference_equal_entries_leave_positions_unverifiable(reference: str) -> None:
    conflicting = _sample("conflict", [_variant(reference, 150, "C"), _variant(reference, 150, "G")])
    reference_equal = _sample("self", [Variant(pos=150, ref=reference[149], seq=reference[149])])
    clean = _sample("clean", [_variant(reference, 200)])
    result = _analyze([conflicting, reference_equal, clean])
    assert {row.position for row in result.variants} == {200}
    assert result.variants[0].carrier_count == 1
    assert result.summary.n_samples_complete_haplotype == 1
    reasons = {mapping.sample_id: mapping.incomplete_reasons for mapping in result.sample_haplotypes}
    for sample_id in ("conflict", "self"):
        assert "HV2:unverifiable_variant_representation" in reasons[sample_id]


def test_misreferenced_insertion_leaves_its_anchor_unverifiable() -> None:
    corrupt = _sample("corrupt", [Variant(pos="200.1", ref="A", seq="C")])
    result = _analyze([corrupt, _sample("reference")])
    assert result.variants == []
    assert result.summary.n_samples_complete_haplotype == 1
    assert "HV2:unverifiable_variant_representation" in result.sample_haplotypes[0].incomplete_reasons


def test_identical_duplicate_variants_count_once(reference: str) -> None:
    variants = [_variant(reference, 100), _variant(reference, 450)]
    duplicated = _sample("duplicated", [*variants, *reversed(variants)])
    result = _analyze([duplicated])
    assert {row.position for row in result.variants} == {100, 450}
    assert all(row.carrier_count == 1 for row in result.variants)
    assert len(result.haplotypes) == 1


def test_unresolved_allele_blocks_only_its_own_position(reference: str) -> None:
    carrier = _sample(
        "carrier",
        [_variant(reference, 104), _variant(reference, 105, "N"), _variant(reference, 106)],
    )
    result = _analyze([carrier, _sample("reference")])
    rows = {row.position: row for row in result.variants}
    assert set(rows) == {104, 106}
    assert rows[104].callable_sample_count == 2
    assert rows[106].callable_sample_count == 2
    assert result.summary.n_samples_complete_haplotype == 1
    assert result.sample_haplotypes[0].sample_id == "carrier"
    assert not result.sample_haplotypes[0].complete_hv_profile


def test_unresolved_insertion_removes_its_anchor() -> None:
    clean = _sample("clean", [Variant(pos="200.1", ref="-", seq="C")])
    murky = _sample("murky", [Variant(pos="200.1", ref="-", seq="c")])
    result = _analyze([clean, murky, _sample("reference")])
    rows = {row.position: row for row in result.variants}
    assert set(rows) == {"200.1"}
    assert rows["200.1"].carrier_count == 1
    assert rows["200.1"].callable_sample_count == 2
    assert result.summary.n_samples_complete_haplotype == 2


def test_iupac_is_an_exact_resolved_allele_not_fractional_c_or_t(reference: str) -> None:
    position = next(pos for pos in range(100, 120) if reference[pos - 1] == "T")
    result = _analyze(
        [
            _sample("Y1", [_variant(reference, position, "Y")]),
            _sample("Y2", [_variant(reference, position, "Y")]),
            _sample("C", [_variant(reference, position, "C")]),
            _sample("T"),
        ]
    )
    assert result.summary.n_samples_complete_haplotype == 4
    assert result.summary.haplotypes.n_distinct_haplotypes == 3
    assert {row.alternate: row.carrier_count for row in result.variants} == {"C": 1, "Y": 2}
    assert {row.alternate: row.population_frequency for row in result.variants} == {"C": 0.25, "Y": 0.5}


@pytest.mark.parametrize("alternate", ["a", "c", "g", "t", "y", "n", "N"])
def test_lowercase_and_n_are_unresolved(reference: str, alternate: str) -> None:
    result = _analyze([_sample("unresolved", [_variant(reference, 100, alternate)]), _sample("ref")])
    assert result.summary.n_samples_complete_haplotype == 1
    assert result.variants == []
    assert not result.sample_haplotypes[1].complete_hv_profile


@pytest.mark.parametrize(
    ("denominator", "field"), [(1000, "frequency_lt_0_001"), (200, "frequency_lt_0_005"), (100, "frequency_lt_0_01")]
)
@pytest.mark.parametrize("offset", [-1, 0, 1])
def test_frequency_thresholds_are_strict(reference: str, denominator: int, field: str, offset: int) -> None:
    samples = [_sample("carrier", [_variant(reference, 100)])]
    samples += _copies(_sample("reference"), denominator + offset - 1)
    threshold = 1 / denominator
    result = _analyze(samples, AnalysisOptions(confidence_level=None, additional_frequency_thresholds=[threshold]))
    (row,) = result.variants
    assert row.population_frequency == 1 / (denominator + offset)
    assert getattr(row, field) is (offset > 0)
    assert row.additional_frequency_flags[str(threshold)] is (offset > 0)
    assert result.summary.variants["all"].additional_frequency_counts[str(threshold)] == int(offset > 0)


def test_duplicate_json_sample_ids_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "approved.json"
    path.write_text('{"duplicate":{"variants":[]},"duplicate":{"variants":[]}}', encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate JSON"):
        load_cohort(path)


def test_duplicate_in_memory_sample_ids_are_rejected() -> None:
    sample = _sample("duplicate")
    with pytest.raises(ValueError, match="Duplicate sample IDs"):
        _analyze([sample, sample])


@pytest.mark.parametrize(
    "groups",
    [{"unknown": []}, {"snps": [], "HV1": []}, {"snps": {}}, {"HV1": "bad"}, {"snps": [None]}, {"HV1": [{"pos": 100}]}],
)
def test_malformed_variant_groups_are_rejected(tmp_path: Path, groups: object) -> None:
    path = tmp_path / "approved.json"
    path.write_text(json.dumps({"sample": {"variants": groups}}), encoding="utf-8")
    with pytest.raises(ValueError, match=r"variant|Variant|canonical"):
        load_cohort(path)


@pytest.mark.parametrize("group", ["snps", "HV2"])
def test_canonical_grouped_exports_load(tmp_path: Path, reference: str, group: str) -> None:
    sample = _sample("sample", [_variant(reference, 100)])
    raw = sample.model_dump(mode="json")
    raw["variants"] = {group: raw["variants"]}
    path = tmp_path / "approved.json"
    path.write_text(json.dumps({sample.sample_id: raw}), encoding="utf-8")
    cohort = load_cohort(path)
    assert cohort.samples[0].variants == sample.variants
    assert _analyze(cohort.samples).summary.n_samples_complete_haplotype == 1


def test_loading_refreshes_changed_export_without_mutating_old_snapshot(tmp_path: Path, reference: str) -> None:
    path = tmp_path / "approved.json"
    _export(path, [_sample("sample")])
    before = load_cohort(path)
    before_result = analyze_population(before, analysis_timestamp=TIMESTAMP)
    _export(path, [_sample("sample", [_variant(reference, 100)])])
    after = load_cohort(path)
    after_result = analyze_population(after, analysis_timestamp=TIMESTAMP)
    assert before.input_sha256 != after.input_sha256
    assert after.input_sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    assert before_result.variants == []
    assert before.samples[0].variants == []
    assert len(after_result.variants) == 1
    assert before_result.haplotypes[0].haplotype_id != after_result.haplotypes[0].haplotype_id


def test_reports_are_deterministic_and_round_trip_full_precision(tmp_path: Path, reference: str) -> None:
    samples = [_sample("carrier", [_variant(reference, 100)]), *_copies(_sample("reference"), 2)]
    options = AnalysisOptions(presentation_precision=2)
    result = _analyze(samples, options)
    reordered = _analyze(list(reversed(samples)), options)
    assert PopulationResult.model_validate_json(result.model_dump_json()) == result
    first = write_reports(result, tmp_path / "first")
    second = write_reports(reordered, tmp_path / "second")
    assert {path.name: path.read_bytes() for path in first} == {path.name: path.read_bytes() for path in second}
    summary = PopulationSummary.model_validate_json((tmp_path / "first" / "mtDNA_population_summary.json").read_text())
    assert summary == result.summary
    assert summary.haplotypes.empirical_random_match_probability == 1 / 3
    for filename, digest in summary.output_sha256.items():
        assert hashlib.sha256((tmp_path / "first" / filename).read_bytes()).hexdigest() == digest
    with (tmp_path / "first" / "mtDNA_variant_population_frequencies.tsv").open(newline="", encoding="utf-8") as handle:
        (row,) = list(csv.DictReader(handle, delimiter="\t"))
    assert float(row["population_frequency"]) == 1 / 3
    assert float(row["frequency_ci_lower"]) == result.variants[0].frequency_ci_lower
    assert float(row["frequency_ci_upper"]) == result.variants[0].frequency_ci_upper
    assert summary.analysis_timestamp == TIMESTAMP.isoformat()


def test_adjacent_intervals_do_not_prove_insertion_anchor_callability() -> None:
    samples = [_sample("reference"), _sample("insertion", [Variant(pos="200.1", ref="-", seq="C")])]
    for sample in samples:
        assert sample.intervals is not None
        start, end = get_settings().regions.REGIONS["HV2"]
        sample.intervals["HV2"] = [[start, 200], [201, end]]
        generated = generate_sequence(sample)
        sample.hv2 = "".join(generated.hv_seqs["HV2"]["seq"])
        assert "N" not in sample.hv2
    result = _analyze(samples)
    assert result.summary.n_samples_complete_haplotype == 0
    assert result.summary.region_haplotypes["HV2"].n_complete_profiles == 0
    assert result.haplotypes == []
    assert result.variants == []
    assert all(mapping.haplotype_id is None for mapping in result.sample_haplotypes)
    assert all(not mapping.region_complete["HV2"] for mapping in result.sample_haplotypes)


class _ReferenceEvidence:
    def annotate(
        self,
        variant: Variant,
        *,
        reference_sha256: str,
        region_definitions: dict[str, list[int]],
    ) -> dict[str, ReferenceVariantAnnotation]:
        assert variant.pos == 100
        assert len(reference_sha256) == 64
        assert set(region_definitions) == {"HV1", "HV2", "HV3"}
        return {
            "unknown": ReferenceVariantAnnotation(observed=None),
            "absent": ReferenceVariantAnnotation(observed=False),
            "present": ReferenceVariantAnnotation(observed=True, reference_frequency=0.8),
        }


def test_only_explicit_reference_false_counts_as_not_observed(reference: str) -> None:
    result = analyze_population(
        PopulationCohort(samples=[_sample("carrier", [_variant(reference, 100)]), _sample("reference")]),
        reference_provider=_ReferenceEvidence(),
        analysis_timestamp=TIMESTAMP,
    )
    assert result.summary.reference_not_observed_counts == {"absent": 1, "present": 0, "unknown": 0}
    assert result.variants[0].reference_annotations["unknown"].observed is None
    assert result.variants[0].population_frequency == 0.5
    assert result.variants[0].reference_annotations["present"].reference_frequency == 0.8


def test_trace_peak_fractions_do_not_replace_carrier_frequency(reference: str) -> None:
    position = next(pos for pos in range(100, 120) if reference[pos - 1] == "T")
    low_c = _variant(reference, position, "Y").model_copy(update={"peaks": [[0, 10, 0, 90]], "quality": [20]})
    high_c = low_c.model_copy(update={"peaks": [[0, 90, 0, 10]], "quality": [40]})
    result = _analyze([_sample("low", [low_c]), _sample("high", [high_c]), _sample("reference")])
    (row,) = result.variants
    assert row.alternate == "Y"
    assert row.carrier_count == 2
    assert row.population_frequency == 2 / 3
    assert sorted(haplotype.count for haplotype in result.haplotypes) == [1, 2]


@pytest.mark.parametrize("position", [True, False])
def test_export_variant_positions_reject_coercion(tmp_path: Path, position: object) -> None:
    path = tmp_path / "approved.json"
    path.write_text(json.dumps({"sample": {"variants": [{"pos": position, "ref": "A", "seq": "C"}]}}), encoding="utf-8")
    with pytest.raises(TypeError, match="positions"):
        load_cohort(path)


@pytest.mark.parametrize(("numeric_literal", "insertion"), [("309.1", "309.1"), ("309.10", "309.10")])
def test_export_numeric_insertion_tokens_decode_losslessly(
    tmp_path: Path, numeric_literal: str, insertion: str
) -> None:
    numeric_export = tmp_path / "numeric.json"
    numeric_export.write_text(
        '{"sample": {"variants": [{"pos": ' + numeric_literal + ', "ref": "-", "seq": "C"}]}}',
        encoding="utf-8",
    )
    string_export = tmp_path / "string.json"
    string_export.write_text(
        json.dumps({"sample": {"variants": [{"pos": insertion, "ref": "-", "seq": "C"}]}}),
        encoding="utf-8",
    )
    numeric_cohort = load_cohort(numeric_export)
    string_cohort = load_cohort(string_export)
    assert [variant.pos for variant in numeric_cohort.samples[0].variants] == [
        variant.pos for variant in string_cohort.samples[0].variants
    ]
    assert numeric_cohort.samples[0].variants[0].pos == insertion


def test_export_zero_insertion_index_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "approved.json"
    path.write_text('{"sample": {"variants": [{"pos": 100.0, "ref": "-", "seq": "C"}]}}', encoding="utf-8")
    samples = load_cohort(path).samples
    with pytest.raises(ValueError, match="Insertion indices must be positive"):
        _analyze(samples)


@pytest.mark.parametrize("endpoint", [True, False, 73.0, 73.5])
@pytest.mark.parametrize("index", [0, 1])
def test_export_interval_endpoints_reject_coercion(tmp_path: Path, endpoint: object, index: int) -> None:
    interval: list[object] = [73, 340]
    interval[index] = endpoint
    path = tmp_path / "approved.json"
    path.write_text(json.dumps({"sample": {"variants": [], "intervals": {"HV2": [interval]}}}), encoding="utf-8")
    with pytest.raises(TypeError, match="endpoints"):
        load_cohort(path)


@pytest.mark.parametrize("region", ["HV1", "HV2", "HV3"])
def test_export_hv_strings_are_neither_required_nor_consumed(tmp_path: Path, region: str) -> None:
    path = tmp_path / "approved.json"
    path.write_text(
        json.dumps({"sample": {"variants": [], "intervals": {region: [list(get_settings().regions.REGIONS[region])]}}}),
        encoding="utf-8",
    )
    cohort = load_cohort(path)
    assert getattr(cohort.samples[0], region.lower()) is None
    result = _analyze(cohort.samples)
    assert result.summary.region_haplotypes[region].n_complete_profiles == 1
    assert result.summary.n_samples_complete_haplotype == 0
