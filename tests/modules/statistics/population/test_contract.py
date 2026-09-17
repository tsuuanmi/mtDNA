"""Public report and strict representation contracts using synthetic exports."""

from __future__ import annotations

import csv
import json
from datetime import UTC, datetime
from typing import TYPE_CHECKING

import pytest

from src.modules.statistics.population import analyze_population, load_cohort, write_reports

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize("group", [[], {}, {"snps": []}, {"HV1": [], "HV2": [], "HV3": []}])
def test_entries_without_intervals_never_become_reference_profiles(tmp_path: Path, group: object) -> None:
    path = tmp_path / "input.json"
    path.write_text(json.dumps({"001": {"variants": group}}), encoding="utf-8")
    result = analyze_population(load_cohort(path))
    assert result.summary.n_samples_analyzed == 1
    assert result.summary.n_samples_complete_haplotype == 0
    assert result.sample_haplotypes[0].sample_id == "001"
    assert result.sample_haplotypes[0].combined_haplotype is None
    assert result.haplotypes == []


def test_empty_reports_keep_headers_nulls_and_concise_executive(tmp_path: Path) -> None:
    path = tmp_path / "input.json"
    path.write_text("{}", encoding="utf-8")
    result = analyze_population(load_cohort(path), analysis_timestamp=datetime(2025, 1, 1, tzinfo=UTC))
    outputs = write_reports(result, tmp_path / "output")
    assert len(outputs) == 7
    for output in outputs:
        if output.suffix == ".tsv":
            with output.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                assert reader.fieldnames
                assert list(reader) == []
    summary = json.loads((tmp_path / "output" / "mtDNA_population_summary.json").read_text(encoding="utf-8"))
    assert summary["haplotypes"]["empirical_random_match_probability"] is None
    executive = (tmp_path / "output" / "mtDNA_population_executive_summary.txt").read_text(encoding="utf-8")
    assert "unavailable" in executive
    assert "matern" in executive.lower()
    assert len(executive.splitlines()) <= 40


def test_report_directory_cannot_overwrite_existing_files(tmp_path: Path) -> None:
    path = tmp_path / "input.json"
    path.write_text("{}", encoding="utf-8")
    result = analyze_population(load_cohort(path))
    with pytest.raises(ValueError, match="new or empty"):
        write_reports(result, tmp_path)
    assert path.read_text(encoding="utf-8") == "{}"


def test_private_metadata_is_not_loaded_or_exported(tmp_path: Path) -> None:
    path = tmp_path / "input.json"
    path.write_text(
        json.dumps({"sample": {"variants": [], "information": {"name": "private-marker"}}}), encoding="utf-8"
    )
    cohort = load_cohort(path)
    assert cohort.samples[0].information is None
    result = analyze_population(cohort)
    paths = write_reports(result, tmp_path / "output")
    assert all("private-marker" not in path.read_text(encoding="utf-8") for path in paths)
