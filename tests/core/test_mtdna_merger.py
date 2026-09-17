"""Tests for the single-pass mtDNA merger."""

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from src.core.models import Sample, Tool, Variant
from src.core.mtdna_merger import MergeResult, MergeValidationError, MtDnaMerger, main, merge_by_regions, merge_samples

SAMPLE_ID = "SAMPLE_001"


def _profile_payload(variants: dict, intervals: dict) -> dict:
    """Return unwrapped per-sample profile (no top-level sample_id key)."""
    return {
        "variants": variants,
        "intervals": intervals,
    }


def _write_json(path: Path, payload: dict) -> None:
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh)


def _base_intervals() -> dict:
    return {
        "HV1": [[16024, 16024]],
        "HV2": [[68, 300]],  # enlarged to allow SNP testing at position 73
        "HV3": [[438, 438]],
    }


def test_merge_writes_merged_output(tmp_path, ref_sequence) -> None:
    """Successful merge of two compatible profiles (identical alleles)."""
    file_a = tmp_path / "a.json"
    file_b = tmp_path / "b.json"

    variants_a = {
        "snps": [],
        "insertions": [{"pos": 68.1, "ref": "-", "seq": "A", "file": "a.json"}],
        "deletions": [],
    }
    variants_b = {
        "snps": [],
        "insertions": [{"pos": 68.1, "ref": "-", "seq": "A", "file": "b.json"}],
        "deletions": [],
    }

    _write_json(file_a, _profile_payload(variants_a, _base_intervals()))
    _write_json(file_b, _profile_payload(variants_b, _base_intervals()))

    merger = MtDnaMerger(ref_path=ref_sequence)
    result = merger.merge([file_a, file_b], SAMPLE_ID, output_dir=tmp_path)

    assert isinstance(result, MergeResult)
    assert result.warnings == []
    assert result.sample.sample_id == SAMPLE_ID
    assert len(result.sample.variants) == 1
    assert result.sample.variants[0].seq == "A"

    per_sample_path = tmp_path / SAMPLE_ID / f"{SAMPLE_ID}.json"
    assert per_sample_path.exists()

    # Verify the per-sample JSON has v4 format keys
    with per_sample_path.open() as f:
        sample_dict = json.load(f)
    assert "variants" in sample_dict
    assert "intervals" in sample_dict
    assert "HV1" in sample_dict
    assert "HV2" in sample_dict
    assert "HV3" in sample_dict


def test_merge_raises_on_conflicting_alleles(tmp_path, ref_sequence) -> None:
    """Explicit conflicting alleles → immediate hard error."""
    file_a = tmp_path / "a.json"
    file_b = tmp_path / "b.json"

    variants_a = {
        "snps": [],
        "insertions": [{"pos": 68.1, "ref": "-", "seq": "A", "file": "a.json"}],
        "deletions": [],
    }
    variants_b = {
        "snps": [],
        "insertions": [{"pos": 68.1, "ref": "-", "seq": "C", "file": "b.json"}],
        "deletions": [],
    }

    _write_json(file_a, _profile_payload(variants_a, _base_intervals()))
    _write_json(file_b, _profile_payload(variants_b, _base_intervals()))

    merger = MtDnaMerger(ref_path=ref_sequence)
    with pytest.raises(MergeValidationError) as excinfo:
        merger.merge([file_a, file_b], SAMPLE_ID, output_dir=tmp_path)

    assert "Conflicting alleles at position 68.1" in str(excinfo.value)
    assert not (tmp_path / SAMPLE_ID / f"{SAMPLE_ID}.json").exists()


def test_merge_raises_on_ref_allele_mismatch(tmp_path, ref_sequence) -> None:
    """Ref allele does not match rCRS → hard error."""
    file_a = tmp_path / "a.json"

    # Position 73: rCRS = A, we deliberately set wrong ref = C
    bad_variants = {
        "snps": [{"pos": 73, "ref": "C", "seq": "G", "file": "a.json"}],
        "insertions": [],
        "deletions": [],
    }
    _write_json(file_a, _profile_payload(bad_variants, _base_intervals()))

    merger = MtDnaMerger(ref_path=ref_sequence)
    with pytest.raises(MergeValidationError) as excinfo:
        merger.merge([file_a], SAMPLE_ID, output_dir=tmp_path)

    assert "ref allele 'C' does not match rCRS" in str(excinfo.value)
    assert not (tmp_path / SAMPLE_ID / f"{SAMPLE_ID}.json").exists()


def test_merge_raises_on_variant_outside_own_interval(tmp_path, ref_sequence) -> None:
    """Variant lies outside the profile's own declared intervals → hard error."""
    file_a = tmp_path / "out_of_interval.json"
    bad_variants = {
        "snps": [{"pos": 1000, "ref": "T", "seq": "G", "file": "bad.json"}],
        "insertions": [],
        "deletions": [],
    }
    _write_json(file_a, _profile_payload(bad_variants, _base_intervals()))

    merger = MtDnaMerger(ref_path=ref_sequence)
    with pytest.raises(MergeValidationError) as excinfo:
        merger.merge([file_a], SAMPLE_ID, output_dir=tmp_path)

    assert "lies outside the profile's declared intervals" in str(excinfo.value)
    assert not (tmp_path / SAMPLE_ID / f"{SAMPLE_ID}.json").exists()


def test_merge_raises_on_validation_warnings(tmp_path, ref_sequence) -> None:
    """Invalid data (e.g. bad base 'N') should raise validation error."""
    file_a = tmp_path / "bad.json"
    bad_variants = {
        "snps": [{"pos": 73, "ref": "A", "seq": "N", "file": "bad.json"}],
        "insertions": [],
        "deletions": [],
    }
    _write_json(file_a, _profile_payload(bad_variants, _base_intervals()))

    merger = MtDnaMerger(ref_path=ref_sequence)
    with pytest.raises(MergeValidationError):
        merger.merge([file_a], SAMPLE_ID, output_dir=tmp_path)

    assert not (tmp_path / SAMPLE_ID / f"{SAMPLE_ID}.json").exists()


def test_cli_exits_nonzero_when_merge_fails(tmp_path, ref_sequence, monkeypatch) -> None:
    """CLI should exit with code 1 on any merge/validation error."""
    file_a = tmp_path / "bad.json"
    bad_variants = {
        "snps": [{"pos": 73, "ref": "A", "seq": "N", "file": "bad.json"}],
        "insertions": [],
        "deletions": [],
    }
    _write_json(file_a, _profile_payload(bad_variants, _base_intervals()))

    monkeypatch.setattr(
        "sys.argv",
        [
            "mtdna_merger",
            "--files",
            str(file_a),
            "--sample-id",
            SAMPLE_ID,
            "--ref-path",
            str(ref_sequence),
            "--output-dir",
            str(tmp_path),
        ],
    )

    with pytest.raises(SystemExit) as exc:
        main()

    assert exc.value.code == 1


# ----------------------------------------------------------------------
# Tests for unified source_tool and region_sources provenance tracking
# ----------------------------------------------------------------------


def test_merge_sets_source_tool_to_unified(tmp_path, ref_sequence) -> None:
    """Merged result must use Tool.UNIFIED regardless of input tools."""

    file_a = tmp_path / "a.json"
    file_b = tmp_path / "b.json"

    variants = {
        "snps": [{"pos": 73, "ref": "A", "seq": "G", "file": "a.json"}],
        "insertions": [],
        "deletions": [],
    }

    payload_a = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "sequencher"},
    }
    payload_b = {
        "variants": {
            "snps": [],
            "insertions": [],
            "deletions": [],
        },
        "intervals": _base_intervals(),
        "information": {"source_tool": "tracy"},
    }

    _write_json(file_a, payload_a)
    _write_json(file_b, payload_b)

    # We need both profiles to have compatible variants/overlap
    # File b has no variants but covers the same position — that's an implicit conflict.
    # Instead, let's have b provide the same variant.
    payload_b["variants"] = variants

    _write_json(file_b, payload_b)

    merger = MtDnaMerger(ref_path=ref_sequence)
    result = merger.merge([file_a, file_b], SAMPLE_ID, output_dir=tmp_path)

    assert result.sample.source_tool == Tool.UNIFIED
    info = result.sample.information or {}
    assert info.get("source_tool") == "unified"


def test_merge_populates_region_sources_with_tool_names(tmp_path, ref_sequence) -> None:
    """Merged result must include region_sources mapping regions to source tool names."""
    file_a = tmp_path / "a.json"
    file_b = tmp_path / "b.json"

    variants = {
        "snps": [{"pos": 73, "ref": "A", "seq": "G", "file": "a.json"}],
        "insertions": [],
        "deletions": [],
    }

    payload_a = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "sequencher"},
    }
    payload_b = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "tracy"},
    }

    _write_json(file_a, payload_a)
    _write_json(file_b, payload_b)

    merger = MtDnaMerger(ref_path=ref_sequence)
    result = merger.merge([file_a, file_b], SAMPLE_ID, output_dir=tmp_path)

    info = result.sample.information or {}
    region_sources = info.get("region_sources", {})
    # Whole-sample merge maps all regions to all contributing tools
    assert "HV1" in region_sources
    assert "HV2" in region_sources
    assert "HV3" in region_sources
    assert "sequencher" in region_sources["HV1"]
    assert "tracy" in region_sources["HV1"]


def test_merge_region_sources_no_duplicate_tools(tmp_path, ref_sequence) -> None:
    """When both files report the same source tool, region_sources should not duplicate it."""

    file_a = tmp_path / "a.json"
    file_b = tmp_path / "b.json"

    variants = {
        "snps": [{"pos": 73, "ref": "A", "seq": "G", "file": "a.json"}],
        "insertions": [],
        "deletions": [],
    }

    payload_a = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "tracy"},
    }
    payload_b = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "tracy"},
    }

    _write_json(file_a, payload_a)
    _write_json(file_b, payload_b)

    merger = MtDnaMerger(ref_path=ref_sequence)
    result = merger.merge([file_a, file_b], SAMPLE_ID, output_dir=tmp_path)

    assert result.sample.source_tool == Tool.UNIFIED
    info = result.sample.information or {}
    region_sources = info.get("region_sources", {})
    # Same tool deduplicated: each region maps to ["tracy"]
    for region in ("HV1", "HV2", "HV3"):
        assert region_sources.get(region) == ["tracy"]


def test_merge_samples_sets_unified_source_tool_and_region_sources(ref_sequence) -> None:
    """merge_samples() must set Tool.UNIFIED and region_sources on the merged Sample."""

    s_a = Sample(
        sample_id="TEST_001",
        variants=[Variant(pos=73, ref="A", seq="G", files=["a.json"], quality=[30])],
        source_tool=Tool.SEQUENCHER,
        intervals={"HV1": [[16024, 16365]], "HV2": [[68, 300]], "HV3": [[438, 576]]},
        batch_id="BATCH_001",
    )
    s_b = Sample(
        sample_id="TEST_001",
        variants=[Variant(pos=73, ref="A", seq="G", files=["b.json"], quality=[25])],
        source_tool=Tool.TRACY,
        intervals={"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
        batch_id="BATCH_001",
    )

    tmp = tempfile.mkdtemp()
    try:
        result = merge_samples([s_a, s_b], ref_path=str(ref_sequence), output_dir=tmp)
        assert result.sample.source_tool == Tool.UNIFIED
        info = result.sample.information or {}
        assert info.get("source_tool") == "unified"
        region_sources = info.get("region_sources", {})
        # Whole-sample merge: all regions map to both tools
        assert "sequencher" in region_sources.get("HV1", [])
        assert "tracy" in region_sources.get("HV1", [])
    finally:
        shutil.rmtree(tmp)


def test_merge_json_output_contains_unified_and_region_sources(tmp_path, ref_sequence) -> None:
    """The merged JSON file must contain source_tool='unified' and region_sources in information."""
    file_a = tmp_path / "a.json"
    file_b = tmp_path / "b.json"

    variants = {
        "snps": [{"pos": 73, "ref": "A", "seq": "G", "file": "a.json"}],
        "insertions": [],
        "deletions": [],
    }

    payload_a = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "sequencher"},
    }
    payload_b = {
        "variants": variants,
        "intervals": _base_intervals(),
        "information": {"source_tool": "tracy"},
    }

    _write_json(file_a, payload_a)
    _write_json(file_b, payload_b)

    merger = MtDnaMerger(ref_path=ref_sequence)
    _ = merger.merge([file_a, file_b], SAMPLE_ID, output_dir=tmp_path)

    per_sample_path = tmp_path / SAMPLE_ID / f"{SAMPLE_ID}.json"
    with per_sample_path.open() as f:
        sample_dict = json.load(f)

    info = sample_dict["information"]
    assert info["source_tool"] == "unified"
    assert "region_sources" in info
    assert "HV1" in info["region_sources"]
    assert "HV2" in info["region_sources"]
    assert "HV3" in info["region_sources"]


def test_merge_by_regions_ms_hv1_seq_hv23(ref_sequence) -> None:
    """HV1 from MS + HV2+HV3 from Sequencher → correct region_sources."""

    ms_hv1 = Sample(
        sample_id="TEST_001",
        variants=[Variant(pos=16223, ref="C", seq="T", files=["ms.json"], quality=[30])],
        source_tool=Tool.MUTATION_SURVEYOR,
        intervals={"HV1": [[16024, 16365]]},
        batch_id="BATCH_001",
    )
    seq_hv23 = Sample(
        sample_id="TEST_001",
        variants=[
            Variant(pos=73, ref="A", seq="G", files=["seq.json"], quality=[25]),
            Variant(pos=455, ref="G", seq="A", files=["seq.json"], quality=[20]),
        ],
        source_tool=Tool.SEQUENCHER,
        intervals={"HV2": [[73, 340]], "HV3": [[438, 576]]},
        batch_id="BATCH_001",
    )

    tmp = tempfile.mkdtemp()
    try:
        result = merge_by_regions(
            region_samples={"HV1": [ms_hv1], "HV2-3": [seq_hv23]},
            sample_id="TEST_001",
            ref_path=str(ref_sequence),
            output_dir=tmp,
        )
        assert result.sample.source_tool == Tool.UNIFIED
        info = result.sample.information or {}
        assert info["source_tool"] == "unified"
        region_sources = info["region_sources"]
        assert region_sources["HV1"] == "mutation_surveyor"
        assert region_sources["HV2"] == "sequencher"
        assert region_sources["HV3"] == "sequencher"
        # Unified JSON has 3 variants total
        assert len(result.sample.variants) == 3
    finally:
        shutil.rmtree(tmp)


def test_merge_by_regions_single_tool(ref_sequence) -> None:
    """Only one tool → region_sources reflects single tool."""

    ms_hv1 = Sample(
        sample_id="TEST_001",
        variants=[Variant(pos=16223, ref="C", seq="T", files=["ms.json"], quality=[30])],
        source_tool=Tool.MUTATION_SURVEYOR,
        intervals={"HV1": [[16024, 16365]]},
        batch_id="BATCH_001",
    )

    tmp = tempfile.mkdtemp()
    try:
        result = merge_by_regions(
            region_samples={"HV1": [ms_hv1]},
            sample_id="TEST_001",
            ref_path=str(ref_sequence),
            output_dir=tmp,
        )
        info = result.sample.information or {}
        region_sources = info["region_sources"]
        assert region_sources["HV1"] == "mutation_surveyor"
    finally:
        shutil.rmtree(tmp)


def test_merge_no_merged_from_in_output(ref_sequence) -> None:
    """Merged output should NOT contain merged_from (replaced by region_sources)."""

    s_a = Sample(
        sample_id="TEST_001",
        variants=[Variant(pos=73, ref="A", seq="G", files=["a.json"], quality=[30])],
        source_tool=Tool.SEQUENCHER,
        intervals={"HV1": [[16024, 16365]], "HV2": [[68, 300]], "HV3": [[438, 576]]},
        batch_id="BATCH_001",
    )

    tmp = tempfile.mkdtemp()
    try:
        result = merge_samples([s_a], ref_path=str(ref_sequence), output_dir=tmp)
        info = result.sample.information or {}
        assert "merged_from" not in info
        assert "region_sources" in info
    finally:
        shutil.rmtree(tmp)
