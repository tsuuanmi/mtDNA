"""Tests for the FASTA validation module.

Covers:
- _merge_con_dict helper
- _build_alignment (extraction ranges, insertions, deletions, SNP, unread)
- _parse_fasta
- validate_sample (PASS/FAIL, alignment construction)
- write_report (output file creation)
- CLI exit codes
"""

import json
from pathlib import Path

import pytest

from src.core.sample import generate_sequence
from src.validation.validation_fasta import (
    FastaValidator,
    _merge_con_dict,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ref_path():
    """Return the path to the reference FASTA file."""
    return Path(__file__).resolve().parent.parent / "ref" / "rCRS.fasta"


@pytest.fixture
def sample_json_dir(tmp_path, sample_data):
    """Create a temporary JSON directory for one sample."""
    json_dir = tmp_path / "json"
    json_dir.mkdir()
    sample_id = "TEST_SAMPLE"
    sample_dir = json_dir / sample_id
    sample_dir.mkdir()
    json_path = sample_dir / f"{sample_id}.json"
    with json_path.open("w") as f:
        json.dump(sample_data, f)
    return json_dir


@pytest.fixture
def sample_data():
    """Return sample variant data for testing (SNPs only)."""
    return {
        "variants": {
            "HV1": [
                {"pos": "16223", "ref": "C", "seq": "T"},
            ],
            "HV2": [
                {"pos": "73", "ref": "A", "seq": "G"},
                {"pos": "263", "ref": "A", "seq": "G"},
            ],
            "HV3": [],
        },
        "intervals": {
            "HV1": [[16024, 16365]],
            "HV2": [[68, 340]],
            "HV3": [[438, 576]],
        },
    }


@pytest.fixture
def sample_data_with_insertion():
    """Return sample variant data including insertions."""
    return {
        "variants": {
            "HV1": [],
            "HV2": [
                {"pos": "73", "ref": "A", "seq": "G"},
                {"pos": "309.1", "ref": "-", "seq": "C"},
                {"pos": "315.1", "ref": "-", "seq": "C"},
            ],
            "HV3": [],
        },
        "intervals": {
            "HV1": [[16024, 16365]],
            "HV2": [[68, 340]],
            "HV3": [[438, 576]],
        },
    }


@pytest.fixture
def sample_data_with_deletion():
    """Return sample variant data including a deletion."""
    return {
        "variants": {
            "HV1": [],
            "HV2": [
                {"pos": "249", "ref": "A", "seq": "-"},
            ],
            "HV3": [],
        },
        "intervals": {
            "HV1": [[16024, 16365]],
            "HV2": [[68, 340]],
            "HV3": [[438, 576]],
        },
    }


@pytest.fixture
def sample_data_multi_ins():
    """Return sample data with a multi-base insertion (10 bases at 217)."""
    return {
        "variants": {
            "HV1": [],
            "HV2": [
                {"pos": "73", "ref": "A", "seq": "G"},
                {"pos": "217.1", "ref": "-", "seq": "T"},
                {"pos": "217.2", "ref": "-", "seq": "T"},
                {"pos": "217.3", "ref": "-", "seq": "A"},
                {"pos": "217.4", "ref": "-", "seq": "A"},
                {"pos": "217.5", "ref": "-", "seq": "T"},
                {"pos": "217.6", "ref": "-", "seq": "T"},
                {"pos": "217.7", "ref": "-", "seq": "G"},
                {"pos": "217.8", "ref": "-", "seq": "A"},
                {"pos": "217.9", "ref": "-", "seq": "T"},
                {"pos": "217.10", "ref": "-", "seq": "T"},
                {"pos": "263", "ref": "A", "seq": "G"},
                {"pos": "309.1", "ref": "-", "seq": "C"},
                {"pos": "309.2", "ref": "-", "seq": "C"},
                {"pos": "315.1", "ref": "-", "seq": "C"},
            ],
            "HV3": [],
        },
        "intervals": {
            "HV1": [[16024, 16365]],
            "HV2": [[73, 340]],
            "HV3": [[438, 576]],
        },
    }


# ---------------------------------------------------------------------------
# _merge_con_dict tests
# ---------------------------------------------------------------------------


class TestMergeConDict:
    """Tests for the _merge_con_dict helper."""

    def test_string_variant_overrides_int_ref(self):
        """String variant key should override integer ref key for SNP positions."""
        con_dict = {16129: "G", "16129": "A"}
        merged = _merge_con_dict(con_dict)
        assert merged[16129] == "A"

    def test_insertion_keys_preserved(self):
        """Insertion string keys like '217.1' should be preserved."""
        con_dict = {217: "A", "217.1": "T"}
        merged = _merge_con_dict(con_dict)
        assert merged[217] == "A"
        assert merged["217.1"] == "T"

    def test_int_keys_unchanged(self):
        """Integer keys without string variants should pass through."""
        con_dict = {100: "A", 200: "C"}
        merged = _merge_con_dict(con_dict)
        assert merged[100] == "A"
        assert merged[200] == "C"

    def test_mixed_keys(self):
        """Test a realistic mix of int keys, SNP string keys, and insertion keys."""
        con_dict = {
            73: "A",
            "73": "G",  # SNP
            263: "A",
            "263": "G",  # SNP
            "309.1": "C",  # Insertion
            "315.1": "C",  # Insertion
        }
        merged = _merge_con_dict(con_dict)
        assert merged[73] == "G"  # SNP overrides
        assert merged[263] == "G"  # SNP overrides
        assert merged["309.1"] == "C"  # Insertion preserved
        assert merged["315.1"] == "C"  # Insertion preserved


# ---------------------------------------------------------------------------
# _build_alignment tests
# ---------------------------------------------------------------------------


class TestBuildAlignment:
    """Tests for the _build_alignment method."""

    @pytest.fixture
    def validator(self, ref_path, tmp_path):
        """Create a minimal validator for alignment testing."""
        return FastaValidator(
            fasta_dir=tmp_path / "fasta",
            json_dir=tmp_path / "json",
            ref_path=ref_path,
            batch_id="TEST",
        )

    def test_alignment_uses_extraction_ranges(self, validator, ref_path, make_sample):
        """Alignment should use extraction_ranges, not REGIONS."""
        # Use HV2 with extraction range [73, 340] (not [68, 340])
        variants = [{"pos": "73", "ref": "A", "seq": "G"}, {"pos": "263", "ref": "A", "seq": "G"}]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        alignment = validator._build_alignment(ref_dict, merged_con, hv_seqs["HV2"]["ranges"])

        # Position 73 should be first (extraction range starts at 73, not 68)
        assert alignment.pos_markers[0] == 73
        # Position 68-72 should NOT appear
        assert 68 not in alignment.pos_markers
        assert 69 not in alignment.pos_markers
        assert 70 not in alignment.pos_markers
        assert 71 not in alignment.pos_markers
        assert 72 not in alignment.pos_markers

    def test_snp_shown_in_alignment(self, validator, ref_path, make_sample):
        """SNP should show ref base at position in diff line."""
        variants = [{"pos": "263", "ref": "A", "seq": "G"}]
        seq_regions = [[68, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        alignment = validator._build_alignment(ref_dict, merged_con, hv_seqs["HV2"]["ranges"])

        # Position 263 should show ref=A, alt=G, diff=A (SNP shows ref base)
        idx_263 = alignment.pos_markers.index(263)
        assert alignment.ref_aligned[idx_263] == "A"
        assert alignment.alt_aligned[idx_263] == "G"
        assert alignment.diff_line[idx_263] == "A"  # SNP shows ref base

    def test_insertion_shows_gap_in_ref(self, validator, ref_path, make_sample):
        """Insertion should show '-' in ref and inserted bases in alt."""
        variants = [
            {"pos": "73", "ref": "A", "seq": "G"},
            {"pos": "309.1", "ref": "-", "seq": "C"},
            {"pos": "315.1", "ref": "-", "seq": "C"},
        ]
        seq_regions = [[68, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        alignment = validator._build_alignment(ref_dict, merged_con, hv_seqs["HV2"]["ranges"])

        # Count insertions (ref='-')
        ins_count = alignment.ref_aligned.count("-")
        assert ins_count == 2  # 309.1 and 315.1

        # Diff line should show '^' for insertions
        assert "^" in alignment.diff_line
        caret_count = alignment.diff_line.count("^")
        assert caret_count == 2

        # Insertion positions should have pos_markers=0
        zero_count = sum(1 for p in alignment.pos_markers if p == 0)
        assert zero_count == 2

    def test_deletion_shows_gap_in_alt(self, validator, ref_path, make_sample):
        """Deletion should show ref base in ref and '-' in alt."""
        variants = [{"pos": "249", "ref": "A", "seq": "-"}]
        seq_regions = [[68, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        alignment = validator._build_alignment(ref_dict, merged_con, hv_seqs["HV2"]["ranges"])

        # Position 249 should show ref=A, alt=-, diff=D
        idx_249 = alignment.pos_markers.index(249)
        assert alignment.ref_aligned[idx_249] == "A"
        assert alignment.alt_aligned[idx_249] == "-"
        assert alignment.diff_line[idx_249] == "D"

    def test_ten_base_insertion(self, validator, ref_path, make_sample):
        """10-base insertion at 217.1-217.10 should produce 10 gaps in ref."""
        variants = [{"pos": "73", "ref": "A", "seq": "G"}]
        variants += [{"pos": f"217.{i}", "ref": "-", "seq": "TTAATTGATT"[i - 1]} for i in range(1, 11)]
        variants += [
            {"pos": "263", "ref": "A", "seq": "G"},
            {"pos": "309.1", "ref": "-", "seq": "C"},
            {"pos": "309.2", "ref": "-", "seq": "C"},
            {"pos": "315.1", "ref": "-", "seq": "C"},
        ]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        alignment = validator._build_alignment(ref_dict, merged_con, hv_seqs["HV2"]["ranges"])

        # 10 insertion gaps at 217 + 2 at 309 + 1 at 315 = 13
        ins_count = alignment.ref_aligned.count("-")
        assert ins_count == 13

        # Verify the 10-base insertion specifically
        idx_217 = alignment.pos_markers.index(217)
        # After position 217, the next 10 columns should be insertion gaps
        for i in range(10):
            assert alignment.ref_aligned[idx_217 + 1 + i] == "-"
            assert alignment.pos_markers[idx_217 + 1 + i] == 0
            assert alignment.diff_line[idx_217 + 1 + i] == "^"

        # The inserted bases should be TTAATTGATT
        inserted = alignment.alt_aligned[idx_217 + 1 : idx_217 + 11]
        assert inserted == "TTAATTGATT"

    def test_ref_alt_aligned_equal_length(self, validator, ref_path, make_sample):
        """ref_aligned and alt_aligned must always be the same length."""
        variants = [
            {"pos": "73", "ref": "A", "seq": "G"},
            {"pos": "249", "ref": "A", "seq": "-"},  # deletion
            {"pos": "309.1", "ref": "-", "seq": "C"},  # insertion
        ]
        seq_regions = [[68, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        for region_name in ["HV1", "HV2", "HV3"]:
            if region_name in hv_seqs:
                alignment = validator._build_alignment(
                    ref_dict,
                    merged_con,
                    hv_seqs[region_name]["ranges"],
                )
                assert len(alignment.ref_aligned) == len(alignment.alt_aligned)
                assert len(alignment.ref_aligned) == len(alignment.diff_line)
                assert len(alignment.ref_aligned) == len(alignment.pos_markers)

    def test_multi_range_extraction(self, validator, ref_path, make_sample):
        """Alignment should handle multiple extraction ranges for a region."""
        variants = [{"pos": "73", "ref": "A", "seq": "G"}]
        # Two separate ranges
        seq_regions = [[73, 200], [250, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        ref_dict = result.ref_dict
        con_dict = result.con_dict
        hv_seqs = result.hv_seqs

        merged_con = _merge_con_dict(con_dict)
        alignment = validator._build_alignment(ref_dict, merged_con, hv_seqs["HV2"]["ranges"])

        # Should cover positions 73-340 (continuous REGIONS range for HV2)
        assert 73 in alignment.pos_markers
        assert 340 in alignment.pos_markers
        # All positions in the REGIONS range are in pos_markers
        # (Gap positions in seq_regions are N-filled but still in the extraction range)
        for p in range(73, 341):
            assert p in alignment.pos_markers


# ---------------------------------------------------------------------------
# _parse_fasta tests
# ---------------------------------------------------------------------------


class TestParseFasta:
    """Tests for FASTA file parsing."""

    @pytest.fixture
    def validator(self, ref_path, tmp_path):
        return FastaValidator(
            fasta_dir=tmp_path / "fasta",
            json_dir=tmp_path / "json",
            ref_path=ref_path,
            batch_id="TEST",
        )

    def test_parse_standard_fasta(self, validator, tmp_path):
        """Parse a standard 3-region FASTA file."""
        fasta_content = """>HV1
ATCGATCG
>HV2|68-340
GCTAGCTA
>HV3
TTTTAAAA
"""
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(fasta_content)

        regions = validator._parse_fasta(fasta_path)
        assert regions["HV1"] == "ATCGATCG"
        assert regions["HV2"] == "GCTAGCTA"
        assert regions["HV3"] == "TTTTAAAA"

    def test_parse_multiline_sequence(self, validator, tmp_path):
        """Parse FASTA with sequence on multiple lines."""
        fasta_content = """>HV1
ATCG
ATCG
>HV2|68-340
GCTA
GCTA
"""
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(fasta_content)

        regions = validator._parse_fasta(fasta_path)
        assert regions["HV1"] == "ATCGATCG"
        assert regions["HV2"] == "GCTAGCTA"


# ---------------------------------------------------------------------------
# validate_sample tests
# ---------------------------------------------------------------------------


class TestValidateSample:
    """Tests for sample validation (PASS/FAIL logic)."""

    def test_pass_sample(self, ref_path, tmp_path, sample_data, make_sample):
        """A correctly generated FASTA should PASS validation."""

        # Create JSON
        json_dir = tmp_path / "json"
        json_dir.mkdir()
        sample_id = "TEST_SAMPLE"
        sample_dir = json_dir / sample_id
        sample_dir.mkdir()
        with (sample_dir / f"{sample_id}.json").open("w") as f:
            json.dump(sample_data, f)

        # Generate FASTA
        variants = []
        for region_vars in sample_data["variants"].values():
            variants.extend(region_vars)
        seq_regions = []
        for value in sample_data["intervals"].values():
            seq_regions.extend(value)

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        # Write FASTA file
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        fasta_path = fasta_dir / f"{sample_id}.fasta"
        with fasta_path.open("w") as f:
            for region in ["HV1", "HV2", "HV3"]:
                seq_str = "".join(hv_seqs[region]["seq"])
                start = hv_seqs[region]["start"]
                end = hv_seqs[region]["end"]
                if region == "HV2":
                    f.write(f">{region}|{start}-{end}\n")
                else:
                    f.write(f">{region}\n")
                # Write sequence in 70-char lines
                f.writelines(seq_str[i : i + 70] + "\n" for i in range(0, len(seq_str), 70))

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        result = validator.validate_sample(sample_id)
        assert result.passed is True

    def test_fail_corrupt_fasta(self, ref_path, tmp_path, sample_data):
        """A corrupted FASTA should FAIL validation."""

        # Create JSON
        json_dir = tmp_path / "json"
        json_dir.mkdir()
        sample_id = "CORRUPT"
        sample_dir = json_dir / sample_id
        sample_dir.mkdir()
        with (sample_dir / f"{sample_id}.json").open("w") as f:
            json.dump(sample_data, f)

        # Write a corrupt FASTA
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        fasta_path = fasta_dir / f"{sample_id}.fasta"
        with fasta_path.open("w") as f:
            f.write(">HV1\nAAAA\n>HV2|68-340\nCCCC\n>HV3\nGGGG\n")

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        result = validator.validate_sample(sample_id)
        assert result.passed is False
        # At least one region should fail
        assert any(not rr.passed for rr in result.regions.values())

    def test_missing_json_returns_error(self, ref_path, tmp_path):
        """Missing JSON file should return a failed SampleResult."""
        json_dir = tmp_path / "json"
        json_dir.mkdir()
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        result = validator.validate_sample("NONEXISTENT")
        assert result.passed is False
        assert "JSON file not found" in result.error

    def test_missing_fasta_returns_error(self, ref_path, tmp_path, sample_data):
        """Missing FASTA file should return a failed SampleResult."""
        json_dir = tmp_path / "json"
        json_dir.mkdir()
        sample_id = "NOFASTA"
        sample_dir = json_dir / sample_id
        sample_dir.mkdir()
        with (sample_dir / f"{sample_id}.json").open("w") as f:
            json.dump(sample_data, f)

        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        result = validator.validate_sample(sample_id)
        assert result.passed is False
        assert "FASTA file not found" in result.error


# ---------------------------------------------------------------------------
# Report writing tests
# ---------------------------------------------------------------------------


class TestWriteReport:
    """Tests for report file writing."""

    def test_report_file_created(self, ref_path, tmp_path, sample_data, make_sample):
        """Report file should be created at the expected path."""

        json_dir = tmp_path / "json"
        json_dir.mkdir()
        sample_id = "REPORT_TEST"
        sample_dir = json_dir / sample_id
        sample_dir.mkdir()
        with (sample_dir / f"{sample_id}.json").open("w") as f:
            json.dump(sample_data, f)

        # Generate and write FASTA
        variants = []
        for region_vars in sample_data["variants"].values():
            variants.extend(region_vars)
        seq_regions = []
        for value in sample_data["intervals"].values():
            seq_regions.extend(value)

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        with (fasta_dir / f"{sample_id}.fasta").open("w") as f:
            for region in ["HV1", "HV2", "HV3"]:
                seq_str = "".join(hv_seqs[region]["seq"])
                start = hv_seqs[region]["start"]
                end = hv_seqs[region]["end"]
                if region == "HV2":
                    f.write(f">{region}|{start}-{end}\n")
                else:
                    f.write(f">{region}\n")
                f.writelines(seq_str[i : i + 70] + "\n" for i in range(0, len(seq_str), 70))

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        validator.validate_all()

        output_dir = tmp_path / "results"
        report_path = validator.write_report(output_dir)

        assert report_path.exists()
        assert report_path.name == "validation_fasta_TEST.txt"
        assert report_path.parent.name == "validation"

        # Check report content
        content = report_path.read_text()
        assert "FASTA Validation Report" in content
        assert "TEST" in content
        assert "Summary" in content


# ---------------------------------------------------------------------------
# Integration test with real data
# ---------------------------------------------------------------------------


class TestMultiBaseInsertionAndAllRegions:
    """Integration tests: 10-base insertion alignment and all-region PASS validation.

    Uses synthetic data (no external files), so these tests never skip.
    Replaces the former TestIntegrationRealData which skipped when
    MS_090426_004 data was unavailable.
    """

    @pytest.fixture
    def setup_multi_ins_sample(self, ref_path, tmp_path, sample_data_multi_ins, make_sample) -> tuple[Path, Path, str]:
        """Create JSON + FASTA for the sample_data_multi_ins fixture."""

        # Create JSON
        json_dir = tmp_path / "json"
        json_dir.mkdir()
        sample_id = "MULTI_INS_SAMPLE"
        sample_dir = json_dir / sample_id
        sample_dir.mkdir()
        with (sample_dir / f"{sample_id}.json").open("w") as f:
            json.dump(sample_data_multi_ins, f)

        # Generate FASTA
        variants = []
        for region_vars in sample_data_multi_ins["variants"].values():
            variants.extend(region_vars)
        seq_regions = []
        for value in sample_data_multi_ins["intervals"].values():
            seq_regions.extend(value)

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        # Write FASTA file
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        fasta_path = fasta_dir / f"{sample_id}.fasta"
        with fasta_path.open("w") as f:
            for region in ["HV1", "HV2", "HV3"]:
                seq_str = "".join(hv_seqs[region]["seq"])
                start = hv_seqs[region]["start"]
                end = hv_seqs[region]["end"]
                if region == "HV2":
                    f.write(f">{region}|{start}-{end}\n")
                else:
                    f.write(f">{region}\n")
                f.writelines(seq_str[i : i + 70] + "\n" for i in range(0, len(seq_str), 70))

        return fasta_dir, json_dir, sample_id

    def test_hv2_ten_base_insertion_alignment(self, ref_path, setup_multi_ins_sample):
        """HV2 has a 10-base insertion at 217: verify alignment shows 10 consecutive gaps + correct sequence."""
        fasta_dir, json_dir, sample_id = setup_multi_ins_sample

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        result = validator.validate_sample(sample_id)
        assert result.passed is True

        hv2 = result.regions["HV2"]

        # 10 insertion gaps at 217 + 2 at 309 + 1 at 315 = 13 total
        ins_count = hv2.ref_aligned.count("-")
        assert ins_count == 13

        # Verify 10 consecutive gaps after position 217
        idx_217 = hv2.pos_markers.index(217)
        gaps_after_217 = 0
        for i in range(idx_217 + 1, len(hv2.pos_markers)):
            if hv2.pos_markers[i] == 0:
                gaps_after_217 += 1
            else:
                break
        assert gaps_after_217 == 10

        # Verify the inserted sequence
        inserted = hv2.alt_aligned[idx_217 + 1 : idx_217 + 11]
        assert inserted == "TTAATTGATT"

    def test_all_regions_pass(self, ref_path, setup_multi_ins_sample):
        """All 3 regions should PASS validation."""
        fasta_dir, json_dir, sample_id = setup_multi_ins_sample

        validator = FastaValidator(
            fasta_dir=fasta_dir,
            json_dir=json_dir,
            ref_path=ref_path,
            batch_id="TEST",
        )
        result = validator.validate_sample(sample_id)

        for region_name, rr in result.regions.items():
            assert rr.passed is True, f"{region_name} failed validation"
