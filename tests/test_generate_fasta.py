"""Tests for the generate_fasta module.

Integration tests exercise the canonical generate_sequence() path (src.core.sample).
FastaGenerator tests cover write_fasta_file and generate_reviewed_fasta.
"""

import json
import logging
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from src.core.sample import (
    _apply_variants_to_dict,
    _extract_region,
    _load_reference_sequence,
    generate_sequence,
)
from src.core.variants import (
    normalize_position,
    pos_base,
)
from src.generate_fasta import FastaGenerator

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ref_path():
    """Return the path to the reference FASTA file."""
    return str(Path(__file__).parent.parent / "ref" / "rCRS.fasta")


@pytest.fixture
def temp_input_dir(tmp_path):
    """Create a temporary input directory for testing."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    return input_dir


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory for testing."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture
def sample_data():
    """Return sample variant data for testing."""
    return {
        "variants": {
            "HV1": [
                {"pos": 16024, "ref": "G", "seq": "A"},
                {"pos": 16223, "ref": "C", "seq": "T"},
            ],
            "HV2": [
                {"pos": 73, "ref": "G", "seq": "A"},
                {"pos": 263, "ref": "A", "seq": "G"},
            ],
            "HV3": [
                {"pos": 490, "ref": "T", "seq": "C"},
                {"pos": 519, "ref": "G", "seq": "A"},
            ],
        },
        "intervals": {
            "HV1": [[16024, 16365]],
            "HV2": [[73, 340]],
            "HV3": [[438, 576]],
        },
    }


@pytest.fixture
def sample_with_complex_variants():
    """Return complex variant data with multiple insertions and deletions."""
    return {
        "variants": {
            "HV1": [
                {"pos": 16189.1, "ref": "-", "seq": "C"},
                {"pos": 16111, "ref": "T", "seq": "-"},
            ],
            "HV2": [
                {"pos": 309.1, "ref": "-", "seq": "C"},
                {"pos": 309.2, "ref": "-", "seq": "C"},
                {"pos": 315.1, "ref": "-", "seq": "T"},
                {"pos": 333, "ref": "G", "seq": "-"},
            ],
            "HV3": [
                {"pos": 555, "ref": "A", "seq": "-"},
            ],
        },
        "intervals": {
            "HV1": [[16024, 16365]],
            "HV2": [[73, 340]],
            "HV3": [[438, 576]],
        },
    }


@pytest.fixture
def create_sample_json(temp_input_dir):
    """Factory fixture to create sample JSON files."""

    def _create_sample(sample_name: str, data: dict) -> tuple[Path, Path]:
        sample_dir = temp_input_dir / sample_name
        sample_dir.mkdir()
        json_path = sample_dir / f"{sample_name}.json"
        with json_path.open("w") as f:
            json.dump(data, f)
        return sample_dir, json_path

    return _create_sample


# ---------------------------------------------------------------------------
# Test: pos_base (redirected from FastaGenerator.normalize_pos)
# ---------------------------------------------------------------------------


class TestPosBase:
    """Tests for pos_base from src.core.variants (replaces FastaGenerator.normalize_pos)."""

    def test_integer_returns_as_is(self):
        """Test that integer positions return as integers."""
        assert pos_base(5) == 5
        assert pos_base(0) == 0
        assert pos_base(-10) == -10

    def test_float_returns_integer_part(self):
        """Test that float positions return integer part."""
        assert pos_base(5.0) == 5
        assert pos_base(309.1) == 309
        assert pos_base(309.9) == 309

    def test_string_position(self):
        """Test that string positions return integer part."""
        assert pos_base("309") == 309
        assert pos_base("309.1") == 309
        assert pos_base("217.10") == 217


# ---------------------------------------------------------------------------
# Integration tests: generate_sequence
# ---------------------------------------------------------------------------


class TestGenSeqFromVariantsRegionsIntegration:
    """Integration tests for generate_sequence from src.core.sample."""

    def test_with_snps(self, ref_path, make_sample):
        """Test SNP application through generate_sequence."""
        variants = [
            {"pos": 16024, "ref": "G", "seq": "A"},
            {"pos": 73, "ref": "G", "seq": "A"},
        ]
        seq_regions = [[16024, 16365], [73, 340], [438, 576]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV1"]["seq"]) == 342
        assert hv_seqs["HV1"]["seq"][0] == "A"  # 16024: G -> A

    def test_with_single_base_insertions(self, ref_path, make_sample):
        """Test insertion handling through generate_sequence."""
        variants = [
            {"pos": 309.1, "ref": "-", "seq": "C"},
            {"pos": 309.2, "ref": "-", "seq": "C"},
            {"pos": 315.1, "ref": "-", "seq": "T"},
        ]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        # HV2: 268 + 3 insertion bases = 271
        assert len(hv_seqs["HV2"]["seq"]) == 271

    def test_with_deletions(self, ref_path, make_sample):
        """Test deletion handling through generate_sequence."""
        variants = [
            {"pos": 263, "ref": "A", "seq": "-"},
        ]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        # HV2: 268 - 1 deletion = 267
        assert len(hv_seqs["HV2"]["seq"]) == 267

    def test_with_mixed_variants(self, ref_path, make_sample):
        """Test combined SNP + insertion + deletion."""
        variants = [
            {"pos": 16024, "ref": "G", "seq": "A"},
            {"pos": 309.1, "ref": "-", "seq": "C"},
            {"pos": 309.2, "ref": "-", "seq": "C"},
            {"pos": 263, "ref": "A", "seq": "-"},
        ]
        seq_regions = [[16024, 16365], [73, 340], [438, 576]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        # HV2: 268 + 2 insertions - 1 deletion = 269
        assert len(hv_seqs["HV2"]["seq"]) == 269
        hv2_str = "".join(hv_seqs["HV2"]["seq"])
        assert "-" not in hv2_str

    def test_multiple_ranges_per_region(self, ref_path, make_sample):
        """Test handling of multiple ranges per region."""
        variants = [
            {"pos": 16024, "ref": "G", "seq": "A"},
        ]
        seq_regions = [[16024, 16365], [73, 340], [438, 576]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV1"]["seq"]) == 342
        assert len(hv_seqs["HV2"]["seq"]) == 268
        assert len(hv_seqs["HV3"]["seq"]) == 139

    def test_del_and_insertion_at_same_position(self, ref_path, make_sample):
        """Test DEL and INSERT at same position."""
        variants = [
            {"pos": 309, "ref": "C", "seq": "-"},
            {"pos": 309.1, "ref": "-", "seq": "C"},
            {"pos": 309.2, "ref": "-", "seq": "D"},
            {"pos": 309.3, "ref": "-", "seq": "B"},
        ]
        seq_regions = [[16024, 16365], [73, 340], [438, 576]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 270
        hv2_str = "".join(hv_seqs["HV2"]["seq"])
        assert "-" not in hv2_str

    def test_polyc_region(self, ref_path, make_sample):
        """Test handling of polyC region (common in mtDNA)."""
        variants = [
            {"pos": 309.1, "ref": "-", "seq": "C"},
            {"pos": 309.2, "ref": "-", "seq": "C"},
            {"pos": 309.3, "ref": "-", "seq": "C"},
            {"pos": 309.4, "ref": "-", "seq": "C"},
        ]
        seq_regions = [[16024, 16365], [73, 340], [438, 576]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        hv2_str = "".join(hv_seqs["HV2"]["seq"])
        assert hv2_str.count("C") > 1

    def test_reference_base_mismatch_warning(self, ref_path, make_sample, caplog):
        """Test that reference base mismatch logs a warning but still applies variant."""
        caplog.set_level(logging.WARNING)

        variants = [
            {"pos": 16024, "ref": "Z", "seq": "A"},  # Wrong ref base
        ]
        seq_regions = [[16024, 16365]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert hv_seqs["HV1"]["seq"][0] == "A"

    def test_n_marking_outside_regions(self, ref_path, make_sample):
        """Test that positions outside specified regions are marked as N."""
        variants = [
            {"pos": 73, "ref": "G", "seq": "A"},
        ]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert hv_seqs["HV2"]["seq"][0] == "A"

    def test_no_regions(self, ref_path, make_sample):
        """Test generate_sequence with no seq_regions (uses default regions)."""
        variants = [
            {"pos": 16024, "ref": "G", "seq": "A"},
        ]

        result = generate_sequence(make_sample(variants, []), ref_path)
        hv_seqs = result.hv_seqs

        assert "HV1" in hv_seqs

    def test_fasta_regions_override(self, ref_path, make_sample):
        """Test that REGIONS produces HV2 starting at 73."""
        variants = []

        result = generate_sequence(make_sample(variants, [[73, 340], [438, 576], [16024, 16365]]), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 268

    def test_default_regions_hv2_starts_at_73(self, ref_path, make_sample):
        """Test that default REGIONS produces HV2 starting at 73."""
        variants = []

        result = generate_sequence(make_sample(variants, [[73, 340], [438, 576], [16024, 16365]]), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 268


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Edge case tests for generate_sequence."""

    def test_multiple_deletions_at_different_positions(self, ref_path, make_sample):
        """Test deletions at various positions (333, 555, 16111)."""
        variants = [
            {"pos": 333, "ref": "G", "seq": "-"},
            {"pos": 555, "ref": "A", "seq": "-"},
            {"pos": 16111, "ref": "T", "seq": "-"},
        ]
        seq_regions = [[73, 340], [438, 576], [16024, 16365]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 267
        assert len(hv_seqs["HV3"]["seq"]) == 138
        assert len(hv_seqs["HV1"]["seq"]) == 341

    def test_mixed_insertions_and_deletions(self, ref_path, make_sample):
        """Test combination of insertions and deletions."""
        variants = [
            {"pos": 309.1, "ref": "-", "seq": "C"},
            {"pos": 309.2, "ref": "-", "seq": "C"},
            {"pos": 315.1, "ref": "-", "seq": "T"},
            {"pos": 525.1, "ref": "-", "seq": "AC"},
            {"pos": 525.2, "ref": "-", "seq": "GT"},
            {"pos": 16189.1, "ref": "-", "seq": "C"},
            {"pos": 333, "ref": "G", "seq": "-"},
            {"pos": 555, "ref": "A", "seq": "-"},
            {"pos": 16111, "ref": "T", "seq": "-"},
        ]
        seq_regions = [[73, 340], [438, 576], [16024, 16365]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 270
        assert len(hv_seqs["HV3"]["seq"]) == 142
        assert len(hv_seqs["HV1"]["seq"]) == 342

        for region in hv_seqs:
            seq_str = "".join(hv_seqs[region]["seq"])
            assert "-" not in seq_str

    def test_insertion_sequence_integrity(self, ref_path, make_sample):
        """Test that inserted sequences are in correct order."""
        variants = [
            {"pos": 309.1, "ref": "-", "seq": "A"},
            {"pos": 309.2, "ref": "-", "seq": "B"},
            {"pos": 309.3, "ref": "-", "seq": "C"},
            {"pos": 315.1, "ref": "-", "seq": "XYZ"},
        ]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        hv2_seq = "".join(hv_seqs["HV2"]["seq"])
        assert len(hv2_seq) == 274

    def test_complex_variants_with_fixture(self, ref_path, make_sample, sample_with_complex_variants):
        """Test using the sample_with_complex_variants fixture."""
        data = sample_with_complex_variants
        all_variants = []
        for region in data["variants"]:
            all_variants.extend(data["variants"][region])

        seq_regions = [
            data["intervals"]["HV2"][0],
            data["intervals"]["HV3"][0],
            data["intervals"]["HV1"][0],
        ]

        result = generate_sequence(make_sample(all_variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 270
        assert len(hv_seqs["HV3"]["seq"]) == 138
        assert len(hv_seqs["HV1"]["seq"]) == 342

        for region in hv_seqs:
            seq_str = "".join(hv_seqs[region]["seq"])
            assert "-" not in seq_str

    def test_multi_base_insertion(self, ref_path, make_sample):
        """Test multi-base insertion (e.g., 'ABCD' at a single position)."""
        variants = [
            {"pos": 309.1, "ref": "-", "seq": "ABCD"},
        ]
        seq_regions = [[73, 340]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV2"]["seq"]) == 272

    def test_snps_with_real_reference(self, ref_path, make_sample):
        """Test with real reference sequence."""
        variants = [
            {"pos": 16024, "ref": "G", "seq": "A"},
            {"pos": 73, "ref": "G", "seq": "A"},
            {"pos": 263, "ref": "A", "seq": "G"},
            {"pos": 490, "ref": "T", "seq": "C"},
        ]
        seq_regions = [[16024, 16365], [73, 340], [438, 576]]

        result = generate_sequence(make_sample(variants, seq_regions), ref_path)
        hv_seqs = result.hv_seqs

        assert len(hv_seqs["HV1"]["seq"]) == 342
        assert len(hv_seqs["HV2"]["seq"]) == 268
        assert len(hv_seqs["HV3"]["seq"]) == 139


# ---------------------------------------------------------------------------
# FastaGenerator tests
# ---------------------------------------------------------------------------


class TestWriteFastaFile:
    """Tests for the write_fasta_file method."""

    def test_write_basic(self, temp_output_dir, ref_path):
        """Test basic FASTA file writing."""
        generator = FastaGenerator(temp_output_dir, temp_output_dir, ref_path)

        region_sequences = {
            "HV1": "ACGTACGTACGT",
            "HV2": "GCTGCTGCTGC",
            "HV3": "TATATATATATA",
        }

        output_path = temp_output_dir / "test.fasta"
        generator.write_fasta_file(output_path, region_sequences)

        assert output_path.exists()

        with output_path.open() as f:
            content = f.read()

        assert ">HV1" in content
        assert ">HV2" in content
        assert ">HV3" in content

    def test_hv2_plain_header(self, temp_output_dir, ref_path):
        """Test that HV2 header is plain (no coordinates)."""
        generator = FastaGenerator(temp_output_dir, temp_output_dir, ref_path)

        region_sequences = {
            "HV1": "ACGT",
            "HV2": "GCTA",
            "HV3": "TATA",
        }

        output_path = temp_output_dir / "test.fasta"
        generator.write_fasta_file(output_path, region_sequences)

        with output_path.open() as f:
            content = f.read()

        assert ">HV2\n" in content

    def test_all_headers_are_plain(self, temp_output_dir, ref_path):
        """Test that all region headers are plain (no coordinates)."""
        generator = FastaGenerator(temp_output_dir, temp_output_dir, ref_path)

        region_sequences = {
            "HV1": "ACGT",
            "HV2": "GCTA",
            "HV3": "TATA",
        }

        output_path = temp_output_dir / "test.fasta"
        generator.write_fasta_file(output_path, region_sequences)

        with output_path.open() as f:
            content = f.read()

        # All region headers should have no pipe character
        for line in content.strip().split("\n"):
            if line.startswith(">"):
                assert "|" not in line, f"Header should not have coordinates: {line}"


class TestGenerateReviewedFasta:
    """Tests for the generate_reviewed_fasta method."""

    def test_generate_from_json(self, temp_input_dir, temp_output_dir, ref_path, sample_data):
        """Test generating FASTA from JSON file and verify output content."""
        sample_name = "sample001"
        sample_dir = temp_input_dir / sample_name
        sample_dir.mkdir()

        json_path = sample_dir / f"{sample_name}.json"
        with json_path.open("w") as f:
            json.dump(sample_data, f)

        generator = FastaGenerator(temp_input_dir, temp_output_dir, ref_path)
        success, failure, _ = generator.generate_reviewed_fasta()

        assert success == 1
        assert failure == 0

        # Verify FASTA file content
        output_file = temp_output_dir / "sample001.fasta"
        assert output_file.exists()
        content = output_file.read_text()
        # All headers should be plain (no coordinates)
        assert ">HV1\n" in content
        assert ">HV2\n" in content
        assert ">HV3\n" in content

    def test_generate_reviewed_fasta_plain_headers(self, temp_input_dir, temp_output_dir, ref_path):
        """Test that all region headers are plain regardless of interval values."""
        sample_data = {
            "variants": {"HV1": [], "HV2": [], "HV3": []},
            "intervals": {
                "HV1": [[16024, 16365]],
                "HV2": [[73, 340]],
                "HV3": [[438, 576]],
            },
        }

        sample_name = "sample002"
        sample_dir = temp_input_dir / sample_name
        sample_dir.mkdir()

        json_path = sample_dir / f"{sample_name}.json"
        with json_path.open("w") as f:
            json.dump(sample_data, f)

        generator = FastaGenerator(temp_input_dir, temp_output_dir, ref_path)
        success, _, _ = generator.generate_reviewed_fasta()

        assert success == 1

        # Verify all headers are plain (no coordinates)
        output_file = temp_output_dir / "sample002.fasta"
        content = output_file.read_text()
        assert ">HV2\n" in content

    def test_generate_reviewed_fasta_missing_json(self, temp_input_dir, temp_output_dir, ref_path):
        """Test handling of missing JSON files."""
        sample_name = "sample003"
        sample_dir = temp_input_dir / sample_name
        sample_dir.mkdir()

        generator = FastaGenerator(temp_input_dir, temp_output_dir, ref_path)
        success, failure, _ = generator.generate_reviewed_fasta()

        assert success == 0
        assert failure == 1


# ---------------------------------------------------------------------------
# Regression test: SNP/deletion must overwrite ref int key (not add str key)
# ---------------------------------------------------------------------------


class TestApplyVariantsDictKeyConsistency:
    """Ensure SNP/deletion variants overwrite the ref_dict int key,
    not add a duplicate string key.

    This was Bug 6: _apply_variants_to_dict stored SNP/deletion variants
    under string keys (e.g., "16129") while ref_dict uses int keys (16129).
    Both persisted in con_dict, causing _extract_region to emit 2 bases per SNP.
    Fix: convert pos to int for SNP/deletion so it overwrites the ref key.
    """

    def test_snp_overwrites_int_key(self, ref_path):
        """SNP at position should overwrite the ref int key, not add a str key."""

        ref_dict = _load_reference_sequence(ref_path)
        variants = [{"pos": "16129", "ref": "G", "seq": "A"}]
        con_dict = _apply_variants_to_dict(ref_dict, variants)

        # Int key should be overwritten to the SNP value
        assert con_dict[16129] == "A"
        # No duplicate string key should exist
        assert "16129" not in con_dict

    def test_deletion_overwrites_int_key(self, ref_path):
        """Deletion should overwrite the ref int key with '-', not add a str key."""

        ref_dict = _load_reference_sequence(ref_path)
        variants = [{"pos": "249", "ref": "A", "seq": "-"}]
        con_dict = _apply_variants_to_dict(ref_dict, variants)

        assert con_dict[249] == "-"
        assert "249" not in con_dict

    def test_insertion_keeps_string_key(self, ref_path):
        """Insertion should keep the string key (e.g., '309.1'), not convert to int."""

        ref_dict = _load_reference_sequence(ref_path)
        variants = [{"pos": "309.1", "ref": "-", "seq": "C"}]
        con_dict = _apply_variants_to_dict(ref_dict, variants)

        assert con_dict["309.1"] == "C"
        # The insertion key is a string, not an int
        # (Note: int 3091 is a real rCRS position, not related to insertion 309.1)

    def test_no_duplicate_keys_after_snp(self, ref_path):

        ref_dict = _load_reference_sequence(ref_path)
        variants = [{"pos": "489", "ref": "T", "seq": "C"}]
        con_dict = _apply_variants_to_dict(ref_dict, variants)

        ref_extracted = _extract_region(ref_dict, 438, 576)
        con_extracted = _extract_region(con_dict, 438, 576)

        # Lengths must match — SNP replaces one base, doesn't add one
        assert len(con_extracted) == len(ref_extracted)


# ---------------------------------------------------------------------------
# Regression tests: normalize_position and type consistency
# ---------------------------------------------------------------------------


class TestNormalizePositionTypes:
    """Ensure normalize_position returns int for base positions,
    str for insertion positions, matching the PerSampleVariant contract.
    """

    def test_int_base_returns_int(self):
        """Base position as int returns int."""
        result = normalize_position(309)
        assert result == 309
        assert isinstance(result, int)

    def test_str_base_returns_int(self):
        """String base position returns int."""
        result = normalize_position("309")
        assert result == 309
        assert isinstance(result, int)

    def test_float_base_returns_int(self):
        """Float base position (e.g., 309.0) returns int."""
        result = normalize_position(309.0)
        assert result == 309
        assert isinstance(result, int)

    def test_str_insertion_returns_str(self):
        """Insertion position as string stays string."""
        result = normalize_position("309.1")
        assert result == "309.1"
        assert isinstance(result, str)

    def test_str_insertion_preserves_precision(self):
        """Multi-digit insertion index preserved: '217.10' stays '217.10'."""
        result = normalize_position("217.10")
        assert result == "217.10"
        assert isinstance(result, str)

    def test_float_insertion_returns_str(self):
        """Float insertion (309.1) returns string preserving decimal."""
        result = normalize_position(309.1)
        assert result == "309.1"
        assert isinstance(result, str)

    def test_comma_stripped_base_returns_int(self):
        """Comma-formatted base position returns int."""
        result = normalize_position("1,234")
        assert result == 1234
        assert isinstance(result, int)

    def test_comma_stripped_insertion_returns_str(self):
        """Comma-formatted insertion position returns string."""
        result = normalize_position("1,234.10")
        assert result == "1234.10"
        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# Test: Skip-existing and --force behavior (single-folder merge)
# ---------------------------------------------------------------------------


class TestSkipExistingForce:
    """Tests for FastaGenerator skip-existing (default) and --force behavior.

    Default: skip samples that already have .fasta files in output dir.
    --force: wipe output dir and regenerate all (same as old unconditional behavior).
    """

    def test_default_skips_existing(self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json):
        """Default mode: skips samples that already have .fasta in output dir."""
        # Create two samples
        create_sample_json("S001", sample_data)
        create_sample_json("S002", sample_data)

        # First run: generate all
        gen1 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        success1, _failure1, skip1 = gen1.generate_reviewed_fasta()
        assert success1 == 2
        assert skip1 == 0

        # Write a marker to S001.fasta to verify it's NOT overwritten on second run
        s001_path = temp_output_dir / "S001.fasta"
        original_content = s001_path.read_text()

        # Second run: same input, should skip existing
        gen2 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        success2, _failure2, skip2 = gen2.generate_reviewed_fasta()
        assert skip2 == 2  # both S001 and S002 already exist
        assert success2 == 0  # nothing new generated

        # Verify original files are untouched
        assert s001_path.read_text() == original_content

    def test_default_generates_all_when_empty(
        self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json
    ):
        """Default mode with empty output dir: generates all samples."""
        create_sample_json("S001", sample_data)
        create_sample_json("S002", sample_data)

        # Remove the output dir so FastaGenerator creates it fresh
        if temp_output_dir.exists():
            shutil.rmtree(str(temp_output_dir))

        gen = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        success, _failure, skipped = gen.generate_reviewed_fasta()
        assert success == 2
        assert skipped == 0
        assert (temp_output_dir / "S001.fasta").exists()
        assert (temp_output_dir / "S002.fasta").exists()

    def test_force_wipes_and_regenerates(
        self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json
    ):
        """--force: wipes output dir and regenerates all samples."""
        create_sample_json("S001", sample_data)
        create_sample_json("S002", sample_data)

        # First run: generate all
        gen1 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        success1, _, _ = gen1.generate_reviewed_fasta()
        assert success1 == 2

        # Add a stale file from a previous pass (simulating a sample no longer in input)
        stale_file = temp_output_dir / "S999.fasta"
        stale_file.write_text(">stale\nACGT\n")

        # Force run: should wipe and regenerate only current samples
        gen2 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=True)
        success2, _failure2, skip2 = gen2.generate_reviewed_fasta()
        assert success2 == 2
        assert skip2 == 0
        assert not stale_file.exists()  # stale file wiped

    def test_force_empty_input(self, temp_input_dir, temp_output_dir, ref_path):
        """--force with empty input dir: creates output dir, generates nothing."""
        gen = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=True)
        success, failure, skipped = gen.generate_reviewed_fasta()
        assert success == 0
        assert failure == 0
        assert skipped == 0

    def test_return_value_includes_skip_count(
        self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json
    ):
        """generate_reviewed_fasta returns (success, failure, skip) tuple."""
        create_sample_json("S001", sample_data)
        create_sample_json("S002", sample_data)
        create_sample_json("S003", sample_data)

        # First run: generate all
        gen1 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        success1, failure1, skip1 = gen1.generate_reviewed_fasta()
        assert (success1, failure1, skip1) == (3, 0, 0)

        # Second run: skip existing
        gen2 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        success2, failure2, skip2 = gen2.generate_reviewed_fasta()
        assert (success2, failure2, skip2) == (0, 0, 3)

    def test_skip_logging(self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json):
        """Skipped samples are logged — verified by behavior (loguru caplog is unreliable)."""

        create_sample_json("S001", sample_data)

        # First run
        gen1 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        gen1.generate_reviewed_fasta()

        # Second run: mock logger to verify skip message
        with patch("src.generate_fasta.logger") as mock_logger:
            gen2 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
            gen2.generate_reviewed_fasta()

            skip_calls = [call for call in mock_logger.info.call_args_list if "Skipping existing" in str(call)]
            assert len(skip_calls) >= 1, f"Expected skip log, got: {mock_logger.info.call_args_list}"

    def test_force_creates_fresh_dir(self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json):
        """--force creates fresh output dir even when it doesn't exist."""
        if temp_output_dir.exists():
            shutil.rmtree(str(temp_output_dir))

        create_sample_json("S001", sample_data)

        gen = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=True)
        success, _, _ = gen.generate_reviewed_fasta()
        assert success == 1
        assert temp_output_dir.exists()

    def test_default_creates_dir_if_missing(self, temp_input_dir, ref_path, sample_data, create_sample_json, tmp_path):
        """Default mode creates output dir if it doesn't exist."""
        new_output = tmp_path / "fasta" / "batch_new"
        # Don't create it — FastaGenerator should create it

        create_sample_json("S001", sample_data)

        gen = FastaGenerator(str(temp_input_dir), str(new_output), ref_path, force=False)
        success, _, _ = gen.generate_reviewed_fasta()
        assert success == 1
        assert new_output.exists()
        assert (new_output / "S001.fasta").exists()

    def test_skip_does_not_overwrite_content(
        self, temp_input_dir, temp_output_dir, ref_path, sample_data, create_sample_json
    ):
        """Default mode: existing .fasta content is preserved, NOT overwritten."""
        create_sample_json("S001", sample_data)

        # First run
        gen1 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        gen1.generate_reviewed_fasta()
        original_content = (temp_output_dir / "S001.fasta").read_text()

        # Second run: should skip S001 (not overwrite)
        gen2 = FastaGenerator(str(temp_input_dir), str(temp_output_dir), ref_path, force=False)
        gen2.generate_reviewed_fasta()
        assert (temp_output_dir / "S001.fasta").read_text() == original_content
