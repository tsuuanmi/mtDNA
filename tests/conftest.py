import json
from pathlib import Path

import pytest

from src.core.models import Sample, Tool, Variant
from src.core.variants import normalize_position


@pytest.fixture
def ref_sequence():
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
            "HV2": [[68, 340]],
            "HV3": [[438, 576]],
        },
    }


@pytest.fixture
def sample_with_insertions():
    """Return sample variant data including insertions."""
    return {
        "variants": {
            "HV1": [{"pos": 16024, "ref": "G", "seq": "A"}],
            "HV2": [
                {"pos": 73, "ref": "G", "seq": "A"},
                {"pos": 309.1, "ref": "-", "seq": "CCCC"},  # Insertion at polyC region
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
def sample_with_deletions():
    """Return sample variant data including deletions."""
    return {
        "variants": {
            "HV1": [{"pos": 16024, "ref": "G", "seq": "A"}],
            "HV2": [
                {"pos": 73, "ref": "G", "seq": "-"},  # Deletion
                {"pos": 263, "ref": "A", "seq": "G"},
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
            "HV2": [[68, 340]],
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


@pytest.fixture
def make_sample():
    """Factory: build a Sample from raw variant dicts + seq_regions for generate_sequence.

    Mirrors the inputs previously passed to gen_seq_from_variants_regions: a flat list of
    {"pos","ref","seq"} dicts and a flat list of [start, end] seq_regions. The seq_regions
    are stored under an opaque interval key; generate_sequence only flattens interval
    values for N-marking, so the key is irrelevant to the output. Region extraction always
    uses the configured REGIONS, matching the former default behaviour.
    """

    def _make(variant_dicts: list[dict], seq_regions: list[list[int]], sample_id: str = "S1") -> Sample:
        variants = [
            Variant(pos=normalize_position(v["pos"]), ref=v["ref"], seq=v["seq"], peaks=None) for v in variant_dicts
        ]
        intervals = {"_": [list(r) for r in seq_regions]} if seq_regions else None
        return Sample(
            sample_id=sample_id,
            source_tool=Tool.UNIFIED,
            variants=variants,
            intervals=intervals,
            hv1=None,
            hv2=None,
            hv3=None,
            no_snps=None,
            no_ins=None,
            no_dels=None,
            no_identicals=None,
            no_unread=None,
            batch_id=None,
            information=None,
        )

    return _make
