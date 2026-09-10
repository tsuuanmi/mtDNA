"""Pytest regression tests for the 72-sample mtDNA pipeline.

Runs the pipeline against archived groundtruth data and compares results
per sample for both JSON variants/intervals and FASTA region sequences.

Set ``MTDNA_TESTS_DIR`` to override the groundtruth directory
(default: ``/mnt/nas/bca/mtDNA/science/tests``).

When the groundtruth directory is not available, the entire module is
skipped with a clear message rather than failing at collection time.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

import pytest
from Bio import SeqIO

from src.config import get_settings
from src.core.sample import _normalize_keyed_sample
from src.core.variants import get_hv_region_for_position, normalize_position, pos_sort_key

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------

SETTINGS = get_settings()
REGIONS = tuple(SETTINGS.regions.REGIONS.keys())
FASTA_EXTENSION = SETTINGS.parameters.FASTA_EXTENSION

TESTS_DIR = Path(os.environ.get("MTDNA_TESTS_DIR", "/mnt/nas/bca/mtDNA/science/tests"))
SAMPLE_LIST_FILE = TESTS_DIR / "data" / "raw" / "20240101_mtDNA_00.txt"
GROUNDTRUTH_JSON_DIR = TESTS_DIR / "results" / "manual_pipeline" / "regenerate" / "20240101_mtDNA_00"
GROUNDTRUTH_FASTA_DIR = TESTS_DIR / "results" / "fasta" / "20240101_mtDNA_00"

BATCH_ID = "20240101_mtDNA_00"

# Skip entire module when groundtruth data is unavailable
if not SAMPLE_LIST_FILE.exists():
    pytest.skip(
        f"Groundtruth sample list not found at {SAMPLE_LIST_FILE}. Set MTDNA_TESTS_DIR to the test data root.",
        allow_module_level=True,
    )

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_sample_ids() -> list[str]:
    """Read sample IDs from the groundtruth test list."""
    return [line.strip() for line in SAMPLE_LIST_FILE.read_text().splitlines() if line.strip()]


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as f:
        return json.load(f)


def _normalize_variant(variant: dict[str, Any]) -> dict[str, Any]:
    """Normalize a variant dict to only pos/ref/seq with normalized position."""
    pos = variant.get("pos")
    return {
        "pos": normalize_position(pos) if pos is not None else None,
        "ref": variant.get("ref"),
        "seq": variant.get("seq"),
    }


def _variant_sort_key(variant: dict[str, Any]) -> tuple[int, int]:
    pos = variant["pos"]
    return pos_sort_key(pos) if pos is not None else (10**9, 10**9)


def _variants_by_region(data: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Bucket variants into per-region lists, supporting both JSON layouts.

    Reuses the canonical ``_normalize_keyed_sample`` flattener (handles both the
    by-type ``{snps, insertions, deletions}`` and by-region ``{HV1, HV2, HV3}``
    layouts), then assigns each flat variant to its HV region by position via
    ``get_hv_region_for_position`` so the two layouts become comparable.
    """
    flat_variants = _normalize_keyed_sample(data).get("variants", [])
    bucket: dict[str, list[dict[str, Any]]] = {region: [] for region in REGIONS}
    for variant in flat_variants:
        region = get_hv_region_for_position(variant.get("pos"))
        if region is not None:
            bucket[region].append(variant)
    return bucket


def _normalize_json(data: dict[str, Any]) -> dict[str, Any]:
    """Normalize JSON to comparable form: sorted variants (pos/ref/seq) + intervals per region."""
    variants = _variants_by_region(data)
    return {
        "variants": {
            region: sorted(
                [_normalize_variant(v) for v in variants.get(region, [])],
                key=_variant_sort_key,
            )
            for region in REGIONS
        },
        "intervals": data.get("intervals", {}),
    }


def _parse_fasta_header(header: str) -> str:
    """Extract region name from FASTA header, ignoring |range suffix."""
    return header.split("|", maxsplit=1)[0]


def _read_fasta_sequences(path: Path) -> dict[str, str]:
    """Read a FASTA file and return {region: sequence} mapping."""
    records: dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        region = _parse_fasta_header(record.description)
        records[region] = str(record.seq).upper()
    return records


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def sample_ids() -> list[str]:
    """Return the list of sample IDs from the groundtruth test set."""
    return _read_sample_ids()


@pytest.fixture(scope="session")
def pipeline_output(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Run the pipeline (steps 3→4→5) once per session and return the output directory."""
    output_dir = tmp_path_factory.mktemp("mtdna_pipeline_output")
    project_root = Path(__file__).resolve().parent.parent

    result = subprocess.run(
        [
            str(project_root / "scripts" / "pipeline.sh"),
            "-s",
            "3,4,5",
            "-d",
            str(TESTS_DIR / "data"),
            "-o",
            str(output_dir),
            "-l",
            str(SAMPLE_LIST_FILE),
            BATCH_ID,
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        pytest.fail(f"Pipeline failed (exit {result.returncode}):\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    return output_dir

    # tmp_path_factory cleans up automatically after the session


# ---------------------------------------------------------------------------
# Parametrized test IDs
# ---------------------------------------------------------------------------


def _parametrize_sample_ids(metafunc: pytest.Metafunc) -> None:
    """Dynamically parametrize sample_id from the groundtruth file."""
    if "sample_id" in metafunc.fixturenames:
        sample_ids = _read_sample_ids()
        metafunc.parametrize("sample_id", sample_ids, ids=sample_ids)


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    _parametrize_sample_ids(metafunc)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_sample_count(sample_ids: list[str]) -> None:
    """Verify the groundtruth sample list is non-empty."""
    assert len(sample_ids) > 0, "Sample list is empty"


def test_json_compatibility(sample_id: str, pipeline_output: Path) -> None:
    """Compare JSON variants and intervals for each sample against groundtruth."""
    expected_path = GROUNDTRUTH_JSON_DIR / sample_id / f"{sample_id}.json"
    actual_path = pipeline_output / "manual_pipeline" / "regenerate" / BATCH_ID / sample_id / f"{sample_id}.json"

    assert expected_path.exists(), f"missing expected JSON for {sample_id}"
    assert actual_path.exists(), f"missing generated JSON for {sample_id}"

    expected = _normalize_json(_read_json(expected_path))
    actual = _normalize_json(_read_json(actual_path))

    assert expected == actual, (
        f"JSON mismatch for {sample_id}\n"
        f"Expected variants:\n{json.dumps(expected['variants'], indent=2)}\n"
        f"Actual variants:\n{json.dumps(actual['variants'], indent=2)}"
    )


def test_fasta_compatibility(sample_id: str, pipeline_output: Path) -> None:
    """Compare FASTA sequences per region for each sample against groundtruth."""
    expected_path = GROUNDTRUTH_FASTA_DIR / f"{sample_id}{FASTA_EXTENSION}"
    actual_path = pipeline_output / "fasta" / BATCH_ID / f"{sample_id}{FASTA_EXTENSION}"

    assert expected_path.exists(), f"missing expected FASTA for {sample_id}"
    assert actual_path.exists(), f"missing generated FASTA for {sample_id}"

    expected_sequences = _read_fasta_sequences(expected_path)
    actual_sequences = _read_fasta_sequences(actual_path)

    for region in REGIONS:
        assert region in expected_sequences, f"missing {region} in expected FASTA"
        assert region in actual_sequences, f"missing {region} in generated FASTA"

        assert expected_sequences[region] == actual_sequences[region], (
            f"FASTA sequence mismatch for {sample_id} region {region}"
        )
