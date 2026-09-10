"""Data-driven tests for ``src.fasta_to_sample`` (single ``FastaToSample`` class).

Tests build synthetic query sequences by applying known edits to rCRS region
slices, write them to a tmp FASTA, call ``call_variants``, and assert the
canonical variants/intervals.  No batch-specific result directories
(``FASTA_DIR`` / ``CANON_DIR``) are used -- only the stable rCRS reference, so
the tests do not depend on a particular batch path.  See ``docs/tools/fasta.md``
for the algorithm and normalization rules under test.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import pytest
from Bio import SeqIO

from src.core.batch import Batch
from src.core.models import Tool
from src.core.sample import generate_sequence
from src.fasta_to_sample import FastaToSample

REPO = Path(__file__).parent.parent
REF_PATH = str(REPO / "ref" / "rCRS.fasta")

# A case builder: (ref, regions) -> (records, expected_variants, expected_intervals).
CaseBuilder = Callable[
    [str, dict[str, list[int]]],
    "tuple[dict[str, str], set[tuple[str, str, str]], dict[str, list[list[int]]]]",
]

# Return type of a case builder.
CaseResult = tuple[dict[str, str], set[tuple[str, str, str]], dict[str, list[list[int]]]]


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


@pytest.fixture
def ref() -> str:
    """The rCRS reference sequence (stable reference, not a batch path)."""
    return str(next(SeqIO.parse(REF_PATH, "fasta")).seq).upper()


@pytest.fixture
def caller() -> FastaToSample:
    """A ``FastaToSample`` caller with the reference loaded once."""
    return FastaToSample(REF_PATH)


@pytest.fixture
def regions(caller: FastaToSample) -> dict[str, list[int]]:
    """Region name -> (start, end) from config (HV1/HV2/HV3)."""
    return caller._regions


def slice_region(ref: str, regions: dict[str, list[int]], region: str) -> str:
    """rCRS slice for one region (1-based inclusive)."""
    start, end = regions[region]
    return ref[start - 1 : end]


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    """Write a synthetic FASTA (id -> sequence)."""
    path.write_text("".join(f">{rid}\n{seq}\n" for rid, seq in records.items()), encoding="utf-8")
    return path


def variants_set(sample) -> set[tuple[str, str, str]]:
    return {(str(v.pos), v.ref, v.seq) for v in sample.variants}


def sub(seq: str, pos: int, region_start: int, base: str) -> str:
    """Substitute the base at 1-based ``pos`` within a region slice."""
    idx = pos - region_start
    return seq[:idx] + base + seq[idx + 1 :]


def ins_after(seq: str, pos: int, region_start: int, bases: str) -> str:
    """Insert ``bases`` after 1-based ``pos`` within a region slice."""
    idx = pos - region_start
    return seq[: idx + 1] + bases + seq[idx + 1 :]


def variants_set_from_json(entry: dict) -> set[tuple[str, str, str]]:
    """Extract (pos, ref, seq) from a serialized sample entry's variants."""
    out: set[tuple[str, str, str]] = set()
    for group in entry.get("variants", {}).values():
        for v in group:
            out.add((str(v["pos"]), v["ref"], v["seq"]))
    return out


# ---------------------------------------------------------------------------
# Parametrized edit -> expected canonical call
# ---------------------------------------------------------------------------


def _case_snp(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    return {"HV2": sub(hv2, 152, regions["HV2"][0], "C")}, {("152", "T", "C")}, {"HV2": [[73, 340]]}


def _case_iupac(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # A real IUPAC ambiguity code is a callable base -> SNP with seq = the code.
    return {"HV2": sub(hv2, 151, regions["HV2"][0], "Y")}, {("151", "C", "Y")}, {"HV2": [[73, 340]]}


def _case_n_unread(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # Query N is unread: no variant, and the position is excluded from the interval.
    return {"HV2": sub(hv2, 151, regions["HV2"][0], "N")}, set(), {"HV2": [[73, 150], [152, 340]]}


def _case_n_internal_split(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # An internal N splits the callable interval into two sub-intervals.
    return {"HV2": sub(hv2, 200, regions["HV2"][0], "N")}, set(), {"HV2": [[73, 199], [201, 340]]}


def _case_n_trailing(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # Replace the trailing 3 bases (338-340) with N: no variant, interval clipped.
    return {"HV2": hv2[: 338 - regions["HV2"][0]] + "NNN"}, set(), {"HV2": [[73, 337]]}


def _case_ins_run1(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # +C in the first polyC run (303-309) lands at 309.1 (split-at-T convention).
    return {"HV2": ins_after(hv2, 309, regions["HV2"][0], "C")}, {("309.1", "-", "C")}, {"HV2": [[73, 340]]}


def _case_ins_run2(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # +C in the second polyC run (311-315) lands at 315.1.
    return {"HV2": ins_after(hv2, 315, regions["HV2"][0], "C")}, {("315.1", "-", "C")}, {"HV2": [[73, 340]]}


def _case_ins_double(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # Two +C after 309 -> 309.1, 309.2 (sequential numbering at one anchor).
    return (
        {"HV2": ins_after(hv2, 309, regions["HV2"][0], "CC")},
        {("309.1", "-", "C"), ("309.2", "-", "C")},
        {"HV2": [[73, 340]]},
    )


def _case_ins_dinuc(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv3 = slice_region(ref, regions, "HV3")
    # +CA in the poly-CA run (514-524) lands at 524.1 + 524.2 (block normalization).
    return (
        {"HV3": ins_after(hv3, 524, regions["HV3"][0], "CA")},
        {("524.1", "-", "A"), ("524.2", "-", "C")},
        {"HV3": [[438, 576]]},
    )


def _case_del_homopolymer(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    # Delete one C at 303 -> right-normalized to the run end (309).
    return {"HV2": sub(hv2, 303, regions["HV2"][0], "")}, {("309", "C", "-")}, {"HV2": [[73, 340]]}


def _case_del_dinuc(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv3 = slice_region(ref, regions, "HV3")
    # Delete a full CA unit (514,515) -> right-normalized to 523,524.
    start = regions["HV3"][0]
    return (
        {"HV3": hv3[: 514 - start] + hv3[516 - start :]},
        {("523", "A", "-"), ("524", "C", "-")},
        {"HV3": [[438, 576]]},
    )


def _case_conv_polyc_a1(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    start = regions["HV2"][0]
    # 309 C>T + 310 T>C -> the polyC C-run shift convention: del 309 + ins 315.1.
    seq = sub(hv2, 309, start, "T")
    seq = sub(seq, 310, start, "C")
    return {"HV2": seq}, {("309", "C", "-"), ("315.1", "-", "C")}, {"HV2": [[73, 340]]}


def _case_conv_polyc_a2(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv2 = slice_region(ref, regions, "HV2")
    start = regions["HV2"][0]
    # 308 C>T + del 310 -> del 308 + del 309 + ins 315.1.
    seq = sub(hv2, 308, start, "T")
    seq = sub(seq, 310, start, "")
    return {"HV2": seq}, {("308", "C", "-"), ("309", "C", "-"), ("315.1", "-", "C")}, {"HV2": [[73, 340]]}


def _case_conv_513_514(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv3 = slice_region(ref, regions, "HV3")
    start = regions["HV3"][0]
    # del 513 G + del 514 C -> SNP 513 G>A + del 523 A + del 524 C.
    return (
        {"HV3": hv3[: 513 - start] + hv3[515 - start :]},
        {("513", "G", "A"), ("523", "A", "-"), ("524", "C", "-")},
        {"HV3": [[438, 576]]},
    )


def _case_del_16189_convention(ref: str, regions: dict[str, list[int]]) -> CaseResult:
    hv1 = slice_region(ref, regions, "HV1")
    # Drop the singleton T at 16189 (between two polyC runs) -> SNP 16189 T>C + del 16193.
    return (
        {"HV1": sub(hv1, 16189, regions["HV1"][0], "")},
        {("16189", "T", "C"), ("16193", "C", "-")},
        {"HV1": [[16024, 16365]]},
    )


CASES: list[tuple[str, CaseBuilder]] = [
    ("snp", _case_snp),
    ("iupac", _case_iupac),
    ("n_unread", _case_n_unread),
    ("n_internal_split", _case_n_internal_split),
    ("n_trailing", _case_n_trailing),
    ("ins_run1", _case_ins_run1),
    ("ins_run2", _case_ins_run2),
    ("ins_double", _case_ins_double),
    ("ins_dinuc", _case_ins_dinuc),
    ("del_homopolymer", _case_del_homopolymer),
    ("del_dinuc", _case_del_dinuc),
    ("del_16189_convention", _case_del_16189_convention),
    ("conv_polyc_a1", _case_conv_polyc_a1),
    ("conv_polyc_a2", _case_conv_polyc_a2),
    ("conv_513_514", _case_conv_513_514),
]


@pytest.mark.parametrize(("name", "build"), CASES, ids=[c[0] for c in CASES])
def test_edit_to_canonical_call(
    name: str,
    build: CaseBuilder,
    caller: FastaToSample,
    ref: str,
    regions: dict[str, list[int]],
    tmp_path: Path,
) -> None:
    """A known edit on an rCRS slice -> the expected canonical variants+intervals."""
    records, expected_variants, expected_intervals = build(ref, regions)
    sample = caller.call_variants(write_fasta(tmp_path / f"{name}.fasta", records))
    assert variants_set(sample) == expected_variants
    assert sample.intervals == expected_intervals


# ---------------------------------------------------------------------------
# Structural / provenance / robustness behavior
# ---------------------------------------------------------------------------


def test_missing_region_ok(caller: FastaToSample, ref: str, regions: dict[str, list[int]], tmp_path: Path) -> None:
    """A FASTA with only HV1 carries only HV1 and does not crash."""
    hv1 = slice_region(ref, regions, "HV1")
    sample = caller.call_variants(write_fasta(tmp_path / "hv1only.fasta", {"HV1": hv1}))
    assert sample.intervals is not None
    assert set(sample.intervals) == {"HV1"}
    assert sample.intervals["HV1"] == [[16024, 16365]]
    assert sample.variants == []


def test_lowercase_query_and_id_robustness(
    caller: FastaToSample, ref: str, regions: dict[str, list[int]], tmp_path: Path
) -> None:
    """A lowercase consensus with a lowercase record id is uppercased and mapped."""
    hv2 = slice_region(ref, regions, "HV2")
    records = {"hv2": sub(hv2, 152, regions["HV2"][0], "C").lower()}
    sample = caller.call_variants(write_fasta(tmp_path / "lower.fasta", records))
    assert variants_set(sample) == {("152", "T", "C")}


def test_source_tool_provenance(caller: FastaToSample, ref: str, regions: dict[str, list[int]], tmp_path: Path) -> None:
    """Source tool and information are tagged as fasta."""
    hv2 = slice_region(ref, regions, "HV2")
    sample = caller.call_variants(write_fasta(tmp_path / "prov.fasta", {"HV2": hv2}))
    assert sample.source_tool == Tool.FASTA
    assert sample.information is not None
    assert sample.information["source_tool"] == "fasta"


def test_insertion_position_is_fractional_str(
    caller: FastaToSample, ref: str, regions: dict[str, list[int]], tmp_path: Path
) -> None:
    """Insertion positions are fractional strings (e.g. '309.1'), not ints."""
    hv2 = slice_region(ref, regions, "HV2")
    seq = ins_after(hv2, 309, regions["HV2"][0], "C")
    sample = caller.call_variants(write_fasta(tmp_path / "frac.fasta", {"HV2": seq}))
    ins_positions = [v.pos for v in sample.variants if v.ref == "-"]
    assert len(ins_positions) >= 1
    assert all(isinstance(p, str) and "." in p for p in ins_positions)


def test_multi_region_round_trip(
    caller: FastaToSample, ref: str, regions: dict[str, list[int]], tmp_path: Path
) -> None:
    """Reconstructed HV sequences from called variants equal the input consensus."""
    hv1 = slice_region(ref, regions, "HV1")
    hv2 = slice_region(ref, regions, "HV2")
    hv3 = slice_region(ref, regions, "HV3")
    records = {"HV1": hv1, "HV2": sub(hv2, 152, regions["HV2"][0], "C"), "HV3": hv3}
    sample = caller.call_variants(write_fasta(tmp_path / "multi.fasta", records))
    rec = generate_sequence(sample, REF_PATH)
    for region, seq in records.items():
        assert "".join(rec.hv_seqs[region]["seq"]) == seq


def test_16189_round_trip(caller: FastaToSample, ref: str, regions: dict[str, list[int]], tmp_path: Path) -> None:
    """The 16189 singleton-deletion convention still round-trips the consensus."""
    hv1 = slice_region(ref, regions, "HV1")
    seq = sub(hv1, 16189, regions["HV1"][0], "")
    sample = caller.call_variants(write_fasta(tmp_path / "t16189.fasta", {"HV1": seq}))
    rec = generate_sequence(sample, REF_PATH)
    assert "".join(rec.hv_seqs["HV1"]["seq"]) == seq


def test_tool_enum_smoke() -> None:
    """Tool.FASTA exists and serializes as 'fasta'."""
    assert Tool.FASTA == "fasta"


# ---------------------------------------------------------------------------
# Batch output layout (synthetic input dir, not a real batch)
# ---------------------------------------------------------------------------


def _synthetic_batch_dir(tmp_path: Path, ref: str, regions: dict[str, list[int]]) -> Path:
    """Two synthetic per-sample FASTAs: one with a SNP, one clean."""
    hv2 = slice_region(ref, regions, "HV2")
    inp = tmp_path / "input"
    inp.mkdir()
    write_fasta(inp / "sampleA.fasta", {"HV2": sub(hv2, 152, regions["HV2"][0], "C")})
    write_fasta(inp / "sampleB.fasta", {"HV2": hv2})
    return inp


def test_process_batch_writes_layout(tmp_path: Path, ref: str, regions: dict[str, list[int]]) -> None:
    """process_batch writes statistic_fullbatch.json + <id>/<id>.json flat layout."""
    inp = _synthetic_batch_dir(tmp_path, ref, regions)
    out_path = FastaToSample(REF_PATH).process_batch(inp, tmp_path / "out", "synbatch")
    assert out_path is not None
    assert out_path.name == "statistic_fullbatch.json"
    assert (out_path.parent / "sampleA" / "sampleA.json").exists()
    assert (out_path.parent / "sampleB" / "sampleB.json").exists()
    data = json.loads(out_path.read_text(encoding="utf-8"))
    assert set(data) == {"sampleA", "sampleB"}
    assert variants_set_from_json(data["sampleA"]) == {("152", "T", "C")}
    assert variants_set_from_json(data["sampleB"]) == set()


def test_process_sample_writes_batch(tmp_path: Path, ref: str, regions: dict[str, list[int]]) -> None:
    """process_sample writes a batch JSON derived from the single file stem."""
    hv2 = slice_region(ref, regions, "HV2")
    fasta = write_fasta(tmp_path / "sampleA.fasta", {"HV2": sub(hv2, 152, regions["HV2"][0], "C")})
    out_path = FastaToSample(REF_PATH).process_sample(fasta, tmp_path / "out", "synbatch")
    assert out_path.name == "statistic_fullbatch.json"
    assert (out_path.parent / "sampleA" / "sampleA.json").exists()
    data = json.loads(out_path.read_text(encoding="utf-8"))
    assert "sampleA" in data


def test_batch_write_layout(tmp_path: Path, ref: str, regions: dict[str, list[int]]) -> None:
    """Batch.write produces the flat json/ layout for a synthetic sample."""
    hv2 = slice_region(ref, regions, "HV2")
    sample = FastaToSample(REF_PATH).call_variants(
        write_fasta(tmp_path / "sampleA.fasta", {"HV2": sub(hv2, 152, regions["HV2"][0], "C")})
    )
    sample.batch_id = "synbatch"
    json_dir = tmp_path / "synbatch" / "json"
    out_path = Batch([sample], REF_PATH).write(json_dir, "synbatch", nest_batch_id=False)
    assert out_path == (json_dir / "statistic_fullbatch.json").resolve()
    assert (json_dir / "sampleA" / "sampleA.json").exists()
