import json
from pathlib import Path

from src.core.models import Sample, Tool, Variant
from src.tools.sequencher.etl import process
from src.tools.sequencher.pipeline import BatchOptions, process_batch
from src.tools.sequencher.quality_control import DUPLICATE_LID_ERROR, SEQUENCHER_TABLE_HEADER_ERROR, run_qc

SEQUENCHER_HEADER_ONLY_TXT = """0 differences between C01_rCRS and the consensus.

LEGEND:
Pos - Offset in C01_rCRS.
Seq - Base in C01_rCRS.
Con - Base in the consensus.
Required Edit - Change needed in C01_rCRS
to match corresponding base in the consensus.

      Pos \t  Seq   \t  Con   \tRequired Edit
"""

SEQUENCHER_VARIANT_TXT = """3 differences between C01_rCRS and the consensus.

LEGEND:
Pos - Offset in C01_rCRS.
Seq - Base in C01_rCRS.
Con - Base in the consensus.
Required Edit - Change needed in C01_rCRS
to match corresponding base in the consensus.

      Pos \t  Seq   \t  Con   \tRequired Edit

       73 \t   A    \t   G    \tChange base
    315.1 \t   :    \t   C    \tInsert base
   16,223 \t   C    \t   T    \tChange base
"""

MALFORMED_NO_HEADER_TXT = """No differences between C01_rCRS and the consensus.

This mock file is not a Sequencher variant-table TXT export.
It intentionally has no Pos/Seq/Con/Required Edit table header.
"""


def _write_txt(input_dir: Path, filename: str, content: str) -> Path:
    txt_file = input_dir / filename
    txt_file.write_text(content, encoding="utf-8")
    return txt_file


def _variant_by_position(sample: Sample, position: float | str) -> Variant:
    return next(variant for variant in sample.variants if str(variant.pos) == str(position))


def test_run_qc_fails_txt_without_variant_table_header(tmp_path) -> None:
    txt_file = _write_txt(tmp_path, "LN_26_BAD.TXT", MALFORMED_NO_HEADER_TXT)

    qc = run_qc(txt_file.name, txt_file)

    assert qc.sample_id == "LN_26_BAD"
    assert not qc.passed
    assert qc.reason == SEQUENCHER_TABLE_HEADER_ERROR


def test_run_qc_passes_txt_with_header_and_zero_variant_rows(tmp_path) -> None:
    txt_file = _write_txt(tmp_path, "LN_26_HEADER_ONLY.TXT", SEQUENCHER_HEADER_ONLY_TXT)

    qc = run_qc(txt_file.name, txt_file)

    assert qc.sample_id == "LN_26_HEADER_ONLY"
    assert qc.passed
    assert qc.reason is None


def test_process_parses_sequencher_txt_variants_and_intervals(tmp_path) -> None:
    txt_file = _write_txt(tmp_path, "LN_26_GOOD.TXT", SEQUENCHER_VARIANT_TXT)

    sample = process("LN_26_GOOD", txt_file, batch_id="TEST_BATCH")

    assert sample.sample_id == "LN_26_GOOD"
    assert sample.source_tool == Tool.SEQUENCHER
    assert sample.batch_id == "TEST_BATCH"
    assert sample.intervals == {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]}
    assert [str(variant.pos) for variant in sample.variants] == ["73", "315.1", "16223"]
    assert _variant_by_position(sample, 73).ref == "A"
    assert _variant_by_position(sample, 73).seq == "G"
    assert _variant_by_position(sample, 315.1).ref == "-"
    assert _variant_by_position(sample, 315.1).seq == "C"
    assert _variant_by_position(sample, 16223).ref == "C"
    assert _variant_by_position(sample, 16223).seq == "T"


def test_process_batch_records_missing_header_as_qc_failure(tmp_path) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    _write_txt(input_dir, "LN_26_BAD.TXT", MALFORMED_NO_HEADER_TXT)

    samples = process_batch(input_dir, output_dir, BatchOptions(batch_id="TEST_BATCH"))

    assert samples == {}
    qc_report = json.loads((output_dir / "preprocess" / "qc_report.json").read_text(encoding="utf-8"))
    assert qc_report["total"] == 1
    assert qc_report["passed"] == 0
    assert qc_report["failed"] == 1
    assert qc_report["samples"] == [
        {
            "sample_id": "LN_26_BAD",
            "filename": "LN_26_BAD.TXT",
            "status": "fail",
            "reason": SEQUENCHER_TABLE_HEADER_ERROR,
        }
    ]


def test_process_batch_rejects_all_files_with_duplicate_lid(tmp_path) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    _write_txt(input_dir, "LN_26_DUP.TXT", SEQUENCHER_VARIANT_TXT)
    _write_txt(input_dir, "LN_26_DUP-73-340 438-573 16024-16365.TXT", SEQUENCHER_VARIANT_TXT)

    samples = process_batch(input_dir, output_dir, BatchOptions(batch_id="TEST_BATCH"))

    assert samples == {}
    qc_report = json.loads((output_dir / "preprocess" / "qc_report.json").read_text(encoding="utf-8"))
    assert qc_report["total"] == 2
    assert qc_report["passed"] == 0
    assert qc_report["failed"] == 2
    assert [entry["status"] for entry in qc_report["samples"]] == ["fail", "fail"]
    assert all(entry["sample_id"] == "LN_26_DUP" for entry in qc_report["samples"])
    assert all(entry["reason"].startswith(DUPLICATE_LID_ERROR) for entry in qc_report["samples"])
    assert not (output_dir / "json" / "LN_26_DUP" / "LN_26_DUP.json").exists()


def test_process_batch_processes_valid_sequencher_txt_and_writes_outputs(tmp_path) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    _write_txt(input_dir, "LN_26_GOOD.TXT", SEQUENCHER_VARIANT_TXT)

    samples = process_batch(input_dir, output_dir, BatchOptions(batch_id="TEST_BATCH"))

    assert list(samples) == ["LN_26_GOOD"]
    sample = samples["LN_26_GOOD"]
    assert [str(variant.pos) for variant in sample.variants] == ["73", "315.1", "16223"]

    qc_report = json.loads((output_dir / "preprocess" / "qc_report.json").read_text(encoding="utf-8"))
    assert qc_report["passed"] == 1
    assert qc_report["failed"] == 0
    assert qc_report["samples"][0]["status"] == "pass"

    sample_json = output_dir / "json" / "LN_26_GOOD" / "LN_26_GOOD.json"
    full_batch_json = output_dir / "json" / "statistic_fullbatch.json"
    assert sample_json.exists()
    assert full_batch_json.exists()
