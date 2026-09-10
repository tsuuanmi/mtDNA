"""Regression tests for FIS-to-Sanger preparation."""

from pathlib import Path

from src.modules.NGS.prepare_fis_comparison import load_sanger_corrections


def test_global_corrections_allow_a_subset_fis_batch(tmp_path: Path):
    corrections = tmp_path / "corrections.tsv"
    corrections.write_text(
        "CE_Sample_ID\tSanger_Sample\tSanger_Batch\tID_Change_Note\n"
        "FIS_1\tSANGER_1\tBATCH\tcorrect\n"
        "FIS_OTHER\tMISSING_TARGET\tBATCH\tother batch\n",
        encoding="utf-8",
    )

    result = load_sanger_corrections(corrections, {"SANGER_1": {"batch": "BATCH"}}, {"FIS_1"})

    assert result == {"FIS_1": "SANGER_1"}
