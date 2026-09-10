"""FIS tool pipeline integration tests."""

from pathlib import Path

import pandas as pd

from src.core.models import Tool
from src.core.sample import load_sample_batch
from src.tools.fis import pipeline
from src.tools.fis.pipeline import process_batch


def _write_workbooks(tmp_path: Path) -> tuple[Path, Path]:
    r1_path, r2_path = tmp_path / "R1.xlsx", tmp_path / "R2.xlsx"
    with pd.ExcelWriter(r1_path) as writer:
        pd.DataFrame({"SampleCode": ["FIS_1", "FIS_fail"], "SampleQC": ["Pass", "QCFail"]}).to_excel(
            writer, sheet_name="T1.BasicStatistics", index=False
        )
    with pd.ExcelWriter(r2_path) as writer:
        pd.DataFrame(
            {
                "SampleCode": ["FIS_1"],
                "Nomenclature correction genotype": ["73G 315.1C"],
                "Nomenclature QC": ["Pass"],
                "Coverage pass rate(%)": [100],
                "Failed position": ["/"],
                "# N in consensus seq": [0],
                "Info": ["/"],
            }
        ).to_excel(writer, sheet_name="T1.MT.Summary.QC", index=False)
        pd.DataFrame(
            {
                "SampleCode": ["FIS_1", "FIS_fail"],
                "Position": ["MT:73", "MT:73"],
                "Marker": ["HVS-II", "HVS-II"],
                "Genotype": ["73g", "73G"],
            }
        ).to_excel(writer, sheet_name="T2.Genotype", index=False)
        pd.DataFrame({"SampleCode": ["FIS_1"], "Marker": ["HVS-II[MT:73-340]"], "Consensus": ["ACGT"]}).to_excel(
            writer, sheet_name="T3.Seq", index=False
        )
    return r1_path, r2_path


def test_default_region_workers_uses_three_quarters_of_available_cores(monkeypatch):
    monkeypatch.setattr(pipeline.os, "cpu_count", lambda: 12)

    assert pipeline.default_region_workers() == 9


def test_pipeline_writes_canonical_and_ngs_artifacts(tmp_path):
    r1_path, r2_path = _write_workbooks(tmp_path)
    output = tmp_path / "output"

    process_batch(r1_path, r2_path, output, Path("ref/rCRS.fasta"), "batch", write_ngs_artifacts=True)

    samples = load_sample_batch(output / "json" / "statistic_fullbatch.json")
    assert [sample.sample_id for sample in samples] == ["FIS_1"]
    assert samples[0].source_tool == Tool.FIS
    assert samples[0].information["fis_metadata"]["FIS_Variants"] == "73g"
    assert (output / "FIS.tsv").is_file()
    assert (output / "FIS_transformed.tsv").is_file()
    assert (output / "FIS.json").is_file()
    assert (output / "FIS.seq.tsv").is_file()
    assert (output / "FIS_QCFail_samples.txt").read_text(encoding="utf-8") == "FIS_fail\n"
    assert (output / "FASTA" / "FIS_1.fasta").read_text(encoding="utf-8") == ">HV2\nGGGGTACGT\n"
    assert (output / "regions" / "HV1" / "HV1_FIS_1.json").is_file()
