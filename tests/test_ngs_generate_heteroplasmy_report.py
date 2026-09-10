"""Tests for the NGS called-heteroplasmy report generator."""

import json
from pathlib import Path

import pandas as pd

from src.modules.NGS.generate_heteroplasmy_report import generate_heteroplasmy_report


def test_generate_heteroplasmy_report_joins_raw_depths_and_filters_samples(tmp_path: Path):
    fis_file = tmp_path / "FIS.tsv"
    fis_file.write_text(
        "sample_id_base\tPosition\tMarker\tGenotype\tAlleleFrequency\tTotalDepth\t"
        "Ref(Depth):Alt(Depth)\tQC_Info\tGenotype(raw)\n"
        "A_run\tMT:204\tHVS-II\t204Y\t0.1\t100\tT(90):C(10)\t-\t204C\n"
        "B_run\tMT:214\tHVS-II\t214R\t0.4\t100\tA(60):G(40)\t-\t214G\n"
        "D_run\tMT:16192\tHVS-I\t16192Y\t0.4\t2433\tC(560):T(1020),CT(853)\t-\t16192T\n",
        encoding="utf-8",
    )
    comparison_file = tmp_path / "variant_level.tsv"
    pd.DataFrame(
        [
            {
                "Sample ID": "A_run",
                "Position": "204",
                "fis Variant": "Y",
                "Variant Flags (fis)": "Heteroplasmy at 204",
                "sanger Variant": "Y",
                "Variant Flags (sanger)": "Heteroplasmy at 204",
            },
            {
                "Sample ID": "B_run",
                "Position": "214",
                "fis Variant": "R",
                "Variant Flags (fis)": "Heteroplasmy at 214",
                "sanger Variant": "None",
                "Variant Flags (sanger)": "",
            },
            {
                "Sample ID": "C_run",
                "Position": "146",
                "fis Variant": "None",
                "Variant Flags (fis)": "",
                "sanger Variant": "Y",
                "Variant Flags (sanger)": "Heteroplasmy at 146",
            },
            {
                "Sample ID": "D_run",
                "Position": "16192",
                "fis Variant": "Y",
                "Variant Flags (fis)": "Heteroplasmy at 16192",
                "sanger Variant": "None",
                "Variant Flags (sanger)": "",
            },
        ]
    ).to_csv(comparison_file, sep="\t", index=False)
    mapping_file = tmp_path / "sample_mapping.tsv"
    mapping_file.write_text(
        "FIS_Sample\tSanger_Sample\nA_run\tSA\nB_run\tSB\nC_run\tSC\nD_run\tSD\n",
        encoding="utf-8",
    )
    correction_file = tmp_path / "correct_mapping_sanger.tsv"
    correction_file.write_text(
        "CE_Sample_ID\tSanger_Sample\tID_Change_Note\nA_run\tSA_corrected\tCorrected ID\n",
        encoding="utf-8",
    )
    sanger_file = tmp_path / "Sanger.json"
    sanger_file.write_text(
        json.dumps(
            {
                "A_run": {
                    "information": {"sanger_sample_id": "SA_corrected"},
                    "variants": [{"pos": 204, "ref": "T", "seq": "Y"}],
                    "variant_flags": {"204|T|Y": ["Heteroplasmy at 204"]},
                },
                "B_run": {"information": {"sanger_sample_id": "SB"}, "variants": [], "variant_flags": {}},
                "C_run": {
                    "information": {"sanger_sample_id": "SC"},
                    "variants": [{"pos": 146, "ref": "T", "seq": "Y"}],
                    "variant_flags": {"146|T|Y": ["Heteroplasmy at 146"]},
                },
                "D_run": {"information": {"sanger_sample_id": "SD"}, "variants": [], "variant_flags": {}},
            }
        ),
        encoding="utf-8",
    )
    sets_file = tmp_path / "sets.tsv"
    sets_file.write_text("A\tbatch\t1\nB\tbatch\t2\nC\tbatch\t3\nD\tbatch\t1\n", encoding="utf-8")
    sample_list_file = tmp_path / "autopass.txt"
    sample_list_file.write_text("A_run\nD_run\n", encoding="utf-8")
    filtered_output_file = tmp_path / "Heteroplasmy_autopass.tsv"

    report = generate_heteroplasmy_report(
        fis_file,
        comparison_file,
        mapping_file,
        sanger_file,
        correction_file,
        sets_file,
        tmp_path / "Heteroplasmy.tsv",
        sample_list_file,
        filtered_output_file,
    )
    filtered_report = pd.read_csv(filtered_output_file, sep="\t", dtype=str, keep_default_na=False)

    assert report["Heteroplasmy_Source"].tolist() == ["Both", "NGS_Only", "Sanger_Only", "NGS_Only"]
    assert report["Group"].tolist() == ["1", "2", "3", "1"]
    assert report.loc[0, "Sanger_Sample"] == "SA_corrected"
    assert report.loc[0, "Sanger_Genotype"] == "204Y"
    assert report.loc[0, "Top_Base"] == "T"
    assert report.loc[0, "Second_Frequency"] == 0.1
    assert report.loc[2, "Top_Base"] == ""
    assert report.loc[3, "FIS_Raw_Ref(Depth):Alt(Depth)"] == "C(560):T(1020),CT(853)"
    assert report.loc[3, "Top_Base"] == ""
    assert filtered_report["FIS_Sample"].tolist() == ["A_run", "D_run"]
    assert report["True_Heteroplasmy (Y/N/U)"].tolist() == ["U", "U", "U", "U"]
