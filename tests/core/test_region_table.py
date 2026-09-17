import json
from pathlib import Path

import pytest

from src.core.region_table import write_region_grouped_table


def test_region_table_groups_authoritative_statistic_fullbatch_variants(tmp_path: Path) -> None:
    statistic_fullbatch = {
        "S1": {
            "no_snps": 2,
            "no_ins": 1,
            "no_dels": 1,
            "no_identicals": 0,
            "no_unread": 0,
            "HV1": "A",
            "HV2": "C",
            "HV3": "G",
            "variants": {
                "snps": [
                    {"pos": 16223, "ref": "C", "seq": "T", "file": ["HV1F"], "quality": []},
                    {"pos": 73, "ref": "A", "seq": "G", "file": ["HV2F"], "quality": []},
                ],
                "insertions": [{"pos": "315.1", "ref": "-", "seq": "C", "file": [], "quality": []}],
                "deletions": [{"pos": 499, "ref": "G", "seq": "-", "file": ["HV3R"], "quality": []}],
            },
            "intervals": {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
            "sample_flags": [],
            "variant_flags": {},
            "information": {"batch_id": "BATCH"},
        },
    }
    statistic_path = tmp_path / "statistic_fullbatch.json"
    table_path = tmp_path / "region_table.json"
    statistic_path.write_text(json.dumps(statistic_fullbatch), encoding="utf-8")

    write_region_grouped_table(statistic_path, table_path)

    table = json.loads(table_path.read_text(encoding="utf-8"))
    final_sample = statistic_fullbatch["S1"]
    table_sample = table["S1"]

    assert {key: value for key, value in table_sample.items() if key != "variants"} == {
        key: value for key, value in final_sample.items() if key != "variants"
    }
    assert table_sample["variants"] == {
        "HV1": [{"pos": 16223, "ref": "C", "seq": "T", "file": ["HV1F"], "quality": []}],
        "HV2": [
            {"pos": 73, "ref": "A", "seq": "G", "file": ["HV2F"], "quality": []},
            {"pos": "315.1", "ref": "-", "seq": "C", "file": [], "quality": []},
        ],
        "HV3": [{"pos": 499, "ref": "G", "seq": "-", "file": ["HV3R"], "quality": []}],
    }


def test_region_table_rejects_out_of_hv_final_variant(tmp_path: Path) -> None:
    statistic_path = tmp_path / "statistic_fullbatch.json"
    statistic_path.write_text(
        json.dumps({"S1": {"variants": {"snps": [{"pos": 1000, "ref": "A", "seq": "G"}]}}}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="outside HV1/HV2/HV3"):
        write_region_grouped_table(statistic_path, tmp_path / "region_table.json")
