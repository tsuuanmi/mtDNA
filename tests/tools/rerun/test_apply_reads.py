"""Thay read chạy lại theo từng vùng."""
from __future__ import annotations

import json

import pytest

from src.tools.rerun.apply_reads import (
    MANIFEST_NAME,
    apply,
    batch_token,
    folder_matches,
    parse_ab1_name,
)

# Thư mục bên Lab mang ngày CHẠY LẠI, không phải ngày batch -- nên nhận
# nhau qua phần định danh còn lại.
BATCH = "20260610_mtDNA_77"
TOKEN = "mtDNA_77"


@pytest.mark.parametrize(
    ("name", "sample", "region", "direction", "date"),
    [
        # Dạng chuẩn
        ("G12_20260610_2555964_HV1F_21.ab1", "2555964", "HV1", "F", "20260610"),
        # Vị trí mồi nằm giữa vùng và chiều -- vẫn là HV1 chiều F
        ("A01_20260715_2555964_HV1-15930F_01.ab1", "2555964", "HV1", "F", "20260715"),
        ("H01_20260731_2554833_HV3-350F_22.ab1", "2554833", "HV3", "F", "20260731"),
        # Chiều đứng trước vị trí mồi
        ("A12_20260823_LN_26_AB6721_HV1R-16410_03.ab1", "LN_26_AB6721", "HV1", "R", "20260823"),
        # Mã mẫu chứa dấu gạch dưới: đếm trường từ trái sẽ cắt nhầm
        ("G11_20250710_257715_2_HV1R_20.ab1", "257715_2", "HV1", "R", "20250710"),
        # Vùng chiếm hai trường
        ("B05_20250830_2522691_HV2_29F_05.ab1", "2522691", "HV2", "F", "20250830"),
    ],
)
def test_parse_ab1_name(name, sample, region, direction, date):
    read = parse_ab1_name(_p(name))
    assert read is not None, name
    assert (read.sample, read.region, read.direction, read.date) == (sample, region, direction, date)


@pytest.mark.parametrize("name", [
    "khong_phai_ab1.txt",
    "B01_20241810_HV1F_242018_04.ab1",  # vùng đứng trước mã mẫu
    "thieu_truong.ab1",
])
def test_parse_ab1_name_tu_choi(name):
    assert parse_ab1_name(_p(name)) is None


def _p(name):
    from pathlib import Path
    return Path("/tmp") / name


def _write(directory, *names):
    directory.mkdir(parents=True, exist_ok=True)
    for n in names:
        (directory / n).write_bytes(b"ab1")


def test_thay_dung_vung_giu_vung_khac(tmp_path):
    """Lab chạy lại một vùng thì chỉ vùng đó bị thay."""
    raw = tmp_path / "raw" / BATCH
    _write(raw,
           "G12_20260610_S1_HV1F_21.ab1", "G12_20260610_S1_HV1R_21.ab1",
           "G12_20260610_S1_HV2F_21.ab1", "G12_20260610_S1_HV3R_21.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "A05_20260715_S1_HV2F_02.ab1")

    apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    con_lai = sorted(p.name for p in raw.glob("*.ab1"))
    assert con_lai == [
        "A05_20260715_S1_HV2F_02.ab1",   # read mới
        "G12_20260610_S1_HV1F_21.ab1",   # ba vùng kia giữ nguyên
        "G12_20260610_S1_HV1R_21.ab1",
        "G12_20260610_S1_HV3R_21.ab1",
    ]


def test_them_khi_goc_chua_co_vung_do(tmp_path):
    raw = tmp_path / "raw" / BATCH
    _write(raw, "G12_20260610_S1_HV3R_21.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "H01_20260731_S1_HV3-350F_22.ab1")

    manifest = apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    assert not manifest["replaced"]
    assert [r["region"] for r in manifest["added"]] == ["HV3F"]
    assert (raw / "G12_20260610_S1_HV3R_21.ab1").exists()


def test_chay_lai_hai_lan_lay_ngay_moi_nhat_theo_tung_vung(tmp_path):
    """Mới nhất tính riêng từng vùng, không phải mới nhất toàn mẫu."""
    raw = tmp_path / "raw" / BATCH
    _write(raw, "G12_20260610_S1_HV2F_21.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "A01_20260715_S1_HV2F_01.ab1",
           "A02_20260715_S1_HV3F_01.ab1")
    _write(tmp_path / "rerun" / f"20260731_{TOKEN}_1 mẫu", "B01_20260731_S1_HV2F_09.ab1")

    manifest = apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    ten = sorted(p.name for p in raw.glob("*.ab1"))
    assert ten == ["A02_20260715_S1_HV3F_01.ab1", "B01_20260731_S1_HV2F_09.ab1"]
    assert manifest["replaced"][0]["new_date"] == "20260731"


def test_ten_khong_doc_duoc_thi_bo_qua_ca_mau(tmp_path):
    """Đè nhầm vùng tệ hơn không đè, nên gặp tên lạ là bỏ nguyên mẫu."""
    raw = tmp_path / "raw" / BATCH
    _write(raw, "B01_20241810_HV1F_S1_04.ab1", "G12_20260610_S1_HV2F_21.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "A05_20260715_S1_HV2F_02.ab1")

    manifest = apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    assert manifest["skipped_samples"] == ["S1"]
    assert not manifest["replaced"] and not manifest["added"]
    assert (raw / "G12_20260610_S1_HV2F_21.ab1").exists()  # không đụng gì


def test_khong_dung_mau_khac(tmp_path):
    raw = tmp_path / "raw" / BATCH
    _write(raw, "G12_20260610_S1_HV2F_21.ab1", "G12_20260610_S2_HV2F_22.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "A05_20260715_S1_HV2F_02.ab1")

    apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    assert (raw / "G12_20260610_S2_HV2F_22.ab1").exists()


def test_dry_run_khong_dung_file(tmp_path):
    raw = tmp_path / "raw" / BATCH
    _write(raw, "G12_20260610_S1_HV2F_21.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "A05_20260715_S1_HV2F_02.ab1")

    manifest = apply(BATCH, {"S1"}, tmp_path / "rerun", raw, dry_run=True)

    assert manifest["replaced"]
    assert sorted(p.name for p in raw.glob("*.ab1")) == ["G12_20260610_S1_HV2F_21.ab1"]
    assert not (raw / MANIFEST_NAME).exists()


def test_ghi_ke_khai_de_buoc_2_biet(tmp_path):
    """Bước 2 nạp lại từ nguồn gốc sẽ mang read cũ về; kê khai để nó cảnh báo."""
    raw = tmp_path / "raw" / BATCH
    _write(raw, "G12_20260610_S1_HV2F_21.ab1")
    _write(tmp_path / "rerun" / f"20260714_{TOKEN}_2 mẫu", "A05_20260715_S1_HV2F_02.ab1")

    apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    ghi = json.loads((raw / MANIFEST_NAME).read_text(encoding="utf-8"))
    assert ghi["batch_id"] == BATCH
    assert ghi["samples"] == ["S1"]
    assert ghi["replaced"][0]["replaced_files"] == ["G12_20260610_S1_HV2F_21.ab1"]


@pytest.mark.parametrize(("batch", "token"), [
    ("20260610_mtDNA_427", "mtDNA_427"),
    ("MS_151225_005", "MS_151225_005"),   # không bắt đầu bằng ngày 8 số
    ("mtDNA_9", "mtDNA_9"),
])
def test_batch_token(batch, token):
    assert batch_token(batch) == token


@pytest.mark.parametrize(("folder", "khop"), [
    ("20260714_mtDNA_32_2 mẫu", True),
    ("20260714_mtDNA_32", True),
    ("20260714_mtDNA_320_1 mẫu", False),   # tiền tố, không được trúng
    ("20260714_mtDNA_328", False),
    ("20260714_mtDNA_3", False),
])
def test_folder_matches_bien(folder, khop):
    assert folder_matches(folder, "mtDNA_32") is khop


def test_bo_qua_thu_muc_cua_batch_khac(tmp_path):
    """Mã mẫu trùng nhau giữa hai batch thì thư mục mới là thứ phân định."""
    raw = tmp_path / "raw" / BATCH
    _write(raw, "G12_20260610_S1_HV2F_21.ab1")
    # cùng mã mẫu, nhưng nằm ở thư mục của batch khác
    _write(tmp_path / "rerun" / "20260714_mtDNA_78_1 mẫu", "A05_20260715_S1_HV2F_02.ab1")

    manifest = apply(BATCH, {"S1"}, tmp_path / "rerun", raw)

    assert not manifest["replaced"] and not manifest["added"]
    assert (raw / "G12_20260610_S1_HV2F_21.ab1").exists()
