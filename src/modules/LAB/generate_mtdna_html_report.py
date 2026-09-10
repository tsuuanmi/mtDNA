#!/usr/bin/env python
"""
Generate HTML tracking report for mtDNA samples.

This script generates an HTML report similar to the STR tracking template,
using concordance results from compare_sequencher.py output.
"""

import sys
from argparse import ArgumentParser
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger
from openpyxl.styles import Alignment, Border, Font, Side

# Default signature names (pre-filled as per template)
DEFAULT_SIGNATURES = {
    "nguoi_thuc_hien": "Triệu Thị Nguyệt",
    "nguoi_kiem_tra": "Đặng Hữu Điện",
    "nguoi_phan_tich_1": "Nguyễn Thị Phương",
    "nguoi_phan_tich_2": "Trần Nhật Tân",
    "nguoi_phe_duyet": "Nguyễn Ngọc Nam",
}


def format_date(dt: datetime) -> str:
    """Format datetime to DD/MM/YYYY display format."""
    return dt.strftime("%d/%m/%Y")


def _parse_date_from_row(row: pd.Series, col: str, batch_suffix: str, date_label: str) -> str | None:
    """Parse a date string from a DataFrame row column."""
    date_str = str(row[col]).strip() if pd.notna(row[col]) else ""  # type: ignore[reportGeneralTypeIssues]
    if not date_str or date_str == "nan":
        return None
    try:
        parsed_date = datetime.strptime(date_str[:10], "%Y-%m-%d")  # noqa: DTZ007
        formatted = format_date(parsed_date)
        logger.info(f"Found {date_label} for {batch_suffix}: {formatted}")
    except ValueError:
        logger.warning(f"Could not parse {date_label}: {date_str}")
        return None
    else:
        return formatted



def load_tracking_tsv(tsv_path: Path, batch_suffix: str) -> tuple[str | None, str | None]:
    """
    Load tracking TSV and find dates for a batch.

    Args:
        tsv_path: Path to tracking TSV file
        batch_suffix: Batch suffix to search for (e.g., 'mtDNA_312')

    Returns:
        Tuple of (sequence_released_date, run1_analysis_released_date) formatted as DD/MM/YYYY,
        or (None, None) if not found
    """
    if not tsv_path or not tsv_path.exists():
        logger.info("No tracking TSV provided or file not found")
        return None, None

    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)

        # Find the batch row
        batch_col = "Batch mtDNA"
        sequence_col = "Sequence Released Date"
        release_col = "Run1 Analysis Results Released Date"

        if batch_col not in df.columns:
            logger.warning(f"Column '{batch_col}' not found in tracking TSV")
            return None, None

        # Search for the batch
        for _, row in df.iterrows():
            batch_name = str(row[batch_col]).strip() if pd.notna(row[batch_col]) else ""  # type: ignore[reportGeneralTypeIssues]
            if batch_name == batch_suffix:
                # Parse sequence released date
                sequence_date_formatted = (
                    _parse_date_from_row(row, sequence_col, batch_suffix, "sequence date")
                    if sequence_col in df.columns
                    else None
                )

                # Parse run1 analysis released date
                release_date_formatted = (
                    _parse_date_from_row(row, release_col, batch_suffix, "release date")
                    if release_col in df.columns
                    else None
                )

                return sequence_date_formatted, release_date_formatted

        logger.info(f"Batch {batch_suffix} not found in tracking TSV")

    except (OSError, pd.errors.ParserError, ValueError, KeyError) as e:
        logger.error(f"Error loading tracking TSV: {e}")
        return None, None
    else:
        return None, None


def load_metadata_excel(excel_path: Path) -> dict[str, str]:
    """
    Load metadata Excel and create LID to Barcode mapping.

    Returns:
        Dictionary mapping LID (as string) to Barcode
    """
    if not excel_path.exists():
        logger.warning(f"Metadata Excel not found: {excel_path}")
        return {}

    try:
        df = pd.read_excel(excel_path, dtype=str)

        # Find LID and Barcode columns (case-insensitive)
        lid_col = None
        barcode_col = None

        for col in df.columns:
            col_lower = col.lower()
            if "lid" in col_lower:
                lid_col = col
            elif "barcode" in col_lower:
                barcode_col = col

        if lid_col is None or barcode_col is None:
            logger.warning(f"Could not find LID or Barcode columns in {excel_path}")
            return {}

        # Create mapping
        mapping = {}
        for _, row in df.iterrows():
            lid = str(row[lid_col]).strip() if pd.notna(row[lid_col]) else ""  # type: ignore[reportGeneralTypeIssues]
            barcode = str(row[barcode_col]).strip() if pd.notna(row[barcode_col]) else ""  # type: ignore[reportGeneralTypeIssues]
            if lid:
                # Remove .0 suffix if present (from numeric conversion)
                lid = lid.removesuffix(".0")
                mapping[lid] = barcode

    except (OSError, ValueError, KeyError, pd.errors.ParserError) as e:
        logger.error(f"Error loading metadata Excel: {e}")
        return {}
    else:
        logger.info(f"Loaded {len(mapping)} LID-Barcode mappings from metadata")
        return mapping



def _find_rerun_columns(df: pd.DataFrame) -> tuple[str | None, str | None, str | None]:
    """Find ID, Issues, and Requirement columns in a rerun Excel DataFrame."""
    id_col = None
    issues_col = None
    requirement_col = None

    for col in df.columns:
        col_lower = col.lower().strip()
        if col_lower == "id":
            id_col = col
        elif col_lower == "issues":
            issues_col = col
        elif col_lower == "requirement":
            requirement_col = col

    return id_col, issues_col, requirement_col



def load_rerun_excel(excel_path: Path) -> dict[str, dict[str, str]]:
    """
    Load rerun Excel and create ID to Issues/Requirement mapping.

    Args:
        excel_path: Path to rerun Excel file with 'Rerun' sheet

    Returns:
        Dictionary mapping ID (as string) to {"issues": str, "requirement": str}
    """
    if not excel_path or not excel_path.exists():
        logger.info("No rerun Excel provided or file not found")
        return {}

    try:
        df = pd.read_excel(excel_path, sheet_name="Rerun", dtype=str)

        # Find ID, Issues, and Requirement columns
        id_col, issues_col, requirement_col = _find_rerun_columns(df)

        if id_col is None:
            logger.warning(f"Could not find ID column in rerun Excel {excel_path}")
            return {}

        # Create mapping
        mapping = {}
        for _, row in df.iterrows():
            sample_id = str(row[id_col]).strip() if pd.notna(row[id_col]) else ""  # type: ignore[reportGeneralTypeIssues]
            if sample_id:
                # Remove .0 suffix if present
                sample_id = sample_id.removesuffix(".0")

                issues = str(row[issues_col]).strip() if issues_col and pd.notna(row[issues_col]) else ""  # type: ignore[reportGeneralTypeIssues]
                requirement = (
                    str(row[requirement_col]).strip() if requirement_col and pd.notna(row[requirement_col]) else ""  # type: ignore[reportGeneralTypeIssues]
                )

                mapping[sample_id] = {"issues": issues, "requirement": requirement}

    except (OSError, ValueError, KeyError, pd.errors.ParserError) as e:
        logger.error(f"Error loading rerun Excel: {e}")
        return {}
    else:
        logger.info(f"Loaded {len(mapping)} rerun entries from Excel")
        return mapping


def extract_date_and_batch(batch_id: str) -> tuple[datetime | None, str]:
    """
    Extract analysis date and batch suffix from batch_id.

    Args:
        batch_id: Batch ID in one of these formats:
            - YYYYMMDD_mtDNA_XX (e.g., 20251120_mtDNA_312)
            - MS_DDMMYY_XXX (e.g., MS_251125_004) - date will be None, must use tracking TSV

    Returns:
        Tuple of (analysis_date as datetime or None, batch_suffix)
    """
    # Check if batch_id starts with MS_ - don't extract date, use tracking TSV
    if batch_id.startswith("MS_"):
        # Use full batch_id as batch_suffix for MS batches
        return None, batch_id

    # Standard format: YYYYMMDD_mtDNA_XX
    # Extract date from first 8 characters (YYYYMMDD)
    date_str = batch_id[:8]
    year = int(date_str[:4])
    month = int(date_str[4:6])
    day = int(date_str[6:8])
    analysis_date = datetime(year, month, day)  # noqa: DTZ001

    # Extract batch suffix (everything after the first underscore)
    batch_suffix = batch_id.split("_", 1)[1] if "_" in batch_id else batch_id

    return analysis_date, batch_suffix


def generate_html_report(
    metadata_mapping: dict[str, str],
    batch_suffix: str,
    analysis_date: datetime | None,
    rerun_mapping: dict[str, dict[str, str]] | None = None,
    tracking_sequence_date: str | None = None,
    tracking_release_date: str | None = None,
) -> tuple[str, list[dict[str, Any]], str]:
    """Generate HTML report content and data rows for Excel export.

    Returns:
        Tuple of (html_content, excel_rows, result_date_str)
    """
    # Use tracking sequence date if provided for "Ngày thực hiện", otherwise use analysis_date
    if tracking_sequence_date:
        analysis_date_str = tracking_sequence_date
    elif analysis_date:
        analysis_date_str = format_date(analysis_date)
    else:
        analysis_date_str = ""  # No date available

    # Use tracking release date if provided for "Ngày trả kết quả", otherwise calculate +2 days
    if tracking_release_date:
        result_date_str = tracking_release_date
    elif analysis_date:
        result_date = analysis_date + timedelta(days=2)
        result_date_str = format_date(result_date)
    else:
        result_date_str = ""  # No date available

    # Initialize rerun_mapping if None
    if rerun_mapping is None:
        rerun_mapping = {}

    # Generate table rows - iterate through all samples from metadata
    rows_html = ""
    excel_rows = []
    for idx, (lid, barcode) in enumerate(metadata_mapping.items(), start=1):
        sample_id = lid

        # Default values
        chi_tiet = ""
        hanh_dong = "Trả KQ"

        # Override with rerun data if available
        if sample_id in rerun_mapping:
            rerun_info = rerun_mapping[sample_id]
            if rerun_info.get("issues"):
                chi_tiet = rerun_info["issues"]
            if rerun_info.get("requirement"):
                hanh_dong = rerun_info["requirement"]

        # Add row for Excel export
        excel_rows.append(
            {
                "STT": idx,
                "LID": sample_id,
                "Barcode": barcode,
                "Ngày thực hiện": analysis_date_str,
                "Ngày trả kết quả": result_date_str,
                "Chi tiết": chi_tiet,
                "Hành động": hanh_dong,
                "Ngày thực hiện lần 2": "",
                "Kết quả lần 2": "",
                "Kết quả cuối": "",
            },
        )

        rows_html += f"""
        <tr>
            <td>{idx}</td>
            <td>{sample_id}</td>
            <td>{barcode}</td>
            <td>{analysis_date_str}</td>
            <td>{result_date_str}</td>
            <td>{chi_tiet}</td>
            <td>{hanh_dong}</td>
            <td></td>
            <td></td>
            <td></td>
        </tr>"""

    # Generate complete HTML
    html_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
body {{
    font-family: Arial, sans-serif;
    margin: 20px;
    font-size: 12pt;
}}
h1 {{
    text-align: center;
    font-size: 16pt;
    font-weight: bold;
    margin-bottom: 10px;
}}
.header-info {{
    display: flex;
    justify-content: space-between;
    margin-bottom: 15px;
    font-size: 11pt;
}}
table {{
    border-collapse: collapse;
    width: 100%;
    font-size: 10pt;
}}
th, td {{
    border: 1px solid #000;
    padding: 5px;
    text-align: center;
}}
th {{
    background-color: #f0f0f0;
    font-weight: bold;
}}
.signature-section {{
    margin-top: 30px;
    font-size: 11pt;
}}
.signature-table {{
    width: 100%;
    border-collapse: collapse;
    margin-top: 10px;
}}
.signature-table td {{
    border: 1px solid #000;
    padding: 8px;
    vertical-align: top;
}}
.signature-label {{
    font-weight: normal;
}}
.signature-name {{
    color: #000000;
    font-weight: normal;
}}
.signature-date {{
    color: #000000;
}}
@media print {{
    body {{
        margin: 5mm;
    }}
    @page {{
        size: A4 landscape;
        margin: 10mm;
    }}
}}
</style>
</head>
<body>

<h1>PHIẾU THEO DÕI SOÁT XÉT KẾT QUẢ XÉT NGHIỆM {batch_suffix}</h1>

<div class="header-info">
    <span>Tên XN: Giải trình tự ADN ty thể</span>
    <span>Phương pháp: XN-QTKT.011.01</span>
</div>

<table>
<thead>
    <tr>
        <th>STT</th>
        <th>LID</th>
        <th>Barcode</th>
        <th>Ngày thực hiện</th>
        <th>Ngày trả kết quả</th>
        <th>Chi tiết</th>
        <th>Hành động</th>
        <th>Ngày thực hiện lần 2</th>
        <th>Kết quả lần 2</th>
        <th>Kết quả cuối</th>
    </tr>
</thead>
<tbody>
{rows_html}
</tbody>
</table>

<p style="margin-top: 15px; margin-bottom: 10px; font-size: 11pt; font-style: italic;">
    Phê duyệt dựa trên hồ sơ mtDNA tạo lập từ kết quả phân tích trình tự AB1
</p>

<div class="signature-section">
    <table class="signature-table">
        <tr>
            <td colspan="3">
                <span class="signature-label">Ngày phân tích Lần 1:</span>
                <span class="signature-date">{result_date_str}</span>
            </td>
            <td colspan="3">
                <span class="signature-label">Ngày phân tích Lần 2:</span>
            </td>
        </tr>
        <tr>
            <td colspan="3" style="height: 60px;">
                <span class="signature-label">Người phân tích 1 (PT1):</span>
                <span class="signature-name">{DEFAULT_SIGNATURES["nguoi_phan_tich_1"]}</span>
            </td>
            <td colspan="3" style="height: 60px;">
                <span class="signature-label">Người phân tích 1 (PT1):</span>
            </td>
        </tr>
        <tr>
            <td colspan="3" style="height: 60px;">
                <span class="signature-label">Người phân tích 2 (PT2):</span>
                <span class="signature-name">{DEFAULT_SIGNATURES["nguoi_phan_tich_2"]}</span>
            </td>
            <td colspan="3" style="height: 60px;">
                <span class="signature-label">Người phân tích 2 (PT2):</span>
            </td>
        </tr>
        <tr>
            <td colspan="2">
                <span class="signature-label">Ngày thực hiện:</span>
                <span class="signature-date">{result_date_str}</span>
            </td>
            <td colspan="2">
                <span class="signature-label">Ngày kiểm tra:</span>
                <span class="signature-date">{result_date_str}</span>
            </td>
            <td colspan="2">
                <span class="signature-label">Ngày phê duyệt:</span>
                <span class="signature-date">{result_date_str}</span>
            </td>
        </tr>
        <tr>
            <td colspan="2" style="height: 60px;">
                <span class="signature-label">Người thực hiện:</span>
                <span class="signature-name">{DEFAULT_SIGNATURES["nguoi_thuc_hien"]}</span>
            </td>
            <td colspan="2" style="height: 60px;">
                <span class="signature-label">Người kiểm tra:</span>
                <span class="signature-name">{DEFAULT_SIGNATURES["nguoi_kiem_tra"]}</span>
            </td>
            <td colspan="2" style="height: 60px;">
                <span class="signature-label">Người phê duyệt:</span>
                <span class="signature-name">{DEFAULT_SIGNATURES["nguoi_phe_duyet"]}</span>
            </td>
        </tr>
    </table>
</div>

</body>
</html>"""

    return html_content, excel_rows, result_date_str


def generate_excel_report(  # noqa: PLR0915
    excel_rows: list[dict[str, Any]],
    output_path: Path,
    batch_suffix: str,
    result_date_str: str,
) -> None:
    """Generate Excel report from data rows with title and signature section."""
    df = pd.DataFrame(excel_rows)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # Write empty sheet first
        df.to_excel(writer, sheet_name=batch_suffix[:31], index=False, startrow=4)

        # Get worksheet
        worksheet = writer.sheets[batch_suffix[:31]]

        # Define styles
        thin_border = Border(
            left=Side(style="thin"),
            right=Side(style="thin"),
            top=Side(style="thin"),
            bottom=Side(style="thin"),
        )
        bold_font = Font(bold=True)
        center_align = Alignment(horizontal="center", vertical="center")

        # Add title (row 1)
        worksheet.merge_cells("A1:J1")
        title_cell = worksheet["A1"]
        title_cell.value = f"PHIẾU THEO DÕI SOÁT XÉT KẾT QUẢ XÉT NGHIỆM {batch_suffix}"
        title_cell.font = Font(bold=True, size=14)
        title_cell.alignment = center_align

        # Add header info (row 2-3)
        worksheet["A3"] = "Tên XN: Giải trình tự ADN ty thể"
        worksheet["F3"] = "Phương pháp: XN-QTKT.011.01"

        # Apply border to data table
        for row in worksheet.iter_rows(min_row=5, max_row=5 + len(df), min_col=1, max_col=10):
            for cell in row:
                cell.border = thin_border
                cell.alignment = center_align

        # Make header row bold
        for cell in worksheet[5]:
            cell.font = bold_font

        # Calculate signature section start row
        sig_start_row = 5 + len(df) + 2

        # Add note before signature section
        note_row = sig_start_row
        worksheet.merge_cells(f"A{note_row}:J{note_row}")
        note_cell = worksheet[f"A{note_row}"]
        note_cell.value = "Phê duyệt dựa trên hồ sơ mtDNA tạo lập từ kết quả phân tích trình tự AB1"
        note_cell.font = Font(italic=True, size=11)
        note_cell.alignment = Alignment(horizontal="left", vertical="center")

        # Adjust signature section to start after note
        sig_start_row = note_row + 1

        # Signature section - Row 1: Ngày phân tích
        worksheet.merge_cells(f"A{sig_start_row}:E{sig_start_row}")
        worksheet[f"A{sig_start_row}"] = f"Ngày phân tích Lần 1: {result_date_str}"
        worksheet[f"A{sig_start_row}"].border = thin_border

        worksheet.merge_cells(f"F{sig_start_row}:J{sig_start_row}")
        worksheet[f"F{sig_start_row}"] = "Ngày phân tích Lần 2:"
        worksheet[f"F{sig_start_row}"].border = thin_border

        # Apply border to merged cells
        for col in range(1, 11):
            worksheet.cell(row=sig_start_row, column=col).border = thin_border

        # Signature section - Row 2: Người phân tích 1
        sig_row2 = sig_start_row + 1
        worksheet.merge_cells(f"A{sig_row2}:E{sig_row2}")
        worksheet[f"A{sig_row2}"] = f"Người phân tích 1 (PT1): {DEFAULT_SIGNATURES['nguoi_phan_tich_1']}"
        worksheet.row_dimensions[sig_row2].height = 40

        worksheet.merge_cells(f"F{sig_row2}:J{sig_row2}")
        worksheet[f"F{sig_row2}"] = "Người phân tích 1 (PT1):"

        for col in range(1, 11):
            worksheet.cell(row=sig_row2, column=col).border = thin_border

        # Signature section - Row 3: Người phân tích 2
        sig_row3 = sig_start_row + 2
        worksheet.merge_cells(f"A{sig_row3}:E{sig_row3}")
        worksheet[f"A{sig_row3}"] = f"Người phân tích 2 (PT2): {DEFAULT_SIGNATURES['nguoi_phan_tich_2']}"
        worksheet.row_dimensions[sig_row3].height = 40

        worksheet.merge_cells(f"F{sig_row3}:J{sig_row3}")
        worksheet[f"F{sig_row3}"] = "Người phân tích 2 (PT2):"

        for col in range(1, 11):
            worksheet.cell(row=sig_row3, column=col).border = thin_border

        # Signature section - Row 4: Ngày thực hiện / kiểm tra / phê duyệt
        sig_row4 = sig_start_row + 3
        worksheet.merge_cells(f"A{sig_row4}:C{sig_row4}")
        worksheet[f"A{sig_row4}"] = f"Ngày thực hiện: {result_date_str}"

        worksheet.merge_cells(f"D{sig_row4}:F{sig_row4}")
        worksheet[f"D{sig_row4}"] = f"Ngày kiểm tra: {result_date_str}"

        worksheet.merge_cells(f"G{sig_row4}:J{sig_row4}")
        worksheet[f"G{sig_row4}"] = f"Ngày phê duyệt: {result_date_str}"

        for col in range(1, 11):
            worksheet.cell(row=sig_row4, column=col).border = thin_border

        # Signature section - Row 5: Người thực hiện / kiểm tra / phê duyệt
        sig_row5 = sig_start_row + 4
        worksheet.merge_cells(f"A{sig_row5}:C{sig_row5}")
        worksheet[f"A{sig_row5}"] = f"Người thực hiện: {DEFAULT_SIGNATURES['nguoi_thuc_hien']}"
        worksheet.row_dimensions[sig_row5].height = 40

        worksheet.merge_cells(f"D{sig_row5}:F{sig_row5}")
        worksheet[f"D{sig_row5}"] = f"Người kiểm tra: {DEFAULT_SIGNATURES['nguoi_kiem_tra']}"

        worksheet.merge_cells(f"G{sig_row5}:J{sig_row5}")
        worksheet[f"G{sig_row5}"] = f"Người phê duyệt: {DEFAULT_SIGNATURES['nguoi_phe_duyet']}"

        for col in range(1, 11):
            worksheet.cell(row=sig_row5, column=col).border = thin_border

        # Adjust column widths
        worksheet.column_dimensions["A"].width = 5  # STT
        worksheet.column_dimensions["B"].width = 12  # LID
        worksheet.column_dimensions["C"].width = 15  # Barcode
        worksheet.column_dimensions["D"].width = 14  # Ngày thực hiện
        worksheet.column_dimensions["E"].width = 16  # Ngày trả kết quả
        worksheet.column_dimensions["F"].width = 12  # Chi tiết
        worksheet.column_dimensions["G"].width = 12  # Hành động
        worksheet.column_dimensions["H"].width = 18  # Ngày thực hiện lần 2
        worksheet.column_dimensions["I"].width = 14  # Kết quả lần 2
        worksheet.column_dimensions["J"].width = 14  # Kết quả cuối

        # Set header and footer for printing
        worksheet.oddHeader.left.text = "Phiếu theo dõi soát xét kết quả xét nghiệm"
        worksheet.oddHeader.right.text = "BM01 XN-QTQL.021.04"
        worksheet.oddFooter.left.text = "Ngày hiệu lực: 10/04/2025"
        worksheet.oddFooter.right.text = "Trang &P/&N"

    logger.info(f"Excel report saved to {output_path}")


def arg_parser() -> ArgumentParser:
    """Create argument parser."""
    parser = ArgumentParser(description="Generate HTML tracking report for mtDNA samples")
    parser.add_argument(
        "-b",
        "--batch-id",
        required=True,
        help="Batch ID (e.g., 20251015_mtDNA_261)",
    )
    parser.add_argument(
        "-m",
        "--metadata-excel",
        required=True,
        help="Path to metadata Excel file with LID and Barcode columns",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path prefix (without extension). Will generate .html and .xlsx files",
    )
    parser.add_argument(
        "-r",
        "--rerun-excel",
        required=False,
        default=None,
        help="Path to rerun Excel file with 'Rerun' sheet containing ID, Issues, Requirement columns",
    )
    parser.add_argument(
        "-t",
        "--tracking-tsv",
        required=False,
        default=None,
        help="Path to tracking TSV file with 'Batch mtDNA' and 'Run1 Analysis Results Released Date' columns",
    )
    return parser


def main() -> None:
    """Main entry point."""
    parser = arg_parser()
    args = parser.parse_args()

    # Validate input files
    metadata_excel_path = Path(args.metadata_excel)
    output_path = Path(args.output)

    # Extract date and batch suffix from batch_id
    try:
        analysis_date, batch_suffix = extract_date_and_batch(args.batch_id)
        if analysis_date:
            logger.info(f"Extracted date: {format_date(analysis_date)}, batch: {batch_suffix}")
        else:
            logger.info(f"No date in batch_id, batch: {batch_suffix} (dates will come from tracking TSV)")
    except (ValueError, IndexError) as e:
        logger.error(f"Invalid batch_id format. Expected YYYYMMDD_mtDNA_XX or MS_XXXXXX_XXX: {e}")
        sys.exit(1)

    # Load metadata for barcode lookup
    logger.info(f"Loading metadata from {metadata_excel_path}")
    metadata_mapping = load_metadata_excel(metadata_excel_path)

    # Load rerun data if provided
    rerun_mapping = {}
    if args.rerun_excel:
        rerun_excel_path = Path(args.rerun_excel)
        logger.info(f"Loading rerun data from {rerun_excel_path}")
        rerun_mapping = load_rerun_excel(rerun_excel_path)

    # Load tracking data if provided
    tracking_sequence_date = None
    tracking_release_date = None
    if args.tracking_tsv:
        tracking_tsv_path = Path(args.tracking_tsv)
        logger.info(f"Loading tracking data from {tracking_tsv_path}")
        tracking_sequence_date, tracking_release_date = load_tracking_tsv(tracking_tsv_path, batch_suffix)

    # Generate HTML report and Excel data
    logger.info(f"Generating reports for batch: {batch_suffix}")
    html_content, excel_rows, result_date_str = generate_html_report(
        metadata_mapping=metadata_mapping,
        batch_suffix=batch_suffix,
        analysis_date=analysis_date,
        rerun_mapping=rerun_mapping,
        tracking_sequence_date=tracking_sequence_date,
        tracking_release_date=tracking_release_date,
    )

    # Write HTML output
    html_output_path = Path(str(output_path) + ".html")
    html_output_path.parent.mkdir(parents=True, exist_ok=True)
    with html_output_path.open("w", encoding="utf-8") as f:
        f.write(html_content)
    logger.info(f"HTML report saved to {html_output_path}")

    # Write Excel output
    excel_output_path = Path(str(output_path) + ".xlsx")
    generate_excel_report(excel_rows, excel_output_path, batch_suffix, result_date_str)


if __name__ == "__main__":
    main()
