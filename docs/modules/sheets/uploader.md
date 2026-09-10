# Google Sheets Uploader

> Upload batch data to Google Sheets for sharing and review.

## Purpose

`src/modules/sheets/uploader.py` provides authentication and upload functionality for pushing mtDNA batch data (from Excel or TSV files) to Google Sheets.

## Public API

### Authentication

| Function | Signature | Description |
|----------|-----------|-------------|
| `authenticate_google_sheets` | `(credentials_file: str) → gspread.Client` | Authenticate with Google Sheets API using service account credentials |

### Data Loading

| Function | Signature | Description |
|----------|-----------|-------------|
| `load_data` | `(file_path: str) → pd.DataFrame` | Load data from `.xlsx` or `.tsv` file, handling NaN/inf values |

### Upload

| Function | Signature | Description |
|----------|-----------|-------------|
| `upload_to_google_sheets` | `(client, sheet_id, worksheet_name, data_df, replace=False) → bool` | Upload DataFrame to Google Sheets. If `replace=True`, clears existing data; otherwise creates a timestamped worksheet |
| `clean_data_for_sheets` | `(values: list) → list` | Clean data for JSON compatibility (replace NaN/inf with empty strings) |

### CLI Entry Point

```bash
python -m src.modules.sheets.uploader \
    --credentials PATH \      # Service account JSON credentials
    --sheet-id ID \           # Google Sheets ID
    --worksheet NAME \        # Worksheet name
    --input PATH \            # Input data file (.xlsx or .tsv)
    [--replace]               # Replace existing worksheet
```

## Key Concepts

### Authentication

Uses Google Service Account credentials (JSON key file) for authentication. The service account must have write access to the target Google Sheet.

### Data Cleaning

- `NaN`, `inf`, and `-inf` values are replaced with empty strings
- All values are converted to Google Sheets-compatible types
- Batch upload is attempted first; cell-by-cell fallback on failure

### Upload Modes

- **Append (default)**: Creates a new worksheet with timestamped name
- **Replace**: Clears existing worksheet and writes new data

## Configuration

- `credentials/` — Google service account JSON key file
- Google Sheets API must be enabled for the service account

## Cross-References

- [generation.md](../../generation.md) — Report generation produces data that may be uploaded
- [pipeline.md](../../pipeline.md) — Pipeline produces batch statistics for upload
