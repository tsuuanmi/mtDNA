# Output Generation

> FASTA file generation, report metadata creation, and PDF report writing.

## Purpose

Three modules handle the final output stages of the pipeline:

- **`src/generate_fasta.py`** — Generate consensus FASTA files from Sample objects
- **`src/generate_reports.py`** — Generate per-sample JSON reports and trigger PDF generation

## Pipeline Position

```
Batch (sequence generation)
         │
         ▼
generate_fasta.py ──► Consensus FASTA files (HV region sequences)
         │
         ▼
generate_reports.py ──► Per-sample JSON reports + PDF generation
```

---

## FASTA Generation (`src/generate_fasta.py`)

### Purpose

Generates consensus FASTA files from each tool's standard Sample JSON output. It loads each per-sample JSON as a `Sample` (via `src.core.sample.load_sample`) and regenerates HV1/HV2/HV3 consensus sequences via `generate_sequence()` (canonical flow). Includes post-generation validation via `FastaValidator`.

### `FastaGenerator`

```python
FastaGenerator(
    input_dir: str,    # Directory with per-sample JSON files
    output_dir: str,   # Output directory for FASTA files
    ref_path: str,     # Reference FASTA path
    force: bool = False,  # If True, wipe output dir first (clean-slate mode)
)
```

#### Methods

| Method | Description |
|--------|-------------|
| `write_fasta_file(output_path, sample_name, region_sequences)` | Write sequences to FASTA format file |
| `generate_reviewed_fasta() → Tuple[int, int, int]` | Generate FASTA for all samples. Returns `(success_count, failure_count, skip_count)`. Default: skip existing `.fasta` files |

#### Behavior

- **Default (force=False)**: Skip samples that already have a `.fasta` file. Existing files are never modified.
- **Force mode (force=True)**: Wipe the output directory first, then regenerate all samples.
- After generation, runs `FastaValidator` to validate output by replaying `generate_sequence()`.

### CLI Entry Point

```bash
python -m src.generate_fasta \
    -i <input_dir> \
    -o <output_dir> \
    -r <reference_path> \
    [--force]
```

---

## Report Generation (`src/generate_reports.py`)

### Purpose

Orchestrates per-sample report generation including barcode extraction, PDF download, and batch statistics. Uses async I/O for PDF report downloads.

### `ReportGenerator`

```python
ReportGenerator()  # No constructor arguments
```

#### Methods

| Method | Description |
|--------|-------------|
| `extract_sample_info(file_path, output_path, barcode_len=None) → List[Tuple[str, str]]` | Extract LID and barcode from Excel or TSV file. Filters out PC/NTC control samples |
| `get_pdf_link(session, lid, barcode, json_path, auth_token) → str` | Get PDF report API link for a sample (async) |
| `download_pdf(session, pdf_link, output_path) → None` | Download PDF report file (async) |
| `get_existing_pdfs(report_dir) → Set[str]` | Get set of already-downloaded PDF filenames |
| `process_sample(session, barcode, auth_token, output_dir) → bool` | Process single sample: get link, download PDF (async) |
| `find_data_file(batch_id, metadata_dir) → Optional[Path]` | Find batch data file in metadata directory |
| `generate_reports(batch_id, data_dir, results_dir, base_dir, report_dir, auth_token) → None` | Full batch report generation pipeline (async). Extracts sample info, downloads PDFs, and generates per-sample reports |

### Report JSON Structure

Each sample report JSON contains:
- Variant data (SNPs, insertions, deletions)
- HV sequences
- Variant statistics (counts)
- Sample metadata (date, time, sample type, location)
- Interval information

---

## Configuration

- `src/config.py` — `REPORT_BASE_URL`, `REPORT_RETRIES`, `REPORT_RETRY_DELAY`
- `src/core/sample.py` — `generate_sequence()` for HV consensus generation
- `src/core/variants.py` — `statistic_variants`, `split_variant_types`
- `src/core/sample.py` — `load_sample()` (per-sample JSON → `Sample`)
- Reference genome: `ref/rCRS.fasta`

## Cross-References

- [blastn.md](tools/blastn.md) — BLASTN analysis engine producing Sample JSON consumed by FASTA generation
- [pipeline.md](pipeline.md) — Uses report generation in automate/manual pipelines
- [core/models.md](core/models.md) — Sample model used by `generate_fasta`
- [core/batch.md](core/batch.md) — Batch populates HV sequences via `generate_sequence()` before FASTA generation
