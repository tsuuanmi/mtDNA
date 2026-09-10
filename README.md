# mtDNA Analysis Pipeline

A Python pipeline for processing and analyzing mitochondrial DNA (mtDNA) sequencing data. Automates steps from raw AB1 file preparation through BLASTN analysis, variant calling with Tracy, automated and manual pipelines, to PDF report generation.

## Web UI

Chọn quy trình trên trình duyệt thay vì sửa tay `.env` trước mỗi lần chạy:

```bash
make ui                  # http://localhost:8765, mở trong LAN
make ui-test             # tự kiểm giao diện + API
```

Chi tiết: [web/README.md](web/README.md).

## Quick Start

```bash
# 1. Install dependencies
bash scripts/setup.sh

# 2. Configure environment
cp .env.example .env
nano .env

# 3. Run pipeline (BLASTN mode — default)
bash scripts/pipeline.sh MS_191125_004

# 3. Run pipeline (Tracy mode)
bash scripts/pipeline.sh -p tracy MS_191125_004
```

## Requirements

- **Python**: 3.12+
- **Package Manager**: uv
- **External Tools**: UGENE (with BLAST+), Tracy, seqtk

See `AGENTS.md` for complete environment setup.

---

## Usage

### Run Full Pipeline

```bash
# BLASTN mode (default)
bash scripts/pipeline.sh MS_191125_004

# Tracy mode
bash scripts/pipeline.sh -p tracy MS_191125_004
```

### Run Specific Steps

Steps: `1=validate, 2=prepare, 3=automate, 4=manual, 5=fasta, 6=report`

Use numbers or names:

```bash
bash scripts/pipeline.sh -s 2,3,4 MS_191125_004
bash scripts/pipeline.sh -s prepare,automate,manual MS_191125_004
```

### Custom Directories

```bash
bash scripts/pipeline.sh -d /custom/data -m /custom/manual -o /custom/results MS_191125_004
```

### Rerun Specific Samples

```bash
bash scripts/pipeline.sh -l /path/to/samples.txt -s 3,4 MS_191125_004
```

### Docker

```bash
# Build
docker build -t mtdna:v1 .

# Run
docker run -it --rm \
  -v /path/to/data:/app/data \
  -v /path/to/results:/app/results \
  -v $(pwd)/.env:/app/.env \
  mtdna:v1
```

---

## Pipeline Workflow

```
raw AB1 files
    ↓
[Step 1: Validate] — Positive Control validation (optional)
    ↓
[Step 2: Prepare] — Organize AB1 files + extract sample list
    ↓
[Step 3: Automate] — Analysis engine (varies by --pipeline mode)
    ├→ blastn: BLASTN → per-sample Tracy decompose → automate pipeline
    └→ tracy:  Integrated Tracy analysis (src/tools/tracy/pipeline.py)
    ↓
[Step 4: Manual] — Manual verification + Sequencher comparison
    ├→ Sequencher ETL (standard Sample JSON w/ HV1-HV3 sequences)
    ├→ regenerate gate: filter variants to sequenced intervals -> final JSON
    └→ compare automated vs manual (reads RAW tool JSON, so the reviewer still
       sees out-of-range/additional variants; the filtered final JSON is used
       only for downstream FASTA/report)
    ↓
[Step 5: FASTA] — Generate FASTA files (from the filtered regenerate/final JSON)
    ↓
[Step 6: Report] — Generate PDF reports (optional, from the filtered final JSON)
```

---

## Testing

```bash
pytest tests/ -v
```

---

## Key Files

| File | Purpose |
|------|---------|
| `scripts/pipeline.sh` | Unified pipeline script (BLASTN + Tracy modes) |
| `src/config.py` | Centralized configuration |
| `src/core/models.py` | Canonical Pydantic data models (Sample, Variant) |
| `src/core/batch.py` | Sequence generation and statistics for batches |
| `src/core/regenerate.py` | Regenerate gate: filter raw tool JSON to sequenced intervals (final JSON for FASTA/report) |
| `src/core/comparison.py` | Pairwise sample comparison (reads raw tool JSON) |
| `src/tools/blastn/pipeline.py` | BLASTN analysis engine |
| `src/tools/tracy/pipeline.py` | Tracy analysis engine |
| `src/tools/sequencher/pipeline.py` | Sequencher ETL + manual comparison |
| `src/tools/mutation_surveyor/` | Mutation Surveyor ETL and step scripts |
| `src/modules/TNLS/verification/` | TNLS verification helpers |
| `src/generate_fasta.py` | FASTA generation |
| `src/generate_reports.py` | PDF report generation |

---

## Documentation

- **AGENTS.md** — Developer/instructions for AI agents
- **docs/** — Additional documentation

---

## License

See project files for details.
