# Pipeline Orchestration

> End-to-end mtDNA analysis pipeline driven by `scripts/pipeline.sh`.

## Purpose

The pipeline is orchestrated by a single shell script, `scripts/pipeline.sh`, which
runs six named steps. Each step delegates to a Python module under `src/`:

| Step | Name | Module | Description |
|------|------|--------|-------------|
| 1 | validate | `src/modules/quality_control/pc_ntc.py` | Positive Control / NTC validation |
| 2 | prepare | (shell + pandas) | Organize AB1 files and extract the sample list from metadata Excel |
| 3 | automate | `src/tools/blastn/pipeline.py` or `src/tools/tracy/pipeline.py` | Run the analysis engine (`--pipeline` mode) |
| 4 | manual | `src/tools/sequencher/pipeline.py` + `src/core/comparison.py` | Sequencher ETL and pairwise comparison vs automated results |
| 5 | fasta | `src/generate_fasta.py` | Generate consensus FASTA files |
| 6 | report | `src/generate_reports.py` | Generate per-sample JSON reports and download PDFs |

The batch runner `scripts/batch_pipeline.sh` wraps `scripts/pipeline.sh` for
multi-batch and rerun workflows.

## Pipeline Position

```
Step 1: validate   — src.modules.quality_control.pc_ntc (PCNTCProcessor)
Step 2: prepare    — shell: copy AB1 files, extract sample list (pandas)
Step 3: automate   — src.tools.{blastn,tracy}.pipeline (process_batch → Batch.write)
    ├── results mirrored into automate regenerate dir
    └── produces statistic_fullbatch.json + per-region JSONs
Step 4: manual     — src.tools.sequencher.pipeline (process_batch)
    ├── results mirrored into manual regenerate dir
    └── src.core.comparison (SampleComparator.compare_batches)
Step 5: fasta      — src.generate_fasta (FastaGenerator)
Step 6: report     — src.generate_reports (ReportGenerator)
```

Two analysis engines are supported for Step 3, selected with `-p`/`--pipeline`:

- **tracy** (default) — `python -m src.tools.tracy.pipeline`
- **blastn** — `python -m src.tools.blastn.pipeline`

Both engines produce the same canonical output: per-sample `Sample` JSON and a
`statistic_fullbatch.json` batch summary. The shell script mirrors each tool's
JSON output into a stable `regenerate` directory so downstream steps (comparison,
FASTA, reports) read a consistent path regardless of the engine.

---

## Step 1 — Validate (`src/modules/quality_control/pc_ntc.py`)

### `PCNTCProcessor`

Validates Positive Control (PC) and No-Template Control (NTC) samples against
expected Sequencher TXT reference profiles. Validation output is written under
`<results-root>/validation/<batch_id>/`. A failed validation logs a warning but
does not abort the main pipeline.

### CLI Entry Point

```bash
python -m src.modules.quality_control.pc_ntc \
    --batch-id BATCH_ID \
    --data-dir DATA_DIR \
    --results-dir RESULTS_DIR \
    --ref-path ref/rCRS.fasta \
    [--pc-reference-dir DIR]
```

---

## Step 2 — Prepare (shell)

Organizes raw AB1 files for the batch into `data/raw/<batch_id>/`, removes control
files (PC/NTC/PCR), copies the metadata Excel file, and extracts the sample ID
list (LID column, excluding controls) into `<raw_dir>.txt`. When `--sample-list`
is provided, the Excel extraction is skipped. No dedicated Python module — the
logic lives inline in `scripts/pipeline.sh` (using `pandas` for Excel parsing).

---

## Step 3 — Automate (`src/tools/<engine>/pipeline.py`)

Runs the selected analysis engine over the prepared AB1 files. Each engine
module exposes `process_batch(...)` and a `main()` CLI entry point, uses
`ProcessPoolExecutor` for parallel per-sample processing, and writes output via
`Batch.write()` plus per-region intermediate JSONs.

### BLASTN engine

```bash
python -m src.tools.blastn.pipeline \
    --samples data/raw/<batch_id>.txt \
    --batch-id BATCH_ID \
    --input-dir data/raw/<batch_id> \
    --output-dir results/tools/blastn/<batch> \
    --ref-path ref/rCRS.fasta
```

### Tracy engine

```bash
python -m src.tools.tracy.pipeline \
    --samples data/raw/<batch_id>.txt \
    --batch-id BATCH_ID \
    --input-dir data/raw/<batch_id> \
    --output-dir results/tools/tracy/<batch_id> \
    --ref-path ref/rCRS.fasta
```

See [tools/blastn.md](tools/blastn.md) and [tools/tracy.md](tools/tracy.md) for
engine internals.

---

## Step 4 — Manual (`src/tools/sequencher/pipeline.py` + `src/core/comparison.py`)

Unzips Sequencher TXT archives for the batch, runs the Sequencher ETL pipeline,
mirrors its JSON output into the manual regenerate directory, then runs the
tool-agnostic pairwise comparison between automated and manual batch JSONs.

### Sequencher ETL

```bash
python -m src.tools.sequencher.pipeline \
    --batch-id BATCH_ID \
    --input-dir MANUAL_DIR \
    --output-dir results/tools/sequencher/<batch_id> \
    --ref-path ref/rCRS.fasta
```

### Pairwise comparison

```bash
python -m src.core.comparison \
    -a automate_regenerate/statistic_fullbatch.json \
    -b manual_regenerate/statistic_fullbatch.json \
    -o results/modules/comparison/<batch_id> \
    --batch-id BATCH_ID
```

`SampleComparator.compare_batches()` produces a `ComparisonResult` exported as
JSON, summary TSV, and variant-level TSV. See
[core/comparison.md](core/comparison.md).

---

## Step 5 — FASTA (`src/generate_fasta.py`)

### `FastaGenerator`

Generates consensus FASTA files from per-sample JSON (loaded as `Sample` via
`src.core.sample.load_sample`) and regenerates HV1/HV2/HV3 consensus sequences
via `generate_sequence()`. Supports `--force` for clean-slate regeneration.
Output is validated with `FastaValidator`.

### CLI Entry Point

```bash
python -m src.generate_fasta \
    -i manual_regenerate_dir \
    -o results/fasta/<batch_id> \
    -r ref/rCRS.fasta \
    [--force]
```

See [generation.md](generation.md).

---

## Step 6 — Report (`src/generate_reports.py`)

### `ReportGenerator`

Extracts sample info (LID + barcode) from the batch metadata Excel/TSV, obtains
an auth token from the GeneStory API, and downloads per-sample PDF reports with
async I/O. Also writes per-sample JSON reports.

### CLI Entry Point

```bash
python -m src.generate_reports \
    -b BATCH_ID \
    -d DATA_DIR \
    -r manual_regenerate_dir \
    --report_dir results/reports/<batch_id>
```

See [generation.md](generation.md).

## Configuration

- `src/config.py` — Region definitions, tool paths, parameters (`get_settings()`)
- `src/core/sample.py` — `generate_sequence()`, `sample_to_dict()`
- `src/core/variants.py` — Position utilities, `statistic_variants()`
- Reference genome: `ref/rCRS.fasta`

## Cross-References

- [ARCHITECTURE.md](ARCHITECTURE.md) — System overview and standard tool pipeline
- [tools/blastn.md](tools/blastn.md) — BLASTN engine
- [tools/tracy.md](tools/tracy.md) — Tracy engine
- [tools/sequencher.md](tools/sequencher.md) — Sequencher ETL
- [core/comparison.md](core/comparison.md) — Pairwise sample comparison
- [generation.md](generation.md) — FASTA and report generation
- [validation/validation_regenerate_results.md](validation/validation_regenerate_results.md) — Result validation
