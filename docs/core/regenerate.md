# Regenerate Gate

> Module: `src/core/regenerate.py`

Produces filtered final JSON for downstream FASTA and report generation.

## Purpose

Tool JSON preserves every called variant, including variants outside sequenced intervals, so comparison can show reviewer-facing additional calls. The regenerate gate creates the downstream authoritative final JSON by filtering each sample to its declared analyzed intervals.

```text
results/tools/<tool>/<batch>/json/   # raw tool JSON, used by comparison
        ↓ regenerate
regenerate/<batch>/                  # filtered final JSON, used by FASTA/report generation
```

Samples with no `intervals` are treated as full-coverage inputs and pass through unchanged.

## Public API

### `filter_to_intervals(sample: Sample) -> Sample`

Returns a sample with variants scoped to `sample.intervals`.

Behavior:

- If `sample.intervals` is missing or empty, returns the sample unchanged.
- Keeps only configured HV region keys present in the sample (`HV1`, `HV2`, `HV3`).
- Uses `filter_sample_by_regions()` for interval-based filtering.
- Falls back to the original sample if filtering produces `None`.

### `regenerate(input_dir, output_dir, ref_path, batch_id) -> int`

Loads raw per-sample JSON files, filters each sample, and writes standard batch output.

| Argument | Description |
|---|---|
| `input_dir` | Raw tool JSON directory containing `{LID}/{LID}.json` files |
| `output_dir` | Destination regenerate/final directory |
| `ref_path` | Reference FASTA path used for sequence recomputation |
| `batch_id` | Informational batch ID passed to `Batch.write()` |

Returns the number of samples written. Returns `0` if the input directory is missing or no samples are loaded.

## Input Layout

`regenerate()` reads per-sample files from the standard tool JSON layout:

```text
<input_dir>/
├── statistic_fullbatch.json   # ignored by regenerate loader
├── SAMPLE_1/SAMPLE_1.json
└── SAMPLE_2/SAMPLE_2.json
```

The per-sample files are the source of truth.

## Output Layout

The output is written flat because `Batch.write(..., nest_batch_id=False)` is used:

```text
<output_dir>/
├── statistic_fullbatch.json
├── SAMPLE_1/SAMPLE_1.json
└── SAMPLE_2/SAMPLE_2.json
```

## CLI

```bash
python -m src.core.regenerate \
  --input results/tools/blastn/<batch>/json \
  --output regenerate/<batch> \
  --reference ref/rCRS.fasta \
  --batch-id <batch>
```

Arguments:

| Flag | Required | Description |
|---|---:|---|
| `-i`, `--input` | Yes | Raw tool JSON input directory |
| `-o`, `--output` | Yes | Regenerate/final output directory |
| `-r`, `--reference` | Yes | Reference FASTA path |
| `--batch-id` | No | Batch identifier, defaults to `unknown` |

## Cross-References

- [sample.md](sample.md) — `filter_sample_by_regions()` and per-sample loading
- [batch.md](batch.md) — final JSON writing
- ../pipeline.md — pipeline orchestration
