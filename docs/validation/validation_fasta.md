# FASTA Validation

> Validate generated FASTA files against the canonical code path.

## Purpose

`src/validation/validation_fasta.py` provides `FastaValidator` for validating that generated FASTA files match the expected output from `generate_sequence()` (via the standard `Sample`). It loads each per-sample JSON as a `Sample` (`src.core.sample.load_sample`), replays the canonical sequence generation path, and compares results position-by-position.

## Pipeline Position

```
Generated FASTA files
         │
         ▼
   FastaValidator.validate_all()
    ├── For each sample:
    │   ├── Load per-sample JSON → standard Sample (load_sample)
    │   ├── Replay generate_sequence()
    │   └── Compare generated vs expected sequences
    └── Report alignment differences
```

## Public API

### Data Classes

| Class | Description |
|-------|-------------|
| `AlignmentData` | Position-level alignment data for a single region (ref_aligned, alt_aligned, diff_line, variants_applied) |
| `RegionResult` | Validation result for one HV region (passed/failed, expected vs actual sequences, alignment diff) |
| `SampleResult` | Validation result for one sample (aggregate of region results) |

### `FastaValidator`

```python
FastaValidator(
    fasta_dir,    # Directory with generated FASTA files
    json_dir,     # Directory with per-sample JSON data
    ref_path,     # Reference FASTA path
    batch_id,     # Batch identifier
)
```

#### Methods

| Method | Description |
|--------|-------------|
| `validate_all() → bool` | Validate all samples; returns True if all pass |
| `validate_sample(sample_id) → SampleResult` | Validate a single sample |
| `write_report(output_path)` | Write validation report to file |

### CLI Entry Point

```bash
python -m src.validation.validation_fasta \
    --fasta-dir PATH \     # FASTA files directory
    --json-dir PATH \      # JSON data directory
    --ref-path PATH \      # Reference FASTA
    --batch-id ID          # Batch identifier
```

## Key Concepts

### Validation Strategy

FastaValidator does **not** compare FASTA files byte-by-byte. Instead, it:

1. Loads the same per-sample JSON used by `generate_fasta.py` as a standard `Sample`
2. Replays `generate_sequence()` to regenerate expected sequences
3. Compares expected vs actual sequences position-by-position
4. Produces alignment diffs showing where they diverge

This catches bugs in sequence generation, variant application, or data corruption.

### Alignment Report

For each region where sequences differ:
- Shows reference alignment with position markers
- Highlights positions where expected ≠ actual
- Lists variants that were applied
- Provides diff markers (mismatches, insertions, deletions)

### `_merge_con_dict()`

Internal helper that merges reference bases (int keys) with variant bases (str keys) to produce a unified position-to-base mapping. Variant bases override reference bases at the same position.

## Cross-References

- [generation.md](../generation.md) — FASTA generation that this validator checks
- [core/sample.md](../core/sample.md) — `generate_sequence()` — the canonical path being validated
- [core/batch.md](../core/batch.md) — Batch writes the JSON data used by validation
