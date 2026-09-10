# Result Regeneration Validation

> Validate regenerated mtDNA results against expected data.

## Purpose

`src/validation/validation_regenerate_results.py` provides `RegenerateResultsValidator` for validating that regenerated results (from manual or automated pipeline) are consistent with the original data. It checks JSON structure, variant integrity, and manual edit compatibility.

## Pipeline Position

```
Regenerated results (manual or automated)
         │
         ▼
   RegenerateResultsValidator
    ├── Validate batch JSON exists and is valid
    ├── Check all sample files exist
    ├── Verify manual edits were applied correctly
    └── Generate validation report (validation.json)
```

## Public API

### `RegenerateResultsValidator`

```python
RegenerateResultsValidator(
    regenerate_dir: Union[str, Path],    # Regenerated results directory
    comparison_dir: Union[str, Path] = None,  # Comparison results directory
    strict: bool = False,                 # Enable strict validation
)
```

#### Methods

| Method | Description |
|--------|-------------|
| `run_validation() → bool` | Run complete validation workflow; returns True if all checks pass |
| `validate_batch_json() → bool` | Validate that batch JSON exists and has correct structure |
| `validate_all_samples()` | Check that all expected sample files exist |
| `validate_manual_compatibility()` | Verify manual edits were applied correctly |
| `generate_report() → Dict` | Generate comprehensive JSON validation report |

### CLI Entry Point

```bash
python -m src.validation.validation_regenerate_results \
    --regenerate-dir PATH \
    [--comparison-dir PATH] \
    [--strict]
```

## Key Concepts

### Validation Steps

1. **Batch JSON Validation** — Verify `statistic_fullbatch.json` exists and contains valid data
2. **Sample Validation** — Check that every sample in the batch has a corresponding JSON file
3. **Manual Compatibility** — Verify that ADD/REMOVE/EDIT operations from the diff file were correctly applied

### Validation Report

The `validation.json` output contains:

```json
{
  "batch_id": "...",
  "regenerate_dir": "...",
  "overall_status": "PASS|FAIL",
  "validation_status": {
    "batch_json_valid": true,
    "all_samples_valid": true,
    "manual_compatibility_valid": true
  },
  "statistics": {
    "total_samples": 10,
    "valid_samples": 10,
    "invalid_samples": 0,
    "missing_json_files": 0
  },
  "samples": {
    "valid": [...],
    "invalid": [...],
    "missing_json_files": [...]
  },
  "issues": {
    "validation_errors": [],
    "compatibility_issues": []
  }
}
```

### Manual Compatibility Check

The `validate_manual_compatibility()` method verifies that:
- **ADD** edits: New variants from the manual diff exist in regenerated results
- **REMOVE** edits: Removed variants are absent from regenerated results
- **EDIT** edits: Modified variants have the correct sequence

### Strict Mode

When `strict=True`, the validator:
- Requires all samples to be present (no missing JSONs)
- Checks that variant sequences exactly match expected values
- Reports any unexpected variants

## Cross-References

- [pipeline.md](../pipeline.md) — Manual pipeline uses this validator after regeneration
- [core/models.md](../core/models.md) — Sample model used for data loading
- [variants.md](../variants.md) — Position utilities for validation
