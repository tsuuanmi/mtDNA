# Batch Processing

> Orchestrates sequence generation and statistics computation for batches of mtDNA samples.

## Purpose

`src/core/batch.py` provides the `Batch` class that takes `Sample` objects (with variants but no HV sequences), computes consensus sequences, calculates statistics, and writes per-sample JSON output. It is the canonical point where `sample_to_dict()` (and thus `generate_sequence()`) runs.

## Pipeline Position

```
ETL Modules produce Sample (variants + intervals, hv1/hv2/hv3 = None)
         │
         ▼
   Batch.write()
    ├── sample_to_dict() → generate_sequence() → hv1, hv2, hv3
    ├── statistic_variants() → no_snps, no_ins, no_dels, ...
    ├── flag_variants() → per-variant flags
    ├── sample_to_dict() → canonical JSON format
    └── Write per-sample JSON
         │
         ▼
   Downstream consumers (Comparison, FASTA Gen, Reports)
```

## Public API

### `Batch`

```python
Batch(
    samples: list[Sample],   # Samples with variants populated, HV sequences = None
    ref_path: str = "ref/rCRS.fasta",  # Reference FASTA path
)
```

#### Methods

| Method | Description |
|--------|-------------|
| `write(output_dir, batch_id)` | Process all samples and write per-sample JSON files |
| `to_json() → dict[str, Any]` | Convert batch data to JSON-serializable dict (for testing/introspection) |

#### Processing Steps

For each sample in the batch:

1. **Serialization** — Call `sample_to_dict()` for each sample, which internally:
   - Calls `generate_sequence()` to compute HV1, HV2, HV3 consensus sequences
   - Calls `statistic_variants()` to count SNPs, insertions, deletions, identicals, unread
   - Calls `flag_variants()` for per-variant quality flags
   - Builds the canonical dict with variants grouped by type, sample/variant flags, and information
2. **Write** — Save `statistic_fullbatch.json` and per-sample JSONs to `{output_dir}/{batch_id}/`

### Error Handling

If `generate_sequence()` or `statistic_variants()` fails for a sample:
- A `computation_error` field is included in the output
- Processing continues for remaining samples
- The sample still appears in output but with incomplete data

## Key Concepts

### HV Sequence Contract

The Batch class is the **sole compute point** for HV sequences. ETL modules must NOT call `generate_sequence()` directly — they produce `Sample` objects with `hv1`/`hv2`/`hv3 = None`, and `sample_to_dict()` (called by Batch) populates them.

### Canonical Output Format

Each sample JSON follows the format produced by `sample_to_dict()`:

```json
{
  "no_snps": 5,
  "no_ins": 1,
  "no_dels": 0,
  "no_identicals": 3,
  "no_unread": 0,
  "HV1": "ACGT...",
  "HV2": "ACGT...",
  "HV3": "ACGT...",
  "variants": {"snps": [...], "insertions": [...], "deletions": [...]},
  "intervals": {"HV1": [[16024, 16365]], ...},
  "sample_flags": [...],
  "variant_flags": {"73|A|G": ["Low quality"], ...},
  "information": {"batch_id": "...", "source_tool": "..."}
}
```

## Configuration

- `src/core/sample.py` — `sample_to_dict()` for serialization (calls `generate_sequence()`, `statistic_variants()`, and `flag_variants()` internally)
- `src/core/sample.py` — `generate_sequence()` for HV consensus sequence generation
- `src/core/variants.py` — `statistic_variants()`, `split_variant_types()`, `variant_to_dict()`
- `src/core/flagging.py` — `flag_variants()` for per-variant flags
- Reference genome: `ref/rCRS.fasta`

## Cross-References

- [ARCHITECTURE.md](../ARCHITECTURE.md) §1.1 — Batch step in the standard pipeline
- [core/models.md](./models.md) — Sample and Variant model definitions
- [core/sample.md](./sample.md) — `sample_to_dict()` serialization
- [variants.md](variants.md) — `statistic_variants()`, `split_variant_types()`, position utilities
- [core/sample.md](./sample.md) — `generate_sequence()` and `sample_to_dict()`
- [core/flagging.md](./flagging.md) — Per-variant flagging
