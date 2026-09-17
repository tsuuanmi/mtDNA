# TNLS/HCLS verification

The TNLS verification module compares known HCLS base samples with TNLS target
samples over mtDNA HV1, HV2, and HV3. It reports only pairs classified as
`CANNOT_EXCLUDE`; this means the available overlapping evidence does not
exclude a common source, not that the samples are proven to be from the same
individual.

For the matching rules, IUPAC handling, result statuses, and family summaries,
read the [matching guide](matching-guide.md).

## Run the module

From the repository root:

```bash
uv run python -m src.modules.TNLS.verification \
  --base-samples path/to/HCLS.json \
  --target-samples path/to/TNLS_merged.json \
  --matching-rule src/modules/TNLS/verification/resources/matching_rule.json \
  --results-dir results/modules/TNLS/verification
```

Optional arguments:

```text
--filter-list path/to/filter.txt       Restrict TNLS samples to selected IDs.
--tnls-metadata path/to/metadata.xlsx  Add identity fields to statistics outputs.
--max-ratio N                          Include 1:1 through 1:N pairs (default: 5).
```

## Inputs

| Input | Purpose |
|---|---|
| HCLS JSON | Base samples, keyed by sample ID, with intervals and variants |
| TNLS merged JSON | Target samples, keyed by sample ID, with intervals and variants |
| `matching_rule.json` | Minimum overlap, selected HV regions, and mismatch thresholds |
| Optional filter list | One TNLS sample ID per line |
| Optional metadata workbook | TNLS identity metadata keyed by `mtDNA LID` |

The default matching rule requires at least 150 overlapping base pairs across
HV1, HV2, and HV3. It permits zero mismatches only, and requires at least one
matched variant before a pair is emitted.

## Outputs

The selected results directory receives:

| File | Contents |
|---|---|
| `output_compare_tnls_hcls.tsv` | Per-pair `CANNOT_EXCLUDE` results |
| `stats_TNLS_HCLS_with_metadata.tsv` | Expanded 1:1 through 1:N matches, with optional metadata |
| `stats_TNLS_HCLS_clean.tsv` | Simplified matched-pair view |
| `stats_TNLS_HCLS_family_summary.tsv` | TNLS profile families and their matched HCLS IDs |

The command creates the results directory when necessary. It does not modify
input JSON files. Generated outputs should be regenerated from source logic or
input data rather than edited manually.

## Python API

```python
from src.modules.TNLS.verification import CompareManager, compare_1n, group_families

matches = compare_1n(base_sample, target_samples, matching_rule)
families = group_families(target_samples, matching_rule["region"], matching_rule["minimum"])
```

`compare_1n` returns a pandas DataFrame containing `CANNOT_EXCLUDE` pairs only.
