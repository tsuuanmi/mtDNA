# 01 — `ms_trace_process.py`

## Purpose

`ms_trace_process.py` is the first script in the MS mtDNA post-processing pipeline.

It reads a Mutation Surveyor custom report, validates batch/sample/trace inventory with QC Layer 0, normalizes MS trace-level variant calls, and writes the trace-level processed output used by downstream HV23/HV1 merge scripts.

This script does **not** use a truthset and does **not** calculate concordance.

## Pipeline position

```text
MS custom report
+
LID manifest
        ↓
01 ms_trace_process.py
        ↓
{BATCH}_MS_trace_processed.xlsx
        ↓
optional: ms_control_QC.py
        ↓
HV23 / HV1 merge scripts
```

## Inputs

```bash
--ms             {BATCH}_anything_df_custom_report.xlsx
--lid-manifest   YYYYMMDD_Danh sách mẫu_{BATCH}.xlsx
--out            {BATCH}_MS_trace_processed.xlsx
--inventory-out  optional explicit Layer 0 report path
```

Example:

```bash
python3 ms_trace_process.py \
  --ms "MS_060526_002_rerun_df_custom_report.xlsx" \
  --lid-manifest "20260522_Danh sách mẫu_MS_060526_002.xlsx" \
  --out "MS_060526_002_MS_trace_processed.xlsx"
```

## Required file naming

### MS custom report

```text
{BATCH}_anything_df_custom_report.xlsx
```

Valid examples:

```text
MS_060526_002_df_custom_report.xlsx
MS_060526_002_rerun_df_custom_report.xlsx
MS_060526_002_v2_df_custom_report.xlsx
```

### LID manifest

```text
YYYYMMDD_Danh sách mẫu_{BATCH}.xlsx
```

The manifest must contain at least:

```text
LID
MẺ CHẠY
```

## QC Layer 0

Layer 0 is a strict batch-level inventory/input integrity gate.

It checks:

- manifest filename format,
- MS custom report filename format,
- batch code consistency,
- blank/dirty/duplicated LIDs in the manifest,
- missing real sample LIDs,
- extra non-control LIDs,
- expected primer inventory per sample,
- unexpected primers,
- primer count mismatch.

Expected production primers are:

```text
HV1F
HV1R
HV2F
HV3R
```

Each real sample LID must have exactly one trace for each expected primer.

PC/NTC controls are allowed in the custom report but are excluded from manifest sample matching. Recognized control formats are:

```text
NTC1, NTC2, ...
PC1_mtDNA_247547, PC2_mtDNA_247547, ...
```

If Layer 0 fails:

```text
- inventory QC report is written
- {BATCH}_MS_trace_processed.xlsx is not created
- exit code = 1
```

## Trace processing logic

After Layer 0 passes, the script performs the original MS trace-level processing logic:

```text
Read all MS header blocks
→ keep real .ab1 rows with Read Start / Read End
→ build Variant Summary
→ build Variant Converted
→ build Variant Converted 2
→ derive Well / LID / Primer / MS Range
→ add preliminary repeat-region flags
→ write {BATCH}_MS_trace_processed.xlsx
```

## Variant columns

### `Variant Summary`

Combines `Variant1`, `Variant2`, ... into one cleaned string.

### `Variant Converted`

Converts raw Mutation Surveyor variant notation into internal variant tokens.

Examples:

```text
73A>G  → 73G
146T>C → 146C
```

### `Variant Converted 2`

Applies repeat-region normalization used by downstream merge scripts, including:

- 303-309 C insertion shifting to 309,
- 311-315 C insertion shifting to 315,
- 452-455 T insertion shifting to 455,
- 568-573 C insertion shifting to 573,
- 248DEL to 249DEL normalization,
- AC repeat 513-525 canonicalization.

## Output columns

```text
Well
LID
Primer
Trace #
Sample Name
# of Mutation
Read Start
Read End
MS Range
Variant Summary
Variant Converted
Variant Converted 2
C_shift_rule_applied
AC_repeat_rule_applied
HV1_PolyC_flag
HV2-3_PolyC_flag
```

## Output

```text
{BATCH}_MS_trace_processed.xlsx
```

This file is the canonical trace-level processed MS output for downstream scripts.

## Design principle

`ms_trace_process.py` handles input integrity and MS trace-level variant normalization only.

Truthset comparison is intentionally deferred to the final profile-vs-truth step.
