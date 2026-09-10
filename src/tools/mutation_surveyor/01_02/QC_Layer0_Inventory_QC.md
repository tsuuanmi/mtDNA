# QC Layer 0 — Inventory / Input Integrity QC

## Purpose

QC Layer 0 is the mandatory first gate in `ms_trace_process.py`.

It verifies that the MS custom report and LID manifest belong to the same batch and that every real sample LID has the expected trace inventory before the pipeline produces `{BATCH}_MS_trace_processed.xlsx`.

If Layer 0 fails, the batch is **not eligible for downstream processing**.

## Position

```text
MS custom report
+
LID manifest
        ↓
QC Layer 0: Inventory / Input Integrity QC
        ↓ PASS only
Trace processing and variant normalization
        ↓
{BATCH}_MS_trace_processed.xlsx
```

## Script

```text
ms_trace_process.py
```

## Required inputs

### MS custom report

```text
{BATCH}_anything_df_custom_report.xlsx
```

### LID manifest

```text
YYYYMMDD_Danh sách mẫu_{BATCH}.xlsx
```

The manifest must contain:

```text
LID
MẺ CHẠY
```

## Expected primers

For each real sample LID:

```text
HV1F = 1
HV1R = 1
HV2F = 1
HV3R = 1
```

Primer total counts for real samples must equal the number of manifest LIDs.

## PC/NTC handling

The LID manifest contains only real sample LIDs.

PC/NTC controls may be present in the MS custom report and are excluded from real-sample manifest matching and primer count checks.

Recognized control formats:

```text
NTC1, NTC2, ...
PC1_mtDNA_247547, PC2_mtDNA_247547, ...
```

Control analytical behavior is handled separately by optional `ms_control_QC.py`.

## Blocking errors

Any of the following fails the whole batch:

- manifest filename format error,
- custom report filename format/batch error,
- batch code mismatch between manifest filename and `MẺ CHẠY`,
- blank LID,
- leading/trailing spaces in LID or batch code,
- duplicate LID in manifest,
- missing real sample LID in MS report,
- extra non-control LID in MS report,
- blank/unparsed LID or primer in MS rows,
- unexpected primer,
- missing primer trace,
- duplicate primer trace,
- primer total count mismatch.

Layer 0 is intentionally strict and does not silently repair dirty input values.

## PASS behavior

If Layer 0 passes:

```text
- inventory QC report is written
- trace processing continues
- {BATCH}_MS_trace_processed.xlsx is created
- exit code = 0
```

## FAIL behavior

If Layer 0 fails:

```text
- inventory QC report is written
- {BATCH}_MS_trace_processed.xlsx is not created
- exit code = 1
```

The user should review AB1 files, the Mutation Surveyor project/export, the custom report, and the LID manifest.

## Inventory report

Default report name is derived from the output path:

```text
{BATCH}_MS_trace_processed_inventory_QC.xlsx
```

Typical sheets:

```text
Inventory_Summary
Issue_Counts
Inventory_Details
Primer_Counts
LID_Primer_Matrix
```

## Design principle

Layer 0 answers:

```text
Is this batch input structurally trustworthy enough to create MS trace-processed data?
```
