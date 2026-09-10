# QC Layer 1 — PC / NTC Control QC

## Purpose

QC Layer 1 is an optional control QC step implemented by `ms_control_QC.py`.

It checks whether batch controls behave as expected after `ms_trace_process.py` has generated `{BATCH}_MS_trace_processed.xlsx`.

It answers:

```text
1. Is the positive control profile correct?
2. Are the negative controls clean?
```

If enabled and failed, the batch is not eligible for auto-pass.

## Position

```text
ms_trace_process.py
        ↓
{BATCH}_MS_trace_processed.xlsx
        ↓ optional
ms_control_QC.py
        ↓
{BATCH}_control_QC.xlsx
```

This step is intended to be controlled by `pipeline.sh`. A control reference file can exist without forcing this step to run.

## Command line usage

```bash
python3 ms_control_QC.py \
  --trace-processed "MS_060526_002_MS_trace_processed.xlsx" \
  --control-ref "mtDNA_PC_reference.xlsx" \
  --out "MS_060526_002_control_QC.xlsx"
```

A deprecated alias `--compare-qc` is accepted for backward compatibility, but new pipeline commands should use `--trace-processed`.

## Required input columns

The trace-processed file must contain:

```text
LID
Primer
Variant Converted 2
```

Expected primers:

```text
HV1F
HV1R
HV2F
HV3R
```

## Control LID formats

### NTC

```text
NTC1
NTC2
NTC3
...
```

### PC

```text
PC1_mtDNA_247547
PC2_mtDNA_247547
...
```

For PC controls, the part after `PC#_` is treated as the reference LID.

Example:

```text
PC1_mtDNA_247547
→ reference LID = mtDNA_247547
```

## PC reference file

The PC reference file is global, not batch-specific.

Minimum structure:

```text
Reference_LID | Expected_Profile
mtDNA_247547  | 73G 263G 309.1C 315.1C ...
```

Accepted key columns:

```text
Reference_LID
LID
Control_LID
Control_ID
PC_LID
```

Accepted profile columns:

```text
Expected_Profile
Profile
Expected_Variants
Variants
```

Expected profile tokens should use the same notation as `Variant Converted 2`.

## Analysis ranges

Control profiles are built only from variants within:

```text
73-340
438-576
16024-16365
```

Variants outside these ranges are ignored for control QC.

## Observed control profile construction

For each control LID, the script:

```text
collects Variant Converted 2 from HV1F + HV1R + HV2F + HV3R
filters tokens to analysis ranges
creates a sorted unique union
```

## PC QC logic

For each PC:

```text
Observed PC profile == expected reference profile
```

Any missing or extra variant fails the PC.

PC replicates using the same reference LID must also be concordant with each other.

## NTC QC logic

For each NTC:

```text
Observed NTC profile must be empty
```

Any in-range variant fails the NTC.

## Control trace inventory

Each PC/NTC control must have exactly one trace for each expected primer:

```text
HV1F
HV1R
HV2F
HV3R
```

Missing, duplicate, or unexpected control primer traces fail the control.

## Optional behavior

Flags:

```bash
--allow-missing-pc
--allow-missing-ntc
```

These can be used for legacy batches where controls are absent or incomplete.

## Output

```text
{BATCH}_control_QC.xlsx
```

Sheets:

```text
Summary
Details
Observed_Profiles
Reference_Profiles
```

## Exit code

```text
PASS → exit code 0
FAIL → exit code 1
```

## Design principle

Layer 1 answers:

```text
Given the standardized trace-level MS variants, do the controls indicate that this batch is analytically trustworthy?
```
