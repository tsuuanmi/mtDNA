# 05 — `ms_final_profiles_merge.py`

## Purpose

`ms_final_profiles_merge.py` combines the already-merged HV2-3 and HV1 outputs into one final MS mtDNA profile per `Batch + LID`.

This step does **not** use a truthset and does **not** calculate concordance. Truthset comparison is performed later by `ms_final_profiles_vs_truth.py`.

---

## Pipeline position

```text
*_HV23_merged.xlsx
+
*_HV1_merged.xlsx
        ↓
ms_final_profiles_merge.py
        ↓
*_final_profiles.xlsx
```

---

## Command-line usage

```bash
python3 ms_final_profiles_merge.py \
  --hv23 "All_15_HV23_merged.xlsx" \
  --hv1 "All_15_HV1_merged.xlsx" \
  --output "All_15_final_profiles.xlsx"
```

Arguments:

| Argument | Required | Meaning |
|---|---:|---|
| `--hv23` | yes | Input workbook from `ms_hv23_merge.py`, usually `*_HV23_merged.xlsx` |
| `--hv1` | yes | Input workbook from `ms_hv1_merge.py`, usually `*_HV1_merged.xlsx` |
| `--output` | yes | Output final profile workbook, usually `*_final_profiles.xlsx` |

---

## Input expectations

### HV2-3 input

The script reads sheet:

```text
LID
```

Important columns used when present:

```text
Batch
Well
LID
Range by HV2-3 Regions Adjusted
Variants by HV2-3 Regions Adjusted
HV2-3 Pair Status
HV3R PHP 341-437 Sites
HV2-3 Flag list
HV2-3 Review Class
```

For compatibility, if an older workbook contains:

```text
HV3R PHP 341_437 Sites
```

it is normalized internally to:

```text
HV3R PHP 341-437 Sites
```

### HV1 input

The script reads sheet:

```text
LID
```

Important columns used when present:

```text
Batch
Well
LID
HV1 Pair Status
Range by HV1 Region Adjusted
Variants by HV1 Region Adjusted
Mismatch by HV1 Region Adjusted
HV1 Flag list
HV1 Review-triggering Flag list
HV1 Review Class
```

---

## Merge key

HV2-3 and HV1 profiles are merged by:

```text
Batch + LID
```

`Well` is **not** used as a merge key because HV2-3 and HV1 traces may not always come from the same well.

After merging:

```text
Well = Well from HV2-3 if available;
       otherwise Well from HV1
```

---

## Final sample profile construction

### `Sample profile - Range`

The final range is the concatenation of adjusted HV2-3 and HV1 ranges:

```text
Range by HV2-3 Regions Adjusted
+
Range by HV1 Region Adjusted
```

Blank parts are skipped.

Example:

```text
73-340 438-576 16024-16365
```

### `Sample profile - Variants`

The final variant profile is built from:

```text
Variants by HV2-3 Regions Adjusted
+
Variants by HV1 Region Adjusted
```

Tokens are:

```text
split
→ deduplicated
→ sorted by numeric position
→ joined back into one string
```

Example:

```text
73G 263G 309.1C 315.1C 16093C 16189C
```

---

## Combined flag logic

The combined flag list is built from:

```text
HV2-3 Flag list
+
HV1 Review-triggering Flag list
```

The script intentionally uses `HV1 Review-triggering Flag list` rather than the full `HV1 Flag list` because HV1 may contain informational-only flags such as `HV1 polyC variant`.

---

## Combined review logic

The script assigns:

```text
Combined Review Class
```

Rule:

```text
Review
    if HV2-3 Review Class = Review
    or HV1 Review Class = Review

Auto - no flag
    otherwise
```

---

## PC / NTC handling

The full output keeps PC/NTC controls in the main `LID` sheet.

The script also creates no-control sheets for production statistics:

```text
LID_no_PC_NTC
Summary_no_PC_NTC
```

Rows are excluded from these sheets if `LID` starts with:

```text
PC
NTC
```

case-insensitively.

---

## Output workbook

The output workbook contains four sheets:

```text
1. LID_no_PC_NTC
2. Summary_no_PC_NTC
3. LID
4. Summary
```

The no-PC/NTC sheets are placed first because they are usually the preferred sheets for reporting production statistics.

---

## Main output columns

```text
Batch
Well
LID
Sample profile - Range
Sample profile - Variants
Range by HV2-3 Regions Adjusted
Variants by HV2-3 Regions Adjusted
HV2-3 Pair Status
HV3R PHP 341-437 Sites
HV2-3 Flag list
HV2-3 Review Class
HV1 Pair Status
Range by HV1 Region Adjusted
Variants by HV1 Region Adjusted
Mismatch by HV1 Region Adjusted
HV1 Flag list
HV1 Review-triggering Flag list
HV1 Review Class
Combined Flag list
Combined Review Class
```

Only columns that exist in the merged input data are retained.

---

## Summary sheets

Both `Summary_no_PC_NTC` and `Summary` include:

```text
Total LID
Combined Review counts
HV2-3 Review counts
HV1 Review counts
Combined flag counts
```

Each summary row includes:

```text
Metric
Value
%
```

Percent is calculated relative to the total number of LIDs in that sheet and rounded to two decimal places.

---

## Design principle

This script answers:

```text
Given the HV2-3 and HV1 regional MS profiles, what is the final MS profile for each Batch + LID?
```

It deliberately does not answer:

```text
Does the final MS profile match the truthset?
```

That question is handled by `ms_final_profiles_vs_truth.py`.
