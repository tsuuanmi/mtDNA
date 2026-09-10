# 06 — `ms_final_profiles_vs_truth.py`

## Purpose

`ms_final_profiles_vs_truth.py` compares final MS mtDNA profiles against a truthset by normalized `LID`.

It evaluates three scopes:

```text
1. Full profile: HV2-3 + HV1
2. HV2-3 adjusted region only
3. HV1 adjusted region only
```

This helps determine whether discordance comes from HV2-3 processing, HV1 processing, both regions, or only the final combined profile.

---

## Pipeline position

```text
*_final_profiles.xlsx
+
truthset.xlsx
        ↓
ms_final_profiles_vs_truth.py
        ↓
*_final_profiles_vs_truth.xlsx
```

---

## Command-line usage

```bash
python3 ms_final_profiles_vs_truth.py \
  --merged "All_15_final_profiles.xlsx" \
  --truth "All_15_Truthset.xlsx" \
  --output "All_15_final_profiles_vs_truth.xlsx"
```

Arguments:

| Argument | Required | Meaning |
|---|---:|---|
| `--merged` | yes | Final MS profile workbook from `ms_final_profiles_merge.py` |
| `--truth` | yes | Truthset workbook |
| `--output` | yes | Output concordance workbook |

---

## MS input workbook

The script preferentially reads sheet:

```text
LID_no_PC_NTC
```

If that sheet is not available, it looks for a sheet name containing `lid`, `pc`, and `ntc`. If no such sheet is found, it falls back to the first sheet.

Required MS columns include:

```text
LID
Sample profile - Range
Sample profile - Variants
Variants by HV2-3 Regions Adjusted
Variants by HV1 Region Adjusted
```

The source MS rows and columns are preserved in the output, and concordance columns are appended.

---

## Truthset input workbook

Required truthset columns:

```text
LID
Sample Profile - Range
Sample Profile - Variants
```

Optional truthset column:

```text
Sample Profile - Flag
```

If the truth range value is:

```text
FULL REGION
```

it is expanded internally to:

```text
73-340 438-576 16024-16365
```

---

## LID matching

MS and truth rows are joined by normalized `LID`.

Normalization includes:

```text
- converting numeric-like Excel values such as 2543704.0 to 2543704
- removing incidental whitespace
- preserving non-numeric IDs
```

This prevents false missing-truth results caused by Excel type differences.

---

## Regions used

### Full profile

```text
73-340
438-576
16024-16365
```

Comparison:

```text
MS Sample profile - Variants
vs
Truth Sample Profile - Variants filtered to full region
```

### HV2-3 region

```text
73-340
438-576
```

Comparison:

```text
MS Variants by HV2-3 Regions Adjusted
vs
Truth Sample Profile - Variants filtered to HV2-3 regions
```

### HV1 region

```text
16024-16365
```

Comparison:

```text
MS Variants by HV1 Region Adjusted
vs
Truth Sample Profile - Variants filtered to HV1 region
```

---

## Truthset region splitting

The truthset usually stores the full profile in one variant column. The script splits this full truth profile by variant position into:

```text
Variants by HV2-3 Region (Truth)
Variants by HV1 Region (Truth)
```

Example:

```text
Truth full profile:
73G 263G 315.1C 16189C 16223T

Truth HV2-3:
73G 263G 315.1C

Truth HV1:
16189C 16223T
```

---

## Variant parsing

Variant strings are split on:

```text
space
semicolon
comma
```

The following values are treated as empty:

```text
OK
NaN
None
Null
-
blank
```

Variant position is parsed from the leading numeric component:

```text
146C       → 146
315.1C     → 315.1
16182.1C   → 16182.1
```

---

## HV2-3 length-variant ignore rule

For full-profile and HV2-3 comparisons, the following length variants are ignored:

```text
309.xC
573.xC
```

Specifically:

```text
309 < position < 310 and allele == C
573 < position < 574 and allele == C
```

This ignore rule is **not** applied to HV1 comparison.

---

## Het handling

Het is treated as a distinct biological state.

Therefore:

```text
146_hetC != 146C
146C_het != 146C
```

Equivalent het spellings are normalized before exact comparison:

```text
146C_het       → 146_hetC
146_hetC       → 146_hetC
16182.1het_C   → 16182.1_hetC
16182.1_het_C  → 16182.1_hetC
```

Examples:

```text
MS:    146C_het
Truth: 146_hetC
Result: Y
```

```text
MS:    146_hetC
Truth: 146C
Result: F
```

`Y2` is not emitted.

---

## Concordance labels

Possible outputs:

```text
Y
Y1
F
```

### `Y`

All comparable variants match exactly after canonicalization.

### `Y1`

At least one non-het IUPAC-compatible allele match exists and no discordance is found.

Example:

```text
MS:    146Y
Truth: 146C
```

because:

```text
Y = C/T
```

### `F`

Any of the following causes discordance:

```text
missing MS variant
extra MS variant
incompatible allele
het vs non-het
variant count/content mismatch at the same position
missing truth row by LID
```

---

## Output columns appended to the MS rows

```text
Full Profile Concordance
Full Profile Discordance reason
HV2-3 Concordance
HV2-3 Discordance reason
HV1 Concordance
HV1 Discordance reason
Discordance Category
Range (Truth)
Variants (Truth)
Flags (Truth)
Range by HV2-3 Region (Truth)
Variants by HV2-3 Region (Truth)
Range by HV1 Region (Truth)
Variants by HV1 Region (Truth)
```

---

## Discordance Category

The script classifies discordance source as:

```text
Both_HV23_and_HV1
HV23_only
HV1_only
Full_profile_only
Concordant
```

Logic:

```text
HV2-3 = F and HV1 = F → Both_HV23_and_HV1
HV2-3 = F and HV1 != F → HV23_only
HV2-3 != F and HV1 = F → HV1_only
HV2-3 != F and HV1 != F and Full = F → Full_profile_only
otherwise → Concordant
```

---

## Output workbook

The output workbook contains three sheets:

```text
1. Main sheet, usually LID_no_PC_NTC
2. Summary
3. Discordance_detail
```

### Main sheet

Preserves the selected MS final profile sheet and appends truth/concordance columns.

### Summary

Reports counts and proportions for:

```text
Full Profile Concordance
HV2-3 Concordance
HV1 Concordance
Discordance Category
```

### Discordance_detail

Stores row-level details for discordant events, including:

```text
LID
Scope
Concordance
Discordance reason
Position
MS variant
Truth variant
Issue
MS Sample profile - Variants
MS HV2-3 Variants Adjusted
MS HV1 Variants Adjusted
Truth profile - Variants
Truth HV2-3 Variants
Truth HV1 Variants
```

---

## Duplicate truth LIDs

If duplicate truthset `LID` values are found:

```text
- the first truth row is used
- duplicate truth records are reported in Discordance_detail
```

---

## Important limitation

HV1 is compared strictly.

No dedicated HV1 polyC equivalence normalization is applied in this script. Therefore, biologically related but differently represented HV1 polyC patterns may still be classified as discordant.
