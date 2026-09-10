# 04 — MS trace-processed output → HV1 merged profiles

## Script

`ms_hv1_merge.py`

Functionally, this is the no-truthset v3 HV1 merge script with the production naming convention.

## Pipeline position

```text
*_MS_trace_processed.xlsx
    ↓
ms_hv1_merge.py
    ↓
*_HV1_merged.xlsx
```

This script is designed to run after `ms_trace_process.py`.

## Purpose

This script processes Mutation Surveyor HV1 trace-level output and produces one merged HV1 profile per LID.

It keeps only HV1F and HV1R traces, trims each trace to the HV1 target region, selects the best HV1F/HV1R pair per LID, merges the selected traces, detects F/R mismatch in overlap, applies HV1-specific flags, and classifies each LID as either auto-pass or review.

The v3 no-truthset version does **not** require truthset columns and does **not** calculate truthset concordance.

## Input

The input is an Excel file produced by the `ms_trace_process.py` step.

Required columns:

```text
LID
Primer
MS Range
Variant Converted 2
```

Common optional columns are preserved in the `Traces` sheet because the full trace-level dataframe is copied to output, but they are not required by this script.

## Output

The script writes an Excel file with three sheets:

```text
Traces
LID
Summary
```

## Output sheet: `Traces`

This sheet contains the filtered HV1 trace-level records. Only traces whose `Primer` is `HV1F` or `HV1R` are retained.

The script appends trace-level adjusted columns:

```text
HV1 Trace Range Adjusted
HV1 Trace Variants Adjusted
```

It also appends LID-level results back onto each trace belonging to the same LID:

```text
HV1 Pair Status
Range by HV1 Region Adjusted
Variants by HV1 Region Adjusted
Mismatch by HV1 Region Adjusted
HV1 Flag list
HV1 Review-triggering Flag list
HV1 Review Class
```

## Output sheet: `LID`

This sheet contains one merged HV1 profile per LID.

Columns written by the script:

```text
Batch
Well
LID
HV1 Pair Status
HV1F Candidate Count
HV1R Candidate Count
HV1F Selected Trace #
HV1R Selected Trace #
HV1F Selected Sample Name
HV1R Selected Sample Name
Range by HV1 Region Adjusted
Variants by HV1 Region Adjusted
Mismatch by HV1 Region Adjusted
HV1 Flag list
HV1 Review-triggering Flag list
HV1 Review Class
```

## Output sheet: `Summary`

This sheet summarizes:

```text
Total LID
Pair Status counts
Review Class counts
Flag counts
Review-triggering Flag counts
```

## Biological regions

### HV1 target region

```text
16024-16365
```

Stored in code as:

```python
HV1_TARGET = [(16024.0, 16365.0)]
```

### HV1 polyC region

```text
16180-16195
```

Stored in code as:

```python
HV1_POLYC = [(16180.0, 16195.0)]
```

Note that this v3 no-truthset script uses `16180-16195`, not the narrower `16180-16193` used in earlier versions.

## High-level workflow

```text
Read input Excel
→ normalize column names
→ validate required columns
→ keep only HV1F/HV1R traces
→ trim each trace to 16024-16365
→ select best HV1F and best HV1R per LID
→ determine pair status
→ merge selected trace ranges without merging adjacent ranges
→ merge selected trace variants within merged adjusted ranges
→ detect HV1F/HV1R mismatch inside true overlap only
→ scan selected raw source traces for het calls
→ build full flag list and review-triggering flag list
→ assign review class
→ write Traces, LID, Summary sheets
```

## Trace filtering

The script keeps only rows where:

```text
Primer = HV1F or HV1R
```

Primer matching is case-insensitive because the script converts `Primer` to uppercase before filtering.

All non-HV1F/HV1R rows are removed before any downstream logic.

## Trace trimming

Each trace is trimmed to the HV1 target region:

```text
16024-16365
```

The input `MS Range` is parsed as one range pair. The trace range is intersected with the HV1 target.

Examples:

```text
MS Range = 15990-16200
→ HV1 Trace Range Adjusted = 16024-16200
```

```text
MS Range = 16190-16400
→ HV1 Trace Range Adjusted = 16190-16365
```

```text
MS Range = 15000-15900
→ HV1 Trace Range Adjusted = blank
```

Variants from `Variant Converted 2` are retained only if their position falls inside the adjusted trace range.

The adjusted variants are written to:

```text
HV1 Trace Variants Adjusted
```

## Variant parsing

Variant tokens are parsed from whitespace- or comma-separated text.

The expected token format is:

```text
position + allele
```

Examples:

```text
16189C
16223T
16182.1het_C
16189DEL
16192Y
```

The parser accepts integer and decimal positions.

`_HET` is normalized to `_het`.

Tokens that do not match the expected pattern are ignored by the parser.

## Trace scoring and selection

For each LID, the script selects one best HV1F trace and one best HV1R trace.

Each primer is evaluated independently:

```text
all HV1F candidates → select one HV1F
all HV1R candidates → select one HV1R
```

The score is calculated as:

```text
coverage length
+ full region bonus
- polyC variant penalty
- PHP penalty
- het penalty
- indel outside polyC penalty
```

In code:

```text
coverage
+ 1000 if one adjusted trace range covers 16024-16365
- 20 × number of polyC variants
- 50 if trace has PHP
- 50 if trace has het
- 80 × number of indels outside polyC
```

Tie-breakers:

```text
1. higher score
2. smaller Trace #
3. Sample Name alphabetically
```

## Pair status

After selecting the best available HV1F and HV1R traces, the script assigns one of four pair statuses:

```text
OK
Missing HV1F
Missing HV1R
Missing HV1F and HV1R
```

`OK` means both selected HV1F and selected HV1R exist.

Any missing primer status is review-triggering.

## Range merge behavior in v3

This is one of the key v3 changes.

The script preserves trace-boundary ranges and does **not** merge adjacent ranges.

The function call is:

```python
normalize_ranges(merged_ranges, merge_adjacent=False)
```

This means only truly overlapping ranges are merged.

### Adjacent ranges are not merged

```text
16024-16189 16190-16365
→ 16024-16189 16190-16365
```

```text
16024-16263 16264-16365
→ 16024-16263 16264-16365
```

These are not collapsed into `16024-16365`.

### Overlapping ranges are merged

```text
16024-16200 16190-16365
→ 16024-16365
```

because the ranges overlap between `16190-16200`.

### Gapped ranges remain separate

```text
16024-16188 16190-16365
→ 16024-16188 16190-16365
```

## Full HV1 region rule

A LID is treated as having full HV1 only if at least one merged adjusted range continuously covers the full HV1 target:

```text
16024-16365
```

This is checked by:

```python
covers_segment(merged_ranges, HV1_TARGET[0])
```

Because adjacent ranges are not merged in v3, the following does **not** count as full HV1:

```text
16024-16189 16190-16365
```

The following also does **not** count as full HV1:

```text
16024-16263 16264-16365
```

These cases receive:

```text
No full HV1 region
```

and are sent to Review.

## Variant merge behavior

The selected HV1F and HV1R adjusted variants are merged into one LID-level variant string.

The script keeps variants only if they fall inside the final merged adjusted ranges.

Duplicate variants are collapsed by:

```text
position + allele uppercased
```

Output variants are sorted by position and allele.

The merged variant string is written to:

```text
Variants by HV1 Region Adjusted
```

## HV1F/HV1R mismatch detection

Mismatch detection is performed only when both selected HV1F and selected HV1R exist.

The script first computes true overlap between the selected adjusted HV1F range and selected adjusted HV1R range.

Adjacent ranges are not treated as overlap.

Example:

```text
HV1F = 16024-16189
HV1R = 16190-16365
→ no overlap
→ no mismatch positions are checked
```

Example:

```text
HV1F = 16024-16200
HV1R = 16190-16365
→ overlap = 16190-16200
→ variants inside 16190-16200 are compared
```

For each variant position in the overlap:

- If one trace has a variant at the position and the other does not, the position is reported as mismatch.
- If both traces have variants but alleles are incompatible, the position is reported as mismatch.
- If alleles are compatible by exact match or IUPAC compatibility, no mismatch is reported.

Mismatch positions are written to:

```text
Mismatch by HV1 Region Adjusted
```

Any mismatch is review-triggering via:

```text
HV1 trace mismatch
```

## Allele compatibility

The script uses IUPAC ambiguity compatibility for mismatch comparison.

Supported IUPAC codes:

```text
W = A/T
R = A/G
M = A/C
K = G/T
Y = C/T
S = C/G
D = A/G/T
H = A/C/T
V = A/C/G
B = C/G/T
N = A/C/G/T
```

Examples considered compatible:

```text
Y vs C
Y vs T
R vs A
R vs G
N vs C
```

Exact allele matches are also compatible.

`_het` is removed during this allele compatibility check, but het is still separately flagged elsewhere.

## Source het scan before adjusted merge

v3 scans the selected raw source traces before relying on the merged adjusted profile.

The function checks `Variant Converted 2` from the selected HV1F and HV1R rows.

If either selected source trace contains `het`, the script adds:

```text
HV1 source het
```

If the source het token is inside the HV1 polyC region `16180-16195`, the script also adds:

```text
HV1 source polyC het
```

Both source het flags are review-triggering.

This protects against cases where het evidence exists in the original selected trace but could be obscured by later adjusted/merged representation.

## Flag system

v3 separates flags into two outputs:

```text
HV1 Flag list
HV1 Review-triggering Flag list
```

### `HV1 Flag list`

This is the full flag list. It includes both informational and review-triggering flags.

### `HV1 Review-triggering Flag list`

This contains only the flags that force manual Review.

The review class is derived only from this field.

## Informational-only flag

The key informational-only flag is:

```text
HV1 polyC variant
```

This means at least one merged adjusted variant falls inside the HV1 polyC region:

```text
16180-16195
```

By itself, `HV1 polyC variant` does **not** trigger Review.

A sample with only this flag is classified as:

```text
Auto - no flag
```

## Review-triggering flags

The following flags trigger Review:

```text
Missing HV1F
Missing HV1R
Missing HV1F and HV1R
No full HV1 region
HV1 trace mismatch
HV1 polyC insertion
HV1 polyC deletion
HV1 polyC PHP
PHP
het
HV1 indel outside polyC
HV1 source het
HV1 source polyC het
```

## PolyC flag logic

The script examines merged adjusted variants inside:

```text
16180-16195
```

If any variant exists in this region:

```text
HV1 polyC variant
```

is added as informational only.

If any polyC token has a decimal position, it is treated as an insertion:

```text
HV1 polyC insertion
```

If any polyC allele contains `DEL`, it is treated as a deletion:

```text
HV1 polyC deletion
```

If any polyC allele is a PHP/IUPAC code from:

```text
R M S K Y W V H B D N
```

then:

```text
HV1 polyC PHP
```

is added.

PolyC insertion, deletion, and PHP are review-triggering.

## General PHP flag

If any merged adjusted variant allele is a PHP/IUPAC code from:

```text
R M S K Y W V H B D N
```

then the script adds:

```text
PHP
```

This is review-triggering.

## General het flag

If the merged adjusted variant string contains `het`, the script adds:

```text
het
```

This is review-triggering.

This is separate from the source-level het scan.

A sample can therefore have:

```text
het
HV1 source het
HV1 source polyC het
```

depending on where het is observed.

## Indel outside polyC flag

For every merged adjusted variant, the script checks whether the token is an insertion or deletion.

Insertion is defined as:

```text
decimal position
```

Deletion is defined as:

```text
allele contains DEL
```

If the insertion/deletion is outside `16180-16195`, the script adds:

```text
HV1 indel outside polyC
```

This is review-triggering.

## Review class

The final review class is based only on `HV1 Review-triggering Flag list`.

```text
blank review-triggering flag list → Auto - no flag
non-blank review-triggering flag list → Review
```

Output values:

```text
Auto - no flag
Review
```

## No truthset / no concordance behavior

This v3 script intentionally removes truthset handling.

It does not require columns such as:

```text
HV1 Range
HV1 Variants
Sample Profiles - Range
Sample Profile - Flag
```

It also does not write:

```text
HV1 Truth Status
Concordance_auto_2 (by HV1 Region Adjusted)
```

The output is intended for production-style processing where truthset comparison is no longer part of this HV1 merge step.

## Excel formatting

The output workbook is formatted with:

```text
freeze panes at A2
auto-filter on each sheet
auto column width capped at 50 characters
```

## Command-line usage

```bash
python3 ms_hv1_merge.py \
  --input input_MS_trace_processed.xlsx \
  --output output_HV1_merged.xlsx
```

## Important examples

### PolyC SNP-only sample

```text
Range by HV1 Region Adjusted = 16024-16365
Variants by HV1 Region Adjusted = 16189C 16223T
HV1 Flag list = HV1 polyC variant
HV1 Review-triggering Flag list = blank
HV1 Review Class = Auto - no flag
```

### Adjacent split HV1 range

```text
Range by HV1 Region Adjusted = 16024-16189 16190-16365
HV1 Flag list = No full HV1 region; HV1 polyC variant
HV1 Review-triggering Flag list = No full HV1 region
HV1 Review Class = Review
```

### PolyC insertion

```text
Variants by HV1 Region Adjusted = 16182.1het_C 16189C
HV1 Flag list = HV1 polyC insertion; HV1 polyC variant; HV1 source het; HV1 source polyC het; het
HV1 Review-triggering Flag list = HV1 polyC insertion; HV1 source het; HV1 source polyC het; het
HV1 Review Class = Review
```

### F/R mismatch in overlap

```text
HV1F overlap variant = 16223T
HV1R overlap variant = blank
Mismatch by HV1 Region Adjusted = 16223
HV1 Flag list includes HV1 trace mismatch
HV1 Review-triggering Flag list includes HV1 trace mismatch
HV1 Review Class = Review
```

## Main v3 behavior changes to remember

Compared with earlier HV1 versions, this v3 no-truthset script is stricter and more production-oriented:

1. Adjacent ranges are preserved and are not merged into one continuous range.
2. Full HV1 requires one continuous range covering `16024-16365`.
3. `HV1 polyC variant` is informational only.
4. PolyC insertion, deletion, and PHP remain review-triggering.
5. Selected raw source traces are scanned for het before relying on adjusted/merged profile.
6. Truthset concordance is removed.
