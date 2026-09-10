# 03 — MS trace-processed output → HV2-3 merged profiles

## Script

`ms_hv23_merge.py`

## Purpose

This script processes Mutation Surveyor (MS) HV2-3 trace-level output and merges selected `HV2F` and `HV3R` traces into one HV2-3 profile per `LID`.

Version 18 is designed to work with the MS-only upstream output from:

```text
ms_trace_process.py
```

It also remains compatible with older truthset-containing MS_trace_processed files if truth columns are present. In no-truth mode, the script leaves truth/concordance/discordance fields blank instead of failing.

## Pipeline position

```text
MS custom report
    ↓
ms_trace_process.py
    ↓
*_MS_trace_processed.xlsx
    ↓
ms_hv23_merge.py
    ↓
*_HV23_merged.xlsx
```

## Command line usage

```bash
python3 ms_hv23_merge.py \
  --input  input_MS_trace_processed.xlsx \
  --output output_HV23_merged.xlsx
```

Arguments:

| Argument | Required | Meaning |
|---|---:|---|
| `--input` | yes | Input Excel file from the MS trace-processing step |
| `--output` | yes | Output Excel file containing `Traces`, `LID`, and `Summary` sheets |

## Input expectations

The script expects a trace-level Excel file with at least these columns:

```text
Well
LID
Primer
Trace #
Sample Name
Read Start
Read End
MS Range
Variant Converted 2
```

The no-truth upstream script also provides these useful QC columns when available:

```text
Variant Summary
Variant Converted
C_shift_rule_applied
AC_repeat_rule_applied
HV1_PolyC_flag
HV2-3_PolyC_flag
```

Older truthset-based files may also contain:

```text
HV2-3 Range
HV2-3 Variants
Sample Profiles - Range
Sample Profile - Flag
Concordance_auto
Concordance_auto_2
Error_type
```

In v18, `HV2-3 Variants` is optional.

## Primer filtering

The script keeps only rows where `Primer` is:

```text
HV2F
HV3R
```

All other primers are removed from this stage.

## Biological target regions

The HV2-3 biological target used by the script is:

```text
73-340
438-576
```

These regions are stored internally as:

```python
HV23_TARGET = [(73.0, 340.0), (438.0, 576.0)]
```

## Adjusted trace ownership

Version 18 uses adjusted primer ownership with a true overlap anchor at position `303`:

```text
HV2F: 73-303
HV3R: 303-340 and 438-576
```

Internally:

```python
HV2F_ADJUSTED = [(73.0, 303.0)]
HV3R_ADJUSTED = [(303.0, 340.0), (438.0, 576.0)]
```

The position `303` is intentionally included in both `HV2F` and `HV3R`. This creates a real one-base overlap, so the script can merge `73-303` and `303-340` into `73-340` only because the ranges truly overlap, not because they are merely adjacent.

## Output workbook

The output Excel workbook contains three sheets:

```text
Traces
LID
Summary
```

## Sheet 1 — `Traces`

This sheet keeps trace-level rows for HV2F/HV3R, with merged LID-level HV2-3 results repeated on each trace row.

Important columns include:

```text
Batch
Well
LID
Primer
Trace #
Sample Name
Read Start
Read End
MS Range
Variant Summary
Variant Converted
Variant Converted 2
HV2-3 Pair Status
HV2-3 Truth Status
HV2-3 Trace Range Adjusted
HV2-3 Trace Variants Adjusted
Range by HV2-3 Regions
Variants by HV2-3 Regions
Mismatch by HV2-3 Regions
Concordance_auto_2 (by HV2-3 Regions)
Range by HV2-3 Regions Adjusted
Variants by HV2-3 Regions Adjusted
Mismatch by HV2-3 Regions Adjusted
Concordance_auto_2 (by HV2-3 Regions Adjusted)
HV1_PolyC_flag
HV2-3_PolyC_flag
```

Only columns present in the input are retained. Missing optional columns are skipped safely.

## Sheet 2 — `LID`

This sheet contains one row per `LID`.

Important columns include:

```text
Batch
Well
LID
HV2-3 Pair Status
HV2-3 Truth Status
Range by HV2-3 Regions
Variants by HV2-3 Regions
Mismatch by HV2-3 Regions
Concordance_auto_2 (by HV2-3 Regions)
Range by HV2-3 Regions Adjusted
Variants by HV2-3 Regions Adjusted
Mismatch by HV2-3 Regions Adjusted
Concordance_auto_2 (by HV2-3 Regions Adjusted)
HV3R PHP 341-437 Count
HV3R PHP 341-437 Sites
HV2-3 Flag list
HV2-3 Review Class
HV2-3 Severity
HV2-3 Discordance reason
HV2-3 Range
HV2-3 Variants
Sample Profiles - Range
Sample Profile - Flag
```

In no-truth mode, truth-related fields such as `HV2-3 Variants`, concordance, and discordance are blank.

## Sheet 3 — `Summary`

This sheet summarizes the LID-level output.

It includes counts for:

```text
Total LID
Review Class
Severity
Pair Status
Adjusted Concordance
Flag frequency
```

In no-truth mode, adjusted concordance will usually be counted under `(blank)`.

## Trace selection logic

For each `LID`, the script selects one representative `HV2F` and one representative `HV3R` trace.

If multiple traces exist for the same primer, selection is deterministic:

1. longest `MS Range` span
2. smallest numeric `Trace #`
3. `Sample Name` alphabetical order

This is handled by `select_trace_pair()`.

Pair status is assigned as:

| Condition | `HV2-3 Pair Status` |
|---|---|
| both HV2F and HV3R selected | `OK` |
| HV2F present, HV3R missing | `Missing HV3R` |
| HV3R present, HV2F missing | `Missing HV2F` |
| both missing | `Missing HV2F and HV3R` |

Any pair status other than `OK` triggers the flag:

```text
Pair status not OK
```

and this is a high-severity condition.

## Base HV2-3 merge

The base view uses the original `MS Range` and `Variant Converted 2` from selected HV2F/HV3R traces.

Steps:

1. Parse each selected trace's `MS Range`.
2. Intersect the range with the biological HV2-3 target:

   ```text
   73-340
   438-576
   ```

3. Merge ranges only when they truly overlap.
4. Do not merge merely adjacent ranges.
5. Merge unique variant tokens from `Variant Converted 2` if they fall inside the final base ranges.
6. Detect mismatch positions between HV2F and HV3R in any overlapping region.

Base output columns:

```text
Range by HV2-3 Regions
Variants by HV2-3 Regions
Mismatch by HV2-3 Regions
Concordance_auto_2 (by HV2-3 Regions)
```

In no-truth mode, base concordance is blank.

## Adjusted HV2-3 merge

The adjusted view is the preferred production view.

Trace-level trimming:

```text
HV2F → keep only 73-303
HV3R → keep only 303-340 and 438-576
```

For each trace, the script writes:

```text
HV2-3 Trace Range Adjusted
HV2-3 Trace Variants Adjusted
```

Then at LID level, it merges the adjusted selected traces into:

```text
Range by HV2-3 Regions Adjusted
Variants by HV2-3 Regions Adjusted
Mismatch by HV2-3 Regions Adjusted
Concordance_auto_2 (by HV2-3 Regions Adjusted)
```

## 303 overlap anchor

Version 18 intentionally creates a true overlap at position `303`:

```text
HV2F adjusted range: 73-303
HV3R adjusted range: 303-340
```

This allows the script to merge the first HV2-3 block into:

```text
73-340
```

because the two component ranges overlap at `303`.

This is different from merging adjacent ranges. The script's `normalize_ranges(..., merge_adjacent=False)` does not merge ranges just because one ends at `302` and the next starts at `303`. The overlap exists because both traces include `303`.

## Mismatch detection

Mismatch detection is performed by `mismatch_positions_between_two_rows()`.

For adjusted mismatch detection, the script compares:

```text
HV2-3 Trace Variants Adjusted
```

between the selected HV2F and HV3R traces, but only inside their overlapping adjusted ranges.

With v18 adjusted ranges, the expected overlap is:

```text
303-303
```

A mismatch is reported if:

1. one trace has a variant at a position and the other trace does not; or
2. both traces have a variant at the same position but the alleles are incompatible.

Alleles are considered compatible if:

- they are identical after `_HET` normalization; or
- one allele is compatible with an IUPAC ambiguity code; or
- two IUPAC ambiguity codes share at least one possible base.

If there is a mismatch at `303`, the output will include:

```text
Mismatch by HV2-3 Regions Adjusted = 303
```

Any adjusted mismatch triggers:

```text
Trace mismatch
```

This sends the LID to:

```text
Review
HIGH
```

## Range display when mismatch exists

The current v18 code builds adjusted range continuity based on range overlap, not mismatch status.

Therefore, if HV2F covers `73-303` and HV3R covers `303-340`, the adjusted range can still display:

```text
73-340 438-576
```

even if the overlap position `303` is listed in:

```text
Mismatch by HV2-3 Regions Adjusted
```

The mismatch is not hidden; it is reported separately and triggers high-severity review.

## Variant parsing and sorting

The script expects variant tokens in the form:

```text
<position><allele>
```

Examples:

```text
73G
146C
309.1C
315.1C
324S
249DEL
```

Parsing is performed by `parse_variant_tokens()`.

Token sorting is numeric by position, then by allele text.

## Duplicate variant handling

During merge, duplicate variants are kept only once.

The uniqueness key is:

```text
(position, allele.upper())
```

This means repeated identical calls from HV2F/HV3R are collapsed into one output token.

## Optional truthset compatibility

Version 18 supports two input modes.

### No-truth mode

If the input does not contain:

```text
HV2-3 Variants
```

then:

```text
HV2-3 Truth Status = blank
Concordance_auto_2 columns = blank
HV2-3 Discordance reason = blank
```

The script still performs:

```text
trace selection
base merge
adjusted merge
mismatch detection
flagging
review classification
severity classification
summary generation
```

### Truth-compatible mode

If the input contains:

```text
HV2-3 Variants
```

then the script uses it as optional truth data.

Truth status:

| Condition | `HV2-3 Truth Status` |
|---|---|
| all non-empty `HV2-3 Variants` values for an LID are consistent | `OK` |
| multiple non-empty truth profiles disagree | `Inconsistent HV2-3 truth` |
| no truth column | blank |

Concordance is calculated only when truth data exists.

## Concordance logic

Concordance is calculated by `concordance_hv23()`.

It compares merged MS variants against `HV2-3 Variants` within the merged range only.

Possible results:

| Value | Meaning |
|---|---|
| `Y` | exact match |
| `Y1` | compatible via IUPAC ambiguity code |
| `Y2` | compatible with heteroplasmy-related compatibility |
| `F` | discordant |
| blank | no comparable truth data or no variants to compare |

The comparison ignores these length-variant insertions:

```text
309.xC
573.xC
```

These are not treated as discriminating discordances.

## Discordance reason

If adjusted concordance is `F`, the script generates:

```text
HV2-3 Discordance reason
```

Possible components:

```text
Missing [...]
Extra [...]
Mismatch [...]
```

Examples:

```text
Missing [204Y]
Extra [310C]
Mismatch [146C!=146Y]
```

The discordance reason uses the same ignored length-variant rule as concordance:

```text
309.xC
573.xC
```

## HV3R source-trace PHP QC: 341-437

Version 18 includes source-trace QC for selected HV3R.

The script checks the selected HV3R trace's original:

```text
Variant Converted 2
```

not only the adjusted merged variants.

It counts PHP/IUPAC calls in:

```text
341-437
```

This window is important because it sits outside the adjusted merged profile but may reveal poor-quality HV3R signal that would otherwise disappear after trimming.

PHP/IUPAC codes are:

```text
R M S K Y W V H B D N
```

Output columns:

```text
HV3R PHP 341-437 Count
HV3R PHP 341-437 Sites
```

Flag rule:

```text
count 0-1 → no HV3R PHP 341-437 flag
count >=2 → add flag HV3R PHP 341-437
```

Severity rule:

```text
count >=2 → MEDIUM, unless a HIGH rule is also present
count >=3 → HIGH
```

## Flag list

The script builds `HV2-3 Flag list` from adjusted merged data and selected-source-trace QC.

Flags include:

```text
No full region
Pair status not OK
Truth status not OK
Trace mismatch
HV3R PHP 341-437
No 315.1C
PHP
het
310 variant
non-C in HV2 polyC
non-C in HV3 polyC
unexpected insertion
unexpected deletion
```

Flags are deduplicated and sorted before output.

## Full-region rule

A LID is considered full-region only if the adjusted merged range covers all required segments:

```text
73-303
303-340
438-576
```

If not, the script adds:

```text
No full region
```

This is a high-severity flag.

Because `73-303` and `303-340` overlap at `303`, a successful full adjusted merge commonly displays as:

```text
73-340 438-576
```

## 315.1C rule

If the adjusted merged variants do not contain:

```text
315.1C
```

then the script adds:

```text
No 315.1C
```

This is treated as a low-severity flag unless other higher-severity flags are present.

## PHP rule in adjusted merged variants

If any adjusted merged variant allele is an IUPAC/PHP code:

```text
R M S K Y W V H B D N
```

then the script adds:

```text
PHP
```

This is a medium-severity flag unless overridden by a high-severity rule.

## Heteroplasmy rule

If adjusted merged variants contain text matching `het` case-insensitively, the script adds:

```text
het
```

This is a medium-severity flag unless overridden by a high-severity rule.

## 310 variant rule

If any adjusted merged variant has integer position `310`, the script adds:

```text
310 variant
```

This is a medium-severity flag unless overridden by a high-severity rule.

## HV2 poly-C rule

The script checks adjusted merged variants in the HV2 C-stretch windows:

```text
303-309
311-315
```

If a variant in these windows has an allele other than `C`, the script adds:

```text
non-C in HV2 polyC
```

This is a medium-severity flag unless overridden by a high-severity rule.

## HV3 poly-C rule

The script checks adjusted merged variants in:

```text
568-573
```

If a variant in this window has an allele other than `C`, the script adds:

```text
non-C in HV3 polyC
```

This is a medium-severity flag unless overridden by a high-severity rule.

## Unexpected insertion rule

Any non-integer variant position is treated as an insertion candidate.

The following insertions are considered common/allowed and do not trigger `unexpected insertion`:

```text
309.xC
315.xC
573.xC
452.xT through 454.xT
```

Any other insertion triggers:

```text
unexpected insertion
```

This is high severity.

## Unexpected deletion rule

Deletion tokens are detected when the allele contains:

```text
DEL
```

The following deletions are considered common/allowed and do not trigger `unexpected deletion`:

```text
249DEL
DEL in 513-525
```

Other deletions trigger:

```text
unexpected deletion
```

This is high severity.

## Review class

The review class is assigned from `HV2-3 Flag list`.

| Condition | `HV2-3 Review Class` |
|---|---|
| flag list is blank | `Auto - no flag` |
| any flag is present | `Review` |

## Severity classification

Severity is assigned by `hv23_severity()`.

Possible values:

```text
AUTO
LOW
MEDIUM
HIGH
```

### AUTO

If the flag list is empty and concordance is one of:

```text
Y
Y1
Y2
```

then severity is:

```text
AUTO
```

In no-truth mode, concordance is blank. If there are no flags, the final fallback also returns `AUTO`.

### HIGH

High severity is assigned for:

```text
HV3R PHP 341-437 Count >= 3
No full region
Pair status not OK
Trace mismatch
unexpected deletion
unexpected insertion
concordance F with discordance reason
concordance F fallback
```

### MEDIUM

Medium severity is assigned for:

```text
HV3R PHP 341-437 Count >= 2
PHP
het
Truth status not OK
310 variant
non-C in HV2 polyC
non-C in HV3 polyC
HV3R PHP 341-437
```

### LOW

Low severity is assigned for:

```text
No 315.1C
```

when no higher-severity rule applies.

## Important behavior in no-truth mode

In no-truth mode, review decisions are driven by:

```text
pair completeness
adjusted range completeness
trace mismatch
variant-content flags
HV3R source-trace PHP QC
insertion/deletion rules
```

They are not dependent on truthset concordance.

This makes v18 suitable as a production-oriented MS-only HV2-3 merge and review-prep step.

## Implementation notes

1. Adjusted range display is based on range overlap, not mismatch status. A mismatch at `303` is reported separately and triggers high-severity review, but the adjusted range may still display as `73-340 438-576`.

2. HV3R source-trace PHP QC in `341-437` uses the current production threshold:

   ```text
   >=2 PHP in 341-437 → add HV3R PHP 341-437 flag / usually MEDIUM
   >=3 PHP in 341-437 → HIGH
   ```

## Output formatting

The script writes the workbook using `openpyxl`.

For each output sheet:

```text
freeze top row
apply autofilter
auto-size columns up to width 50
```

The script also prints a terminal summary after saving:

```text
Saved: <output path>
Summary:
  Total LID: ...
  Review Class - ...
  Severity - ...
  Adjusted Concordance - ...
```
