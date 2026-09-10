# ms_trace_process.py — Trace Processing and Variant Conversion Logic

## Purpose

`ms_trace_process.py` is the first processing step in the MS mtDNA post-processing pipeline.

It reads a Mutation Surveyor custom report, extracts real `.ab1` trace rows, standardizes MS variant calls, applies repeat-region normalization, performs mandatory Layer 0 inventory QC, and writes a trace-level processed output file.

The script does **not** use a truthset and does **not** calculate concordance. Truthset comparison is performed later by `ms_final_profiles_vs_truth.py`.

---

## Main inputs

```bash
python3 ms_trace_process.py \
  --ms "{BATCH}_anything_df_custom_report.xlsx" \
  --lid-manifest "YYYYMMDD_Danh sách mẫu_{BATCH}.xlsx" \
  --out "{BATCH}_MS_trace_processed.xlsx"
```

Optional:

```bash
--inventory-out "{BATCH}_inventory_QC.xlsx"
```

---

## Main output

```text
{BATCH}_MS_trace_processed.xlsx
```

This is the standardized trace-level output used by downstream scripts:

```text
ms_control_QC.py
ms_hv23_merge.py
ms_hv1_merge.py
```

---

## High-level workflow

```text
MS custom report
    ↓
Read all repeated MS report blocks
    ↓
Keep real .ab1 trace rows only
    ↓
Build Variant Summary
    ↓
Build Variant Converted
    ↓
Build Variant Converted 2
    ↓
Derive Well / LID / Primer / MS Range
    ↓
Layer 0 Inventory QC
    ↓
Add HV1 / HV2-3 repeat-region QC flags
    ↓
Write {BATCH}_MS_trace_processed.xlsx
```

Important note: in the current implementation, `Variant Summary`, `Variant Converted`, and `Variant Converted 2` are created before Layer 0 is executed. If Layer 0 fails, the script exits before writing the production `MS_trace_processed` output.

---

# 1. Reading the Mutation Surveyor custom report

Mutation Surveyor custom reports can contain multiple repeated table blocks.

The script scans the first 500 rows and detects header rows containing both:

```text
Trace #
Sample Name
```

For each detected block, it reads rows until the next header block.

Only rows where `Sample Name` contains:

```text
.ab1
```

are retained.

Duplicate rows are removed with `drop_duplicates()`.

---

# 2. Cleaning real trace rows

After reading all MS blocks, the script keeps only rows that satisfy both conditions:

```text
Sample Name contains .ab1
Read Start is not blank
Read End is not blank
```

This removes non-trace rows, summary rows, and incomplete trace rows.

---

# 3. Variant Summary

## Purpose

`Variant Summary` combines all MS variant columns into one clean string.

The script detects columns named:

```text
Variant1
Variant2
Variant3
...
```

These columns are sorted numerically and joined into a space-separated list.

## Values ignored

The following are ignored:

```text
blank cells
NaN
n.a.
```

## Cleaning rules

After joining the tokens, the script:

```text
removes c. prefix
removes commas
collapses repeated whitespace
converts nan / None strings to blank
```

## Example

Input MS columns:

```text
Variant1 = c.73A>G$148
Variant2 = c.263A>G$104
Variant3 = 309.1C$16
```

Output:

```text
Variant Summary = 73A>G$148 263A>G$104 309.1C$16
```

---

# 4. Variant Converted

## Purpose

`Variant Converted` is the first standardized representation of MS variant calls.

Each token from `Variant Summary` is processed independently.

## Token-level cleaning

For each token, the script removes:

```text
MS signal suffixes such as $148
commas
c. prefix
```

## Plain tokens

Tokens already in internal form are kept as-is.

Examples:

```text
309.1C  → 309.1C
521DEL  → 521DEL
```

## Substitution tokens

Tokens in this format:

```text
position REF>ALT
```

are converted to:

```text
position ALT
```

Examples:

```text
73A>G   → 73G
263A>G  → 263G
146T>C  → 146C
```

## Multi-base ALT and IUPAC conversion

If the ALT contains more than one base, the unordered ALT set is converted to an IUPAC ambiguity code when possible.

Examples:

```text
A/G → R
C/T → Y
A/C → M
A/T → W
G/T → K
C/G → S
```

Example conversions:

```text
73A>AG  → 73R
73A>CT  → 73Y
73A>AC  → 73M
```

If the ALT combination is not in the IUPAC map, the ALT string is retained.

---

# 5. Variant Converted 2

## Purpose

`Variant Converted 2` is the main normalized variant representation used by downstream scripts.

It starts from `Variant Converted` and applies row-level, context-aware normalization for known repeat/stutter regions.

The function also creates two debug columns:

```text
C_shift_rule_applied
AC_repeat_rule_applied
```

---

## 5.1 C-stretch normalization: 303–309

C insertions in positions `303–309` are shifted to position `309`.

Rule:

```text
303-309 → 309
```

This rule is applied only if there is **no non-C variant** anywhere in `303–309`.

Example:

```text
303.1C → 309.1C
304.1C → 309.1C
```

If a non-C variant is present in the same region, the shift is not applied.

When applied, the debug column records:

```text
303-309_to_309
```

---

## 5.2 C-stretch normalization: 311–315

C insertions in positions `311–315` are shifted to position `315`.

Rule:

```text
311-315 → 315
```

This rule is applied only if there is **no non-C variant** anywhere in `311–315`.

Example:

```text
311.1C → 315.1C
```

When applied, the debug column records:

```text
311-315_to_315
```

---

## 5.3 T-stretch normalization: 452–455

T insertions in positions `452–455` are shifted to position `455`.

Rule:

```text
452-455 → 455
```

This rule is applied only if there is **no non-T variant** anywhere in `452–455`.

Example:

```text
452.1T → 455.1T
```

When applied, the debug column records:

```text
452-455_to_455T
```

---

## 5.4 C-stretch normalization: 568–573

C insertions in positions `568–573` are shifted to position `573`.

Rule:

```text
568-573 → 573
```

This rule is applied only if there is **no non-C variant** anywhere in `568–573`.

Example:

```text
568.1C → 573.1C
```

When applied, the debug column records:

```text
568-573_to_573
```

---

## 5.5 248 deletion normalization

The script applies a specific deletion correction:

```text
248DEL → 249DEL
```

This is applied only if `248DEL` is present and `249DEL` is not already present.

When applied, the debug column records:

```text
248DEL_to_249DEL
```

---

## 5.6 AC repeat 513–525 normalization

The script contains a dictionary/consensus-based normalization module for the HVIII AC repeat region:

```text
513–525
```

The local rCRS reference used by the script is:

```text
513 G
514 C
515 A
516 C
517 A
518 C
519 A
520 C
521 A
522 C
523 A
524 C
525 C
```

Combined sequence:

```text
GCACACACACACC
```

## Supported token forms

Inside `513–525`, the module supports:

```text
SNP tokens: 513A, 523G, etc.
DEL tokens: 513DEL, 524DEL, etc.
```

## Unsupported token forms

The module does not normalize unsupported local forms such as:

```text
insertions, e.g. 523.1C
heteroplasmy-like unsupported forms, e.g. 523het_C
other ambiguous/complex tokens
```

If unsupported local tokens are found within `513–525`, the AC repeat region is left unchanged.

## Canonicalization process

For supported local tokens, the script:

```text
1. Builds a local consensus sequence from the observed SNP/DEL tokens.
2. Enumerates possible alignments against the rCRS 513–525 reference.
3. Scores each alignment.
4. Selects the best canonical representation.
5. Converts the best alignment back to variant tokens.
```

## Alignment scoring preference

The best alignment is selected by minimizing, in order:

```text
1. AC motif break count
2. Number of contiguous deletion blocks
3. Deletion placement score favoring right alignment
4. Number of substitutions
5. Substitution placement score favoring right alignment
```

If canonicalization changes the local token representation, the debug column records:

```text
dictionary_consensus_canonicalized
```

---

# 6. Derived columns

After variant conversion, the script derives several trace-level fields.

## Well

```text
Well = first 3 characters of Sample Name
```

## LID

The LID is extracted from `Sample Name` using this pattern:

```regex
^[^_]+_\d{8}_(.*?)(?=_HV)
```

This assumes a sample name structure like:

```text
<well_or_prefix>_<date>_<LID>_HV...
```

## Primer

The primer is extracted from `Sample Name` using:

```text
HV1F, HV1R, HV2F, HV2R, HV3F, HV3R
```

The result is converted to uppercase.

Current production Layer 0 expects only:

```text
HV1F
HV1R
HV2F
HV3R
```

## MS Range

`MS Range` is created from:

```text
Read Start-Read End
```

Trailing `.0` is removed.

Example:

```text
73.0 and 302.0 → 73-302
```

---

# 7. Layer 0 Inventory QC

Layer 0 is a strict batch-level input integrity gate.

It checks:

```text
custom report filename
manifest filename
batch code consistency
manifest LID formatting
manifest duplicate/blank LID
missing real sample LID
extra non-control LID
expected primer inventory per LID
primer total count
unexpected primer
```

The required custom report filename format is:

```text
{BATCH}_anything_df_custom_report.xlsx
```

The required manifest filename format is:

```text
YYYYMMDD_Danh sách mẫu_{BATCH}.xlsx
```

The manifest must contain:

```text
LID
MẺ CHẠY
```

PC/NTC control LIDs are allowed in the MS report but excluded from manifest matching.

Recognized control formats:

```text
NTC1, NTC2, ...
PC1_mtDNA_247547, PC2_mtDNA_247547, ...
```

If Layer 0 fails:

```text
- inventory QC report is written
- production MS_trace_processed output is not generated
- script exits with code 1
```

---

# 8. Repeat-region QC flags

After Layer 0 passes, the script adds preliminary repeat-region QC flags.

These flags are based on `Variant Converted`, not `Variant Converted 2`, so they reflect the original pre-normalized repeat signal.

---

## 8.1 HV1_PolyC_flag

Applied only to HV1 traces.

The script counts variant tokens in:

```text
16180–16193
```

Thresholds:

| Number of hits | Flag |
|---:|---|
| 0 | blank |
| 1 | LOW |
| 2–4 | SUSPECT |
| ≥5 | BAD |

---

## 8.2 HV2-3_PolyC_flag

Applied only to HV2/HV3 traces.

The script counts C-insertion tokens matching:

```regex
^\d+\.\d+C$
```

within:

```text
303–315
```

Thresholds:

| Number of hits | Flag |
|---:|---|
| 0 | blank |
| 1 | LOW |
| 2–4 | SUSPECT |
| ≥5 | BAD |

---

# 9. Final output columns

The final `{BATCH}_MS_trace_processed.xlsx` output keeps these columns:

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

---

# 10. What this script does not do

`ms_trace_process.py` does not:

```text
- compare with truthset
- calculate concordance
- create final sample profiles
- merge HV23/HV1 regions
- determine final auto-pass status
- evaluate PC/NTC control profiles
```

PC/NTC control QC is handled by:

```text
ms_control_QC.py
```

Truthset comparison is handled later by:

```text
ms_final_profiles_vs_truth.py
```

---

# 11. Summary

`ms_trace_process.py` converts raw Mutation Surveyor custom report rows into a standardized trace-level representation.

Its key responsibilities are:

```text
1. Read and clean MS trace rows.
2. Build Variant Summary from Variant1/Variant2/... columns.
3. Convert raw MS variant notation into Variant Converted.
4. Normalize repeat-region representation into Variant Converted 2.
5. Derive Well, LID, Primer, and MS Range.
6. Enforce Layer 0 inventory QC.
7. Add preliminary HV1/HV2-3 repeat-region flags.
8. Write {BATCH}_MS_trace_processed.xlsx for downstream processing.
```

The output is trace-level and MS-only. It is the foundation for control QC, HV23 merge, HV1 merge, final profile generation, and final comparison against truthset.
