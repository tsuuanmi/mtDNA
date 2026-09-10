# TNLS/HCLS Verification Guide

This guide explains how the TNLS verification code compares HCLS and TNLS
mtDNA variant data, including uncertain bases, insertion positions, result
statuses, family summaries, and the FAM-001 example.

## 1. What is being compared?

- **HCLS** is the known or base sample.
- **TNLS** is the target sample being checked against HCLS.
- The comparison uses the configured mtDNA hypervariable regions: **HV1, HV2,
  and HV3**.
- A variant is compared only where the HCLS and TNLS sequenced intervals
  overlap.
- The default minimum overlap is **150 base pairs**.

The matcher compares one HCLS sample against every selected TNLS sample. The
main pair output contains only `CANNOT_EXCLUDE` rows. Internally, pairs can
also be classified as `INCONCLUSIVE`, `EXCLUSION`, or
`MR_MINIMUM_NOT_REACH`.

## 2. Input and output files

The command-line comparison receives:

| Input | Purpose |
|---|---|
| HCLS JSON | Base samples, usually `statistic_fullbatch.json` or equivalent |
| TNLS merged JSON | Target samples and their variants |
| `matching_rule.json` | Overlap, region, and mismatch thresholds |
| Optional filter list | Restricts TNLS samples to selected IDs |
| Optional metadata workbook | Adds TNLS identity fields to statistics outputs |

The comparison writes these files under the selected results directory:

| File | Contents |
|---|---|
| `output_compare_tnls_hcls.tsv` | Per-pair `CANNOT_EXCLUDE` matches |
| `stats_TNLS_HCLS_with_metadata.tsv` | Expanded 1:1 through 1:N matched pairs |
| `stats_TNLS_HCLS_clean.tsv` | Simplified matched-pair view |
| `stats_TNLS_HCLS_family_summary.tsv` | TNLS profile families and their matched HCLS IDs |

Because non-matching pairs are filtered out of the main pair output, an absent
pair does not by itself show whether the reason was insufficient overlap, one
mismatch, or multiple mismatches. Use the Stage 2/Stage 3 logic below or
intermediate logs when investigating a missing pair.

## 3. Matching configuration

The current `resources/matching_rule.json` is:

| Setting | Current value | Meaning |
|---|---:|---|
| `minimum` | `150` | Minimum overlapping base pairs required |
| `region` | `HV1`, `HV2`, `HV3` | Regions included in the comparison |
| `conclusion_threshold` | `0` | `CANNOT_EXCLUDE` allows zero mismatches only |
| `exclusion_threshold` | `2` | Two or more mismatches cause `EXCLUSION` |

The thresholds are read when `CompareManager` is created. The code does not
change the input JSON files.

## 4. The three comparison stages

### Stage 1: Sequenced interval overlap

For each HCLS/TNLS pair, the matcher calculates the intersection of their
sequenced intervals. Intervals are inclusive: `[start, end]` contains
`end - start + 1` base pairs.

- Overlap of **150 bp or more**: continue to Stage 2.
- Overlap below **150 bp**: classify internally as `MR_MINIMUM_NOT_REACH`.

A variant outside the common overlap is not compared. This is important when
one sample has a wider region or when the two samples have a coverage gap.

### Stage 2: Position-only fast filter

Stage 2 compares only the set of variant positions in the common overlap. It
ignores the actual base values for speed.

A position mismatch occurs when a normalized position exists on only one side.
The mismatch count is symmetric:

```text
positions only in HCLS + positions only in TNLS
```

If the count is **2 or more**, the pair is immediately classified as
`EXCLUSION`; Stage 3 is not run for that pair.

If the count is **0 or 1**, the pair continues to Stage 3.

#### Position normalization

Variant positions are normalized before Stage 2 and Stage 3. Therefore these
values are treated as the same insertion position:

```text
HCLS: "315.1"   (JSON string)
TNLS: 315.1     (JSON number)
```

This prevents the same position from being counted twice as a false mismatch.
The original JSON files are not modified.

Position comparison does not use the variant type. For example, a SNP and an
insertion at the same numeric position can cancel in Stage 2 because their
positions are equal. Stage 3 still compares typed keys (`snp_`, `ins_`, and
`del_`) and can classify the different variant types as a mismatch.

### Stage 3: Base and variant comparison

Stage 3 compares typed variant keys inside the common overlap:

- `snp_<position>` for SNPs
- `ins_<position>` for insertions
- `del_<position>` for deletions

For a two-sided variant, the reference allele must match exactly. The sequence
alleles are then checked according to the rules in the next section.

For a one-sided variant, the variant is a mismatch unless it is an HCLS SNP
with `N`, which is ignored.

## 5. Base and IUPAC rules

### What does `N` mean?

`N` means the base is unresolved. It may be A, C, G, or T.

The current rule is intentionally different for the two sides:

- **HCLS SNP `N`: ignored.** It is neither a match nor a mismatch.
- **TNLS SNP `N`: not automatically ignored.** If the same typed HCLS SNP
  exists, it is compared as an IUPAC code. If TNLS has an `N`-only variant that
  HCLS does not have, it is a one-sided mismatch.

The HCLS `N` rule applies even when TNLS has no variant at that position. This
means an HCLS-only unresolved call does not reject the pair.

The ignore rule applies to an HCLS **SNP** call. It does not make a different
TNLS variant type at the same position disappear. For example, an HCLS
`snp_100` with `N` and a TNLS `ins_100` are still not the same typed variant;
the TNLS insertion remains a mismatch.

### IUPAC compatibility

The matcher expands both SNP calls into their possible concrete bases. Two
SNP calls match when their possible-base sets have a non-empty intersection.
The comparison is case-insensitive for IUPAC letters.

| Code | Possible bases |
|---|---|
| `A` | A |
| `C` | C |
| `G` | G |
| `T` | T |
| `R` | A or G |
| `Y` | C or T |
| `S` | C or G |
| `W` | A or T |
| `K` | G or T |
| `M` | A or C |
| `B` | C, G, or T |
| `D` | A, G, or T |
| `H` | A, C, or T |
| `V` | A, C, or G |
| `N` | A, C, G, or T |

Examples:

| HCLS SNP | TNLS SNP | Compatible? | Reason |
|---|---|---|---|
| `A` | `A` | Yes | Same possible base |
| `R` | `G` | Yes | `G` is included in `R` |
| `R` | `M` | Yes | `{A,G}` intersects `{A,C}` at `A` |
| `R` | `Y` | No | `{A,G}` and `{C,T}` do not intersect |
| `N` | `A` | HCLS `N` is ignored | HCLS unresolved call is not evidence or a mismatch |
| `N` | `N` | HCLS `N` is ignored | The HCLS call is skipped before allele comparison |
| `?` | `A` | No | `?` is not a valid IUPAC code |

Invalid codes expand to an empty set and are incompatible. Invalid codes are
not treated as matching just because the same invalid character appears on
both sides.

### Insertions and deletions

IUPAC expansion applies only to SNP sequence alleles. Insertions and deletions
must match exactly after position normalization:

| HCLS | TNLS | Result |
|---|---|---|
| `ins_315.1: C` | `ins_315.1: C` | Match |
| `ins_315.1: C` | `ins_315.1: G` | Mismatch |
| `del_523: -` | `del_523: -` | Match |
| `del_523: -` | No `del_523` | One-sided mismatch |
| `ins_315.1` | `del_315.1` | Mismatch; variant types differ |

## 6. Result classification

After Stage 3, the internal result is determined by the base-level mismatch
count and the number of matched variants:

| Condition | Internal result | Included in main pair output? |
|---|---|---|
| Overlap below 150 bp | `MR_MINIMUM_NOT_REACH` | No |
| Stage 2 position mismatches >= 2 | `EXCLUSION` | No |
| Stage 3 mismatches >= 2 | `EXCLUSION` | No |
| Stage 3 mismatches = 1 | `INCONCLUSIVE` | No |
| Stage 3 mismatches = 0 and matched variants = 0 | `INCONCLUSIVE` | No |
| No variants are available in the common overlap | `INCONCLUSIVE` | No |
| Stage 3 mismatches = 0 and matched variants > 0 | `CANNOT_EXCLUDE` | Yes |

`CANNOT_EXCLUDE` does not prove that the samples came from the same person.
It means the available overlapping mtDNA evidence does not exclude that
possibility under the configured rules.

An ignored HCLS `N` does not count as matched evidence. A pair with only
ignored HCLS `N` calls and no other matched variant remains `INCONCLUSIVE`.

## 7. FAM-001 example

FAM-001 contains TNLS samples `2525154` and `2525366`. Their representative
profile includes:

```text
HV1: 16172C 16223T 16311C
HV2: 73G 146C 263G 315.1C
HV3: 489C 523DEL 524DEL
```

HCLS `NCC5279` includes:

```text
73G 146C 16037A>N 16172C 16223T 315.1C
```

The important facts are:

1. HCLS stores `315.1` as a string while TNLS stores it as a number. Position
   normalization removes that false mismatch.
2. HCLS has `16037A>N`, while FAM-001 has no variant at position `16037`.
   Because this is an unresolved HCLS SNP call, it is ignored.
3. The resolved shared variants provide positive evidence. For the current
   coverage overlap, the matched displayed variants are:

   ```text
   73G 146C 263G 315.1C 16172C 16223T
   ```

4. Both `2525154` and `2525366` therefore produce `CANNOT_EXCLUDE` against
   `NCC5279` after the code fix.
5. The family summary aggregates those per-pair matches, so `NCC5279` appears
   in `Matched_HCLS_IDs` for FAM-001 after the pipeline is rerun.

The other HCLS samples in the inspected batch still have at least two real
position mismatches with FAM-001 and remain `EXCLUSION`. Variants outside the
common sequenced intervals, such as a TNLS variant in a region that the HCLS
sample did not sequence, do not contribute to that pair's comparison.

## 8. Family grouping and family summary

Family grouping is a separate TNLS-only step. It groups TNLS samples with
compatible profiles over their common coverage before attaching HCLS matches.

### Family grouping rules

- The same configured regions and 150 bp minimum are used.
- Profiles use typed positions, so SNP, insertion, and deletion positions are
  distinct profile entries.
- SNP profiles use the same IUPAC possible-base intersection rule.
- Insertions and deletions require exact sequence equality.
- Exact coverage/profile duplicates are collapsed first.
- Compatible links are transitive. If A matches B and B matches C, all three
  can belong to one family even if A and C were not directly compared.
- A family ID such as `FAM-001` is generated sequentially and is separate from
  any source barcode family ID.
- The family representative is the member with the broadest total configured
  region coverage; sample ID breaks ties.

HCLS `N` is not part of family grouping because family grouping compares TNLS
samples only. HCLS `N` is handled by the HCLS/TNLS pair matcher.

### Family summary columns

`stats_TNLS_HCLS_family_summary.tsv` contains:

| Column | Meaning |
|---|---|
| `Family_ID` | Generated family ID, such as `FAM-001` |
| `Member_IDs` | Comma-and-space-separated TNLS members |
| `Profile_Frequency` | Family size / total selected TNLS sample count, formatted as a fraction such as `15/38721` |
| `Profile_Variants` | Complete representative profile, combining the space-joined HV1, HV2, and HV3 variants |
| `Matched_HCLS_IDs` | HCLS IDs with at least one `CANNOT_EXCLUDE` member pair |
| `HCLS_matched_variants` | Distinct matched variants for each matched HCLS |
| `HCLS_matched_variant_sites` | Distinct matched site count, aligned with matched HCLS IDs |
| `HCLS_matched_variant_bp` | Matched base count, aligned with matched HCLS IDs |
| `HCLS_ratios` | HCLS-to-TNLS ratio, such as `1:2` |
| `HV1_variants` | Full representative profile in HV1 |
| `HV2_variants` | Full representative profile in HV2 |
| `HV3_variants` | Full representative profile in HV3 |

The `HV*_variants` columns are display-only. They show the representative's
full regional profile and are not themselves used to decide whether a family
matches HCLS.

A family has no `Matched_HCLS_IDs` when no member produces a
`CANNOT_EXCLUDE` pair. This can happen because every pair has insufficient
overlap, one mismatch, or at least two mismatches.

## 9. Rerunning after a code change

Generated result files are not changed automatically when source code changes.
Rerun the verification pipeline to regenerate the pair and family outputs.
Use the verification script from the TNLS verification module, for example:

```bash
cd src/modules/TNLS/verification
bash verification.sh --run-step0
```

Use `--skip-filter` when the existing TNLS filtered input should be reused. The
pipeline may also accept both flags when HCLS regeneration and the existing
filter behavior are both required.

After rerunning, inspect:

1. `output_compare_tnls_hcls.tsv` for exact HCLS/TNLS pairs.
2. `stats_TNLS_HCLS_clean.tsv` for the simplified matched view.
3. `stats_TNLS_HCLS_family_summary.tsv` for family-level HCLS matches.

Do not edit generated TSV or JSON files manually; fix the source logic or
input data and regenerate the outputs.
