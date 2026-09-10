# TNLS Verification Pipeline

Compares **HCLS** (known/base) samples against **TNLS** (target/unknown) samples using mtDNA
hypervariable region (HV1, HV2, HV3) matching. Produces a list of `CANNOT_EXCLUDE` pairs.

---

## Quick Start

```bash
# Full pipeline (Step 0 regenerates HCLS JSON from raw data, then Steps 1-4)
bash verification.sh --run-step0

# Use existing tnls_filtered.json — skips Step 2
bash verification.sh --skip-filter

# Both flags together
bash verification.sh --run-step0 --skip-filter
```

Results are written to `results/`.

---

## Folder Structure

```
verification/
├── resources/
│   └── matching_rule.json       # Matching thresholds
├── data/
│   ├── base/                   # HCLS base samples
│   │   └── HCLS_*/statistic_fullbatch.json
│   └── target/                 # TNLS target samples
│       ├── filter.txt          # Optional LID filter list
│       └── TNLS_metadata.xlsx  # Identity metadata (optional)
│       #  NOTE: merged_statistics.json is read from the mnt mount directly.
│       #  Use --skip-filter to run without filter.txt.
├── results/                    # All outputs (TSV, TXT)
│
├── verification.sh              # ← Run this: end-to-end pipeline
├── export_sample_lists.py       # Export sample lists (Step 1)
├── compare_TNLS_HCLS.py        # Main comparison with inline LID filter (Step 2)
├── compare_TNLS_HCLS_stats.py   # Summary statistics (Step 3)
├── compare_1n_class.py         # CompareManager class (logic only)
└── utils.py                   # Shared utilities
```

---

## Matching Rules (`resources/matching_rule.json`)

| Key | Value | Meaning |
|-----|-------|---------|
| `minimum` | 150 | Minimum overlap (bp) between base and target to proceed |
| `region` | HV1, HV2, HV3 | mtDNA hypervariable regions to use |
| `conclusion_threshold` | 0 | Max mismatches for `CANNOT_EXCLUDE`; a pair must also have at least one matched variant in the overlap |
| `exclusion_threshold` | 2 | Min mismatches for `EXCLUSION` |

---

## Pipeline (3 Steps)

```
Step 1: export_sample_lists.py    Export sample IDs from all JSON files → <json>.txt
Step 2: compare_TNLS_HCLS.py      Compare HCLS vs TNLS → output_compare_tnls_hcls.tsv
                               (inline LID filter applied here if filter.txt is present)
Step 3: compare_TNLS_HCLS_stats.py  Summary statistics + breakdowns
```

---

## Results

### Conclusion Values

| Result | Meaning |
|--------|---------|
| `CANNOT_EXCLUDE` | Samples may be from the same individual (≤ 0 mismatches and at least one matched variant in the overlap) |
| `INCONCLUSIVE` | Results between thresholds, including no-mismatch pairs with zero matched variants |
| `EXCLUSION` | Samples are from different individuals (≥ 2 mismatches or insufficient overlap) |
| `MR_MINIMUM_NOT_REACH` | Overlap region too short (< minimum bp, default 150) |

### Output Files

| File | Description |
|------|-------------|
| `results/modules/TNLS/verification/output_compare_tnls_hcls.tsv` | Per-pair: overlap, matched/mismatched variants & bases |
| `results/modules/TNLS/verification/stats_TNLS_HCLS_summary.txt` | Overall summary + ratio + variant breakdown |
| `results/modules/TNLS/verification/stats_TNLS_HCLS_1to1_1to2_with_tnls_metadata.tsv` | 1:1 & 1:2 pairs with identity metadata |
| `results/modules/TNLS/verification/stats_TNLS_HCLS_1to1_1to2_clean.tsv` | Clean view: HCLS, status, TNLS, overlap, variants, identity |
| `results/modules/TNLS/verification/stats_TNLS_HCLS_family_summary.tsv` | Family-level summary: mtDNA-profile families with per-HCLS match stats (ratio, exclusivity, matched TNLS, matched variant sites/bp) |

---

## Family Grouping by mtDNA Profile

After the per-pair comparison, `compare_TNLS_HCLS.py` uses
`family_profiles.py` to group TNLS target samples into **mtDNA-profile families**
and writes `stats_TNLS_HCLS_family_summary.tsv`.

- **Grouping:** individuals with an identical mtDNA profile over their common
  overlapping interval (≥ the matching rule `minimum`, default 150 bp) belong to
  the same family. SNP calls are IUPAC-tolerant (a hetero/homo pair such as `A`
  vs `R` is compatible); insertions/deletions must match exactly. Grouping is
  transitive (connected components via union-find), so a family may contain samples
  that do not all pairwise overlap. Exact coverage/profile duplicates are collapsed
  first; remaining comparisons are indexed by coverage, typed variant positions,
  and overlap-restricted canonical profiles instead of enumerating every sample pair.
- **Family ID:** `FAM-001`, `FAM-002`, … (sequential; distinct from the barcode
  `FAMxxxx` families).
- **Matching:** aggregated from the existing per-pair `CANNOT_EXCLUDE` output —
  a family matches an HCLS remain if any member is `CANNOT_EXCLUDE` against it. The
  mtDNA matching engine is unchanged.
- **Columns** (in order): `Family_ID`, `Member_IDs`, `Profile_Frequency`,
  `Profile_Variants`, `Matched_HCLS_IDs`,
  `HCLS_matched_variants`, `HCLS_matched_variant_sites`,
  `HCLS_matched_variant_bp`, `HCLS_ratios`, `HV1_variants`,
  `HV2_variants`, `HV3_variants`. The per-HCLS columns are comma-and-space-joined
  aligned with `Matched_HCLS_IDs`, empty for families with no matched HCLS.
  `HCLS_matched_variants` is the distinct union of matched variants across the
  HCLS's matched TNLS pairs in the simplified display form (e.g.
  `16093C 16189C 309.1C`), so a 1:2 match against two identical TNLS reports the
  variants once (5, not 10). `Profile_Frequency` is the family size and total
  selected TNLS sample count formatted as a fraction (for example, `15/38721`).
  `Profile_Variants` is the representative's complete profile,
  combining the space-joined HV1, HV2, and HV3 variants. Member / matched-HCLS
  counts and the per-HCLS matched-TNLS list are omitted as derivable from
  `Member_IDs` / `Matched_HCLS_IDs` (the exact HCLS->TNLS pairs live in the per-pair
  `output_compare_tnls_hcls.tsv`). The `HV*_variants` columns contain the full
  variant set (SNPs + insertions + deletions) from the representative member with
  the broadest interval coverage, space-joined in the same simplified form as
  `HCLS_matched_variants`; they are reference only and are **not** used for
  matching. This merges the former `stats_TNLS_HCLS_matching.tsv` (per-HCLS
  matched variants, ratio, exclusivity, and distinct matched variant totals)
  into the family summary.

---

## Dependencies

- Python 3.12+
- `pandas`
- `loguru`

---

## Notes

- Scripts use `Path(__file__).parent` — safe to run from any working directory
- `results/` is always wiped before each run
- Only `CANNOT_EXCLUDE` pairs appear in the main comparison output
- SNP comparison uses IUPAC possible-base compatibility; unresolved HCLS `N` calls are ignored and indels match exactly
- Numeric and string forms of insertion positions (for example `315.1` and `"315.1"`) are normalized before comparison
- See [`docs/modules/TNLS/verification.md`](../../../../docs/modules/TNLS/verification.md) for examples and the FAM-001 case
