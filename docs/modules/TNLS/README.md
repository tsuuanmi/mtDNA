# TNLS Module

> Tìm Nguồn Lành Sáng (Kinship/Ancestry) analysis utilities for mtDNA data.

## Purpose

The TNLS module (`src/modules/TNLS/`) provides utilities for processing mtDNA variant data for kinship analysis and metadata management, including barcode processing, family grouping, and kinship classification.

## Module Files

| File | Purpose |
|------|---------|
| `extract.py` | Extract soldier names from metadata (`SoldierNameExtractor`) |
| `filter_family_variants.py` | Filter and normalize family variants for comparison |
| `fix_barcode_merge.py` | Fix and merge barcode data with family information |
| `kinship_analysis.py` | Kinship classification using likelihood ratios |
| `map_cccd.py` | Map CCCD (national ID) to sample data |
| `merge_family_info.py` | Merge family information into sample data |
| `merge_str_mtdna_ids.py` | Merge STR and mtDNA sample identifiers |
| `upload_mapping_to_gsheets.py` | Upload mapping data to Google Sheets |

### Verification Submodule

| File | Purpose |
|------|---------|
| `verification/compare_1n_class.py` | Compare one-generation classification results (`Compare_Manager`) |
| `verification/compare_TNLS_HCLS.py` | Compare TNLS results with HCLS classification + mtDNA-profile family summary |
| `verification/utils.py` | Shared verification utilities (JSON loading, variant comparison) |
| [Verification matching guide](verification.md) | Simple explanation of HCLS/TNLS matching, `N`/IUPAC calls, `315.1`, and the FAM-001 example |

## Key Concepts

### Kinship Analysis

`kinship_analysis.py` uses likelihood ratio calculations to classify relationships:
- Parent-child, sibling, and other familial relationships
- Probabilistic classification based on mtDNA haplotype matching

### Family Variant Filtering

`filter_family_variants.py` normalizes variants across family members:
- IUPAC-aware equivalence checking
- Common variant identification within families
- Unique variant detection for individual identification

### Family Grouping by mtDNA Profile

`compare_TNLS_HCLS.py` uses the indexed `family_profiles.py` implementation to
group TNLS target samples by identical mtDNA profile (IUPAC-tolerant, over the
common ≥150 bp overlap; transitive via union-find). Exact profiles are collapsed
before overlap-restricted canonical profiles are compared, avoiding exhaustive
sample-pair enumeration. The verification writes a single family-level summary
of HCLS (remains) matches
(`stats_TNLS_HCLS_family_summary.tsv`) that also carries the per-HCLS match
statistics (matched variants, ratio, exclusivity, matched variant sites/bp) — merging the former
`stats_TNLS_HCLS_matching.tsv`. Family IDs use the `FAM-NNN`
format, distinct from the barcode `FAMxxxx` families. The variant columns are
reference only and are not part of the matching criteria.

### Barcode Processing

`fix_barcode_merge.py` handles:
- Barcode validation and formatting
- Vietnamese name normalization
- Family ID creation from barcode data

### Blind Copy Processing

The `blind_copy.py` module provides functions for:
- Loading metadata and mapping files
- Parsing FASTA sequences and splitting by HV region
- Processing samples with variant data
- Creating JSON output files for each sample

### Variant Integration

`add_variants_to_samples.py` enriches the sample mapping with:
- Per-sample variant information grouped by HV region
- Variant statistics (SNP counts, insertion counts, deletion counts)
- Region-specific variant formatting
- Analyzed intervals (`Analyzed_Range` column), reusing the shared TNLS interval formatter from `add_variants_to_metadata.py`

## Cross-References

- [core/models.md](../../core/models.md) — Sample and Variant data models
- [modules/sheets/uploader.md](../sheets/uploader.md) — Google Sheets upload
- [variants.md](../../variants.md) — Variant utilities
