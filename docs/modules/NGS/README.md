# NGS Module

> Next-Generation Sequencing analysis utilities for mtDNA data.

## Purpose

The NGS module (`src/modules/NGS/`) provides utilities for processing Next-Generation Sequencing mtDNA data, including FIS (Forensic Identity System) analysis, variant frequency calculation, and comparison with Sanger sequencing data.

## Module Files

| File | Purpose |
|------|---------|
| `FIS_merge_fasta.py` | Merge FASTA sequences from multiple Excel files |
| `FIS_merge_variants.py` | Merge variant data from multiple Excel files |
| `calculate_variant_freq.py` | Calculate variant frequency from JSON data |
| `compare_fis_fasta_direct.py` | Compare FIS FASTA data directly with pipeline results |
| `compare_fis_sanger_direct.py` | Legacy standalone FIS-to-Sanger comparator (retained, not used by the NGS runner) |
| `prepare_fis_comparison.py` | Map and filter canonical FIS Samples for shared comparison |
| `prepare_sanger_comparison.py` | Select and normalize mapped Sanger samples for shared comparison |
| `enrich_fis_sanger_comparison.py` | Join FIS metadata onto shared comparison results |
| `split_fis_sanger_comparison_sets.py` | Split a FIS-to-Sanger report into Set 1–3 TSV/XLSX pairs |
| `generate_heteroplasmy_report.py` | Write a batch-local called-heteroplasmy TSV |
| `convert_variants_to_tsv.py` | Convert JSON variant data to TSV format |
| `src/tools/fis/` | Reusable FIS Excel-to-canonical-Sample tool; invoked by the NGS runner |
| `merge_comparison.py` | Merge comparison results across multiple batches |
| `merge_json_files.py` | Merge multiple JSON result files |
| `variant_frequency_analyzer.py` | Analyze variant frequency across samples |

## Key Concepts

### FIS Analysis

The reusable `src.tools.fis` pipeline provides FIS-specific conversion:
- FIS Excel conversion to canonical Samples plus raw TSV, JSON, and FASTA artifacts
- Per-sample nomenclature, coverage, and QC metadata
- Stable tool and artifact outputs consumed by the staged comparison modules

The FIS tool writes stable NGS artifact filenames in its output directory, regardless of the input Excel filename. Compatibility FASTA files retain region order (`HV1`, `HV2`, `HV3`), prepend `GGGGT` to an HVS-II sequence starting at position 73, and exclude sequences with more than 90% `N` content:
- `FIS.tsv`
- `FIS_transformed.tsv`
- `FIS.json`
- `FIS.seq.tsv`
- `FIS_QCFail_samples.txt`
- `FASTA/`

The reusable [`FIS tool`](../../tools/fis.md) writes canonical `json/statistic_fullbatch.json` for comparison. Root-level `FIS.json` remains a raw per-position NGS artifact and `FIS_transformed.tsv` remains a report artifact. Canonical FIS variants and FIS flag reasons are calculated from `FIS_Nomenclature correction genotype`; blank, `-`, and `/` nomenclature values fall back to the raw `FIS_Variants` audit profile.

The NGS runner uses explicit comparison boundaries:

1. `src.tools.fis.pipeline` writes canonical `json/statistic_fullbatch.json` plus the established root-level FIS artifacts. `prepare_fis_comparison.py` consumes the canonical JSON, applies documented CE-to-Sanger ID changes from `results/modules/NGS/correct_mapping_sanger.tsv`, then uses the existing `results/modules/NGS/mapping.tsv` aliases as fallback; it writes comparison-scoped `comparison/FIS.json` and records every resolved FIS-to-Sanger match in `comparison/sample_mapping.tsv`.
2. `prepare_sanger_comparison.py` selects only mapped Sanger records and writes canonical `comparison/Sanger.json` using the same FIS comparison keys.
3. `src.core.comparison` writes intermediate standard JSON, TSV, Excel, and variant-level TSV outputs in the same `comparison/` folder using the centralized exact comparison rules.
4. `enrich_fis_sanger_comparison.py` joins those results with canonical FIS metadata and unmatched FIS samples into the official final sheets, `comparison/FIS_sanger_comparison_direct.tsv` and `comparison/FIS_sanger_comparison_direct.xlsx`.
5. `split_fis_sanger_comparison_sets.py` uses `results/modules/NGS/sets.tsv` to write the Set 1–3 TSV/XLSX report pairs in `comparison/`.
6. `generate_heteroplasmy_report.py` writes `Heteroplasmy.tsv` at the batch root from FIS calls in the variant-level comparison report and Sanger calls in the correction-resolved `comparison/Sanger.json`. Its `FIS_Sample` and `Sanger_Sample` columns distinguish the source FIS ID from the authoritative resolved Sanger ID. It validates authoritative IDs from `correct_mapping_sanger.tsv` against the prepared Sanger records, joins biallelic FIS depth data where available, appends complete matched raw FIS fields (including multiallelic depth text), assigns Group through the same longest-prefix `sets.tsv` lookup, and marks unreviewed truth as `U`. It also accepts a one-sample-ID-per-line TXT to write a separately filtered TSV. It does not infer uncalled minor-depth candidates from unavailable per-position data.

The final enriched sheet—not the intermediate core sheet—is the user-facing comparison output. It keeps only `FIS_Sample` as its sample identifier, omits replicate-specific columns, and joins `Nomenclature correction genotype`, `Nomenclature QC`, `Coverage pass rate(%)`, `Failed position`, and `# N in consensus seq` from `FIS_transformed.tsv`. `Concordant` remains the direct core result; `FIS_Correct` is the FIS-review refinement that accepts reference-covering IUPAC calls and 16193 indels when they are the only FIS-side difference. When an alias resolves through `mapping.tsv` but no Sanger record exists, the manifest preserves the canonical ID and the final sheet displays `<ID> (Not found)`. The NGS runner removes prior generated FIS outputs before each run so stale results are not retained.

### Variant Comparison

The active FIS-to-Sanger path filters both prepared profiles to the Sanger analyzed intervals, then delegates concordance and unique-variant decisions to `src/core/comparison.py`. Its exact matching and special-position behavior are therefore identical to other core comparisons. Uninformative FIS `N` calls are removed at the preparation boundary. `compare_fis_sanger_direct.py` retains the older standalone IUPAC-aware behavior only for legacy use.

The FIS preparation stage accepts `-t/--mapping-file` with these columns:
- `mtDNA`: canonical mtDNA ID used by Sanger results
- `STR`: corresponding STR/HID ID

Both IDs are treated as aliases. Run suffixes in FIS sample names are resolved using the longest matching alias.

Before exact/alias lookup, stage 2 requires `results/modules/NGS/correct_mapping_sanger.tsv`. It applies rows with a nonempty `ID_Change_Note` that match the current FIS batch, using `CE_Sample_ID`, `Sanger_Sample`, and `Sanger_Batch`. Corrections match the complete FIS sample ID, take precedence over normal mapping, and are validated against the merged Sanger JSON target ID and canonical batch (the source batch may carry an added descriptive suffix). The additive `Mapping_Source` field in `comparison/sample_mapping.tsv` records `correction`, `exact`, `alias`, or `unresolved`.

`split_fis_sanger_comparison_sets.py` accepts a headerless `sets.tsv` with `Sample_ID`, `Sanger_Batch`, and Set (`1`, `2`, or `3`) columns. It assigns each `FIS_Sample` through the longest matching complete sample-ID prefix; `Sanger_Batch` remains source provenance because an authoritative ID correction can target a different Sanger batch. It then writes `<comparison-stem>_set_1`, `_set_2`, and `_set_3` as TSV/XLSX pairs. Rows absent from `sets.tsv` are excluded without modifying the source report.

## Cross-References

- [core/comparison.md](../../core/comparison.md) — Core pairwise comparison logic
- [core/models.md](../../core/models.md) — Sample and Variant data models
