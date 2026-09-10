## Unreleased
### Added
- **tools/tracy**: Add default-on per-trace noise QC and masking for Tracy AB1 analysis. The pipeline evaluates overlapping peak-metric windows after `tracy decompose`, writes `preprocess/<sample_id>/qc_report.json` with ranges and exclusions, and removes variants plus the same spans from the originating trace's intervals before trace merging. A clean mate trace can retain coverage in the final LID profile. Set `MTDNA_TRACY_NOISE_MASK_ENABLED=false` to retain calls and intervals while continuing to write QC reports; the window thresholds and support rule are configurable with `MTDNA_TRACY_NOISE_*` variables.
- **fis**: Add reusable `src/tools/fis/` Excel-to-canonical-Sample tool pipeline with standard batch and region JSON output.
- **ngs**: Add `split_fis_sanger_comparison_sets.py` to split a FIS-to-Sanger TSV into Set 1–3 TSV/XLSX pairs using the headerless `sets.tsv` sample definitions and longest sample-ID-prefix matching; run it automatically after final FIS-to-Sanger enrichment.
- **ngs**: Add full and Set-1-autopass-filtered batch-local called-heteroplasmy TSVs generated from FIS depth, variant-level FIS calls, and correction-resolved prepared Sanger calls; preserve matched raw FIS fields for multiallelic calls, label source/resolved IDs as `FIS_Sample`/`Sanger_Sample`, and validate authoritative Sanger ID corrections.
- **statistics**: Add `src/modules/statistics/variant_frequencies.py` to calculate population-style variant frequencies from `merged_statistics.json` and write frequency-sorted variant/count/frequency TSV output.
- **statistics**: `src/modules/statistics/find_matching_profiles.py` now accepts a TXT file of sample IDs and writes matching-profile results as TSV columns `Sample ID`, `Intervals`, `Profile`, and `Matched Sample ID`.
- **tnls**: `src/modules/TNLS/verification/compare_TNLS_HCLS.py` now groups TNLS target samples into mtDNA-profile families (IUPAC-tolerant SNP equivalence over the common >=150 bp overlap; transitive union-find grouping) with sequential `FAM-NNN` IDs, and writes a new `stats_TNLS_HCLS_family_summary.tsv` (Family_ID, Num_Members, Member_IDs, Num_Matched_HCLS, Matched_HCLS_IDs, HV1/HV2/HV3_variants from the broadest-coverage representative). Family matches are aggregated from the existing per-pair `CANNOT_EXCLUDE` output; the mtDNA matching engine is unchanged and all existing outputs are preserved. Variant columns are reference only and are not part of matching.
- **fasta**: Add FASTA->Sample JSON variant caller (`Tool.FASTA`, single `FastaToSample` class in `src/fasta_to_sample.py`) with canonical rightmost-in-tandem-repeat indel normalization. Insertions are right-normalized per base for homopolymer (period-1) runs and as a block (re-ordered to the run phase) for dinucleotide (period-2) runs: polyC (HV2 303-315) `+C` uses the split-at-T convention (run1 -> `309.x`, run2 -> `315.x`; e.g. 2557844 -> `309.1` + `315.1`) and poly-CA (HV3 514-524) `+CA` lands at `524.1` + `524.2` (e.g. 2557921), both round-tripping the consensus FASTA faithfully and matching the regenerate/sequencher/tracy pipelines. A real IUPAC ambiguity code (`Y`/`R`/...) in the query is called as a SNP (seq = the code); a query `N` is unread (no variant). Site-specific nomenclature conventions that no alignment can infer (equivalent placements that all round-trip but differ in nomenclature) are encoded as an explicit, data-driven `REGION_CONVENTIONS` table applied after biological normalization — e.g. `16189 T>C` + `del 16193` (e.g. 2557851, never a `T` deletion); a polyC run shift `309 C>T` + `310 T>C` (or `308 C>T` + `del 310`) is named `del 309` + `ins 315.1` (e.g. MS_200426_005), not adjacent SNPs; and a `513 G`+`514 C` deletion is named `513 G>A` + `del 523` + `del 524` (the flank base is a SNP, the length change moves to the poly-CA run end), never a `G` deletion; add entries as data when a new convention is agreed.  Per-region `intervals` default to the config region span minus any query `N` positions (the *callable* range; an internal `N` splits sub-intervals), so `N` is preserved in the reconstructed HV strings and the interval matches the curated/regenerate convention (e.g. 2557858 HV3 `NNN` tail -> `[438,573]`). A SNP inside a tandem-repeat run that also contains a right-normalized deletion is re-anchored to the canonical rightmost-gap coordinate (shifted left by the deletion length) so the variant list round-trips the consensus (e.g. an `A>G` in the poly-CA run with a concurrent `CA` deletion lands at `519`, not `521`). Query bases are uppercased and record ids normalized so mixed/lowercase external consensus FASTAs import correctly.
- **comparison**: Add opt-in legacy summary output mode to `src/core/comparison.py` (`--legacy-format` CLI flag / `legacy_format` param on `compare_batches` and `compare_batch_files`). When enabled, the summary TSV, Excel, and JSON outputs reproduce the legacy 13-column layout (hardcoded `Pipeline`/`Sequencher` labels, single `Variants Count (Sequencher)`, merged `Flagged (Pipeline)`/`Flagged (Sequencher)` columns) as in `Batch_MS_280326_001.tsv`; the variant-level TSV is unchanged and the default 16-column layout is byte-identical when the flag is absent
- **sheets**: `src/modules/sheets/merger.py` adds an optional `-b/--batches` filter (comma-separated batch IDs) to `arg_parser`, `_load`, and `merge_excel_files`. When set, only the per-batch subfolders under `input_dir` whose name matches a listed batch ID are merged; `None`/omitted keeps the existing merge-everything behavior. `p.parent.name` is matched against the batch IDs since comparison output lives under `results/modules/comparison/<BATCH>/<BATCH>.xlsx`
- **pipeline**: `scripts/batch_pipeline.sh` merge step now passes the run's `BATCHES` array to the merger via `-b`, so `all_batches_comparison` consolidates only the batches processed in that run instead of every historical batch subfolder already present under `results/modules/comparison`
- **pipeline**: `scripts/batch_pipeline.sh --rerun` now refreshes `data/raw/${BATCH}.txt` from the metadata Excel `LID` column without running the full prepare step, builds rerun lists by parsing Sequencher ZIP member sample IDs with `parse_filename_ranges` (for example `LN_25_AA1630-73-340 ...` matches `LN_25_AA1630`), logs how many unique ZIP samples were found in the raw TXT, and writes `${BATCH}_zip_not_in_txt.txt` only for Sequencher ZIP TXT entries whose parsed sample ID is absent from the raw batch TXT.
- **metadata**: Add `src/modules/statistics/metadata_to_raw_txt.py` to generate one raw sample-ID TXT per metadata Excel workbook from the `LID` column, excluding `PC`/`NTC` controls; existing TXT files are compared against Excel and rewritten only when missing or mismatched, then verified.
- **tnls**: `stats_TNLS_HCLS_clean.tsv` (and `stats_TNLS_HCLS_with_metadata.tsv`) now include an `overlap_intervals` column with the actual `[start, end]` overlap regions (sourced from the per-pair `list_overlap`), so the overlapping base-pair ranges are written alongside the existing `overlap_bp` count. A `matched_variants` column was also added, listing the matched variants per pair in the simplified display form (e.g. `16093C 16189C 309.1C 249DEL`) via `src.core.variants.format_variants_simplified` (rebuilt from the typed per-pair `variant_match` keys; IUPAC codes preserved, insertions/deletions rendered as `<pos>DEL`/`<pos><seq>`). `HCLS_variants` and `TNLS_variants` columns were added showing each sample's complete variant profile (all regions, not only the overlap) in the same simplified form, built from the base/target sample JSONs via `format_variants_simplified`
### Changed
- **statistics**: Move the blinded test-data helpers to `src/modules/statistics/`; run them with `scripts/modules/statistics/test.sh`. The editable blind-ID mapping is now the ignored local file `data/modules/statistics/sample_mapping.tsv`, and generated test data is written to `temp/test`.
- **statistics**: Move and rename the data/results merger to `src/modules/statistics/merge.py`; its runner is `scripts/modules/statistics/run.sh`.
- **ngs**: Add `FIS_Correct` to the final FIS-to-Sanger report, preserving the core `Concordant` result while accepting reference-covering IUPAC calls and 16193 indels as FIS-review-correct when no other difference remains.
- **fis**: Write independent per-sample region JSON artifacts concurrently with 75% of the available CPU cores (minimum one worker).
- **ngs**: Apply documented global CE-to-Sanger ID corrections from `correct_mapping_sanger.tsv` before stage-2 comparison preparation; each full FIS sample ID takes precedence over exact/alias lookup and is validated against its target Sanger ID and canonical batch.
- **ngs**: Split the active FIS-to-Sanger workflow into analyzer JSON metadata, canonical FIS preparation, filtered Sanger preparation, shared core comparison, and final FIS enrichment boundaries; keep all prepared and standard comparison artifacts in one `comparison/` folder and retain the standalone direct comparator as legacy code.
- **ngs**: Add nomenclature genotype/QC, coverage pass rate, failed positions, and consensus-N count from `FIS_transformed.tsv` to the official enriched FIS-to-Sanger TSV and Excel sheets.
- **core/models**: Add `fis` and `sanger` tool provenance values for canonical comparison batches.
- **ngs**: Write FIS-to-Sanger comparison results as both TSV and Excel, and omit `Sample_ID` when it duplicates `FIS_Sample` for every result row.
- **ngs**: Replace separate exact/enhanced FIS-to-Sanger results with one direct comparison using only exact or same-position IUPAC-compatible alleles; remove alignment exceptions and `_Enhanced` columns while retaining `N` and established always-concordant position filtering.
- **ngs**: Format nomenclature-corrected FIS genotypes like `FIS_Variants` and use them for direct Sanger comparison.
- **ngs**: Prefix all FIS-derived columns in `FIS_transformed.tsv` with `FIS_` for clear FIS-to-Sanger comparisons.
- **ngs**: Add R2 T1 nomenclature and coverage QC fields to `FIS_transformed.tsv` and use the actual `Coverage pass rate(%)` column for QC flagging.
- **ngs**: Clean generated FIS outputs before each NGS runner execution to prevent stale results.
- **ngs**: Standardize NGS output filenames with the `FIS_` prefix, including QC-fail and FIS-to-Sanger comparison results.
- **ngs**: Resolve FIS-to-Sanger sample IDs through the compact MT-to-STR/HID mapping format, including run-suffixed FIS names.
- **ngs**: Give FIS analyzer outputs stable `FIS`-based filenames instead of deriving them from the input Excel filename.
- **tools/blastn**: Remove the unused per-trace BLASTN `-outfmt 3` preprocessing alignment; the AB1 path now runs only the concatenated BLASTN alignment consumed by ETL, reducing one BLASTN subprocess per trace.
- **tnls**: TNLS/HCLS verification now requires at least one matched variant in the overlap before reporting a pair as `CANNOT_EXCLUDE`; no-mismatch pairs with zero matched variants are treated as `INCONCLUSIVE` and are omitted from match outputs.
- **tnls**: Format family-summary `Member_IDs` with a space after each comma (for example, `sample1, sample2`).
- **tnls**: Add `Profile_Frequency` to the family summary as a fraction such as `15/38721`, plus `Profile_Variants` for the complete representative profile.
- **tnls**: `scripts/modules/TNLS/verification.sh` now writes generated HCLS sample output under `results/modules/HCLS/<batch>/`, reads `HCLS_SAMPLES` from that batch-specific HCLS results directory, reads TNLS metadata from `data/TNLS/TNLS_metadata.xlsx`, cleans only TNLS/HCLS verification outputs, and logs progress through Loguru.
- **tnls**: Rename the matched-variant count columns back to their original names: `matched_variant_sites`->`matched_variants` (count of matched variant loci) and `matched_variant_bp`->`matched_bases` (total matched variant sequence characters) in `output_compare_tnls_hcls.tsv`, `stats_TNLS_HCLS_with_metadata.tsv`, and `stats_TNLS_HCLS_clean.tsv`. The redundant per-pair `matched_variants` string column is dropped (it was identical to the simplified `variant_match`); the matched-variant list now lives in `variant_match`, and `with_metadata`/`clean` carry it as a `variant_match` column. The family summary keeps its prefixed `HCLS_matched_variant_sites`/`HCLS_matched_variant_bp` count columns alongside the `HCLS_matched_variants` union string
- **tnls**: Merge `stats_TNLS_HCLS_matching.tsv` into `stats_TNLS_HCLS_family_summary.tsv`. Column order is now `Family_ID`, `Member_IDs`, `Matched_HCLS_IDs`, `HCLS_matched_variants`, `HCLS_matched_variant_sites`, `HCLS_matched_variant_bp`, `HCLS_ratios`, `HV1_variants`, `HV2_variants`, `HV3_variants`. The per-HCLS columns are comma-joined and aligned with `Matched_HCLS_IDs`, empty for families with no matched HCLS; `HCLS_matched_variants` is the distinct union of matched variants across the HCLS's matched TNLS pairs in the simplified display form (e.g. `16093C 16189C 309.1C`), so a 1:2 match against two identical TNLS reports the variants once. Redundant columns are dropped: the per-HCLS matched-TNLS list (the matched TNLS are the family `Member_IDs`; exact HCLS->TNLS pairs stay in the per-pair TSV), the derivable member/matched-HCLS counts (`Num_Members`, `Num_Matched_HCLS`), and the per-HCLS `HCLS_exclusive` flag (TNLS-level exclusivity, only meaningful for 1:1 pairs). `stats_TNLS_HCLS_matching.tsv` is no longer written. The `HV1/HV2/HV3_variants` reference columns are now space-joined (e.g. `16093C 16189C 309.1C`) instead of comma-joined, matching the other simplified variant columns. The shared per-HCLS aggregation is extracted into `_hcls_match_summary` so the two views stay consistent
- **tnls**: `output_compare_tnls_hcls.tsv` now stores `variant_match` and `variant_mismatch` in the simplified display form (matching `matched_variants`/`HCLS_variants`/`TNLS_variants`) instead of Python dict literals. `variant_match` is the space-separated matched variants (e.g. `16093C 16189C 16189C`); `variant_mismatch` is `<pos><base>><target>` per mismatched site (`DEL` for a deletion allele, `·` for a one-sided/absent allele), e.g. `16172C>A 16223T>·`, or `None` when there are no mismatches. Downstream stats now parse the simplified `matched_variants` column via `src.core.variants.parse_variant_position_allele` instead of `ast.literal_eval`-ing the dict literal
- **modules**: Move all `src/modules/PYQG/` logic into `src/modules/TNLS/`
- **modules**: Deduplicate `src/modules/TNLS/add_variants_to_pyqg.py` against `src/core/variants.py`.
- **modules**: `src/modules/TNLS/add_variants_to_pyqg.py` CLI cleanup
- **modules**: `src/modules/TNLS/add_variants_to_metadata.py` CLI cleanup
- **modules**: `src/modules/TNLS/blind_copy.py` always splits FASTA
- **modules**: `src/modules/TNLS/add_variants_to_samples.py` now emits
- **modules**: `src/modules/TNLS/add_variants_to_samples.py` — drop the `is_position_in_region` wrapper and call `src.core.variants.is_position_in_intervals` directly in `get_variants_for_region`
 an `Analyzed_Range` column with the sample's analyzed intervals (from the JSON `intervals` field), reusing the shared `get_intervals_for_sample` helper from `src/modules/TNLS/add_variants_to_metadata.py` instead of duplicating the interval formatter. Output columns are now `Sample_ID, Barcode, Variants, HV1, HV2, HV3, Analyzed_Range`
 files by HV region — removed the `--split-fasta`/`--no-split-fasta` CLI flags and the `process_samples(..., split_fasta=...)` parameter; the non-splitting copy branch is now dead code removed. No behavior change (splitting was already the default)
 — drop the `-e/--output_excel` flag; the Excel output is now always written alongside the TSV in the output directory (same stem, `.xlsx` suffix). Long flags renamed to hyphenated form (`--input-file`, `--json-file`, `--output-file`, `--stats-file`) to match the repo convention; short flags (`-i`, `-j`, `-o`, `-s`) and `args.*` attribute access are unchanged

- **modules**: Rename the moved TNLS scripts
- **pipeline**: `scripts/modules/TNLS/prepare.sh` now writes all TNLS prepare outputs (blinded copies, metrics JSON, and the `TNLS_with_variants_<YYYYMMDD>.tsv`/`.xlsx` annotation) under `results/modules/TNLS/` instead of `data/TNLS/`, matching the repo `results/modules/<module>/` convention
 to reflect their logic and drop the "PYQG" label — `src/modules/TNLS/PYQG.py`→`blind_copy.py`, `add_variants_to_pyqg.py`→`add_variants_to_samples.py`, `PYQG.tsv`→`sample_mapping.tsv`. Internal `load_pyqg_data`/`process_pyqg_data`/`pyqg_df` renamed to `load_sample_data`/`process_sample_data`/`sample_df`. `scripts/modules/TNLS/prepare.sh` now invokes `src.modules.TNLS.blind_copy` and `src.modules.TNLS.add_variants_to_samples`, writes outputs under `data/TNLS/`, and names the variant-annotated output `TNLS_with_variants_<YYYYMMDD>.tsv` (Excel written alongside as `.xlsx`); default mapping path is `sample_mapping.tsv`
 — drop the `-e/--output_excel` flag; the Excel output is now always written alongside the TSV in the output directory (same stem, `.xlsx` suffix). Long flags renamed to hyphenated form (`--input-file`, `--json-file`, `--output-file`) to match the repo convention; short flags (`-i`, `-j`, `-o`) and `args.*` attribute access are unchanged
 `format_variants_standard` now delegates to the canonical `format_variants_simplified` (preserving the PYQG `"No variants"` empty sentinel) and `is_position_in_region` delegates to `is_position_in_intervals`; the duplicated local formatting/interval logic is removed. No output change
 (`PYQG.py`, `add_variants_to_pyqg.py`, `PYQG.tsv`). `scripts/modules/TNLS/prepare.sh` now invokes `python -m src.modules.TNLS.PYQG` / `src.modules.TNLS.add_variants_to_pyqg` and reads `src/modules/TNLS/PYQG.tsv`; `docs/modules/PYQG/` folded into `docs/modules/TNLS/README.md` and `docs/ARCHITECTURE.md`. No public API or behavior change
- **sheets**: Move `src/core/uploader.py` and `src/core/merger.py` to `src/modules/sheets/` (new `sheets` module folder under `src/modules/`). No active Python code imports them; only `scripts/batch_pipeline.sh` invokes them via CLI
- **sheets**: `src/modules/sheets/merger.py` `_load` now recursively globs per-batch `*.xlsx` files under the input dir (excluding `all_batches*` merged outputs and `~$` temp files) instead of the old flat `*_sequencher_comparison.xlsx` glob, matching the new `results/modules/comparison/<BATCH>/<BATCH>.xlsx` layout
- **pipeline**: `scripts/batch_pipeline.sh` Google Sheets upload/merge paths updated — individual batch comparison files now read from `${RESULTS_DIR}/modules/comparison/${BATCH}/${BATCH}.xlsx` (was `${RESULTS_DIR}/modules/sequencher/${BATCH}_sequencher_comparison.xlsx`); merge output now `${RESULTS_DIR}/modules/comparison/all_batches_comparison`; uploader/merger invoked as `python src/modules/sheets/uploader.py`/`merger.py`
- **pipeline**: Introduce a regenerate gate (`src/core/regenerate.py`) as the boundary between raw tool output and the final JSON consumed by downstream tasks. `scripts/pipeline.sh` now filters each tool's raw JSON (`results/tools/<tool>/json`) into the regenerate (final) directory via `run_regenerate` (replacing the old `copy_to_regenerate` plain copy). Variants are scoped to the sequenced intervals (`sample.intervals`), so out-of-range/additional variants (e.g. an insertion called inside a coverage gap like 303-315) are excluded from the final JSON and from FASTA/report output. The comparison step now reads the RAW tool JSON directly (not the filtered regenerate copy), so the reviewer still sees out-of-range/additional variants in the comparison. FASTA (`generate_fasta`) and report (`generate_reports`) generation continue to read the regenerate (final) directory, which is now filtered. `generate_sequence` (`src/core/sample.py`) is unchanged — it faithfully applies the variants it is given, so each output is self-consistent: the raw tool JSON reflects all called variants (variant list and HV sequence agree), while the filtered regenerate/final JSON reflects only the in-range variants. Filtering lives solely at the regenerate gate, the single source of truth
- **pipeline**: `scripts/pipeline.sh` no longer runs `src.manual_pipeline.py` in Step 4 (its comparison/regenerate outputs were unused downstream). Step 4 now only extracts manual data, runs the Sequencher ETL, and runs `src.core.comparison`. Step 6 (`src.generate_reports -r`) now sources per-sample JSONs from `FASTA_INPUT_DIR` (the Sequencher `json/` output, same input as Step 5 FASTA) instead of the former `MANUAL_REGENERATE_DIR`
### Removed
- **tnls**: Remove unused `src/modules/TNLS/prepare.py`.
- **ngs**: Remove the legacy `src.modules.NGS.fis_analyzer` API and positional CLI; use `python -m src.tools.fis.pipeline` for FIS workbook processing.
- **ngs**: Remove `FIS_Replicate`, `Replicate`, and `Replicate_Concordance` generation and comparison logic from FIS analysis and FIS-to-Sanger outputs.
- **pipeline**: Move unused legacy root modules `src/automate_pipeline.py`, `src/compare_automate_pipeline.py`, and `src/write_report.py` to `src/backup/`. No active code imports or invokes them (only already-backed-up `src/backup/process_pc_samples.py`/`src/backup/blastn.py` referenced them); `src/manual_pipeline.py` is retained in `src/` for standalone use but is no longer invoked by `scripts/pipeline.sh` (see Changed below)
### Fixed
- **tools/tracy**: Use shared `src.core.polyc` directional policy to exclude called HV2F/HV3F candidates at `>=304`, HV2R/HV3R candidates at `<=315`, and 16189-C-triggered HV1 artifacts before same-LID consensus, while retaining the trace and recording each removal with trace, quality, and centralized reason.
- **comparison**: Treat an insertion immediately after an analyzed interval endpoint as outside coverage (for example, `16193.1C` is excluded from `16024-16193`), preventing false discordance in FIS/Sanger comparisons.
- **fis**: Build canonical FIS variants and flag reasons from the nomenclature-corrected profile with rCRS-backed references; fall back to raw calls when the corrected profile is absent or uses the FIS `-`/`/` sentinels. Restore valid, ordered compatibility FASTA output, HV2 prefixing, and sequence QC.
- **ngs**: Preserve a canonical Sanger ID resolved from `mapping.tsv` even when the current Sanger JSON has no matching record, so suffix-bearing FIS IDs such as `243311_HID_45` report `243311 (Not found)` instead of losing the mapping.
- **pipeline**: Make batch reruns require and refresh from the exact canonical metadata workbook (`DATA_DIR/metadata/<batch>.xlsx`) before selecting samples, failing clearly instead of continuing with a stale raw sample list; Step 3 continues to use the supplied rerun sample list while reusing AB1 files prepared by the normal run.
- **tnls**: Speed up `merge_data.py` unchanged reruns by caching FASTA batch scans, include FASTA files from `results/archive`, correctly reusing regenerate scan entries, reusing fresh merged statistics output, and removing the unused PDF/report status check from batch statistics.
- **sequencher**: Reject every TXT for a duplicated Sequencher LID during pre-ETL QC, preventing same-LID files from being merged or overwriting per-region outputs.
- **sequencher**: Fail pre-ETL QC for TXT files missing the required Sequencher variant table header (`Pos`, `Seq`, `Con`, `Required Edit`), so non-variant/report-style TXT inputs are skipped instead of producing misleading manual outputs.
- **fasta**: Fix the post-generation FASTA validation report ignoring the pipeline `--output-dir`. `src/generate_fasta.py` now accepts `--results-dir` (priority: explicit arg > `RESULTS_DIR` env > `results`) and `scripts/pipeline.sh` passes `--results-dir "${RESULTS_DIR}"`, so the report lands under the same results root as the rest of the pipeline outputs (`<results-dir>/validation/validation_fasta_<batch>.txt`) instead of always falling back to the `.env` `RESULTS_DIR` default
- **sheets**: Fix merged comparison output losing data across batches while preserving legacy headers. `src/modules/sheets/merger.py` now normalizes tool-parameterized column names (e.g. `Variants (tracy)`/`Variants (sequencher)`) or generic `Variants (A)`/`Variants (B)` names to `Pipeline`/`Sequencher` before concatenating, so batches using different automate tools (tracy vs blastn) align without changing merged headers to `(A)`/`(B)`. `_VARIANT_COLUMNS` uses the legacy `Pipeline`/`Sequencher` names so variant-separator normalization applies to the merged output
- **tools/blastn**: Fix missing per-region intermediate JSON output. `process_batch` now writes region JSONs from the **combined** per-LID samples (after `_build_combined_samples`) using the data-driven `REGION_TO_KEYS` groups (`HV1`, `HV2-3`), mirroring Tracy/Sequencher's `_write_full_region_jsons`. Previously it iterated the uncombined per-amplicon `processed` items and derived region groups via `region_groups_from_flags(sample.sample_flags)` — but BLASTN never sets auto-pass flags (`Autopass HV1` / `Autopass HV2 and HV3`, used only by Mutation Surveyor), so `region_groups_from_flags` always returned `{}` and `regions/` was left empty. Each LID now produces one `HV1_{LID}.json` and one `HV2-3_{LID}.json` (in `regions/HV1/` and `regions/HV2-3/`) with merged forward/reverse primer data, plus range-QC warnings via `validate_region_group_intervals`
- **tools/blastn**: Close the deferred `PostProcess` parity gaps so `src/tools/blastn/etl.py` now reproduces the legacy `Analysis → PostProcess.post_process()` final variant set, not just the consensus set. Ported into `src/tools/blastn/transform.py`: `_handle_polyc_removal` (post-hoc removal of low-quality variants in/after polyC regions, keeping 309.1C/309.2C and 16189T>C), `_apply_region_validation_rules` (remove non-native-region low-quality variants), and `_polyc_warning_flags` (emit polyC/after-polyC/16189T>C warnings as `variant_flags` strings, since the standardized `Variant` model has no metadata fields and `core/flagging` is the sole flagging authority). `_filter_positions` now applies the legacy `is_special_position` 309.x dynamic rule. `etl.process` runs the removal passes in the legacy `post_process` order (handle_polyC → filter → conversion → region_validation → warnings). `_compute_consensus_variants` and `_build_variant_objects` now use the `file` key consistently so `Variant.files` propagate (previously dropped by a `files`/`file` key mismatch) and quality is built as `list[int]` to match `Variant.quality`. Final-output `(pos,ref,seq)` parity is locked by new `tests/test_blastn_etl_parity.py` cases (T11-T16, P3) against real fixtures. The new module also fixes a pre-existing legacy crash: `Analysis.concensus_variant` stored scalar quality for single-alignment variants, so `PostProcess.apply_region_validation_rules` raised `max(int)` on ~6% of real samples (12/187); the new module always stores quality as a list and `_max_quality` handles both, so those samples now produce output instead of being skipped
- **tools/blastn**: Port legacy BLASTn variant-calling algorithm logic into the standardized ETL so `src/tools/blastn/etl.py` produces the same consensus variant set as legacy `Analysis.concensus_variant()`. `_call_row_variants` now calls SNPs, insertions (`.1` positions, `ref == "-"`) and deletions (`seq == "-"`) with insertion-aware position arithmetic; `_right_norm_indel` (ported from `Analysis.right_norm_indel`) groups and right-normalizes consecutive indels; and `_compute_consensus_variants` applies the legacy strand-aware polyC `flagged_for_filtering` consensus filter (remove heterozygous single-strand and polyC-single-strand-flagged variants; retain homozygous-both-strand and non-polyC heterozygous variants for `Heteroplasmy` flagging). Consensus-level parity (`pos, ref, seq`) is locked by `tests/test_blastn_etl_parity.py` against real fixtures. Out-of-scope `PostProcess` passes (`handle_polyC_regions` post-hoc removal, `apply_region_validation_rules`, `warning_variants`) remain deferred follow-ups
- **blastn**: Fix the legacy `src/analysis.py` path (`scripts/pipeline.sh -p blastn` → `python -m src.blastn`) crash that the modular `src/tools/blastn` ETL had already worked around. `Analysis.concensus_variant` stored a scalar `int` in `variant["quality"]` for single-alignment variants, so `PostProcess.apply_region_validation_rules`/`handle_polyC_regions` raised `TypeError: 'int' object is not iterable` at `max(quality)` for non-native-region variants with one trace; the per-sample exception skipped `<LID>.json`, which then surfaced as a `FileNotFoundError` in `WriteBatchReport.write_files`. `quality` is now always stored as a `list[int]` to match the `Variant.quality` contract, so affected samples produce output instead of crashing the batch
- **tnls**: Fix `stats_TNLS_HCLS_matching.tsv` double-counting matched variants for 1:N HCLS matches. `matched_variant_sites`/`matched_variant_bp` are now the distinct union of matched variant keys across all matched TNLS pairs (sourced from the per-pair `variant_match`), so an HCLS matching two TNLS with identical variants reports 5, not 10. Falls back to per-pair sums only when `variant_match` is unavailable
- **pipeline**: Fix Positive Control validation (Step 1) failing with `Source directory not found: <LAB_DATA_DIR>/<batch>` in isolated/local runs (`-d/--data-dir`). `src/modules/quality_control/pc_ntc.py` no longer references `LAB_DATA_DIR`; it sources its AB1 files directly from the `--data-dir` value (the batch directory containing the HV subdirectories), so it works with relative or absolute paths. `scripts/pipeline.sh` now passes the user-provided data dir for an isolated run and `LAB_DATA_DIR/<batch_id>` for a normal NAS run. The copied-control working area (`pc_raw_dir`) moved under `<results>/validation/<batch_id>/raw` so it never writes into the AB1 source, and the copy `find` excludes the destination as a safety net
### Added
- **tools/blastn**: `process_batch` now writes `Batch.write()` flat into the `json/` subdirectory (`nest_batch_id=False`) to match the Tracy/Sequencher reference layout, avoiding a duplicated `<batch_id>` path segment; output lands directly under `results/tools/blastn/<batch_id>/json/` as `statistic_fullbatch.json` plus per-LID `<LID>/<LID>.json`, with per-region JSONs in `regions/` and preprocessing artifacts in `preprocess/`
- **flagging**: Flag the 309 deletion (`309DEL`) as `Deletion at 309` (Level 3). A deletion at the base position 309 is now flagged even though 309 is a special polyC position whose other non-insertion variants remain skipped; this catches the real variant event produced by the Tracy `309C>T` + `310T>C` conversion. Deletions at the other special positions (455, 463, 573) and at sub-positions such as `309.1` are still skipped.
- **tools/blastn**: Add FASTA-input batch pipeline `src/tools/blastn/fasta_pipeline.py` (`process_batch_fasta` + `FastaBatchOptions`, CLI `python -m src.tools.blastn.fasta_pipeline -i -o -r [--batch-id]`) for pre-assembled FASTA files (one sample per file, ID = filename stem). It runs a single BLASTN alignment per FASTA (outfmt 7, same columns as the AB1 pipeline) then reuses the shared ETL (`etl.process`), combine helpers, `Batch.write()`, and `write_region_jsons()`, producing the same `json/` + `regions/` + `preprocess/` output layout as the AB1 pipeline. FASTA carries no per-base quality, so no `.concatenate.quality` is produced and the ETL uses empty quality scores (matching the legacy FASTA path).
### Fixed
- **tnls/verification**: Normalize variant positions, ignore unresolved HCLS `N` SNP calls, and compare other SNP ambiguity codes by IUPAC compatibility while keeping indel matching exact.
- **tnls**: Include the batch ID when logging sample IDs missing from merged JSON statistics, making cross-batch FASTA/statistics discrepancies traceable.
- **tnls**: Report merged TSV/JSON discrepancies using unique sample IDs and distinguish duplicate TSV rows from JSON-only samples.
- **tnls**: Preserve sample identifiers as strings when adding variants, preventing valid numeric IDs from being reported as missing from `merged.json`; region statistics now count missing data correctly.
- **tools**: Fix cross-region sample-level flags (e.g. "Complex: 459 deletion with 16192 issue") missing from combined per-LID samples. `_combine_same_tool_samples` in Tracy, Sequencher, and BLASTN pipelines unioned per-region `sample_flags` without re-running `SampleFlagger.analyze()` on the merged variant set, so combined conditions spanning HV1 and HV2-3 (459DEL in HV2-3 + 16192T in HV1) were never detected when regions were processed separately. Mutation Surveyor ETL had a related gap: it called only per-variant `flag_variants()` and never ran `SampleFlagger.analyze()`, so all core sample-level flags were missing. All four tools now re-run `SampleFlagger.analyze()` on the full variant set so cross-region and sample-level flags are consistently detected. Additionally, `sample_to_dict` now recomputes core sample-level flags against each JSON's own variant set via `recompute_sample_flags`, so per-region JSONs no longer inherit stale cross-region flags from the combined sample (e.g. an HV1-only JSON no longer carries the 459DEL+16192T Complex flag when it lacks 459DEL). This replaces the narrower `recompute_no_315_1_flag` call and subsumes the "No 315.1 variant" absence-flag recompute
- **tools/tracy**: Stop duplicating the batch ID under the JSON output — `Batch.write()` now writes flat into `results/tools/tracy/<batch_id>/json/` so `statistic_fullbatch.json` and per-LID `<LID>/<LID>.json` sit directly under `json/` with no nested `<batch_id>` folder; `scripts/pipeline.sh` points `AUTOMATE_REGENERATE_DIR` at `${TRACY_DIR}/json` accordingly

- **tools/tracy**: Close behavioral gaps so the modular Tracy package reproduces legacy `src/tracy.py` variant/interval output — `_combine_same_tool_samples` now deduplicates by `(pos, ref, seq)` and aggregates per-file `files`/`peaks`/`quality` (matching legacy `merge_variants` `groupby`), re-merges overlapping region intervals, `_build_variant_objects` populates `Variant.files` from the source trace, and `_process_sample_json_files` filters Tracy decompose JSONs by the `HV` primer tag in the filename (also excluding any merged result JSON). Parity is locked by `tests/test_tracy_parity.py` against a real batch
- **tools/tracy**: Close remaining `statistic_fullbatch.json` gaps against the legacy regenerate output — `variant_to_dict` now emits the legacy key order (`peaks` before `quality`) and coerces peak heights/quality to `int` (matching legacy `int(q)` and integer peak data); `sample_to_dict` deduplicates conflicting calls at the same position keeping the last by `(pos, ref, seq)` (mirroring legacy `regenerate_result` → `get_variants_from_json`), resolving forward/reverse strand disagreements (e.g. position 489); `_process_sample_json_files` preserves filesystem (`os.listdir`) order instead of sorting so merged per-file `files`/`peaks`/`quality` provenance matches legacy `Tracy.process_sample()`. New `sample_flags`/`variant_flags`/`information` fields are retained

- **flagging**: Suppress the "No 315.1 variant" flag when position 315 (HV2 region) is not covered by the sample's intervals; `SampleFlagger` now accepts an optional `intervals` argument (checked via `is_position_in_intervals`) so HV1-only / partial-coverage samples are no longer wrongly flagged
- **sequencher**: Constrain flagging intervals to the regions actually sequenced per the filename region prefix (`HV1_` / `HV2-3_`), so a region-prefixed "FULL REGION" file no longer reports coverage for unsequenced regions (which previously caused "No 315.1 variant" to fire for HV1-only samples)
- **flagging**: Add `deduplicate_sample_flags()` and apply it across Sequencher, Tracy, BLASTn, and Mutation Surveyor ETL so per-variant flags (e.g. "Deletion at 16193") are kept only in `variant_flags` and removed from `sample_flags`; the aggregated "16180-16193 region (pos1, pos2)" summary is also dropped from `sample_flags` when the per-variant "16180-16193 region" flag is present (it is consolidated into the same form in the Variant Flags column)
- **comparison**: `ToolResult.flagged` now reflects any flag (`sample_flags` or `variant_flags`) so the comparison "Flag" column stays "Yes" when a sample has only per-variant flags; the "Sample Flags" column now reflects sample-level flags only
- **sequencher**: Retain the "Range - Possibly wrong" flag when `validate_analysis_ranges` rejects ranges that fail validation rules (out-of-bounds, unsorted, or cross-region); previously the flag was discarded and the sample silently fell back to "FULL REGION" with no flag. Also deduplicate the flag when both a filename format error and a range-validation error produce the same "Range - Possibly wrong" string
- **sequencher**: `parse_filename_ranges` now robustly extracts the sample ID (LID) and all interval pairs when the LID-to-range separator is a space, dot, or underscore instead of the rule-#19 hyphen (e.g. `LN_26_AA4522 73-340 438-573 16024-16365.TXT`); previously the first range start was glued to the LID (yielding a wrong `sample_id` and a dropped first pair), so the sample vanished from the comparison. Such filenames are still flagged "Range - Possibly wrong" so they display like other rule-#19 violations

- **flagging**: The "Consecutive indels (5+)" flag now requires truly consecutive indels of the same type: deletions at consecutive integer base positions (e.g. 100DEL 101DEL 102DEL 103DEL 104DEL) or insertions with consecutive indices at the same base position (e.g. 100.1C 100.2C 100.3C 100.4C 100.5C); previously any N+ indels merely adjacent in the sorted variant list triggered the flag
- **flagging**: Rework the "No 315.1 variant" absence flag (position 315, HV2 polyC stretch 303-315). Unlike per-variant flags, this is an *absence* flag (it reports a variant missing from the list), so it is now derived from each output JSON's own variant set at write time in `sample_to_dict` (`src/core/sample.py`) via the new `recompute_no_315_1_flag()` (`src/core/flagging.py`), instead of being computed per region file during ETL and union-merged. This makes every output JSON self-consistent: the combined 3-region JSON is flagged against its full variant set, while each per-region JSON (HV1 / HV2-3) is flagged against only that region's variants — so an HV1-only region JSON is no longer wrongly flagged. Gating is region-level (not exact-position): the flag fires when 315.1C is absent AND the region containing position 315 (HV2) was analyzed (non-empty coverage intervals), so a full-region sample that left a gap over 303-315 (a data problem) is still flagged (e.g. LN_25_AA3838 manual range `73-302 316-340 438-576 16024-16365`) while an HV1-only / HV2-not-sequenced sample is left alone. This also fixes false positives where a per-file read covered 315 but did not call 315.1C left a stale flag even though another read detected 315.1C (e.g. LN_25_AA7617 showed "No 315.1 variant" while its variant set contained 315.1C); the FIS/NGS path (intervals=None) keeps full-coverage legacy behavior
- **flagging**: Remove the redundant "Variant at position 310/460" per-variant/sample reason. It was produced by a `pos in ["310", "460"]` check that only matched *string* positions (the FIS path, where `convert_genotypes_to_variants` keeps positions as strings) and never matched the normalized *int* positions used by the Tracy/Sequencher/BLASTn ETL, so the same variant reported the extra reason in one path and not the other. The consistently-produced `"Has 310 variant"` / `"Has 460 variant"` reason (already emitted for both int and string positions via `has_variant_at`) now also drives the complex/technical classification (`COMPLEX_PATTERNS` and the 460 complex check), so 310/460 flagging and classification are identical across all ETL paths; flag levels (analyze and classify) are unchanged. The dead, never-called `VariantAnalyzer.sorted_variants` property and `_sorted` cache are also removed

- **comparison**: The "Variant Flags (A)/(B)" summary column no longer prefixes each flag with its position (e.g. now "Has 521A-G variant; Heteroplasmy at 16129" instead of "521: Has 521A-G variant; 16129: Heteroplasmy at 16129"); each reason already embeds its position, so the prefix was redundant

### Changed
- **tools/sequencher**: Adopt the standard `OutputPaths` folder layout (`json/`, `regions/`) like Tracy/BLASTn — batch JSONs now write flat to `results/tools/sequencher/<batch_id>/json/` and per-region intermediate JSONs to `.../regions/`; `--json-dir` still overrides the region directory for the unify workflow
- **core/batch**: `Batch.write()` gains a keyword-only `nest_batch_id` flag (default `True`, backward compatible); pass `nest_batch_id=False` when the output dir is already batch-specific to avoid a duplicated `<batch_id>` path segment
- **scripts/pipeline.sh**: Sequencher is invoked with a batch-specific `--output-dir` (`${SEQUENCHER_RESULTS_DIR}/${BATCH_ID}`) and the comparison reads `${SEQUENCHER_RESULTS_DIR}/${BATCH_ID}/json/statistic_fullbatch.json`
- **scripts/pipeline.sh**: Wipe the Sequencher batch results folder before the ETL run (matching the BLASTn/Tracy steps) so the output keeps the clean `OutputPaths` layout (`json/` + `regions/`) with no stale legacy files at the batch root
- **scripts/pipeline.sh**: Step 5 (FASTA) now reads the Sequencher standard Sample JSON output (`${RESULTS_DIR}/tools/sequencher/<batch_id>/json`) instead of the manual regenerate directory. `src/generate_fasta` loads each per-sample JSON as a standard `Sample` and always regenerates HV1/HV2/HV3 consensus sequences via `generate_sequence` (canonical flow); `src/validation/validation_fasta` replays `generate_sequence` for integrity checking. `src/core/sample.py` adds `load_sample()` (single per-sample JSON -> Sample) and `_normalize_keyed_sample()` now also flattens region-grouped variants (`{HV1,HV2,HV3}`)
- **sequencher/utils**: `_split_lid_and_ranges` now tries the standard format `LID-a1-b1 a2-b2 ... an-bn` first (split on first hyphen, validate against `VALID_INTERVAL_PATTERN`), falling back to the robust overlapping scan only for non-standard separators. Fixes LID truncation when trailing LID digits formed a valid range pair (e.g. `LN_26_AA0016-73-...` was parsed as `LN_26_AA`). `LID_RANGE_PAIR_RE` lookbehind tightened from `(?<!\d)` to `(?<![A-Za-z0-9])` for the fallback path
- **comparison**: Comparison summary, variant-level, and Excel column headers now name the actual source tool per batch (e.g. `Variants (sequencher)`, `Variant Flags (mutation_surveyor)`) instead of the generic `(A)`/`(B)`; a batch without tool attribution displays `Unknown`
- **core/models**: Add `Tool.UNKNOWN` sentinel; `load_sample_batch()` and `read_region_json()` now mark a missing or unrecognized `source_tool` as `Tool.UNKNOWN` (displayed as `Unknown`) instead of defaulting to `sequencher`
- **core/flagging**: Change SNP variant flag format from `Has {pos}{ref}-SNP variant` to `Has {pos}{ref}-{seq} variant` (e.g. `Has 519A-G variant` instead of `Has 521A-SNP variant`); change insertion flag format from `Insertion at {pos}` to `Insertion {seq} at {pos}` (e.g. `Insertion Y at 573.1` instead of `Insertion at 573.1`)

- **tools**: Move Tracy, BLASTn, Mutation Surveyor, and Sequencher tool modules from `src/modules/` to `src/tools/` — update all imports and references

- **docs**: Sync all documentation with current codebase — update API signatures, type annotations, and behavioral descriptions to match actual code
- **docs/core/batch.md**: Fix Batch processing steps to reflect that `sample_to_dict()` (not Batch directly) calls `generate_sequence()`, `statistic_variants()`, and `flag_variants()`; update output path to include `batch_id` subdirectory
- **docs/core/sample.md**: Update `sample_to_dict()` to reference `generate_sequence()` instead of `gen_seq_from_variants_regions()`; add `flag_variants()` call in serialization steps; fix default `ref_path` parameter
- **docs/core/mtdna_merger.md**: Add `build_merged_sample()` and `merge_profiles()` methods to public API; add `_load_and_flatten_profiles()` and `_collect_flags_and_info()` private methods
- **docs/core/models.md**: Update type annotations from `Union`/`Optional`/`List`/`Dict` to modern `|` / `list` / `dict` syntax; fix `validate_position` signature to `float | str`; fix `Tool` to `StrEnum`; fix `peaks` type to `list[list[float | None]] | None`
- **docs/core/comparison.md**: Fix `find_unique_variants` signature to use keyword-only `ignore_special`; update type annotations
- **docs/core/flagging.md**: Update type annotations to modern Python syntax
- **docs/core/region.md**: Update type annotations to modern Python syntax
- **docs/generation.md**: Fix `FastaGenerator` description to accurately state it uses `gen_seq_from_variants_regions()` (not `generate_sequence()`); fix `WriteReport` constructor parameter names (`processed_blastn_path`, `sub_outdir`); fix `ReportGenerator.get_pdf_link()` signature to include `session`, `lid`, `json_path` params; fix `generate_reports()` signature
- **docs/processing.md**: Update `PostProcess` processing step numbering to include `apply_region_validation_rules()` (step 3); fix type annotation formatting
- **docs/ARCHITECTURE.md**: Update Batch step 5 to reference `generate_sequence()` via `sample_to_dict()` instead of direct `gen_seq_from_variants_regions()` call; update ETL contract; update remaining callers table
- **docs/variants.md**: Update Sequence Generation Contract to reference `generate_sequence()` as the canonical API; clarify `gen_seq_from_variants_regions()` is legacy; fix type annotations
- **docs/config.md**: Fix environment variable examples from double-underscore `MTDNA_TRACY__TRIM` to single-underscore `MTDNA_TRACY_TRIM` matching actual `env_prefix`; add `TracySettings` and updated `Settings` sections
- **docs/tools/blastn.md**: Add `process_ab1_plots()` method to `SampleProcessor` documentation
- **docs/tools/unify.md**: Create new documentation for `src/tools/unify.py` batch unification module

### Changed

- **tools/blastn**: Decouple inner BLASTN threading from sample-level parallelism to fix CPU oversubscription. `preprocess_sample` no longer passes `-num_threads <cpu_count>` to each BLASTN subprocess; it now uses `get_settings().blastn.blastn_threads` (new setting, default 1). Previously every per-sample worker process (ProcessPoolExecutor, `max_workers` = cpu_count) launched BLASTN with `-num_threads cpu_count`, creating `max_workers` x `cpu_count` threads (e.g. 256 x 256 on a 256-core host). Each BLASTN subprocess is now single-threaded by default so N samples run on N cores without thrashing, matching the Tracy sample-parallel pattern. Raise `MTDNA_BLASTN_BLASTN_THREADS` only for single-sample runs

- **tools/blastn**: Split variant-transform logic out of `etl.py` into a new `src/tools/blastn/transform.py` (indel right-normalization, YAML conversion rules, polyC metadata handling). `etl.py` imports and re-exports the moved helpers so the public ETL entry point and the `tests/test_blastn_etl_parity.py` import surface are unchanged; behavior is identical (pure move)
- **docs/core/models.md**: Updated type annotations from `Union`/`Optional`/`List`/`Dict` to modern Python `|` / `list` / `dict` syntax; fixed `Tool` to `StrEnum`; updated `validate_position` signature to `float | str`; updated `peaks` field type to `list[list[float | None]] | None`
- **docs/core/batch.md**: Fixed `Batch` constructor signature to match code (`samples`, `ref_path`); updated `write()` signature to include `output_dir` and `batch_id` params; updated `to_json()` return type to `dict[str, Any]`
- **docs/core/flagging.md**: Updated type annotations to modern Python syntax; fixed `VariantAnalyzer`, `SampleFlagger`, and `flag_variants` signatures
- **docs/core/comparison.md**: Updated type annotations to modern Python syntax; added `ALWAYS_CONCORDANT_POSITIONS` constant reference; fixed `find_unique_variants` signature to use keyword-only `ignore_special`
- **docs/core/mtdna_merger.md**: Fixed class name from `mtDNAMerger` to `MtDnaMerger` (matching code); updated type annotations; fixed `merge` method signature
- **docs/core/sample.md**: Updated `sample_to_dict` to reference `generate_sequence()` instead of `gen_seq_from_variants_regions()`; fixed type annotations; fixed `filter_sample_by_regions` signature
- **docs/core/region.md**: Updated type annotations to modern Python syntax
- **docs/core/uploader.md**: Updated type annotations to modern Python syntax
- **docs/config.md**: Added `TracySettings` section with all Tracy parameters; added `tracy` field to `Settings` class; added Tracy environment variable overrides
- **docs/variants.md**: Updated `Position` type to `int | str`; fixed `validate_position`, `pos_base`, `pos_sort_key` signatures; updated `generate_sequence` return type; fixed `statistic_variants` and `split_variant_types` return types
- **docs/processing.md**: Updated type annotations to modern Python syntax
- **docs/generation.md**: Updated `FastaGenerator` description to reference `generate_sequence()`; added `generate_sequence` to cross-references
- **docs/TOOLS_STANDARD.md**: Updated Tracy conformance status from "Not conforming" to "Conforms"; updated Mutation Surveyor status; clarified BLASTn description
- **docs/ARCHITECTURE.md**: Updated Tracy conformance status; updated Mutation Surveyor status; noted `generate_sequence()` in Batch table
- **docs/tools/sequencher.md**: Updated type annotations
- **docs/tools/tracy.md**: Updated type annotations
- **docs/tools/mutation_surveyor.md**: Updated `__init__.py` description from "Empty — package marker" to "Public API exports"; updated type annotations
- **docs/modules/mutation_surveyor.md**: Updated `__init__.py` description; updated type annotations
### Fixed
- **mutation_surveyor**: Route per-region JSON output to shared `${ANALYSIS_DIR}/json_temp/${BATCH_ID}/` subdirectories (`mutation_surveyor/` and `sequencher/`) instead of `${OUT_DIR}` — compute `JSON_TEMP_DIR` from `OUT_DIR` so MS and Sequencher region JSONs land in the canonical `analysis/json_temp/${BATCH_ID}/` structure expected by unify step 9
- **mutation_surveyor**: Fix `find_custom_report()` to filter by BATCH_ID — when IN_DIR contains custom reports from multiple batches, the function now prefixes its glob with `${BATCH_ID}_` so each batch picks its own report instead of always selecting the alphabetically-first file (which caused all non-MS_080426_002 batches to fail QC Layer 0)

### Changed
- **region/unify**: Flatten JSON directory structure — remove redundant `JSON/` nesting from per-region intermediate output and redundant `batch_id/` nesting from merge output. Region JSONs now live directly under the tool directory (`mutation_surveyor/{LID}/` instead of `mutation_surveyor/JSON/{LID}/`), and unified output lives directly under `merge/` instead of `merge/{BATCH_ID}/`

### Changed
- **mutation_surveyor**: Route all JSON outputs (MS region JSON, Sequencher region JSON, unified JSON) through `${JSON_TEMP_DIR}/${BATCH_ID}/` subdirectories (`mutation_surveyor/`, `sequencher/`, `merge/`) instead of `OUT_DIR` and `RESULTS_DIR/tools/` — add `JSON_TEMP_DIR` to `.env`/`.env.example`, validate it in shell script, and update Steps 8–10 paths
- **sort_review_data**: Reorganize AB1 trace folder structure from `{category}/{HV1,HV2_HV3}/AB1/` to `{category}/AB1/{HV1F,HV1R,HV2F,HV3R}/` — add `PRIMER_TYPES` and `REGION_TO_PRIMERS` constants, refactor `_collect_trace_filenames()` to group by individual primer type, and rename `_copy_region_traces()` to `_copy_primer_traces()`
### Changed
- **tracy**: Add method-mapping comments to all module files documenting the correspondence between original src/tracy.py.bak methods and refactored src/modules/tracy/ functions — __init__.py, etl.py, pipeline.py, preprocessing.py, transforms.py, and utils.py now include explicit mapping comments and intentional-divergence documentation
- **tracy**: Fix behavioral divergences from original tracy.py.bak — add position 250 to _POSITION_TRANSFORMS dispatch (converts 248/250 deletions to 249 ref=A matching original), add position 523 transform (converts 523A>C to 523A-), add _transform_pos_514 handling for 513G-+514C- removal with 523/524 deletion creation, 513/514 deletion conversion to 523/524, and 514/515 deletion conversion to 523/524, fix _transform_pos_524 to move 524T to position 522 (matching original), restore _apply_polyc_primer_removal HV1 polyC removal with has_16189_T_C guard (matching original), restore pos 73 ref=A seq=G auto-pass in validate_peak_quality (matching original), restore deletion counting in _call_variant_at_index matching original's _call_variant behavior
- **tracy**: Document all intentional divergences from original tracy.py.bak — config via get_settings() (standard pattern), _normalize_ref1pos() bug fix, Variant model wrapping bug fix, Sample/Batch/flag_variants standard integrations, pandas→dict merge replacement, shutil.which() binary validation, and insertion tracking simplification

### Fixed
### Fixed
- **sequencher**: Remove empty batch directory after cleanup — `rmdir` the `${SEQUENCHER_DIR}/${BATCH_ID}` directory after deleting `txt/` so that stale empty stubs like `Sequencher_temp/MS_210426_003` are not left behind
- **mutation_surveyor**: Route Step 10 comparison output to `${OUT_DIR}/comparison` instead of relative `results/modules/comparison/${BATCH_ID}`
- **mutation_surveyor**: Route all JSON temp output to `${OUT_DIR}/json_temp/` subdirectories instead of `${JSON_TEMP_DIR}/${BATCH_ID}/` — remove `JSON_TEMP_DIR` from `.env`/`.env.example` and the validation block since `OUT_DIR` already encodes the batch ID
- **tracy**: Fix `AttributeError: 'Settings' object has no attribute 'variant_regions'` — change `get_settings().variant_regions` to `get_settings().regions.REGIONS` to match the `GenomicRegionsSettings` structure
- **tracy**: Fix `TypeError: 'int' object is not subscriptable` on `data["ref1pos"][index]` — Tracy >= 0.7.8 outputs `ref1pos` as a single integer (starting position) instead of an array; add `_normalize_ref1pos()` to convert integer to position array, use alignment index as variant position in `_call_variant_at_index`, and handle both int/array formats in `_create_variant`
- **tracy**: Fix `Variant` model validation errors where `peaks` (flat list) and `quality` (scalar) were rejected — wrap single-source peaks/quality in outer lists to match `list[list[float | None]]` and `list[float]` types
- **lint**: Add `PLR0915` to ruff ignore list for `process()` function

### Fixed
- **tracy**: Restore sample list from TXT file — add `BatchOptions.samples_path` field and `_read_sample_ids()` helper; `process_batch()` now always reads sample IDs from a TXT file matching the original `src/tracy.py.bak` behavior, no AB1 glob fallback
- **tracy**: Restore parallel sample processing — add `_process_single_sample()` combining decompose + ETL per sample, and use `ProcessPoolExecutor` with `max_workers` (default: `os.cpu_count()`) in `process_batch()`, matching the original `Tracy.process_samples()` pattern where each sample's full pipeline runs in parallel
- **tracy**: Make `--samples` CLI argument required in `src/tools/tracy.py` and wire it to `BatchOptions.samples_path`


### Fixed
- **tracy**: Fix basedpyright `reportUnusedFunction` errors — rename `_apply_all_transformations` → `apply_all_transformations` and `_update_variant_positions` → `update_variant_positions` since these are cross-module API functions, not internal helpers

### Changed
- **config**: Add TracySettings Pydantic settings class (trim, pratio, maxindel, quality_threshold, min_peak_value, heteroplasmy_threshold with MTDNA_TRACY_ prefix) to src.config; remove DEFAULT_* constants from src.modules.tracy.utils; resolve defaults via get_settings().tracy; remove dead TRACY_* env vars from .env.example and .env; fix TRACY binary path env var to MTDNA_TOOLS_TRACY
- **tracy**: Remove all noqa suppressions from Tracy module — refactor complex functions to satisfy lint rules: extract per-position transform helpers with dispatch table in `_apply_position_specific_transforms` (C901), extract `_apply_polyc_insertion`/`_apply_polyc_primer_removal` from `_apply_hv2_polyc_transforms` (C901), extract `_transform_reverse_strand`/`_redetect_heteroplasmy` from `_transform_variant_position` (C901), group optional params into `DecomposeConfig`/`_ProcessConfig` NamedTuples (PLR0913), use dict comprehensions in `detect_variant_conditions` (PERF403), validate binary path with `shutil.which` in `decompose_sample` (S603, added to global ruff ignore as false-positive-prone)

### Fixed
- **tracy**: Fix all 69 ruff lint errors across `src/modules/tracy` and `src/tools/tracy.py` — extract mtDNA position magic numbers into named constants (`POS_73`–`POS_16193`, `NUM_PEAK_CHANNELS`), reduce function argument counts with NamedTuple parameter objects (`_CallParams`, `_CreateParams`, `_TransformParams`, `_FilterConfig`, `_PositionParams`), break up complex functions (`_apply_all_transformations`, `process_batch`, `detect_variant_conditions`, `validate_peak_quality`, `_update_variant_positions`) into smaller focused helpers, fix unused noqa directives, collapse nested if/try-else blocks, replace dict key iteration with `.items()`, use `.values()` where only values are needed, fix import sorting, and add noqa for stable API signatures (`process()`, `decompose_sample()`)

### Changed
- **tracy**: Refactor from 1095-line class monolith (`src/tracy.py`) to standard module structure (`src/modules/tracy/`) following TOOLS_STANDARD.md — pure ETL (`etl.py`), batch orchestration (`pipeline.py`), AB1→decompose preprocessing (`preprocessing.py`), and utility helpers (`utils.py`) with thin CLI wrapper (`src/tools/tracy.py`)

### Removed
- **tracy**: Remove `src/tracy.py` monolith; use `src/modules/tracy/` and `src/tools/tracy.py` instead
- **tools/blastn**: Remove legacy `src/modules/fasta_to_json.py` (`FastaProcessor`) — replaced by the FASTA-input pipeline in `src/tools/blastn/fasta_pipeline.py`. The legacy `src/blastn.py`, `src/preprocessing.py`, `src/analysis.py`, and `src/postprocessing.py` it depended on are moved to `src/backup/` (no longer imported by any active module); `docs/modules/fasta_to_json.md` is removed in favor of the FASTA-input section in `docs/tools/blastn.md`
- **tracy**: Remove legacy `regenerate_result()` output path; Tracy now uses `Batch.write()` and `write_region_jsons()` for standard output
### Changed
- **tracy**: Deslop newly refactored module — remove duplicate CLI entry point from `pipeline.py` (kept only in `src/tools/tracy.py`), remove unused `ref_seq`/`insertions` params from `_call_variant_at_index`, remove unused `insertions_count` param from `_create_variant`, remove unused `insertions` param from `_apply_all_transformations`, collapse identical three-branch peak extraction into single call, move inline `merge_intervals` import to top-level in `etl.py`, move deferred imports to top-level in `src/tools/tracy.py`, move `detect_primer_type` import to top-level in `pipeline.py`, remove dead `filter_sample_by_regions` and `concurrent processing` comments, remove unused `sample_id` param from `_write_region_jsons_for_sample`, fix E501 docstring in `_compute_trim_bounds`
### Changed
- **logging**: Remove per-sample/per-region debug-only log lines across core and tool modules — region JSON writes, per-sample merge/load messages, auto-pass filter details, per-sample success messages, MtDnaMerger init, and per-region JSON summaries removed entirely; INFO logs now focus on batch-level summaries only
- **sequencher**: Fix misleading "N/M samples succeeded" message — now reports "N unique LIDs from M TXT files" to clarify that TXT file count ≠ unique sample count, and adds a WARNING when TXT files fail processing

### Fixed
- **sequencher**: Fix `parse_filename_ranges` to detect and strip `HV1_`/`HV2-3_` region prefixes before extracting sample ID and ranges from Sequencher filenames, fixing incorrect sample_id extraction and range parsing for region-prefixed files like `HV2-3_LN_26_AB5425-73-302 316-340 438-576.TXT`
- **sequencher**: Fix TypeError in `_combine_same_tool_samples` by replacing `normalize_position(v.pos)` with `pos_sort_key(v.pos)` for safe sorting of mixed int/str variant positions
- **sequencher**: Add Rule 6 to `validate_analysis_ranges`: region-prefixed files (HV1/HV2-3) must have ranges matching their declared region; HV2-3 files with HV1 ranges or HV1 files with HV2/HV3 ranges are now flagged as invalid

### Added
- **docs**: Standard Tool Module Template specification (docs/TOOLS_STANDARD.md) defining mandatory files (etl.py, pipeline.py, __init__.py), optional files (utils.py, preprocessing.py, flagging.py), CLI wrapper conventions, ETL/pipeline contracts, documentation requirements, and module communication rules for all tool modules
- **docs**: Architecture section 12 referencing the Standard Tool Module Template with current conformance status for Sequencher, Mutation Surveyor, Tracy, and BLASTn

### Changed
- **mutation_surveyor**: Refactor shell script to use associative array step registry (`STEPS_ENABLED`, `STEP_ALIASES`, `STEP_NUMBERS`) and `step_enabled()` helper instead of 10 boolean flags and verbose `if/else/case` patterns
- **sequencher**: Add cleanup of TXT and JSON directories in `sequencher.sh` before generating new results to avoid stale data from previous runs
- **mutation_surveyor**: Add log file output to `logs/{batch}.log` (matching `pipeline.sh` pattern) — `log_message()` now `tee -a` to `LOG_FILE`; all `python`/`bash` commands append `>> "$LOG_FILE" 2>&1`


### Removed
- **unify**: Removed legacy whole-sample merge path (`_process_legacy_samples`, `discover_samples`, `load_per_sample_json`) — both tool pipelines now produce per-region JSONs exclusively, making the fallback to `merge_samples()` unreachable
- **sequencher**: Removed `comparator.py` — legacy Sequencher comparison module fully superseded by `src/modules/sequencher/etl.py` (TXT → Sample ETL) and `src/core/comparison.py` (tool-agnostic pairwise comparison). Scrubbed all "legacy comparator" references from docs, comments, and test docstrings.
- **sequencher**: Removed `test_date` and `test_time` parameters from `process_batch()` (pipeline.py) and `process()` (etl.py) — unused in ETL
- **sequencher**: Removed `process_sample()` wrapper from `pipeline.py`; `process_batch()` now calls `process()` directly via ProcessPoolExecutor
- **sequencher**: Reduced `process_batch()` complexity by extracting `_write_single_region_jsons()`, `_write_full_region_jsons()`, `_write_region_jsons_for_sample()`, and `_build_combined_samples()` helpers
- **sequencher**: Removed unused `output_path`/`ref_path` from `process_sample()`, unused `cast` import, and alias `etl_process` in `pipeline.py`
- **sequencher**: Removed `--test_date` and `--test_time` CLI arguments from `src/tools/sequencher.py`
- **sequencher**: Removed `TIME` variable and `--test_time` argument from `scripts/tools/sequencher.sh`
- **sequencher**: Removed `extract_date_from_batch` import and auto-extraction logic from `pipeline.py`

### Changed
- **sequencher**: Remove unused `output_path` and `ref_path` from `process_sample()` in `pipeline.py`
- **sequencher**: Fix 7 ruff errors in `etl.py` — replace `open()` with `Path.open()`, extract `_resolve_seq_regions()` to reduce `process()` complexity, remove unused `_output_path`/`_ref_path` params, rename shadowed loop variable, extract magic number to `_MIN_VARIANT_PARTS` constant

### Added
- **sequencher**: Handle multiple TXT files per LID (two-TXT workflow) — `process_batch()` now accumulates Samples per LID, writes per-region JSONs per-TXT immediately, and combines same-LID Samples for `statistic_fullbatch.json`
- **sequencher**: Add `_combine_same_tool_samples()` helper to assemble region Samples from the same tool into one combined Sample for batch-level output
- **core**: Region-level intermediate JSONs + merge provenance — `src/core/region.py` provides per-region JSON I/O (`write_region_jsons`, `read_region_json`) and range QC (`validate_region_group_intervals`)
- **core**: `merge_by_regions()` in `mtdna_merger.py` — per-region merge with `region_sources` provenance tracking, replaces `merged_from`
- **unify**: `unify_batch()` discovers per-region JSONs from both MS and Sequencher, merges by region, falls back to legacy whole-sample merge

### Changed
- **sequencher**: `process_batch()` uses `Dict[str, List[Sample]]` accumulation instead of `Dict[str, Sample]` to support multiple TXT files per LID
- **mutation_surveyor**: Shell script Step 8 delegates Sequencher extraction/ETL to `sequencher.sh` when available, with inline fallback
- **mutation_surveyor**: Shell script Step 8/9 comments updated to describe per-region JSON flow and `region_sources` provenance

### Fixed
- **mutation_surveyor**: Fix false "Batch code mismatch" error when manifest MẺ CHẠY column uses spaces instead of underscores (e.g. "mtDNA 321" vs "mtDNA_321")

### Changed
- **mutation_surveyor**: Strip amino acid notation (`,p.XXX`) from variant cells during TXT→XLSX conversion (e.g. `c.310T>TC,p.Leu104LeuPro,$7` → `c.310T>TC,$7`)
- **mutation_surveyor**: Normalise TXT→XLSX output filename to `{BATCH}_{anything}_df_custom_report.xlsx` regardless of input suffix (`_CR_mutationSurveyor_CRSetting`, `_default_custom_report`, `_df_custom_report`, or `_custom_report`)
- **mutation_surveyor**: Update shell script step 1 to search for all custom report suffixes (`_CR_mutationSurveyor_CRSetting`, `_default_custom_report`, `_df_custom_report`, `_custom_report`) in priority order for both XLSX and TXT discovery; resolve converted output by glob instead of stem substitution to handle filename normalisation
- **mutation_surveyor**: Shell script falls back to module-bundled `mtDNA_PC_reference.xlsx` when no batch-local copy is found in the input directory
- **mutation_surveyor**: Move control_QC.xlsx output to control_QC/ subdirectory for organizational consistency with other pipeline output subdirectories
- **sort_review_data**: Renamed `autopass`/`not_autopass_one` to `no_triggering_flag`/`triggering_flag` in classification, folder names, enum values, and variable names to align with existing "triggering" terminology in the codebase


### Added

- **clean_ab1_files**: Add metadata validation — for each batch, read sample list from TSV/XLSX metadata files and verify every expected sample has exactly 4 AB1 files (one per primer). Report OUTLIER (≠4 files), ORPHAN (not in metadata), and OK batches. Add `--check-only` flag to run validation without cleaning.

### Changed

- **mutation_surveyor**: Remove `NOT_AUTOPASS_BOTH` category from `sort_review_data` — samples failing both regions are now classified as `NOT_AUTOPASS_ONE` with traces in both regional subfolders
- **mutation_surveyor**: Copy AB1 files for autopass samples into `autopass/HV1/AB1/` and `autopass/HV2_HV3/AB1/` (previously only `autopass_lids.txt` was written)
- **mutation_surveyor**: Remove `not_autopass_both_lids` field from `SortResult` model
- **mutation_surveyor**: Update `ReviewCategory` enum to two values: `AUTOPASS` and `NOT_AUTOPASS_ONE`


## [Unreleased] — 2026-06-10 — Remove columns.py; inline column name constants

### Removed
- **mutation_surveyor**: Deleted `columns.py` — only 2 active consumers remained (etl.py, sort_review_data.py), both now use inline string literals
- **mutation_surveyor**: Removed `from src.modules.mutation_surveyor.columns import ...` from `etl.py` and `sort_review_data.py`

### Changed
- **mutation_surveyor**: Replaced column name constants (`HV1_REVIEW_CLASS`, `HV23_REVIEW_CLASS`, `PROFILE_FILES`, `PROFILE_QUALITY`, `PROFILE_RANGE`, `PROFILE_VARIANTS`) with inline string literals in `etl.py`
- **mutation_surveyor**: Replaced column name constants (`HV1_REVIEW_CLASS`, `HV23_REVIEW_CLASS`) with inline string literals in `sort_review_data.py`

## [Unreleased] — 2026-06-10 — Remove unused code from mutation_surveyor etl and pipeline

### Removed
- **mutation_surveyor**: Removed `build_samples()` from `etl.py` — not imported by any active module; batch processing uses `pipeline.process_batch()` instead
- **mutation_surveyor**: Removed `process_sample()` from `pipeline.py` — not imported by any active module; single-LID processing uses `etl.process()` directly
- **mutation_surveyor**: Removed `from src.core.batch import Batch` from `etl.py` — only used by removed `build_samples()`
- **mutation_surveyor**: Removed unused `typing.Any` import from `pipeline.py` (auto-fixed by ruff)

### Changed
- **mutation_surveyor**: Updated module docstring in `etl.py` to reflect removal of `Batch.write()` reference

## [Unreleased] — 2026-06-10 — Remove unused code from mutation_surveyor utils and flagging

### Removed
- **mutation_surveyor**: Removed 8 unused functions from `flagging.py` — `flags_to_text`, `source_het_flags`, `build_hv1_flags`, `build_hv23_flags`, `hv23_severity`, `hv23_review_class`, `classify_review`, `combined_review` — not imported by any active pipeline module
- **mutation_surveyor**: Removed `VARIANT_CONVERTED_2` import and entire `utils` import block from `flagging.py` — only used by removed functions
- **mutation_surveyor**: Removed ~25 unused functions from `utils.py` — `normalize_text`, `format_num`, `_strip_m_prefix`, `parse_range_pair`, `parse_multi_ranges`, `intersect_ranges`, `overlap_pair`, `normalize_ranges`, `format_ranges`, `in_any_range`, `covers_segment`, `token_sort_key`, `merge_unique_tokens`, `merge_quality_tokens`, `first_nonempty`, `consistent_nonempty`, `alleles_compatible`, `compatibility_label`, `is_ignored_hv23_length_variant`, `parse_local_ac_token`, `split_local_ac_tokens`, `build_local_ac_consensus`, `enumerate_ac_alignments`, `_ac_motif_break_count`, `_deletion_positions`, `_substitution_positions`, `_contiguous_blocks`, `score_ac_alignment`, `_ac_steps_to_tokens`, `canonicalize_local_ac_consensus`, `normalize_ac_repeat_tokens`, `mismatch_positions_between_two_rows`, `is_common_hv23_insertion`, `is_common_hv23_deletion`, `parse_control_lid` — not imported by any active pipeline module
- **mutation_surveyor**: Removed unused constants from `utils.py` — `IUPAC`, `PHP_CODES`, `HV1_TARGET`, `HV1_POLYC`, `HV1_POLYC_MERGE`, `HV23_TARGET`, `HV2F_ADJUSTED`, `HV3R_ADJUSTED`, `RCRS_513_525`, `LOCAL_POSITIONS_513_525`, `REF_SEQ_513_525`, `REF_AC_CORE_513_525`, `HV2_POLYC_REGIONS`, `HV3_POLYC_REGIONS`, `COMMON_HV23_INSERTIONS`, `COMMON_HV23_DELETIONS`, `EXPECTED_PRIMERS` (partial — kept list), `_CONTROL_PC_PREFIX_RE` — not imported by any active pipeline module
- **mutation_surveyor**: Removed unused imports from `utils.py` — `collections.abc.Sequence`, `typing.Union`, `src.utils.merge_intervals`, `src.variants.is_position_in_intervals`, `src.variants.pos_base`, `src.variants.pos_sort_key`

### Changed
- **mutation_surveyor**: Modernized type annotations in `utils.py` — replaced `Union[int, str]` with `int | str` (Python 3.12+ syntax)
- **mutation_surveyor**: Updated stale module docstring in `utils.py` — replaced reference to non-existent modules with actual consumers

## [Unreleased] — 2026-06-10 — Deslop mutation_surveyor utils and flagging

### Changed
- **mutation_surveyor**: Modernized type annotations in `utils.py` — replaced `Union[int, str]` with `int | str` (Python 3.12+ syntax) in `in_any_range()` and `is_ignored_hv23_length_variant()`
- **mutation_surveyor**: Removed unused `from typing import Union` import from `utils.py`
- **mutation_surveyor**: Updated stale module docstring in `utils.py` — replaced reference to non-existent `quality_control`, `transform_hv1`, `transform_hv23` modules with actual consumers (`flagging`, `ETL`, `pipeline`, `sort_review_data`)

## [Unreleased] — 2026-06-10 — Add step 10 (compare) to mutation_surveyor.sh

### Added
- **mutation_surveyor**: Added Step 10 (`compare`) to `scripts/tools/mutation_surveyor.sh` — compares regenerate vs unified batch results using `src/core/comparison.py`, outputs to `results/modules/comparison/{BATCH_ID}/`, hard-fails when either batch JSON is missing

### Changed
- **mutation_surveyor**: Renumbered pipeline steps from 1–9 to 1–10 (total 10 steps, with new Step 10 = compare) with `10|compare` step selector in `-s/--steps`

## [Unreleased] — 2026-06-10 — Wire sort_review_data as step in pipeline.py

### Changed
- **mutation_surveyor**: `pipeline.py` `main()` now calls `sort_review_data()` after `process_batch()` — AB1 sort by autopass/review classification is now part of step 8 (ETL + auto-pass filter + AB1 sort) instead of a separate step 9
- **mutation_surveyor**: Added `--raw-dir` optional CLI arg to `pipeline.py` — defaults to `Settings.directories.data / "raw" / {batch_id}` when not provided
- **mutation_surveyor**: `pipeline.py` now imports `sort_review_data` from `sort_review_data.py` and `get_settings` from `src.config`
- **mutation_surveyor**: Removed `main()` and `__main__` block from `sort_review_data.py` — the only entry point is now `pipeline.py`

### Removed
- **mutation_surveyor**: Removed step 9 (AB1 sort) from `scripts/tools/mutation_surveyor.sh` — sort is now part of step 8 via `pipeline.py`
- **mutation_surveyor**: Removed standalone CLI `python -m src.modules.mutation_surveyor.sort_review_data` — replaced by `pipeline.py` orchestration

## [Unreleased] — 2026-06-10 — Wire mtdna_merger as Step 10 (unify) in mutation_surveyor.sh

### Added
- **unify**: New `src/tools/unify.py` — batch unification script that discovers per-sample JSONs from Mutation Surveyor and Sequencher directories, merges pairs via `mtdna_merger.merge_samples()`, passes through single-source (Sequencher-only) as unified, and writes per-sample JSONs + `statistic_fullbatch.json` via `Batch`
- **mutation_surveyor**: Added Step 10 (`unify`) to `scripts/tools/mutation_surveyor.sh` — merges MS + Sequencher results into `results/tools/unified/{BATCH_ID}/`, skips with warning when Sequencher dir is missing, requires MS JSON dir (step 8) as prerequisite

## [Unreleased] — 2026-06-10 — Wire TXT→XLSX conversion as Step 1 in mutation_surveyor.sh

### Changed
- **mutation_surveyor**: Added Step 1 (convert) to `scripts/tools/mutation_surveyor.sh` — checks `IN_DIR` for existing `.xlsx` custom report; if only `.txt` is found, converts it to `.xlsx` using `src/modules/mutation_surveyor/txt_to_excel.py` before proceeding
- **mutation_surveyor**: Renumbered pipeline steps from 1–8 to 1–9 (total 9 steps, with new Step 1 = convert) with `1|convert` step selector in `-s/--steps`

## [Unreleased] — 2026-06-09 — Rewrite mutation_surveyor.sh with step-based CLI
## [Unreleased] — 2026-06-10 — Fix AB1 sort: missing traces and output directory structure

### Fixed
- **mutation_surveyor**: `sort_review_data.py` now reads Traces sheets from BOTH `HV23_merged.xlsx` and `HV1_merged.xlsx` instead of only one — previously only HV2-3 trace filenames were loaded, causing HV1 AB1 files to be skipped entirely during copy
- **mutation_surveyor**: `sort_review_data.py` now copies AB1 files into nested `AB1/` subdirectories matching the documented output structure: `not_autopass_one/HV1/AB1/`, `not_autopass_one/HV2_HV3/AB1/`, `not_autopass_both/AB1/`
- **docs**: Updated `docs/tools/mutation_surveyor.md` §3.8 table to reflect the `AB1/` subdirectory nesting


### Changed
- **mutation_surveyor**: Replaced flat positional-arg script with step-based CLI matching `pipeline.sh` pattern — added `-s/--steps` (comma-separated numbers or names), `-h/--help`, `-f/--force` options; each of 8 steps now gated by a flag with dependency checks and `[Step N/8 - name]` log prefixes; steps: 1=trace, 2=control, 3=hv23, 4=hv1, 5=profiles, 6=truth, 7=etl, 8=sort

## [Unreleased] — 2026-06-09 — Rewrite mutation_surveyor.sh as simple pipeline wiring script

### Changed
- **mutation_surveyor**: Replaced `scripts/tools/mutation_surveyor.sh` with a fresh minimal script that wires the 6 Python pipeline scripts (`ms_trace_process`, `ms_control_QC`, `ms_hv23_merge`, `ms_hv1_merge`, `ms_final_profiles_merge`, `ms_final_profiles_vs_truth`) in dependency order with `BATCH_ID IN_DIR OUT_DIR` CLI args

### Removed
- **mutation_surveyor**: Removed `.env` sourcing, auto-discovery helpers, multi-batch loop, AB1 sorting, JSON filtering, Loguru logging, and `init_batch_dirs` from `scripts/tools/mutation_surveyor.sh`

## [Unreleased] — 2026-06-09 — Accept both TXT and Excel in Mutation Surveyor shell script input discovery

### Changed
- **mutation_surveyor**: `find_mutation_surveyor_input()` in `scripts/tools/mutation_surveyor.sh` now matches custom report files by name pattern regardless of extension (`.*`), accepting both `.txt` and `.xlsx` inputs — aligning with the Python pipeline's `find_custom_report()` which already handles both formats
- **mutation_surveyor**: `find_lid_manifest()` fallback search in `scripts/tools/mutation_surveyor.sh` now matches `_custom_report.*` and `_CR_mutationSurveyor_CRSetting.*` patterns instead of only `*.xlsx`

## [Unreleased] — 2026-06-09 — Wire Layer 1 Control QC into Mutation Surveyor shell script

## [Unreleased] — 2026-06-09 — Align sort_review_data docs with folder structure

### Changed
- **mutation_surveyor**: Fixed module docstring in `sort_review_data.py` — changed `HV2-3/` to `HV2_HV3/` to match the actual output directory name created by the code
- **docs**: Updated `docs/tools/mutation_surveyor.md` §3.8 — changed output path from `results/tools/mutation_surveyor/{BATCH_ID}/data/` to reference §3.7.1 (`results/analysis/mutation_surveyor/{BATCH}/`), changed CLI argument from `--data-dir` to `--raw-dir` to match actual code, and fixed note reference from `--data-dir` to `--raw-dir`
- **docs**: Added `sort_review_data.py` to File Map and new §10b section in `docs/modules/mutation_surveyor.md`
- **docs**: Added `sort_review_data.py` to Mutation Surveyor module table in `docs/ARCHITECTURE.md` §11

### Changed
- **mutation_surveyor**: Added Step 2b (Layer 1 Control QC) to `scripts/tools/mutation_surveyor.sh` — `control_qc.py` is now invoked after the core pipeline produces the trace-processed xlsx, before concordance validation
- **mutation_surveyor**: Added `CONTROL_REF` config variable defaulting to `${BASE_DIR}/ref/PC_control_reference.xlsx`; Layer 1 QC is skipped with a warning when the reference file is not configured or not found
- **docs**: Updated `docs/tools/mutation_surveyor.md` to note Layer 1 Control QC is wired as Step 2b

## [Unreleased] — 2026-06-09 — Document proposed Mutation Surveyor folder structure

### Changed
- **docs**: Updated Section 3.7 of `docs/tools/mutation_surveyor.md` — replaced the old two-location output directory structure with the proposed target structure using `results/analysis/mutation_surveyor/{BATCH}/` for reviewer-facing intermediates and `results/tools/mutation_surveyor/{BATCH}/` for downstream-facing JSON output. Added tree diagrams, naming convention notes, ZIP file documentation, and CLI argument table

## [Unreleased] — 2026-06-08 — Fix AC repeat normalization losing DEL/insertion/het tokens

### Fixed
- **mutation_surveyor**: Fixed `normalize_ac_repeat_tokens()` silently dropping DEL, insertion, and het tokens in the AC repeat region (positions 513-525). Three bugs were fixed:
  1. DEL tokens (e.g., `523DEL`, `524DEL` from `523_524delAC` expansion) produced a malformed consensus string (`DEL` is 3 chars, not 1), causing alignment failure and all local tokens to be lost
  2. Insertion tokens (e.g., `515.1C` from `515dupC`) were skipped in consensus building but not included in canonical output, causing silent data loss
  3. Het tokens (e.g., `523A_het`) had the same malformed-consensus issue as DEL tokens
  4. Contagion effect: when any DEL or het token was present, ALL local tokens (including valid SNPs) were also lost due to failed alignment
- **mutation_surveyor**: `split_local_ac_tokens()` now returns 4 categories instead of 3: `local_snp` (single-nucleotide substitutions for canonicalization), `local_passthrough` (DEL, insertion, het tokens that pass through unchanged), `local_unsupported`, and `nonlocal_tokens`

### Changed
- **mutation_surveyor**: `normalize_ac_repeat_tokens()` now classifies local tokens as SNP (canonicalize) vs passthrough (preserve unchanged). Only single-nucleotide SNP tokens at integer positions with allele in ACGT go through the consensus alignment; all other tokens in the 513-525 range pass through unchanged

## [Unreleased] — 2026-06-08 — Remove unused Sample I/O functions

### Removed
- **sample**: Removed `load_unwrapped_sample_json()`, `write_unwrapped_sample_json()`, and `discover_sample_jsons()` from `src/core/sample.py` — these were moved from the deleted `filter_regions.py` but have no callers in the current codebase

## [Unreleased] — 2026-06-08 — Remove merge_unified.py; Use mtdna_merger Directly in Shell Script

### Changed
- **shell**: Step 6 in `mutation_surveyor.sh` now calls `python -m src.core.mtdna_merger` per-sample instead of `python -m src.modules.mutation_surveyor.merge_unified` — merge logic lives in the shell loop directly, output goes to `results/tools/unified/{BATCH}/`
- **shell**: Single-source samples (only MS or only Sequencher JSON) are now skipped with a warning instead of passed through

### Removed
- **mutation_surveyor**: Removed `merge_unified.py` — its two-dir discovery + pass-through logic is replaced by per-sample checks in the shell script using `mtdna_merger` directly

### Fixed
- **shell**: Fixed `mutation_surveyor.sh` Step 2 variable collision — `MUTATION_SURVEYOR_INPUT_DIR` was used for both the config search path and the computed result, causing `find_mutation_surveyor_input` to always fail and the script to exit at Step 2. Config path renamed to `MUTATION_SURVEYOR_INPUT_ROOTS`; computed result uses `MS_INPUT_DIR`
- **shell**: Added `|| log_message` error handling to the `bash scripts/pipeline.sh -p "tracy" -s "4"` call so a Tracy/manual pipeline failure no longer aborts the entire MS pipeline
- **mtdna_merger**: Fixed `_load_profile()` to accept unwrapped per-sample JSON format directly — previously it required a wrapped format with a top-level sample_id key, which does not match the standard `{LID}/{LID}.json` layout used by all pipeline tools
- **mtdna_merger**: Removed duplicate `{sample_id}_mtdna_merged.json` (wrapped) output from `merge()` and `merge_samples()` — now writes only `{sample_id}/{sample_id}.json` (unwrapped) to match the standard per-sample layout used by all other pipeline tools
- **mtdna_merger**: Removed implicit overlap conflict check — when only one profile reports a variant at a position, the explicit call is now trusted instead of raising a hard error. Tools have different detection sensitivities (e.g., one calls deletions, another doesn't), so implicit silence is not treated as a conflicting allele
- **tests**: Updated `test_mtdna_merger` payloads to use unwrapped per-sample JSON format matching actual pipeline output

## [Unreleased] — 2026-06-08 — Integrate filter_regions into pipeline module

### Removed
- **mutation_surveyor**: Removed standalone `filter_regions.py` module — `AUTOPASS_FLAG_TO_REGIONS` and `autopass_regions_from_flags()` are now defined directly in `pipeline.py`

### Changed
- **mutation_surveyor**: `pipeline.py` now defines `AUTOPASS_FLAG_TO_REGIONS` and `autopass_regions_from_flags()` inline instead of importing from `filter_regions.py`
- **docs**: Updated `docs/ARCHITECTURE.md` — removed `filter_regions.py` row from MS module table, added auto-pass mapping note to `pipeline.py` row
- **docs**: Updated `docs/core/mtdna_merger.md` — removed stale `filter_regions.py` CLI reference, updated Step 5b description

## [Unreleased] — 2026-06-07 — Remove TSV Output from Mutation Surveyor

### Removed
- **mutation_surveyor**: Removed `ResultBase.to_tsv()` method from `models.py` — pipeline stages no longer write TSV files alongside xlsx
- **mutation_surveyor**: Removed TSV file writing from `pipeline.py` — `process_batch()` no longer calls `.to_tsv()` on intermediate results
- **mutation_surveyor**: Removed TSV file writing from `_write_region_subset()` in `merge_profiles.py` — region subset files no longer produce companion TSV files
- **mutation_surveyor**: `sort_review_data.py` now reads xlsx directly instead of falling back to TSV files — removed TSV file path construction and `pd.read_csv(sep="\t")` branches

### Changed
- **mutation_surveyor**: Updated docstrings across `models.py`, `pipeline.py`, `merge_profiles.py`, `flagging.py`, and `sort_review_data.py` to reflect xlsx-only output
- **mutation_surveyor**: Updated `docs/tools/mutation_surveyor.md` to remove all TSV output references


## [Unreleased] — 2026-06-07 — Add Region-Specific Subset Outputs for Merge Profiles

### Added
- **mutation_surveyor**: Added `HV1_SUBSET_COLUMNS` and `HV23_SUBSET_COLUMNS` column group constants to `columns.py` for region-specific subset output
- **mutation_surveyor**: Added `_build_region_summary()` helper to `merge_profiles.py` — builds summary with region-specific review class and flag counts
- **mutation_surveyor**: Added `_write_region_subset()` to `merge_profiles.py` — writes a region-specific xlsx + tsv with 4 sheets mirroring `final_profiles.xlsx`
- **mutation_surveyor**: `run_merge_profiles()` now writes two additional region-specific subset files when `output_path` is provided: `{stem}_HV1.xlsx` (HV1 columns + Sample profile + Combined review) and `{stem}_HV23.xlsx` (HV2-3 columns + Sample profile + Combined review)
- **mutation_surveyor**: Added `hv1_subset_path` and `hv23_subset_path` fields to `LIDMergeResult` model to track subset output paths

### Changed
- **mutation_surveyor**: Updated `docs/tools/mutation_surveyor.md` — merge_profiles section now documents HV1 and HV23 subset outputs

## [Unreleased] — 2026-06-06 — Fix M: Prefix Range Parsing in Mutation Surveyor Pipeline

### Fixed
- **mutation_surveyor**: Fixed `parse_range_pair()` and `parse_multi_ranges()` in `utils.py` to handle `M:` prefixed position strings (e.g. `M:16018-M:16193`) from CR format Excel files — previously only plain numeric ranges (e.g. `73-340`) were supported, causing all adjusted range/variant columns to be empty and all samples to be classified as "Review" with 0 autopass



## [Unreleased] — 2026-06-06 — Remove Truthset Excel Path from Validation

### Changed
- **mutation_surveyor**: Renamed `run_validation_from_batch()` → `run_validation()` with `batch_json_path` parameter — truthset Excel path removed, batch JSON is now the sole validation source
- **mutation_surveyor**: `--batch-json` is now a required CLI argument (was mutually exclusive with `--truth`)
- **mutation_surveyor**: Updated shell script Step 3 — removed `find_truthset()` helper, `TRUTHSET_ROOTS`, and `--truth` fallback branch; validation now requires batch JSON only
- **mutation_surveyor**: Updated docstrings to reference "batch JSON" / "reference DataFrame" instead of "truthset"

### Removed
- **mutation_surveyor**: Removed `run_validation(merged_path, truth_path, ...)` (truthset Excel path) from `validation.py`
- **mutation_surveyor**: Removed `_load_truth()` from `validation.py` — only loaded truthset Excel workbooks
- **mutation_surveyor**: Removed `--truth` CLI argument from `validation.py` main()
- **mutation_surveyor**: Removed `find_truthset()` helper and `TRUTHSET_ROOTS` from shell script

## [Unreleased] — 2026-06-06 — Clean Up mutation_surveyor.sh Shell Script

### Changed
- **mutation_surveyor**: Replaced `log_message` python -c Loguru calls with plain `printf` + timestamp — removes Python process spawn per log line
- **mutation_surveyor**: Extracted `find_ms_input()`, `find_lid_manifest()`, `find_truthset()` helper functions — deduplicates inline search loops from Steps 1, 3, and 5
- **mutation_surveyor**: Converted `MS_INPUT_ROOTS`, `TRUTHSET_ROOTS`, `LID_MANIFEST_ROOTS` from plain strings to proper bash arrays — fixes `"${var[@]}"` usage
- **mutation_surveyor**: Renamed shadowed `RESULTS_DIR` in Step 6 to `BATCH_RESULTS_DIR` — was overwriting the global `RESULTS_DIR` from `.env`
- **mutation_surveyor**: Renumbered steps 3→2, 5→3, 6→4 after removing dead Step 2 (commented-out pipeline.sh) and Step 4 (commented-out Control QC)
- **mutation_surveyor**: Updated header comment to reflect current step numbers

### Removed
- **mutation_surveyor**: Removed ~60 commented-out batch IDs from `BATCHES` array
- **mutation_surveyor**: Removed commented-out Step 2 (`bash scripts/pipeline.sh -p "tracy" -s "2"`) and Step 4 (Layer 1 Control QC) blocks
- **mutation_surveyor**: Removed unused `COMPARISON_DIR` variable

## [Unreleased] — 2026-06-06 — Remove Hardcoded Filename Format Checks from Inventory QC

### Changed
- **mutation_surveyor**: Removed hardcoded filename format validation from `inventory_qc.py` — no longer requires manifest filenames to match `YYYYMMDD_Danh sách mẫu_{BATCH}.xlsx` or custom report filenames to match `{BATCH}_anything_df_custom_report.xlsx`; batch code is now derived solely from the manifest `MẺ CHẠY` column
- **mutation_surveyor**: Removed `_batch_from_manifest_filename()`, `_manifest_filename_matches_required_pattern()`, and `_custom_report_filename_matches_batch()` helper functions from `inventory_qc.py` (no longer needed)
- **mutation_surveyor**: Added `_CR_mutationSurveyor_CRSetting` pattern support to `find_custom_report()` in `mutation_surveyor.py` — recognizes both `_custom_report` and `_CR_mutationSurveyor_CRSetting` file naming conventions

## [Unreleased] — 2026-06-06 — Add Direct Batch JSON Validation Path & Multi-Batch Path Discovery

### Added
- **mutation_surveyor**: Added `run_validation_from_batch()` to `validation.py` — validates merged profiles directly against a tool batch JSON (e.g., Sequencher `statistic_fullbatch.json`) instead of requiring a manually-recreated truthset Excel file
- **mutation_surveyor**: Added `_samples_to_truth_df()` helper to `validation.py` — converts `Sample` objects into truthset-compatible DataFrame format
- **mutation_surveyor**: Added `_run_validation_core()` to `validation.py` — extracted shared comparison logic from `run_validation()` for reuse by both truthset and batch JSON paths
- **mutation_surveyor**: Added `--batch-json` CLI flag to `validation.py` — mutually exclusive with `--truth`, selects batch JSON validation mode
- **mutation_surveyor**: Updated `scripts/tools/mutation_surveyor.sh` — uses dynamic path discovery (`MS_INPUT_ROOTS`, `TRUTHSET_ROOTS`, `LID_MANIFEST_ROOTS`) instead of hardcoded directories; `SEQUCHER_RESULTS_DIR` uses strict `${RESULTS_DIR}/manual_pipeline/sequencher/{BATCH}/statistic_fullbatch.json` path; step 5 prefers `--batch-json` when available, falls back to `--truth`

## [Unreleased] — 2026-06-05 — Align Mutation Surveyor Output with Reference Scripts

### Changed
- **mutation_surveyor**: Changed column name `Combined Flags` → `Combined Flag list` in `columns.py` and all consumers to match reference `ms_final_profiles_merge.py` output
- **mutation_surveyor**: Aligned `merge_profiles.py` output sheets with reference: now writes `LID_no_PC_NTC`, `Summary_no_PC_NTC`, `LID`, `Summary` (was `LID`, `PC_NTC`, `Summary`)
- **mutation_surveyor**: Removed extra columns from `merge_profiles.py` output that were not in the reference: `Sample profile - Quality`, `Sample profile - Files`, `HV2-3 Quality Adjusted`, `HV1 Quality Adjusted`
- **mutation_surveyor**: Removed `merge_quality_text()` and `merge_files_text()` helper functions from `merge_profiles.py` (no longer needed)
- **mutation_surveyor**: Changed `_pick_sheet()` default in `validation.py` from `"LID"` to `"LID_no_PC_NTC"` to match reference `ms_final_profiles_vs_truth.py` behavior
- **mutation_surveyor**: Added `_style_validation_xlsx()` to `validation.py` — applies Excel styling (freeze panes, auto-filter, header formatting, column auto-width) matching reference output
- **mutation_surveyor**: Updated `LIDMergeResult` docstring in `models.py` to reflect new 4-sheet output structure

### Added
- **mutation_surveyor**: Added `lid_manifest_path` parameter to `pipeline.py::process_batch()` — when provided, runs Layer 0 inventory QC as a blocking gate before producing trace-processed output, matching reference `ms_trace_process.py` behavior
- **mutation_surveyor**: Added `InventoryQCError` exception class to `pipeline.py` for Layer 0 gate failures
- **mutation_surveyor**: Added `main()` CLI entry point to `merge_files.py` — now invocable as `python -m src.modules.mutation_surveyor.merge_files`

## [Unreleased] — 2026-06-05 — Separate Excel/TSV Output to Review Directory

### Changed
- **mutation_surveyor**: Added `excels_dir` parameter (required) to `process_batch()` in `pipeline.py` — all xlsx/tsv files are now written to `excels_dir`
- **mutation_surveyor**: Added `--excels-dir` CLI argument (required) to `src/tools/mutation_surveyor.py` — specifies the directory for xlsx/tsv outputs
- **mutation_surveyor**: Added `excels_dir` parameter (required) to `sort_review_data()` in `sort_review_data.py` — reads xlsx/tsv from excels_dir
- **mutation_surveyor**: Updated `scripts/tools/mutation_surveyor.sh` — added `MUTATION_SURVEYOR_EXCELS` variable pointing to `review/{batch_id}/excels/`, all xlsx/tsv outputs (core pipeline, inventory QC, control QC, validation) now go to excels directory
- **mutation_surveyor**: Shell script now passes `--excels-dir` to core pipeline and `--excels-dir` to sort_review_data CLI
- **docs**: Updated `docs/ARCHITECTURE.md` with output directory separation documentation
- **docs**: Updated `docs/tools/mutation_surveyor.md` with excels_dir in pipeline flow diagram and new output directory structure section

## [Unreleased] — 2026-06-05 — Remove auto_format_worksheet and Align Output Filenames

### Changed
- **mutation_surveyor**: Removed `auto_format_worksheet()` from `utils.py` and all callers (`models.py`, `merge_files.py`, `merge_profiles.py`, `control_qc.py`, `inventory_qc.py`) — output xlsx files no longer apply freeze panes, auto-filter, or auto-width formatting
- **mutation_surveyor**: Aligned all output filenames to match reference script naming convention:
  - `compare_QC.xlsx` → `{BATCH}_MS_trace_processed.xlsx`
  - `HV2-3.xlsx` → `{BATCH}_HV23_merged.xlsx`
  - `HV1.xlsx` → `{BATCH}_HV1_merged.xlsx`
  - `HV23_HV1_final_profiles_v5.xlsx` → `{BATCH}_final_profiles.xlsx`
  - `concordance.xlsx` → `{BATCH}_final_profiles_vs_truth.xlsx`
  - `{BATCH}_inventory_QC.xlsx` → `{BATCH}_MS_trace_processed_inventory_QC.xlsx`
- **mutation_surveyor**: Updated `default_stem` values in `models.py` to match reference naming: `MS_trace_processed`, `HV1_merged`, `HV23_merged`, `final_profiles`, `final_profiles_vs_truth`
- **mutation_surveyor**: Updated `KNOWN_SUFFIXES` in `merge_files.py` to reference naming only — removed backward compatibility entries (`_HV2-3`, `_HV1`, `_HV23_HV1_final_profiles_v5`, `_compare_QC`, `_df_compare_QC`)
- **mutation_surveyor**: Removed backward compatibility suffix `_df_compare_QC` from batch ID derivation in `pipeline.py`, `quality_control.py`, and `mutation_surveyor.py`
- **mutation_surveyor**: Removed backward compatibility suffix `_HV23_HV1_final_profiles_v5` from ETL batch ID derivation in `etl.py`
- **mutation_surveyor**: Updated `mutation_surveyor.sh` shell script filenames to match reference convention
- **mutation_surveyor**: Updated `sort_review_data.py` to use new filename patterns with `batch_id` prefix
- **mutation_surveyor**: Updated CLI help texts in `control_qc.py` and `validation.py` to reference new filenames
- **mutation_surveyor**: Updated docstrings in `transform.py` to reference `MS_trace_processed` instead of `compare_QC`
- **docs**: Updated `docs/tools/mutation_surveyor.md` with new filename conventions

# Changelog

All notable changes to the **mtdna-raw** project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — 2026-06-05 — Fix Batch Code Extraction, TXT Handling, and Shell Script Logging

### Fixed
- **mutation_surveyor**: Shell script (`mutation_surveyor.sh`) now uses `log_message` instead of bare `echo` calls for all pipeline step messages
- **mutation_surveyor**: `_batch_from_manifest_filename()` in `inventory_qc.py` now correctly extracts the batch ID from manifest filenames that include a `_metadata` suffix (e.g. `20260515_Danh sÃ¡ch máº«u_MS_210426_003_metadata.xlsx` -> `MS_210426_003` instead of `MS_210426_003_metadata`), fixing the "Batch code mismatch" false positive
- **mutation_surveyor**: Inventory QC (`inventory_qc.py`) now converts TXT custom reports to XLSX before processing — prevents `openpyxl.InvalidFileException` when the input is a `.txt` file
- **mutation_surveyor**: Shell script (`mutation_surveyor.sh`) now prioritizes `.xlsx` files over `.txt` files when searching for custom reports

## [Unreleased] — 2026-06-04 — Mutation Surveyor Slop Cleanup

### Changed
- **mutation_surveyor**: Removed duplicate `_clean_str()` function in `control_qc.py` — replaced with `safe_str()` from `utils.py`
- **mutation_surveyor**: Moved `_parse_control_lid()` from `control_qc.py` to `utils.py` as `parse_control_lid()` — consolidates control LID parsing alongside `is_control_lid()`
- **mutation_surveyor**: Replaced inline `_CONTROL_PC_RE`/`_CONTROL_NTC_RE` regex usage in `control_qc.py` with `is_control_lid()` and `parse_control_lid()` from `utils.py`
- **mutation_surveyor**: Fixed transitive `Position` import — `transform.py` and `validation.py` now import `Position` from `src.core.models` directly instead of through `utils.py`

## [Unreleased] — 2026-06-04 — AB1 Review Sorting

### Added
- **mutation_surveyor**: Add `sort_review_data` module — classifies LIDs by autopass/review status per region and copies AB1 trace files into categorized folders (`autopass/`, `not_autopass_one/HV1/`, `not_autopass_one/HV2-3/`, `not_autopass_both/`) under `results/tools/mutation_surveyor/{BATCH_ID}/data/`
- **mutation_surveyor**: Add `--data-dir` CLI flag to `mutation_surveyor.py` for AB1 sorting integration
- **mutation_surveyor**: Add Step 7 (sort review data) to `pipeline.process_batch()` — runs after validation when `data_dir` is provided

## [Unreleased] — 2026-06-03 — Docs Sync

### Changed
- **docs**: Synced all documentation files with current source code APIs, method signatures, and behavior
- **docs**: Updated `docs/analysis.md` — added missing methods: `analyze_polyc_context`, `is_after_polyc_region`, `is_before_polyc_region`, `load_quality_scores`, `call_variants`, `right_norm_indel`
- **docs**: Updated `docs/tracy.md` — added missing methods: `create_variant`, `convert_variant`, `update_variant`, `update_variant_data`; expanded private methods documentation
- **docs**: Updated `docs/processing.md` — expanded `PostProcess` API with all methods: `load_filter_positions`, `load_no_filter_positions`, `is_special_position`, `filter_special_positions`, `get_regions_containing_position`, `apply_region_validation_rules`, `is_in_polyc_region`, `is_after_polyc_region`, `is_polyc_stretch`, `filter_poliC`, `parse_variant_rule`, `check_variant_exist`, `_check_rule_conditions`, `_apply_rule_changes`; updated `PCProcessor` method names to match source
- **docs**: Updated `docs/pipeline.md` — added missing methods to `MtDNAPipeline` (`setup_directories`, `run_pipeline`), `ManualPipelineProcessor` (`process_manual_results`, `process_sample_pipeline_results`, `diff_manual_pipeline`, `format_diff`, `filter_variants_by_intervals`), `VariantProcessor` (`analyze_position_reverse`, `process_blastn_variants`)
- **docs**: Updated `docs/utils.md` — added missing functions: `get_start_pos`, `reverse_complement`, `read_fastq`, `preprocess_seq`, `check_sequence_get_start_pos`, `get_highlight_positions`, `save_plt`, `get_annotation_pos`, `iupac_to_bases`, `get_blastn_path`, `check_output_path`, `count_index`, `check_folder_exist`, `annotate_variants`, `annotate_positions`, `plot_chromatograph`, `annotate_peaks`, `highlight_peak`, `annote_mutation`, `get_highlight_position`, `get_align_pos`, `get_list_positions`
- **docs**: Updated `docs/generation.md` — added async methods to `ReportGenerator` (`get_pdf_link`, `download_pdf`, `process_sample`, `generate_reports`)
- **docs**: Updated `docs/core/batch.md` — added `to_json()` method
- **docs**: Updated `docs/core/comparison.md` — added `format_intervals`, `compare_batches`, `summary_row`, and module-level functions (`create_variant_level_data`, `write_comparison_json`, `write_comparison_tsv`, `write_variant_level_tsv`, `write_comparison_excel`, `compare_batch_files`)
- **docs**: Updated `docs/core/flagging.md` — added `VariantAnalyzer` methods: `check_combined_conditions`, `has_variant_at`, `has_523_524_ac_deletion`, `has_513a_with_523_524_deletion`; added `SampleFlagger` private methods
- **docs**: Updated `docs/core/models.md` — synced all model fields with source code; added `Sequence`, `validate_position` documentation
- **docs**: Updated `docs/core/sample.md` — synced `sample_to_dict` and `load_sample_batch` with source
- **docs**: Updated `docs/core/mtdna_merger.md` — added all private methods and `merge_samples` function
- **docs**: Updated `docs/core/merger.md` — synced with source code
- **docs**: Updated `docs/variants.md` — added `generate_sequence` function documentation
- **docs**: Updated `docs/config.md` — verified all settings fields match source
- **docs**: Updated `docs/blastn.md` — verified API accuracy
- **docs**: Updated `docs/ARCHITECTURE.md` — added source module index table, updated date

### Added
- **docs**: New `docs/modules/LAB/generate_mtdna_html_report.md` — LAB HTML report generation
- **docs**: New `docs/modules/NGS/README.md` — NGS module overview
- **docs**: New `docs/modules/PYQG/README.md` — PYQG module overview
- **docs**: New `docs/modules/TNLS/README.md` — TNLS module overview
- **docs**: New `docs/modules/statistics/README.md` — Statistics module overview
- **docs**: New `docs/modules/fasta_to_json.md` — FASTA to JSON converter
- **docs**: New `docs/modules/mutation_surveyor_inventory_qc.md` — MS Inventory QC documentation
- **docs**: New `docs/modules/mutation_surveyor_control_qc.md` — MS Control QC documentation
- **docs**: New `docs/modules/mutation_surveyor_txt_to_excel.md` — MS TXT to Excel converter
- **docs**: New `docs/modules/sequencher_comparator.md` — Legacy Sequencher comparator


- **mutation_surveyor**: Centralized `is_control_lid()` and `EXPECTED_PRIMERS` in `utils.py`; removed duplicate implementations from `inventory_qc.py`, `control_qc.py`, and `merge_profiles.py`
- **mutation_surveyor**: Replaced `os.path.basename()` with `Path.name` in `inventory_qc.py`; removed unused `os` import
- **mutation_surveyor**: New `inventory_qc.py` module with `run_inventory_qc()` — Layer 0 batch integrity gate validating MS custom report against LID manifest (PLAN.md Task A1)
- **mutation_surveyor**: New `InventoryQCResult` Pydantic model with 5-sheet report output (Inventory_Summary, Issue_Counts, Inventory_Details, Primer_Counts, LID_Primer_Matrix)
- **mutation_surveyor**: `pipeline.py:process_batch()` now accepts `lid_manifest_path` and `control_ref_path` parameters; inventory QC runs as a blocking gate before QC output is written
- **mutation_surveyor**: Inventory QC checks: manifest filename validation, batch code consistency, blank/duplicate LIDs, missing/extra sample LIDs, primer counts, per-LID primer matrix

- **mutation_surveyor**: New `control_qc.py` module with `run_control_qc()` — Layer 1 PC/NTC Control QC validating PC profiles against reference file and NTC profiles for emptiness (PLAN.md Task A2)
- **mutation_surveyor**: New `ControlQCResult` Pydantic model with 4-sheet report output (Summary, Details, Observed_Profiles, Reference_Profiles)
- **mutation_surveyor**: `pipeline.py:process_batch()` accepts `control_ref_path` parameter for Layer 1 Control QC integration (non-blocking, optional)

- **mutation_surveyor**: `build_hv23_flags()` now includes: "No full region", "Trace mismatch", "No 315.1C", "310 variant", "non-C in HV2 polyC", "non-C in HV3 polyC", "unexpected insertion", "unexpected deletion", "HV3R PHP 341-437" flags (PLAN.md Task B1)
- **mutation_surveyor**: `build_hv1_flags()` now includes: "HV1 trace mismatch", "HV1 polyC insertion", "HV1 polyC deletion", "HV1 polyC PHP", "HV1 indel outside polyC" flags (PLAN.md Task B2)
- **mutation_surveyor**: New constants in `utils.py`: `HV2_POLYC_REGIONS`, `HV3_POLYC_REGIONS`, `COMMON_HV23_INSERTIONS`, `COMMON_HV23_DELETIONS`, plus `is_common_hv23_insertion()` and `is_common_hv23_deletion()` helpers (PLAN.md Tasks B1-B2)

- **mutation_surveyor**: New HV1 trace selection detail columns: `HV1F Selected Trace #`, `HV1R Selected Trace #`, `HV1F Selected Sample Name`, `HV1R Selected Sample Name` in LID output (PLAN.md Task B7)
- **mutation_surveyor**: New HV2-3 trace selection detail columns: `HV2F Candidate Count`, `HV3R Candidate Count`, `HV2F Selected Trace #`, `HV3R Selected Trace #`, `HV2F Selected Sample Name`, `HV3R Selected Sample Name` in LID output (PLAN.md Task B8)
- **mutation_surveyor**: Fixed `_profile_to_string()` in `control_qc.py` — restored `sorted(set(...))` dedup+sort behavior matching the reference script; the plain `". ".join(tokens)` could produce duplicate tokens and inconsistent ordering causing spurious PC profile mismatches
- **mutation_surveyor**: Fixed `read_control_reference()` in `control_qc.py` — restored duplicate reference ID detection that raises `ValueError` on duplicate LIDs; the silent overwrite masked data quality issues
- **mutation_surveyor**: Documented removal of `Concordance_auto_2` columns from HV2-3 transform output as an intentional breaking change; downstream consumers expecting these columns must be updated
- **mutation_surveyor**: New column constants in `columns.py` for trace selection details (PLAN.md Tasks B7-B8)

- **mutation_surveyor**: Three-scope concordance in `validation.py`: Full Profile, HV2-3, and HV1 concordance with per-scope discordance reason columns (PLAN.md Task B5)
- **mutation_surveyor**: New `_compare_ms_truth_region()` and `_classify_discordance_category()` helpers for region-scoped concordance comparison
- **mutation_surveyor**: `Discordance Category` column: Concordant, Full_profile_only, HV23_only, HV1_only, Both_HV23_and_HV1
- **mutation_surveyor**: Truth region columns: `Range by HV2-3 Region (Truth)`, `Variants by HV2-3 Region (Truth)`, `Range by HV1 Region (Truth)`, `Variants by HV1 Region (Truth)`
- **mutation_surveyor**: New column constants in `columns.py` for three-scope concordance output

- **mutation_surveyor**: `merge_profiles.py` now writes PC/NTC sheet in output xlsx alongside LID and Summary sheets (PLAN.md Task C2)

- **mutation_surveyor**: New `merge_files.py` utility with `merge_files()` function for concatenating multiple xlsx files with Batch and Source_File columns (PLAN.md Task D1)
- **mutation_surveyor**: CLI `src/tools/mutation_surveyor.py` now accepts `--lid-manifest` and `--control-ref` args for Layer 0/1 QC integration (PLAN.md Task D3)

- **docs**: `docs/tools/mutation_surveyor.md` updated with Layer 0/1 QC gates in pipeline flow (PLAN.md Task D4)
- **docs**: `docs/ARCHITECTURE.md` updated with Mutation Surveyor module table including inventory_qc, control_qc, merge_files, three-scope concordance (PLAN.md Task D4)

### Added
- **mutation_surveyor**: Explicit `"het"` flag in `build_hv1_flags()` and `build_hv23_flags()` when het variants are present (review-triggering, sample-level) — aligns with reference scripts
- **mutation_surveyor**: `source_het_flags()` detects `"HV1 source het"` and `"HV1 source polyC het"` from source trace variant data — aligns with reference scripts
- **validation**: `normalize_het_variant()` now normalizes `het_DEL` (multi-letter alleles) to `_hetDEL` canonical form — previously `16189het_DEL` was unchanged
- **validation**: `normalize_het_variant()` now normalizes `heteroplasmy` keyword to `het` before pattern matching — aligns with reference scripts
- **validation**: Y2 concordance label restored for compatible alleles with het on either side — previously het vs non-het was always F
- **flagging**: `flag_variants()` in core now detects `"Heteroplasmy at {pos}"` for het deletions (`seq="DEL"`) and het insertions (lowercase `seq`)

- **mutation_surveyor**: `HV2F_ADJUSTED` constant fixed from `[(73, 302)]` to `[(73, 303)]` — position 303 must belong to both HV2F and HV3R for correct overlap-based merging (PLAN.md Task A3)
- **mutation_surveyor**: Added `HV1_POLYC_MERGE = [(16180, 16195)]` constant for merge-stage polyC detection; `HV1_POLYC = [(16180, 16193)]` is kept for QC flag usage (PLAN.md Task A4)
- **mutation_surveyor**: `build_hv1_flags()` and `source_het_flags()` in flagging.py now use `HV1_POLYC_MERGE` for polyC variant detection
- **mutation_surveyor**: `_compute_trace_metrics()` calls in transform.py now pass `polyc_region=HV1_POLYC_MERGE` for HV1 trace scoring

- **mutation_surveyor**: `hv23_severity()` consolidated to single function in `flagging.py` returning AUTO/LOW/MEDIUM/HIGH with concordance-aware priority logic; removed duplicate `_hv23_severity()` from `transform.py` (PLAN.md Task B6)

- **mutation_surveyor**: Column names aligned with reference scripts: `Input Range`→`MS Range`, `C-shift Applied`→`C_shift_rule_applied`, `AC-repeat Applied`→`AC_repeat_rule_applied`, `HV1 PolyC Flag`→`HV1_PolyC_flag`, `HV23 PolyC Flag`→`HV2-3_PolyC_flag`, `HV1 Range Adj`→`Range by HV1 Region Adjusted`, `HV1 Variants Adj`→`Variants by HV1 Region Adjusted`, `HV23 Pair Status`→`HV2-3 Pair Status`, `HV23 Range Raw`→`Range by HV2-3 Regions`, `HV23 Flags`→`HV2-3 Flag list`, `Profile Range`→`Sample profile - Range`, and more (38 column constants renamed in `columns.py`) (PLAN.md Task C1)

### Changed
- **utils**: `parse_variant_tokens()` now normalizes `_HET` to `_het` in allele strings — aligns with reference scripts
- **mutation_surveyor**: `build_hv1_flags()` now accepts `hv1f_row` and `hv1r_row` parameters for source het detection
- **mutation_surveyor**: CLI now accepts `-i`/`--input-dir` (directory) instead of `-m`/`--ms` (file path); dynamically searches for the custom report matching the batch ID
- **mutation_surveyor**: `find_custom_report()` searches input directory for XLSX/TXT custom report files by batch ID with priority ordering
- **mutation_surveyor**: Automatic TXT→XLSX conversion step — if the found custom report is a `.txt` file, it is converted to `.xlsx` before pipeline processing
- **mutation_surveyor**: `--truth` now accepts a directory path and auto-searches for matching truthset files by batch ID
- **mutation_surveyor**: Shell script `scripts/tools/mutation_surveyor.sh` updated to use directory-based input with `--batch-id` argument and conditional comparison step
- **txt_to_excel**: Renamed `Txt_to_Excel.py` → `txt_to_excel.py` for PEP 8 module naming compliance
- **txt_to_excel**: Added `is_txt_custom_report()` helper to check whether a path looks like an MS TXT custom report
- **txt_to_excel**: Added `_derive_xlsx_path()` helper for clearer output path derivation
- **txt_to_excel**: Added logging via `loguru` in `convert_txt_to_excel()` for conversion progress tracking

### Changed
- **models**: `Variant.quality` type changed from `Optional[float]` to `List[float]` (default `[]`) to support per-trace quality scores matching Tracy format
- **mutation_surveyor**: Extract `$N` quality scores from variant tokens into `VARIANT_QUALITY` and `VARIANT_QUALITY_2` columns parallel to `VARIANT_CONVERTED` / `VARIANT_CONVERTED_2`
- **mutation_surveyor**: Thread quality and files through transform (`HV1_QUALITY_ADJ`, `HV23_QUALITY_ADJ`, `HV1_FILES`, `HV23_FILES`) and merge (`PROFILE_QUALITY`, `PROFILE_FILES`) stages
- **mutation_surveyor**: Populate `Variant.quality` and `Variant.files` from `PROFILE_QUALITY` and `PROFILE_FILES` in ETL
- **sample**: `_parse_variant_dict()` now normalizes quality: scalar→[scalar], list→list, None→[]
- **mtdna_merger**: `_flat_variant_dict()` preserves quality; `_merge_profiles()` unions quality lists for matching positions
- **variants**: `variant_to_dict()` always includes quality as a list
- **sequencher/etl**: Quality stored as list instead of scalar mean
- **postprocessing**: Quality comparison updated to always use list format

### Added
- **docs**: Created comprehensive human-readable documentation for all `src/` modules (excluding `src/modules/`):
  - `docs/config.md` — Configuration module (Pydantic Settings, tool paths, regions, parameters)
  - `docs/utils.md` — Core utilities (interval parsing, plotting, sequence handling, JSON generation)
  - `docs/variants.md` — Variant and position domain utilities (position contract, sequence generation, statistics)
  - `docs/analysis.md` — BLASTN consensus analysis (variant calling, indel normalization)
  - `docs/tracy.md` — Tracy integration (AB1 trace decomposition, variant calling, merging)
  - `docs/blastn.md` — BLASTN pipeline orchestration (AB1→variants→reports)
  - `docs/processing.md` — Pre/post-processing and positive controls (preprocessing, postprocessing, PC samples)
  - `docs/pipeline.md` — Pipeline orchestration (automate, manual, comparison pipelines)
  - `docs/generation.md` — Output generation (FASTA, reports, batch statistics)
  - `docs/core/models.md` — Core data models and Sample I/O (Position, Variant, Sample, ComparisonResult)
  - `docs/core/batch.md` — Batch processing (sequence generation, statistics, JSON serialization)
  - `docs/core/comparison.md` — Pairwise sample comparison (concordance, unique variants, export)
  - `docs/core/merger.md` — Basic JSON profile merging
  - `docs/core/uploader.md` — Google Sheets upload
  - `docs/core/sample.md` — Sample serialization/deserialization (sample_to_dict, load_sample_batch)
  - `docs/validation/validation_fasta.md` — FASTA validation against canonical code path
  - `docs/validation/validation_regenerate_results.md` — Regenerated results validation


### Fixed
- **batch**: Added `encoding="utf-8"` to JSON file writes in `Batch.write()` to ensure consistent output across platforms.
- **comparison**: Removed Excel column-width auto-sizing from `write_comparison_excel()` — openpyxl default widths are sufficient and the manual sizing was unnecessary complexity.
- **comparison**: Added warning log when `is_concordant()` receives a sample with no intervals — concordance may be unreliable without analyzed ranges.
- **merger**: Replaced `glob.glob` with `Path.glob` in `_load()` for type-safe file discovery. Removed Excel column-width auto-sizing from `_write_excel()` — openpyxl default widths are sufficient and the manual sizing was unnecessary complexity.
- **mtdna_merger**: Added validation that rCRS reference dict is non-empty after loading — raises `MergeValidationError` on empty/malformed reference files instead of silently producing incorrect results.
- **sample**: Added `computation_error` field to `sample_to_dict()` output when sequence generation or statistics fail — downstream consumers can now detect incomplete data instead of silently operating on defaults.
- **uploader**: Replaced broad `except Exception` blocks with specific exception types (`OSError`, `gspread.exceptions.APIError`, `ValueError`, `pd.errors.EmptyDataError`). Authentication failure now raises `RuntimeError` instead of returning `None`. Cell-by-cell fallback no longer silently writes empty strings on failure.

### Changed
- **core**: Extracted Sample I/O into `src/core/sample.py`. `sample_to_dict()`, `load_sample_batch()`, `_normalize_keyed_sample()`, and `_parse_variant_dict()` moved from `batch.py` and `comparison.py` into the new module. `Batch` class now delegates to `sample_to_dict()` for serialization; `comparison.py` imports `load_sample_batch` from `sample.py`. No behavior changes.
- **core**: Removed misleading V3/V4 variant format labels. The grouped variant format (`{snps, insertions, deletions}`) is the standard on-disk format; flat `List[Variant]` is the internal Sample model representation. The unreachable V4 flat variant branch in `mtdna_merger._flat_variant_dict()` was removed as dead code.

### Removed
- **core**: Removed unreachable V4 flat variant format code from `mtdna_merger._flat_variant_dict()` and misleading V3/V4 labels from `comparison._normalize_keyed_sample()` and `batch.py` docstrings. The grouped variant format ({snps, insertions, deletions}) is the standard on-disk format; flat List[Variant] is the internal Sample model representation.
- **sequencher_to_json**: Removed `src/modules/sequencher_to_json.py` — dead code with zero callers. The Sequencher ETL module (`src/modules/sequencher/`) produces Sample objects via the standard pipeline, making this manual JSON builder unnecessary.

### Changed
- **mtdna_merger**: Refactored to use `Sample` objects internally instead of raw dicts. The `merge()` method now returns a `MergeResult` (with `.sample` and `.warnings` attributes) instead of a raw dict. Output is serialized via `Batch.sample_to_dict()` for consistent format. Added `merge_samples()` convenience function for merging `Sample` objects directly.
- **mtdna_merger**: Output now includes per-sample JSON (`{sample_id}/{sample_id}.json`) alongside the existing `{sample_id}_mtdna_merged.json`, matching the `Batch.write()` output convention.

### Added
- **mtdna_merger**: Added `MergeResult` class with `.sample` (Sample) and `.warnings` (List[str]) attributes, replacing the previous raw dict return type.
- **mtdna_merger**: Added `merge_samples()` function for merging `Sample` objects directly without writing/reading JSON files.


### Added
- **mtdna_merger**: Added "Different Positions, Same Coverage Range" example to `docs/core/mtdna_merger.md` — documents that when two profiles cover the same range with variants at different positions, each triggers an implicit conflict (overlap consistency check).
- **mtdna_merger**: Created `docs/core/mtdna_merger.md` — conceptual reference documentation for the mtDNA profile merger, covering validation rules, merge flow, error handling, pipeline role, CLI usage, and future work linking to ROADMAP.
- **docs**: Added cross-references from `docs/ROADMAP.md` §4.1 and `docs/ARCHITECTURE.md` §2.5 to the new `docs/core/mtdna_merger.md`.
- **docs**: Created `docs/ROADMAP.md` — MS AutoRun + MS1 target workflow covering auto-pass criteria, three non-auto-pass scenarios, review packages, region-based merge via `mtdna_merger.py`, configurable tool priority, het handling, and TXT-to-Excel preprocessing.
- **docs**: Added cross-reference from `docs/ARCHITECTURE.md` §11 (Regenerate Gate) to `docs/ROADMAP.md`.
- **docs**: Added cross-references from `docs/tools/mutation_surveyor.md` §9 and §10 to `docs/ROADMAP.md`.
- **docs**: Added "Future Extensions" section to `docs/core/flagging.md` linking to `docs/ROADMAP.md`.

- **docs**: Created `docs/modules/mutation_surveyor.md` — comprehensive source code reference for all files in `src/modules/mutation_surveyor/`, covering public APIs, data models, column constants, flag classification, review logic, AC-repeat canonicalization, ETL token mapping, concordance validation, and pipeline orchestration.
- **docs**: Added cross-reference from `docs/tools/mutation_surveyor.md` to the new module-level doc.

### Fixed

- **docs/mutation_surveyor**: Corrected outdated references from `processs.py` to `etl.py` and from `processs` to `process` in `docs/tools/mutation_surveyor.md`.

### Changed

- **batch**: `information` dict in `statistic_fullbatch.json` now includes `source_tool` and `batch_id` (moved from top-level). Legacy `information["batch"]` key replaced by `batch_id`. Top-level `source_tool` and `batch_id` removed — both fields now live exclusively inside `information`. Deduplicated inline at serialization.
- **mutation_surveyor/etl**: Removed `information["batch"]` assignment; `batch_id` now propagated to `information` via `_enrich_information()` at serialization time instead of duplicating it in the ETL layer.
- **comparison**: `load_sample_batch()` now reads `source_tool` and `batch_id` from `information` dict when not present at top level, for backward compatibility with v4+ JSON.

### Fixed

- **mutation_surveyor/flagging**: `classify_review()` now uses `_NON_REVIEW_TRIGGER_NAMES` (which includes Autopass flags) instead of `MS_INFORMATIONAL_FLAG_NAMES` (which only included polyC). Autopass flags no longer incorrectly trigger Review classification.

### Changed

- **mutation_surveyor/transform**: `_review_trigger_flags_hv1()` now uses `_NON_REVIEW_TRIGGER_NAMES` for consistency with `MSFlag.is_review_trigger` and `classify_review()`.

### Removed

- **mutation_surveyor/flagging**: Removed `MS_INFORMATIONAL_FLAG_NAMES` (unused after classify_review fix) and `_INFORMATIONAL_FLAG_NAMES` legacy alias. Use `_NON_REVIEW_TRIGGER_NAMES` for review-trigger logic.
- **mutation_surveyor/utils**: Removed unused `is_common_hv23_insertion()` and `is_common_hv23_deletion()` (zero consumers after unexpected insertion/deletion flag removal).

### Changed


- **flagging**: Unified variant-level het flagging across tools. MS pipeline now normalizes `_het` tokens to IUPAC codes in `Variant.seq` (e.g. `309C_het` → `seq="Y"`), so `flag_variants()` from `src/core/flagging.py` detects `"Heteroplasmy at {pos}"` automatically. Removed duplicate `"het"` MSFlag generation from `build_hv1_flags()` and `build_hv23_flags()`. Removed `het_flags` tracking in `etl.py`.

- **flagging**: Replaced MS `"unexpected insertion"` and `"unexpected deletion"` flags with core's `"Insertion at {pos}"` and `"Deletion at {pos}"` (which include exception rules for 315.1, polyC C insertions, 249 deletion, and AC 523-524 pair). Removed flag generation from `build_hv23_flags()` and unused `is_common_hv23_insertion`/`is_common_hv23_deletion` imports.

### Removed

- **mutation_surveyor/flagging**: `"het"` MSFlag generation removed (het now detected via IUPAC codes by core `flag_variants()`).
- **mutation_surveyor/flagging**: `"unexpected insertion"` and `"unexpected deletion"` MSFlag generation removed (now handled by core `flag_variants()` with exception rules).
- **mutation_surveyor/etl**: `het_flags` list removed from ETL loop (het signal now in `Variant.seq` IUPAC codes).

### Changed

- **Consolidated transform modules**: `transform_hv1.py` and `transform_hv23.py` merged into a single `transform.py` with a `Transform` class. Shared trace metric computation deduplicated via `_compute_trace_metrics` returning `TraceMetrics` (Pydantic model). Each region has its own `run_hv1()` / `run_hv23()` method; region-specific logic stays in private methods. Pipeline calls updated to use `Transform(input_path, ...).run_hv1()` / `.run_hv23()`.

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added

- **Region status flags in sample_flags**: Each `Sample` now carries explicit autopass/review flags per HV region in `sample_flags`: `Autopass HV1`, `Autopass HV2 and HV3`, `Review HV1`, `Review HV2 and HV3`. Downstream tasks can check these flags for per-region dataflow routing without parsing combined review class strings.

### Removed

- **`is_auto_pass()`** from `mutation_surveyor/flagging.py` — replaced by direct `sample_flags` checks.
- **Pipeline Step 7** (auto-pass/review classification + `review_manifest.json`) from `mutation_surveyor/pipeline.py` — same information now available via `sample_flags` in `statistic_fullbatch.json`.
- **`information["review_class"]`** from `Sample` objects — replaced by region status flags in `sample_flags`.

### Changed
- **`mutation_surveyor/pipeline.py` and `src/core/batch.py` logging cleanup**: Use `logger.success` for completion messages with full paths. `logger.info` for step-start only. Removed f-strings in `batch.py` warnings, fixed `processs` typo, human-readable step labels (1/6–6/6).

- **Deslop mutation_surveyor module**: Removed 12 dead constants (IUPAC_AMBIGUITY_CODES, FLAGGING_REGION_START/END, AC_REPEAT_REGION_START/END, INSERTION_C_EXCEPTION_POSITIONS, MS_HV1_TARGET, MS_HV1_POLYC, MS_HV23_TARGET, MS_HV2F_ADJUSTED, MS_HV3R_ADJUSTED, MS_HV23_REQUIRED_SEGMENTS) that duplicated src/core/flagging.py or utils.py constants with zero consumers. Replaced MS_PHP_CODES with the already-imported PHP_CODES from utils.py. Removed duplicate _derive_batch_from_path definitions from pipeline.py and quality_control.py; inlined the 3-line filename-suffix-stripping logic at each call site instead of creating a shared utility for a trivial operation. Cleaned stale .pyc files for deleted modules.

- **Centralized Position type to src/core/models.py**: Moved `Position = Union[int, str]` type alias to `src.core.models` as the canonical source of truth, with `validate_position()` boundary validator that rejects `float` positions. Removed duplicate `Position` definitions from `src/variants.py` (now re-exported from models) and `src/modules/mutation_surveyor/utils.py` (now imported from models). All modules should import `Position` from `src.core.models` going forward.
- **Eliminated MSFlagger wrapper class**: Deleted the `MSFlagger` static-method-only class from `mutation_surveyor/flagging.py`. Its methods (`build_hv1_flags`, `build_hv23_flags`, `hv23_severity`, `hv23_review_class`, `classify_review`, `combined_review`, `is_auto_pass`) are now module-level functions. Also deleted 4 `_ms_*` private wrapper functions that just delegated to `utils.py`. All call sites updated to use direct imports and function calls.
- **Renamed Flag → MSFlag in mutation_surveyor module**: The `Flag` Pydantic model in `mutation_surveyor/flagging.py` is now `MSFlag` to disambiguate from `src.core.flagging.VariantFlag`. This makes the two different flag systems (per-reason MS flags vs per-variant core flags) explicitly distinct in code and imports.


### Fixed
- **Replaced 2 skipped integration tests with real unit tests**: `TestIntegrationRealData` tests (`test_ln_26_ab2227_hv2_insertion`, `test_ln_26_ab2227_all_regions_pass`) relied on unavailable `MS_090426_004` data and always skipped. Replaced with `TestMultiBaseInsertionAndAllRegions` using synthetic data via `sample_data_multi_ins` fixture, testing the same 10-base insertion alignment logic and all-region PASS validation without external file dependencies.

### Fixed
- **Fix TypeError in mutation_surveyor validation when sorting mixed Position types**: `_compare_variant_sets` used bare `sorted()` on a set containing both `int` base positions and `str` insertion positions (e.g., "309.1"), which raises `TypeError: '<' not supported between instances of 'str' and 'int'` in Python 3. Added `key=pos_sort_key` to the `sorted()` call, using the canonical position sort key that handles mixed int/str positions correctly.

### Fixed
- **Suppress Pydantic v2 deprecation warnings for model construction**: Explicitly pass all optional fields in `Sample()` and `Variant()` call sites. For `Sample`: `intervals`, `sample_flags`, `variant_flags`, `information`, `hv1`–`no_unread` in all three call sites (`mutation_surveyor/etl.py` lines 216 and 301, `sequencher/etl.py` line 198). For `Variant`: `quality`, `files`, `peaks` in all three call sites (`mutation_surveyor/etl.py` lines 88, 104, 108). Pydantic v2 warns when optional fields with defaults are omitted from construction; will become an error in a future version.

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Mutation Surveyor pipeline stages renamed** to purpose-driven names: `compare_qc` → `quality_control`, `hv23_merge` → `transform_hv23`, `hv1_merge` → `transform_hv1`, `lid_merge` → `merge_profiles`, `etl` → `processs`, `concordance` → `validation`. Module filenames, function names, and imports updated accordingly.
- **Mutation Surveyor documentation updated** to reflect parallel pipeline structure (HV1 and HV2+3 are independent branches after quality_control), optional `validation` as a downstream task (not terminal step), and `processs` as the terminal step producing `Sample` objects for GSTmtDNAv1 downstream consumption (comparison, FASTA generation, multi-tool merging).

### Fixed
- **Flag column format in `comparison.py`**: Changed Flag column from "Y"/"N" to "Yes"/"No" and Sample Flags from "Y"/"N" to "Yes"/"No", matching legacy `comparator.py` output format.
- **Variant Flags consolidation in `comparison.py`**: `_format_variant_flags()` now consolidates per-position "16180-16193 region" entries into a single "16180-16193 region (pos1, pos2, ...)" entry, matching the legacy `comparator.py` format. Other per-variant flags (Deletion at X, Insertion at X, etc.) remain per-position.

### Removed
- **"Heteroplasmy outside polyC" flag removed from `flagging.py`**: Removed sample-level and per-variant "Heteroplasmy outside polyC" flagging from `SampleFlagger.analyze()` and `flag_variants()`, as well as its level classification and technical-issues check. The helper method `get_heteroplasmies_outside_polyc()` is retained for potential future use.

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **FULL REGION display in `comparison.py`**: `format_intervals` now returns "FULL REGION" when all predefined HV regions are present, matching `comparator.py` behavior instead of listing individual ranges like "73-340 438-576 16024-16365".
- **Missing sample-level flags in Sequencher ETL**: `src/modules/sequencher/etl.py` now also runs `SampleFlagger.analyze()` to capture sample-level flags like "No 315.1 variant" and "Consecutive indels (5+)" that were missing from `Sample Flags` columns in comparison output. Previously only range flags and per-variant flags were included.
- **Filesystem validation in `comparison.py` CLI**: `_validate_input_file()` checks that input paths exist and contain valid JSON before processing. Invalid inputs produce clear Loguru error messages and `sys.exit(1)` instead of raw Python tracebacks.
- **Sample Flags columns in variant-level TSV**: `VARIANT_TSV_FIELDS` now includes `Sample Flags (A)` and `Sample Flags (B)` columns, surfacing sample-level flag reasons (e.g. "No 315.1 variant", "Consecutive indels") alongside per-position variant flags.
- **Sequencher ETL step in `pipeline.sh`**: Added `src/tools/sequencher.py` call in Step 4 to produce v4 `statistic_fullbatch.json` (with `sample_flags` + `variant_flags`) before running comparison.
- **Updated comparison CLI call in `pipeline.sh`**: Replaced legacy `comparator.py` call with `python -m src.core.comparison` using two-JSON-input pattern (`-a` automate, `-b` sequencher), matching `test.sh` Step 5.

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Per-sample JSON output** in `Batch.write()`: Each sample now also written to `output_dir/batch_id/{sample_id}/sample.json` with the same schema as its entry in `statistic_fullbatch.json`. Tool-agnostic — all pipelines (Sequencher, Tracy, blastn, mutation_surveyor) benefit automatically.
- **`generate_sequence` function** in `src/variants.py`: New function accepting a `Sample` model and returning a structured `Sequence` result object. Derives `seq_regions` from `sample.intervals`, `regions` from config defaults, and `ref_path` from config when not provided. Replaces the multi-param tuple-returning `gen_seq_from_variants_regions` for new callers.
- **`batch.py` migrated** to use `generate_sequence` instead of `gen_seq_from_variants_regions`.

### Fixed
- **Sequencher comparator sample ID corruption**: `compare_results()` now iterates over pipeline sample IDs (source of truth) instead of manual TXT filenames. Added `find_manual_file()` to match clean sample IDs to manual filenames using non-alphanumeric separator detection. This fixes corrupted sample names like `LN_26_AA4522 73`, `LN_26_AA4523.73`, `LN_26_AA4524_73` that appeared when Sequencher filenames embedded range positions with inconsistent delimiters (space, dot, underscore). Pipeline-only rows are now shown for samples without a matching manual file.
- **`statistic_variants` insertion-count bug**: `ref=="-"` branch now correctly increments `no_ins` instead of `no_dels`.
- **`Sequence` model dict types**: Changed `ref_dict` and `con_dict` from `Dict[int, str]` to `Dict[Union[int, str], str]` to support insertion position keys like `"309.1"` and `"315.1"`.

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **MS flagging unified in `src/core/flagging.py`**: Added `MSFlagger` class with `build_hv1_flags()`, `build_hv23_flags()`, `classify_review()`, `combined_review()`, and `is_auto_pass()`. Refactored `hv1_merge.py` and `hv23_merge.py` to delegate flag building to `MSFlagger` instead of inline logic. Refactored `lid_merge.py` to use `MSFlagger.combined_review()`. Added auto-pass/review routing in `process_batch()` with review manifest JSON output. Added MS-specific constants (`MS_INFORMATIONAL_FLAG_NAMES`, region boundaries, `MS_PHP_CODES`). All flags now flow through `src/core/flagging.py` as the single source of truth.
- **Batch output directory**: `Batch.write()` now takes `output_dir` and `batch_id` args, writing to `{output_dir}/{batch_id}/statistic_fullbatch.json` instead of a flat path. Sequencher pipeline updated accordingly; test.sh comparison path updated accordingly.
- **docs/tools/mutation_surveyor.md**: Complete rewrite (previously `mutation_suveyor.md` with typo filename). Now reflects the actual implemented Python pipeline: documents MS custom report xlsx as input (not AB1 files), 6-stage pipeline with module/function references, src/core integration (Sample, Variant, Flag, FlagLevel, Batch, flag_variants), structured Flag flow, ETL → Batch.write() → JSON path, CLI/API usage, intermediate file naming, and review classification rules. External workflow steps (AB1 processing, manual review, GSTmtDNAv1) moved to out-of-scope section.

### Removed
- Removed `check_type_variant` from `src/variants.py` (dead code; zero callers; classification logic duplicated by `split_variant_types`)
- Moved `format_variants_list` from `src/variants.py` to `src/modules/NGS/utils.py` (NGS-only utility); updated imports in `compare_fis_sanger_direct.py` and `compare_fis_fasta_direct.py`

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- Restored `flatten_variant_dict`, `variant_to_simple_dict`, `variants_to_grouped`, and `_quality_to_list` in `src/variants.py` (uncommitted functions referenced by callers but missing from committed HEAD)

### Fixed
- **TNLS verification module** (`src/modules/TNLS/verification/`):
  - `utils.py`: Corrected IUPAC ambiguous base expansions (R, Y, K, M, S, W were inverted); replaced `list` with `set`; removed dead debug `print()`; cleaned `compare_lists()` dead variable
  - `compare_TNLS_HCLS.py`: Removed `if processed_samples:` guard that skipped the first HCLS sample from comparison
  - `compare_TNLS_HCLS_stats.py`: Replaced hardcoded `CDI_20260301_...` base filename with `--base-samples` CLI arg
  - `compare_1n_class.py`: Replaced all trivial getter methods with direct `self.attr` access; inlined `get_base_region()`, `get_target_region()`, `extract_hv_regions()` into the stage methods; renamed private helpers to `_get_*` convention; `get_matching_rule()` removed as dead code; removed commented-out `convert_to_list_of_dicts()` stub
  - `README.md`: Corrected `minimum` value from 300 to 150 in table and `MR_MINIMUM_NOT_REACH` description; documented `--run-step0` and `--skip-filter` flags
  - `verification.sh`: All paths are now `SCRIPT_DIR`-relative (not hardcoded absolute paths); added `--run-step0` flag to regenerate HCLS base JSON from raw data; consolidated `HCLS_BATCH` variable; Step 1 runs only on `--run-step0`; moved `rm -rf results/*` before Step 1 so it doesn't wipe pre-existing data needed for `--skip-filter`; input validation prints clear error messages; wired `--base-samples` argument to Step 4

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **TNLS verification module** (`src/modules/TNLS/verification/`):
  - `docs/verification.md`: Comprehensive documentation covering module purpose, file inventory, 3-stage matching logic, matching rules reference, CLI usage, output guide, data schemas, IUPAC table, HV region boundaries, and known limitations
- `src/generate_variant_summary.py`: new script that aggregates `all_batches_sequencher_comparison.tsv` into a per-variant summary TSV (`results/variant_summary_report.tsv`). Each row = one unique variant; columns include `Total_in_Sequencher`, `Pipeline_Correct`, `Pipeline_Accuracy_Pct`, `Total_in_Pipeline`, `Pipeline_Extra`, and `Skip` (Yes if accuracy >= 95%). Supports `--input`, `--output`, `--min-accuracy` CLI flags.
- **Filename Validation for Comparator**: Added validation function for Sequencher filename ranges
  - Validates against 4 rules: monotonic, valid positions, single region, valid range
  - Added "Filename Flag" column to comparator output
  - Invalid filenames show "Flag: Range – Possibly wrong" and hide variant results
  - New tests in `tests/test_comparator.py` with 32 test cases
- **Smart Caching for merge_data.py**: Added persistent cache for directory scans
  - Cache stored in `{RESULTS_DIR}/merged/scan_cache.json`
  - Uses mtime-based invalidation to detect modified batches
  - Only re-scans batches that have changed since last run
  - New `--force-refresh` flag to bypass cache

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Improved regenerate directory scanning**: Now updates cache with sample IDs for faster subsequent runs
- **Consolidated duplicate `sequencher_to_json.py`**:
  - Removed `src/sequencher_to_json.py`; canonical version is now `src/modules/sequencher_to_json.py`
  - Added `__init__.py` to all `src/modules/` subdirectories (`LAB`, `NGS`, `PYQG`, `TNLS`, `sequencher`, `TNLS/verification`) to make them proper Python packages
  - Fixed `sys.path` manipulation in `src/modules/sequencher_to_json.py` to use `Path(__file__).resolve().parents[N]` instead of `os.getcwd()`, so imports work from any working directory
  - Updated shell scripts (`src/modules/TNLS/verification/verification.sh`, `scripts/test.sh`) to use the new path
  - Fixed pre-existing unused-variable lint issue in `src/modules/LAB/generate_mtdna_html_report.py`
  - Added **Modules & Packages** section to `CLAUDE.md` documenting the two import styles (`from src.*` vs bare `from config`) and the correct `sys.path` pattern

### Fixed
- **TNLS verification dict-key collision** (`src/modules/TNLS/verification/utils.py`):
  - `combined_variants()` now uses typed dict keys (`snp_<pos>`, `ins_<pos>`, `del_<pos>`) to prevent silent overwrites when different variant types (SNP/insertion/deletion) share the same genomic position
  - `is_in_overlap_range()` now parses the typed key's numeric suffix and is guarded against malformed keys
  - Previously, a SNP and deletion at the same integer position would silently collide (later overwrites earlier), producing false `CANNOT_EXCLUDE` conclusions
  - Added 8 regression tests in `tests/verification/test_utils.py`
- **generate_reports.py metadata lookup**: Now uses the `-d` data_dir argument for metadata file resolution instead of always falling back to `DATA_DIR` env var
  - When a custom data directory is provided via `-d`, its `metadata/` subdirectory is used for Excel/TSV lookup
  - Falls back to `DATA_DIR` env var, then `base_dir/data` in order when `-d` is not provided
- **generate_reports.py JSON lookup**: Added separate `--results_dir` argument for locating sample JSON files
  - Previously the single `-d` argument was used for both metadata (needs `data/metadata/`) and JSON (needs `results/.../regenerate/{BATCH_ID}/`)
  - This caused JSON lookups to fail when using a custom data dir because JSON files live in RESULTS_DIR, not DATA_DIR
  - `-d` / `--data_dir` → metadata manifest (optional, falls back to DATA_DIR)
  - `-r` / `--results_dir` → sample JSON files (required)
  - `--report_dir` → PDF output (unchanged)
- **pipeline.sh step 6**: Updated to pass `--results_dir ${MANUAL_REGENERATE_DIR}` in addition to `-d ${DATA_DIR}`

## [0.6.5] - 2026-03-17

### Fixed
- **Reference Mismatch Warnings**: Removed unnecessary validation that compared upstream `ref` field to rCRS
  - The upstream `ref` field could differ from rCRS due to strand orientation or different reference sources
  - Now always uses rCRS reference base from the loaded reference dictionary
  - Removed misleading warning from documentation

## [0.6.4] - 2026-03-16

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Multiple Ranges Support**: Added support for multiple ranges per region in FASTA generation
  - Regions can now have multiple separate ranges (e.g., HV2: [[79, 302], [322, 340]])
  - Gaps between ranges are filled with 'N' characters
  - FASTA headers show all ranges (e.g., `>HV2|79-302,322-340`)
  - New test case `test_multiple_ranges_per_region` added

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Updated `gen_seq_from_variants_regions`**: Modified to handle list of ranges per region
- **Updated `write_fasta_file`**: Now accepts `region_ranges` parameter for multiple ranges in header

## [0.6.3] - 2026-03-13

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Documentation**: Added cross-references between pipeline.md and analysis.md for better navigation

## [0.6.2] - 2026-03-13

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Complex Variant Tests**: Added comprehensive tests in `tests/dev/test_generate_fasta_v2.py` for complex insertion/deletion scenarios
  - `TestComplexVariants` class with 5 new test methods
  - Tests for multiple insertions at different positions (309.1, 309.2, 315.1, 525.x, 16189.1)
  - Tests for multiple deletions at various positions (333, 555, 16111)
  - Tests for mixed variants (SNPs + insertions + deletions)
  - Tests for insertion sequence integrity and ordering
- **Test Fixture**: Added `sample_with_complex_variants` fixture in `tests/conftest.py`

## [0.6.1] - 2026-03-13

### Fixed
- **Bug 2 - Reference Base Validation**: Added validation in `apply_variants_to_dict()` to check that variant's ref matches reference genome, logs warning on mismatch
- **Bug 3 - N-Marking for Non-Regions**: Added `mark_non_region_as_n()` method to mark positions outside specified regions as 'N'
- **Bug 4 - Wrong Range in con_str Generation**: Fixed range calculation to properly handle fractional insertion keys (filter to integer keys only)
- **Bug 5 - Multi-Base Insertions**: Modified `extract_region()` to iterate through each character in multi-base insertions

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **New Tests**: Added comprehensive tests for bug fixes in `tests/dev/test_generate_fasta_v2.py`
  - Tests for reference base validation
  - Tests for N-marking outside regions
  - Tests for con_str generation with fractional keys
  - Tests for multi-base insertions

## [0.6.0] - 2026-03-13

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Simplified FASTA Generation**: Rewrote `src/dev/generate_fasta_v2.py` using dictionary-based approach
  - Replaced list-based consensus with dictionary using position keys
  - SNPs replace values, deletions marked with '-', insertions use fractional keys (e.g., 309.1)
  - Region extraction now uses simple range query + sort (no complex offset tracking)
  - Much cleaner code: ~370 lines → ~250 lines

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **New Tests**: Added comprehensive test suite in `tests/dev/test_generate_fasta_v2.py`
  - Tests for dictionary-based variant application
  - Tests for region extraction with insertions
  - Edge case tests for boundary regions, non-contiguous insertions

## [0.5.0] - 2026-03-12

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Integrated Batch Pipeline**: Added `--rerun` flag support to `batch_pipeline.sh`
  - Normal mode: `bash scripts/batch_pipeline.sh` - runs all samples in batch
  - Rerun mode: `bash scripts/batch_pipeline.sh --rerun` - filters to samples with manual data
  - ZIP extraction and sample list generation for rerun mode
  - Skips individual batch uploads in rerun mode (only uploads merged results)
  - Added `--help` flag for usage information

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Batch Pipeline Consolidation**: Merged `batch_pipeline_rerun.sh` functionality into `batch_pipeline.sh`
  - Deleted redundant `scripts/batch_pipeline_rerun.sh`
  - Single unified script with mode detection via `--rerun` flag
  - Hardcoded batch lists for both normal and rerun modes

### Documentation
- **Updated docs/pipeline.md**: Added Batch Pipeline section documenting `batch_pipeline.sh`
- **Updated pipeline steps**: Renumbered from 0-7 to 1-6 to match actual pipeline.sh
  - Step 1: Validate
  - Step 2: Prepare
  - Step 3: Automate (combined BLASTN + Tracy + merge/regenerate)
  - Step 4: Manual
  - Step 5: FASTA
  - Step 6: Report
- **Consolidated Pipeline Documentation**: Merged `docs/tracy.md` into `docs/pipeline.md`
  - Added documentation for `-p` flag (blastn vs tracy mode)
  - Included Tracy-specific parameters and workflow differences
  - Updated step table to reflect unified pipeline
  - Deleted redundant `docs/tracy.md` file

## [0.4.1] - 2026-03-11

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Pipeline Consolidation**: Removed legacy pipeline scripts and consolidated into single `pipeline.sh`
  - Removed: `pipeline.sh`, `pipeline_tracy.sh`, `pipeline_sample.sh`
  - Renamed: `pipeline_unified.sh` → `pipeline.sh`
  - Updated all documentation references to use `pipeline.sh`
  - `batch_pipeline.sh` now uses `pipeline.sh -p tracy` for Tracy mode

### Fixed
- **Isolated Run Support**: `-d` option now requires `-l` for proper isolated runs without LAB_DATA_DIR access

## [0.4.0] - 2026-03-10

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Unified Pipeline Script**: New `scripts/pipeline_unified.sh` consolidating all pipeline functionality into one script
  - `-p`/`--pipeline` flag to select analysis engine: `blastn` (default) or `tracy`
  - Simplified 6-step system: `validate`, `prepare`, `automate`, `manual`, `fasta`, `report`
  - Steps selectable by number (`-s 2,3,4`) or name (`-s prepare,automate,manual`)
  - `-l`/`--sample-list` for custom sample lists (works with both pipeline modes)
  - Directory existence validation for custom paths (`-d`, `-m`, `-o`)
  - Pipeline mode and batch ID logged at startup for traceability

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **README.md**: Updated usage documentation to reflect `pipeline_unified.sh` as the primary entry point

## [0.3.0] - 2026-03-06

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Serena Integration**: Added initial instructions and project configuration for Serena AI agent
- **CLAUDE.md Updates**: Added Pre-Work Checklist, Code Quality section, and Changelog guidelines
- **Documentation**: Renamed all docs/*.md files (removed `_documentation` suffix)
- **prepare.py to TNLS.sh**: Added `src/modules/TNLS/prepare.py` to TNLS.sh workflow

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Module Reorganization**:
  - Moved `generate_fasta_v2.py` to `src/dev/` (development version)
  - Renamed `src/modules/compare_sequencher/` to `src/modules/sequencher/`
  - Renamed `compare_sequencher.py` → `comparator.py`
  - Renamed `merge_comparisons.py` → `merger.py`
  - Renamed `upload_to_gsheets.py` → `uploader.py`
  - Updated all script references to use new module paths
- **Logging Migration**: Replaced 55 `print()` statements with Loguru in:
  - `src/modules/NGS/variant_frequency_analyzer.py`
  - `src/modules/NGS/calculate_variant_freq.py`
  - `src/modules/TNLS/merge_family_info.py`
- **Test Reorganization**: Moved test files to mirror source structure
  - Moved `tests/test_generate_fasta_v2.py` to `tests/dev/` to match `src/dev/` location
  - Updated import path to `from src.dev.generate_fasta_v2 import ...`

### Removed
- **Unused NGS Modules** (6 files):
  - `FIS_single_excel_variants.py`
  - `compare_ngs_variants.py`
  - `compare_variants_simple.py`
  - `convert_variants_to_json.py`
  - `location_variant_frequency_analyzer.py`
  - `process_variant_frequency.py`

### Fixed
- **Linting Errors**: Fixed all ruff linting errors (F401, F821, F841) across the codebase
  - Added `__all__` to `src/utils/__init__.py` to fix unused import warnings
  - Fixed undefined `mapping` variable in `src/modules/PYQG/PYQG.py`
- **Utils Module Renaming**: Renamed utils files to remove `_utils` suffix
  - `command_utils.py` → `command.py`
  - `io_utils.py` → `io.py`
  - `plotting_utils.py` → `plotting.py`
  - `sequence_utils.py` → `sequence.py`
  - `variant_utils.py` → `variant.py`
  - Updated all imports across the codebase
- **Import Fix**: Fixed old import in `src/utils/plotting.py`
  - Changed `from .variant_utils import count_index` → `from .variant import count_index`
- **Performance Fix**: Optimized `generate_fasta.py` insertion logic
  - Changed character-by-character list insertion to slice assignment (O(n*m) → O(n))
  - This significantly improves performance for samples with long insertions

## [0.2.0] - 2026-03-05

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **Code Quality**: Integrated `ruff` linter and formatter for Python code
- **Testing**: Added `pytest` framework with test configuration in `pyproject.toml`
- **Shell Linting**: Ran `shellcheck` across all shell scripts and fixed warnings
- **Ruff Configuration**: Added `[tool.ruff]` and `[tool.ruff.lint]` sections in `pyproject.toml`
  - Line length set to 120
  - Target Python 3.12
  - Rules: `E` (errors), `F` (pyflakes), `I` (isort)
  - Ignored: `E501` (line too long), `E402` (module-level import order)

### Fixed
- Removed duplicate `get_highlight_position` function in `src/utils.py`
- Fixed undefined variable `mapping` reference in `src/modules/PYQG/PYQG.py` log statement
- Removed unused variable assignments across multiple modules (`F841` fixes)
- Fixed unquoted variables and command substitutions in shell scripts
- Added missing `#!/bin/bash` shebangs to scripts that lacked them
- Fixed unsafe `cd` commands without `|| exit` fallback in `scripts/setup.sh`
- Replaced `$?` exit code checks with direct `if cmd; then` patterns

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- Auto-formatted all Python files with `ruff format`
- Reorganized imports with `isort` rules via `ruff`
- Updated `pyproject.toml` dependencies to use `uv add` with latest versions
- Updated batch pipeline rerun list with current batch IDs

## [0.1.0] - Initial Release

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- Core pipeline scripts (`pipeline.sh`, `pipeline_tracy.sh`, `pipeline_sample.sh`)
- Batch processing (`batch_pipeline.sh`, `batch_pipeline_rerun.sh`)
- Python analysis modules: `analysis.py`, `blastn.py`, `tracy.py`, `utils.py`
- FASTA generation and report generation
- Manual vs automated pipeline comparison
- Sequencher comparison module
- NGS comparison modules (FIS analyzer, variant comparison)
- TNLS modules (kinship analysis, family variant filtering)
- PYQG blinding module
- Setup script with UGENE and Tracy tool downloads
- Docker support via `scripts/docker.sh`

### Changed
- **Migrated MS module from float to canonical Position type (Union[int, str])**: All position types in `src/modules/mutation_surveyor/` now use `int` for base positions (309, 16569) and `str` for insertions ("309.1", "217.10"), matching `src/variants.py`. Region constants (`HV1_TARGET`, `HV23_TARGET`, etc.) changed from `tuple[float, float]` to `tuple[int, int]`. Functions `in_any_range`, `normalize_ranges`, and `token_sort_key` now delegate to `src/variants.is_position_in_intervals`, `src/utils.merge_intervals`, and `src/variants.pos_sort_key` respectively, eliminating duplicated logic. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validate_concordance.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validate_concordance.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`. Additional float remnants fixed: `quality_control.py` `build_variant_converted_2` and `polyc_flag` now use `normalize_position()` + `pos_base()` instead of `float()` for position parsing and range comparison; `validation.py` `_variant_position` now returns canonical `Position` type via `normalize_position()` instead of `float`; `validation.py` profile variant dicts use `dict[Position, str]` instead of `dict[float, str]`; `transform_hv23.py` `trim_trace_adjusted` uses `dict[tuple[Position, str], str]` and `_safe_numeric_sort` uses `normalize_position()` + `pos_sort_key()` for canonical Position-aware sorting; `transform_hv23.py` PHP site sorting uses `pos_sort_key(normalize_position(s))` instead of `int(float(s))`.
- **Deduplicated `mismatch_positions_between_two_rows`** in Mutation Surveyor module: moved generalized version to `src/modules/mutation_surveyor/utils.py`, removed duplicate definitions from `transform_hv1.py` and `transform_hv23.py`. HV1 now calls the shared function with `variant_key="HV1 Trace Variants Adjusted"` and `range_ranges1=range_ranges2=HV1_TARGET`.
- **Sequencher module deslop**: Removed dead code from `src/modules/sequencher/utils.py` — deleted `SEQUENCHER_HEADER_ROWS`, `ALWAYS_CONCORDANT_POSITIONS`, `IUPAC_AMBIGUITY_CODES`, `FLAGGING_REGION_START/END`, `MAX_COLUMN_WIDTH` (canonical locations: `src/core/comparison.py` and `src/core/flagging.py`), `safe_extract_sample_id` (duplicate of `src/manual_pipeline.py`), `ranges_to_string`, and `parse_manual_txt` (only used by excluded comparator.py). Removed unused `import pandas as pd` and `import normalize_position` from `src.variants`.

- **Fixed stale column name references after rename**: `validation.py` `_find_column()` still looked for old names (`"Sample profile - Range"`, `"Sample profile - Variants"`) instead of new constants (`PROFILE_RANGE`, `PROFILE_VARIANTS`). `etl.py` hardcoded `"Sample profile - Variants"` to read variant data, producing empty variant lists and zero concordance. Both files now use `columns.py` constants with old names kept as fallback candidates in `_find_column()`.
### Added
- **`docs/tools/sequencher.md`**: New tool-specific documentation for the Sequencher module covering public API, file structure, integration points, and filename format conventions.

- **Renamed etl.process_batch → etl.build_samples**: Eliminated name collision between `mutation_surveyor/etl.py:process_batch` (lid_merge→Samples adapter) and `mutation_surveyor/pipeline.py:process_batch` (full pipeline orchestrator). The ETL function is now `build_samples`, matching its actual responsibility. Updated `pipeline.py` import and call site. Updated module docstring to clarify entry points.


## [Unreleased] — 2026-06-03 — Docs Sync

### Changed
- **docs**: Synced all documentation files with current source code APIs, method signatures, and behavior
- **docs**: Updated `docs/analysis.md` — added missing methods: `analyze_polyc_context`, `is_after_polyc_region`, `is_before_polyc_region`, `load_quality_scores`, `call_variants`, `right_norm_indel`
- **docs**: Updated `docs/tracy.md` — added missing methods: `create_variant`, `convert_variant`, `update_variant`, `update_variant_data`; expanded private methods documentation
- **docs**: Updated `docs/processing.md` — expanded `PostProcess` API with all methods: `load_filter_positions`, `load_no_filter_positions`, `is_special_position`, `filter_special_positions`, `get_regions_containing_position`, `apply_region_validation_rules`, `is_in_polyc_region`, `is_after_polyc_region`, `is_polyc_stretch`, `filter_poliC`, `parse_variant_rule`, `check_variant_exist`, `_check_rule_conditions`, `_apply_rule_changes`; updated `PCProcessor` method names to match source
- **docs**: Updated `docs/pipeline.md` — added missing methods to `MtDNAPipeline` (`setup_directories`, `run_pipeline`), `ManualPipelineProcessor` (`process_manual_results`, `process_sample_pipeline_results`, `diff_manual_pipeline`, `format_diff`, `filter_variants_by_intervals`), `VariantProcessor` (`analyze_position_reverse`, `process_blastn_variants`)
- **docs**: Updated `docs/utils.md` — added missing functions: `get_start_pos`, `reverse_complement`, `read_fastq`, `preprocess_seq`, `check_sequence_get_start_pos`, `get_highlight_positions`, `save_plt`, `get_annotation_pos`, `iupac_to_bases`, `get_blastn_path`, `check_output_path`, `count_index`, `check_folder_exist`, `annotate_variants`, `annotate_positions`, `plot_chromatograph`, `annotate_peaks`, `highlight_peak`, `annote_mutation`, `get_highlight_position`, `get_align_pos`, `get_list_positions`
- **docs**: Updated `docs/generation.md` — added async methods to `ReportGenerator` (`get_pdf_link`, `download_pdf`, `process_sample`, `generate_reports`)
- **docs**: Updated `docs/core/batch.md` — added `to_json()` method
- **docs**: Updated `docs/core/comparison.md` — added `format_intervals`, `compare_batches`, `summary_row`, and module-level functions (`create_variant_level_data`, `write_comparison_json`, `write_comparison_tsv`, `write_variant_level_tsv`, `write_comparison_excel`, `compare_batch_files`)
- **docs**: Updated `docs/core/flagging.md` — added `VariantAnalyzer` methods: `check_combined_conditions`, `has_variant_at`, `has_523_524_ac_deletion`, `has_513a_with_523_524_deletion`; added `SampleFlagger` private methods
- **docs**: Updated `docs/core/models.md` — synced all model fields with source code; added `Sequence`, `validate_position` documentation
- **docs**: Updated `docs/core/sample.md` — synced `sample_to_dict` and `load_sample_batch` with source
- **docs**: Updated `docs/core/mtdna_merger.md` — added all private methods and `merge_samples` function
- **docs**: Updated `docs/core/merger.md` — synced with source code
- **docs**: Updated `docs/variants.md` — added `generate_sequence` function documentation
- **docs**: Updated `docs/config.md` — verified all settings fields match source
- **docs**: Updated `docs/blastn.md` — verified API accuracy
- **docs**: Updated `docs/ARCHITECTURE.md` — added source module index table, updated date

### Added
- **docs**: New `docs/modules/LAB/generate_mtdna_html_report.md` — LAB HTML report generation
- **docs**: New `docs/modules/NGS/README.md` — NGS module overview
- **docs**: New `docs/modules/PYQG/README.md` — PYQG module overview
- **docs**: New `docs/modules/TNLS/README.md` — TNLS module overview
- **docs**: New `docs/modules/statistics/README.md` — Statistics module overview
- **docs**: New `docs/modules/fasta_to_json.md` — FASTA to JSON converter
- **docs**: New `docs/modules/mutation_surveyor_inventory_qc.md` — MS Inventory QC documentation
- **docs**: New `docs/modules/mutation_surveyor_control_qc.md` — MS Control QC documentation
- **docs**: New `docs/modules/mutation_surveyor_txt_to_excel.md` — MS TXT to Excel converter
- **docs**: New `docs/modules/sequencher_comparator.md` — Legacy Sequencher comparator
 - 2026-05-29

### Changed

- **Column naming standardization in Mutation Surveyor module** — All DataFrame output columns now use concise, consistent names defined as constants in `src/modules/mutation_surveyor/columns.py`. Key renames:
  - Profile prefix: `Sample Profiles - Range` / `Sample Profile - Range` / `Sample profile - Range` → `Profile Range`
  - Profile variants: `Sample profile - Variants` → `Profile Variants`
  - Profile flags: `Sample Profile - Flag` / `Sample Profile - MSFlag` → `Profile Flags`
  - Flag lists: `HV1 Flag list` → `HV1 Flags`, `HV2-3 Flag list` → `HV23 Flags`, `Combined Flag list` → `Combined Flags`, `HV1 Review-triggering Flag list` → `HV1 Review Flags`
  - Region columns: `Range by HV2-3 Regions Adjusted` → `HV23 Range Adj` (and similar for Variants/Mismatch, HV1)
  - Separator convention: `HV1_PolyC_flag` → `HV1 PolyC Flag`, `HV2-3_PolyC_flag` → `HV23 PolyC Flag`, `Error_type` → `Error Type`, `C_shift_rule_applied` → `C-shift Applied`, `AC_repeat_rule_applied` → `AC-repeat Applied`
  - Misc: `# of Mutation` → `Mutation Count`, `MS Range` → `Input Range`, `HV3R PHP 341_437 Count` → `HV3R PHP Count`, `HV1F Selected Sample Name` → `HV1F Sample`, etc.
- Added `src/modules/mutation_surveyor/columns.py` with all column name constants and canonical column order lists (`QC_COLUMNS`, `HV23_TRACES_COLUMNS`, `HV23_LID_COLUMNS`, `MERGE_LID_COLUMNS`).

### Changed
- **docs**: Synced documentation with source code to fix outdated API references:
  - **variants.md**: Added `generate_sequence`, `check_type_variant`, `format_variants_list`, `parse_variant_position_allele`; moved `get_intervals_from_json` to utils.md where it actually lives
  - **utils.md**: Added `get_intervals_from_json` function (was incorrectly listed in variants.md)
  - **core/models.md**: Added missing `Sample` fields (`no_snps`, `no_ins`, `no_dels`, `no_identicals`, `no_unread`, `variant_flags`); added `ToolResult` and `Sequence` models; fixed `ComparisonResult` schema (now uses `ToolResult` with `unique_variants_a/b` instead of old `sample_id_a/b` with `variant_comparison`)
  - **core/comparison.md**: Added `compare_batches`, `summary_row`, `create_variant_level_data`, `write_comparison_json`, `write_variant_level_tsv`, `compare_batch_files`; fixed `ComparisonResult` output format; changed export methods from class methods to module-level functions
  - **core/flagging.md**: Added `VariantFlag` dataclass, `VariantAnalyzer` class, `SampleFlagger` methods, `flag_variants()` function signatures
  - **core/merger.md**: Replaced incorrect `Merger` class description with actual module-level `merge_excel_files` function; fixed merge strategy description (Excel file merging, not JSON)
  - **core/mtdna_merger.md**: Added `merge_samples()` convenience function and helper functions; fixed `MergeResult.sample` field (was `profile: Dict`, now `sample: Sample`)
  - **generation.md**: Replaced `generate_fasta()` function with `FastaGenerator` class; replaced `generate_sample_reports()` function with `ReportGenerator` class

## [Unreleased] — 2026-06-04 — Extract QC/Validation from MS Pipeline

### Changed
- **mutation_surveyor**: Extracted Layer 0 (inventory_qc), Layer 1 (control_qc), concordance validation, and AB1 sorting from `process_batch()` into standalone `python -m` CLI entry points
- **mutation_surveyor**: Core pipeline (`process_batch`) now contains only: quality_control → transform → merge → ETL → Batch
- **mutation_surveyor**: Removed `--lid-manifest`, `--control-ref`, `--truth`, `--data-dir` CLI flags from `mutation_surveyor.py` tool entry point
- **mutation_surveyor**: Removed `truth_path`, `lid_manifest_path`, `control_ref_path`, `data_dir`, `max_workers` params from `process_batch()`
- **mutation_surveyor**: Updated step numbering from 1-7 to 1-5 in `process_batch()`
- **mutation_surveyor**: Restructured `scripts/tools/mutation_surveyor.sh` to orchestrate all 5 steps independently with proper error handling

### Added
- **mutation_surveyor**: Added `main()` + `if __name__ == "__main__"` CLI to `inventory_qc.py` — standalone Layer 0 blocking gate via `python -m src.modules.mutation_surveyor.inventory_qc`
- **mutation_surveyor**: Added `main()` + `if __name__ == "__main__"` CLI to `control_qc.py` — standalone Layer 1 control QC via `python -m src.modules.mutation_surveyor.control_qc`
- **mutation_surveyor**: Added `main()` + `if __name__ == "__main__"` CLI to `validation.py` — standalone concordance validation via `python -m src.modules.mutation_surveyor.validation`
- **mutation_surveyor**: Added `main()` + `if __name__ == "__main__"` CLI to `sort_review_data.py` — standalone AB1 sorting via `python -m src.modules.mutation_surveyor.sort_review_data`
- **mutation_surveyor**: Added `--output-dir` parameter to `sort_review_data()` for separate output directory

## [Unreleased] — 2026-06-07 — Clean Up 100%-Confidence Dead Code

### Removed
- **analysis**: Removed unused `end_ref` parameter from `call_variants()` — parameter was accepted but never used in the function body
- **mutation_surveyor**: Removed unused `**kwargs` from `etl.process()`, `etl.build_samples()`, `pipeline.process_sample()`, and `pipeline.process_batch()` — forward-compatibility kwargs were never forwarded or used by callers
- **sequencher**: Removed unused `**kwargs` from `etl.process()` — forward-compatibility kwargs were never forwarded or used by callers
- **TNLS**: Removed unreachable duplicate `return` statement in `filter_family_variants.py` — second `return " ".join(canonical_variants)` was dead code after the first return

### Changed
- **analysis**: Updated `call_variants()` call site to remove the `forward_ref_end` argument that corresponded to the removed `end_ref` parameter

## [Unreleased] — 2026-06-08 — Deslop Filter Regions and Merge Unified Modules

### Removed
- **mutation_surveyor**: Removed `--autopass-json` CLI argument and autopass_map JSON file writing from `filter_regions.py` — inter-step communication no longer uses file serialization
- **mutation_surveyor**: Removed `autopass_map` return value from `run_filter_ms()` — function now returns `None`
- **sequencher**: Removed `--autopass-json` CLI argument and autopass_map JSON file reading from `filter_regions.py` — Sequencher filter now always reads auto-pass regions from sample_flags
- **sequencher**: Removed `autopass_map` parameter from `run_filter_sequencher()` — simplified to single path using `autopass_regions_from_flags()`
- **mutation_surveyor**: Removed `import sys` from `filter_regions.py` (unused)
- **mutation_surveyor**: Removed duplicated `ALL_REGIONS` constant from `filter_regions.py` — now imported from `src.core.sample`
- **shell**: Removed `AUTOPASS_JSON` variable and `--autopass-json` args from `mutation_surveyor.sh`

### Changed
- **core**: Moved `load_unwrapped_sample_json()`, `write_unwrapped_sample_json()`, and `discover_sample_jsons()` from `filter_regions.py` to `src/core/sample.py` — these are generic Sample I/O utilities, not filter-specific
- **sequencher**: Moved `complementary_regions()` from MS `filter_regions.py` to Sequencher `filter_regions.py` — function is only used by Sequencher filtering
- **mutation_surveyor**: `merge_unified.py` now imports I/O helpers from `src.core.sample` instead of `filter_regions`
- **mutation_surveyor**: `merge_unified.py` — moved `import tempfile` to module level; removed redundant `tool_sources` tracking and inner try/raise in `_process_sample`
- **sequencher**: `filter_regions.py` now imports `ALL_REGIONS` and I/O helpers from `src.core.sample` instead of MS `filter_regions`

### Added
- **core**: Added `discover_sample_jsons()` to `src/core/sample.py` — shared helper for `<dir>/<LID>/<LID>.json` directory discovery (eliminates 3x duplication)

## [Unreleased] — 2026-06-08 — Simplify merge_unified to Tool-Agnostic Two-Dir Merge

### Changed
- **mutation_surveyor**: `merge_unified.py` now takes `--dir-a` / `--dir-b` instead of `--ms-dir` / `--sequencher-dir` — module no longer cares which tool produced each JSON
- **mutation_surveyor**: `_process_sample()` now takes a `sample_jsons: Dict[str, Path]` instead of separate `ms_path`/`seq_path` — single pass-through path for all single-dir samples
- **mutation_surveyor**: Added `--output-dir` CLI argument (defaults to `--dir-a`) and `output_dir` parameter to `run_merge_unified()`
- **shell**: Updated `mutation_surveyor.sh` Step 7 to use `--dir-a`/`--dir-b` instead of `--ms-dir`/`--sequencher-dir`

### Removed
- **mutation_surveyor**: Removed tool-specific branching ("Sequencher-only copy to MS dir" vs "MS-only keep as-is") — single pass-through path writes to `output_dir` regardless of source

## [Unreleased] — 2026-06-08 — Inline Auto-Pass Region Filter into ETL Pipeline

### Changed
- **mutation_surveyor**: Auto-pass region filtering now happens inline in `pipeline.process_sample()` and `etl.build_samples()` after ETL — no longer a separate CLI step
- **mutation_surveyor**: `process_sample()` returns `status: "skipped"` for samples with no auto-pass regions (previously these would have been written and then deleted by a separate step)
- **mutation_surveyor**: `build_samples()` now applies `autopass_regions_from_flags()` + `filter_sample_by_regions()` per LID before adding to batch — samples without auto-pass regions are excluded before `Batch.write()`
- **mutation_surveyor**: `etl.build_samples()` per-LID loop now continues on ETL failure (previously added partial sample to results)
- **shell**: Removed Step 5 (MS filter CLI) from `mutation_surveyor.sh` — filtering is now inline; renumbered Step 6→5, Step 7→6

### Removed
- **mutation_surveyor**: Removed `run_filter_ms()` and `main()` CLI from `filter_regions.py` — filtering is now inline, module only provides `autopass_regions_from_flags()`
- **mutation_surveyor**: Removed `--ms-dir`, `--batch-id`, `--ref-path` CLI arguments from `filter_regions.py`

## [Unreleased] — 2026-06-08 — Move Auto-Pass Filter from etl.py to pipeline.py

### Changed
- **mutation_surveyor**: `process_batch()` now runs ETL per LID + auto-pass filter + Batch.write() directly — no longer delegates to `etl.build_samples()` for the batch write path
- **mutation_surveyor**: `etl.build_samples()` is unchanged — remains a pure ETL entry point that writes all samples without filtering

### Removed
- **mutation_surveyor**: Removed auto-pass filtering from `etl.py` — filtering lives only in `pipeline.process_sample()` and `pipeline.process_batch()`

## [Unreleased] — 2026-06-09 — Align mutation_surveyor.sh output directories with §3.7

### Changed
- **mutation_surveyor**: Renamed `MUTATION_SURVEYOR_REVIEW_DIR` → `MUTATION_SURVEYOR_ANALYSIS_DIR` in `scripts/tools/mutation_surveyor.sh` to match §3.7 naming (reviewer-facing analysis directory)
- **mutation_surveyor**: Replaced hardcoded NAS output paths with `${RESULTS_DIR}`-relative paths — `MUTATION_SURVEYOR_RESULTS_DIR` now uses `${RESULTS_DIR}/tools/mutation_surveyor`, `MUTATION_SURVEYOR_ANALYSIS_DIR` uses `${RESULTS_DIR}/analysis/mutation_surveyor`, and `SEQUENCHER_RESULTS_DIR` uses `${RESULTS_DIR}/tools/sequencher` per §3.7
- **mutation_surveyor**: Added `init_batch_dirs()` function to `scripts/tools/mutation_surveyor.sh` — creates the full §3.7 output directory structure for each batch upfront (input/, autopass/, not_autopass_one/HV1/AB1/, not_autopass_one/HV2_HV3/AB1/, not_autopass_both/AB1/) and copies MS software source files into input/
- **mutation_surveyor**: Removed `/excels` subdirectory from `MUTATION_SURVEYOR_EXCELS` path — xlsx intermediates now go directly under `{BATCH}/` per §3.7 instead of `{BATCH}/excels/`
- **mutation_surveyor**: Updated sort_review_data output directory to use `MUTATION_SURVEYOR_ANALYSIS_DIR` — review artifacts (autopass/not_autopass directories) co-located with xlsx intermediates per §3.7
- **mutation_surveyor**: Updated shell script header comments to document §3.7 output directory structure
- **mutation_surveyor**: Changed directory name `HV2-3` → `HV2_HV3` in `sort_review_data.py` for consistency with §3.7 tree diagram
- **docs**: Fixed §3.8 table inconsistency — changed `not_autopass_one/HV2-3/` → `not_autopass_one/HV2_HV3/` to match §3.7 directory tree

## [Unreleased] — 2026-06-09 — Wire steps 7 (ETL) and 8 (AB1 sort) into mutation_surveyor.sh

### Added
- **mutation_surveyor**: Added Step 17 (ETL + auto-pass region filter) to `scripts/tools/mutation_surveyor.sh` — invokes `mutation_surveyor_dev/pipeline.py` with `--final-profiles` and `--output-dir`
- **mutation_surveyor**: Added Step 18 (AB1 sort by autopass/review) to `scripts/tools/mutation_surveyor.sh` — invokes `mutation_surveyor_dev/sort_review_data.py` with `--batch-id`, `--results-dir`, `--raw-dir`, `--excels-dir`
- **mutation_surveyor**: Added `main()` CLI entry point to `src/modules/mutation_surveyor_dev/pipeline.py` with `--final-profiles`, `--output-dir`, `--batch-id`, `--ref-path` arguments
- **mutation_surveyor**: Added `RAW_DIR` as 4th positional argument to `scripts/tools/mutation_surveyor.sh` for raw AB1 file directory

### Changed
- **mutation_surveyor**: Refactored `process_batch()` in `src/modules/mutation_surveyor_dev/pipeline.py` — removed broken references to deleted steps (`ms_path`, `lid_manifest_path`, `lid_path`, `lid_result.combined_flags`); now accepts `final_profiles` path directly
- **mutation_surveyor**: Updated module docstring in `pipeline.py` to reflect ETL-only scope (steps 01–06 handled by shell script)
- **mutation_surveyor**: Updated `scripts/tools/mutation_surveyor.sh` usage/help text to document all 4 positional arguments and 8 steps
- **mutation_surveyor_dev**: Fixed cross-module imports — changed `from src.modules.mutation_surveyor.{columns,flagging,models,utils,...}` to `from src.modules.mutation_surveyor_dev.{...}` in all _dev modules and `etl.py` to reflect module restructuring

## [Unreleased] — 2026-06-10 — Remove iupac_flags and legacy fallback from ETL

### Removed
- **etl**: Removed `iupac_flags` tracking (`php:{token}` flag strings) from `process()` — PHP detection is handled by `MSFlag` objects from `build_hv1_flags()` / `build_hv23_flags()` with `FlagLevel.VARIANT`
- **etl**: Removed `_parse_flags_from_row()` function (semicolon-parsing legacy fallback) — structured `MSFlag` routing via `combined_flags` is now the sole path; when `combined_flags` is absent for a LID, flags default to empty lists
- **etl**: Removed `COMBINED_FLAGS` column import — no longer used in ETL

## [Unreleased] — 2026-06-10 — Deslop txt_to_excel.py
### Removed
- **mutation_surveyor**: Removed `--sheet-name` and `--encoding` CLI arguments from `txt_to_excel.py` — not used by any caller, library API still accepts them as keyword args

### Changed
- **mutation_surveyor**: Inlined `_has_trace_header()` and `_derive_xlsx_path()` into their single call sites in `txt_to_excel.py`
- **mutation_surveyor**: Trimmed `__all__` in `txt_to_excel.py` to remove `DEFAULT_SHEET_NAME` and `main()` — not consumed externally

### Removed
- **mutation_surveyor**: Removed `_has_trace_header()` and `_derive_xlsx_path()` private helpers from `txt_to_excel.py` — single-call-site abstractions that added indirection without clarity
- **mutation_surveyor**: Removed `main()` and `DEFAULT_SHEET_NAME` from `__all__` in `txt_to_excel.py` — CLI entry points and internal constants are not library API

## [Unreleased] — 2026-06-10 — Sync docs/tools/mutation_surveyor.md with current source

### Changed
- **docs**: Updated `docs/tools/mutation_surveyor.md` §2 pipeline flow to reflect `01_02/`, `03_04/`, `05_06/` subdirectory structure instead of deleted top-level modules (`quality_control.py`, `transform_hv1.py`, `transform_hv23.py`, `merge_profiles.py`, `columns.py`)
- **docs**: Replaced §3 pipeline stages with subdirectory module map (Steps 01–06) and top-level module map (Steps 07–09)
- **docs**: Fixed §4 flagging architecture — replaced `MSFlagger` with `MSFlag`/`FlagLevel`, updated review classification to match current `flagging.py`
- **docs**: Fixed §5 data models — removed non-existent result models (`CompareQCResult`, `HV23MergeResult`, etc.), replaced with actual `MSFlag`, `ReviewCategory`, `LIDClassification`, `SortResult`
- **docs**: Fixed §6 CLI usage — updated args to match current `src/tools/mutation_surveyor.py` and `pipeline.py`
- **docs**: Fixed §7 intermediate file naming — replaced old stage names with subdirectory script paths
- **docs**: Fixed §9.1 txt_to_excel CLI — removed `--sheet-name` and `--encoding` flags
- **docs**: Fixed §10 review flag categories — updated source module references to use subdirectory names
- **docs**: Removed `Flag` and `FlagLevel` from §1.2 core models table (now MS-module-only)

### Changed
- **core**: Fix 127 ruff lint errors across src/core/ modules — replace `open()` with `Path.open()`, combine nested `if` statements, use `list.extend` instead of loop+append, add type annotations, replace blind `Exception` catches with specific exceptions, annotate mutable class attributes with `ClassVar`, fix naming violations with `noqa` comments for domain-specific names (rCRS, mtDNA), replace `os.path.splitext` with `Path.suffix`, replace `Union[X, Y]` with `X | Y`, fix boolean positional arguments by making them keyword-only, rename unused variables with `_` prefix, move validation raises out of try blocks, and add `noqa` comments for acceptable complexity/magic-value violations
- **core**: Fix remaining noqa comments — move deferred imports to top level (no circular deps), extract mtDNA position magic values to named constants (`POS_249`, `POS_459`, `POS_513`, etc.), promote `_RANGE_FLAGS`/`_MANUAL_TOOLS`/`TOOL_PRIORITY` to module-level constants, rename `mtDNAMerger` → `MtDnaMerger`, rename `_merge_profiles` → `merge_profiles` (public), rename `_ref_path` → `ref_path` (public), rename `_load_rCRS_reference` → `_load_rcrs_reference`, rename `rCRS`/`rCRS_base` local vars → `rcrs`/`rcrs_base`, add return type to `arg_parser()`
- **mtdna_merger**: Add explicit `hv1`, `hv2`, `hv3`, `no_snps`, `no_ins`, `no_dels`, `no_identicals`, `no_unread` parameters (set to `None`) to all three `Sample()` constructions in `merge()`, `merge_samples()`, and `merge_by_regions()`

### Fixed
- **core/flagging**: Remove masking fallback slop — 4 `try/except (ValueError, AttributeError): pass` blocks around `pos_base()` calls were silently swallowing errors; replaced with direct calls since `pos_base()` handles both int and str positions
- **core/flagging**: Eliminate `pass`-with-`else` antipattern — 4 inverted condition branches (`if cond: pass; else: action`) converted to `if not cond: action` for readability
- **core/flagging**: Remove stale `PLR0912` and `PLR0915` from noqa directives after branch reduction
- **core/uploader**: Fix 4 basedpyright errors — add `# type: ignore[reportArgumentType]` for oauth2client stub gaps, use ternary for df assignment to eliminate unbound variable, pass `replace=` as keyword argument
- **tests/mtdna_merger**: Fix broken import — `mtDNAMerger` renamed to `MtDnaMerger`, add missing `main` import

### Changed
- **core/flagging**: Combine nested `if` conditions per SIM102 — deletion and AC-repeat-SNP branches flattened

### Fixed
- **sequencher**: Fix `TypeError: process_batch() got an unexpected keyword argument 'ref_path'` — CLI now passes `ref_path`, `batch_id`, and `json_dir` via `BatchOptions` instead of as separate kwargs to `process_batch()`

### Fixed
- **sequencher**: Fix doubled `JSON/JSON` and doubled batch ID in output paths — OUT_DIR now uses `$RESULTS_DIR/tools/sequencher` (no batch ID), and JSON_DIR no longer appends `/JSON` since `write_region_jsons()` already adds `JSON/` internally
- **mutation_surveyor**: Replace hardcoded Sequencher output path with `$RESULTS_DIR/tools/sequencher`; fix `SEQ_JSON_DIR` to not double-append `/JSON`

### Fixed
- **sequencher**: Separate region JSON output (stays in `$SEQUENCHER_DIR/$BATCH_ID`) from batch JSON output (goes to `$RESULTS_DIR/tools/sequencher/$BATCH_ID`); fix doubled `JSON/JSON` path and doubled batch ID in output
- **mutation_surveyor**: Pass `--json-dir` to sequencher.sh pointing to analysis directory; fix unify step `SEQ_JSON_DIR` to find region JSONs in the analysis directory

### Changed
- **tracy**: Remove duplicate `_TransformParams` NamedTuple from `etl.py` — was identical to `_CreateParams`; all callers now use `_CreateParams` directly, eliminating redundant `_CreateParams(ref_seq=..., ...)` construction inside `_apply_all_transformations`

- **tracy**: Remove redundant CLI args (trim, pratio, maxindel, quality_threshold, min_peak_value, heteroplasmy_threshold) from src/tools/tracy.py — these params are already configurable via TracySettings (env vars MTDNA_TRACY_*) and should not be duplicated as CLI args
- **tracy**: Remove Tracy param fields from BatchOptions — read directly from get_settings().tracy at usage sites instead of duplicating defaults via field(default_factory=...)
- **tracy**: Replace `_tracy = get_settings().tracy` pattern with direct `get_settings().tracy` access across etl.py, preprocessing.py, utils.py, and pipeline.py
- **tracy**: Split src/modules/tracy/etl.py (941→477 lines) by extracting variant transformation logic (position-specific transforms, polyC transforms, strand transforms, _CreateParams, _PositionParams, _extract_peak_data, _create_variant) into new src/modules/tracy/transforms.py module

### Fixed
- **pipeline**: Fix tracy CLI invocation — use `--input-dir`/`--output-dir` with hyphens (standard convention) instead of underscores in tracy.py argparse definitions; fix `--input-dir` to point to batch-specific `${RAW_DIR}` instead of parent `${DATA_DIR}/raw`; fix `--output-dir` to use `${AUTOMATE_TRACY_DIR}` instead of `${AUTOMATE_DIR}`

- **tracy**: Add per-sample progress logging in `process_batch()` — log `[{idx}/{total}] Decomposing sample {id}` and `[{idx}/{total}] Processing ETL for sample {id}` so the pipeline no longer appears silent between batch start and completion
