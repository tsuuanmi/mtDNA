# Mutation Surveyor Pipeline

> **Module**: `src/tools/mutation_surveyor/`
> **Entry point**: `src/tools/mutation_surveyor/pipeline.py`

## Purpose

The Mutation Surveyor pipeline takes the **custom report exported by the
Mutation Surveyor software** (an `.xlsx`, or a `.txt` that is first converted
to `.xlsx`) and carries it through a multi-stage analysis to produce
standardized `Sample` objects, per-sample and batch JSON, and an optional
concordance check against a truthset workbook.

It does **not** process AB1 trace files directly. AB1 files are analyzed by
the external Mutation Surveyor software, which exports the custom report that
this pipeline consumes.

The pipeline is one of the tool entry points in the **GSTmtDNAv1** workflow.
Every `Sample` it produces sets `source_tool = Tool.MUTATION_SURVEYOR`, so the
shared downstream tasks — batch aggregation, pairwise comparison
(`SampleComparator`), multi-tool merging (`mtdna_merger`), batch unification
(`unify`), and FASTA generation (`generate_fasta`) — can consume its output
through the shared `Sample` / `Batch` / `Variant` models in `src/core`.

The main outputs are:

- Per-sample JSON files (`{LID}.json`) carrying variants, HV sequences, flags,
  intervals, and statistics.
- A batch JSON (`statistic_fullbatch.json`) aggregating all samples.
- Intermediate xlsx workbooks at each pipeline stage (written to `excels_dir`).
- An optional concordance report (`{BATCH}_final_profiles_vs_truth.xlsx`),
  produced by the standalone Step 06 script when a truthset workbook is
  supplied.

## Two execution layers

The pipeline has two distinct layers, and it helps to keep them separate:

1. **Standalone shell-script steps (01–06)** live in the subdirectories
   `01_02/`, `03_04/`, `05_06/`. They are invoked as individual scripts by
   `scripts/tools/mutation_surveyor.sh` and are *not* imported by the
   top-level Python pipeline. Each subdirectory has its own README describing
   its internal logic. They turn the exported custom report into merged,
   QC-flagged per-LID profiles (`{BATCH}_final_profiles.xlsx`).
2. **The Python pipeline (steps 08–09)** lives in the top-level
   `src/tools/mutation_surveyor/` modules and is driven by `pipeline.py`.
   Step 08 runs the ETL that converts `final_profiles.xlsx` into `Sample`
   objects plus JSON, and applies the auto-pass region filter. Step 09 sorts
   AB1 traces into review packages by classification.

`txt_to_excel.py` (Step 1) is a small standalone converter used by the shell
script when the input report is a `.txt` rather than an `.xlsx`.

## Stage flow

```text
MS custom report (.xlsx or .txt)
   │  .txt → convert (Step 1: txt_to_excel.py)
   ▼
┌────────────────────────────────────────────────────────────────┐
│  Steps 01–02: trace processing + control QC                    │
│  01_02/ms_trace_process.py → 01_02/ms_control_QC.py            │
│  → {BATCH}_MS_trace_processed.xlsx (+ {BATCH}_control_QC.xlsx)  │
└───────────────────────────────┬────────────────────────────────┘
                                │
              ┌─────────────────┴─────────────────┐
              ▼                                   ▼
  ┌────────────────────────┐          ┌────────────────────────┐
  │ Step 03: HV2-3 merge   │          │ Step 04: HV1 merge     │
  │ 03_04/ms_hv23_merge.py │          │ 03_04/ms_hv1_merge.py  │
  │ → {BATCH}_HV23_merged  │          │ → {BATCH}_HV1_merged   │
  └────────────┬───────────┘          └────────────┬───────────┘
               └──────────────┬────────────────────┘
                              ▼
  ┌────────────────────────────────────────────────────────────────┐
  │  Step 05: final profiles merge                                  │
  │  05_06/ms_final_profiles_merge.py                               │
  │  → {BATCH}_final_profiles.xlsx                                  │
  └───────────────────────────────┬────────────────────────────────┘
                                  │ (optional)
                                  ▼
  ┌────────────────────────────────────────────────────────────────┐
  │  Step 06: final profiles vs truth                               │
  │  05_06/ms_final_profiles_vs_truth.py                            │
  │  → {BATCH}_final_profiles_vs_truth.xlsx                         │
  └───────────────────────────────┬────────────────────────────────┘
                                  ▼
  ┌────────────────────────────────────────────────────────────────┐
  │  Step 08: ETL + auto-pass region filter                         │
  │  pipeline.process_batch() → etl.process()                       │
  │  → per-sample JSON + statistic_fullbatch.json                   │
  └───────────────────────────────┬────────────────────────────────┘
                                  ▼
  ┌────────────────────────────────────────────────────────────────┐
  │  Step 09: AB1 sort by classification                            │
  │  sort_review_data.sort_review_data()                            │
  │  → categorized AB1 folders for manual review                    │
  └────────────────────────────────────────────────────────────────┘
```

The stages, in prose:

1. **TXT→XLSX (Step 1, `txt_to_excel.py`)** — If the input report is a `.txt`,
   it is converted to `.xlsx`. Amino-acid notation (`,p.XXX`) is stripped so
   only nucleotide-level information and quality scores remain, and the output
   filename is normalized to the `{BATCH}_{anything}_df_custom_report.xlsx`
   convention regardless of the input suffix. If an `.xlsx` report already
   exists it is used directly.

2. **Trace processing + control QC (Steps 01–02, `01_02/`)** —
   `ms_trace_process.py` parses the multi-block MS custom report, cleans
   non-real trace rows, produces the `Variant Summary` / `Variant Converted` /
   `Variant Converted 2` sheets, and emits the QC-ready polyC indicators
   `HV1_PolyC_flag` and `HV2-3_PolyC_flag`. Its output is
   `{BATCH}_MS_trace_processed.xlsx`. The optional
   `ms_control_QC.py` (Step 02b) validates PC and NTC controls at the profile
   level (PC observed profile must exactly match the expected profile; NTC
   must be empty; each control must have exactly one trace per expected
   primer `HV1F/HV1R/HV2F/HV3R`) and writes `{BATCH}_control_QC.xlsx`,
   exiting with code 1 if any control fails.

3. **HV2-3 / HV1 merge (Steps 03–04, `03_04/`)** — These two are parallel
   branches that both consume the trace-processed workbook independently.
   `ms_hv23_merge.py` builds a **base** HV2-3 view (`Range`/`Variants`/
   `Mismatch by HV2-3 Regions`, duplicate variants kept once, out-of-range
   variants removed) and an **adjusted** view (HV2F trimmed to 73–303; HV3R
   trimmed to 303–340 and 438–576), then computes the `HV2-3 Flag list`,
   `HV2-3 Pair Status`, `HV2-3 Review Class`, `HV2-3 Truth Status` (blank when
   truth columns are absent), `HV2-3 Severity`, and the `HV3R PHP 341-437` QC.
   `ms_hv1_merge.py` keeps HV1F/HV1R, trims each trace to HV1 `16024-16365`,
   selects one primary HV1F + HV1R per LID, merges them into one HV1 profile,
   preserves trace-boundary ranges (no adjacent-range merging), detects strict
   mismatch in overlap, flags same-position multi-calls, scans source traces
   for het calls, and splits flags into a full `HV1 Flag list` and a
   `HV1 Review-triggering Flag list`. It does not require truthset columns.
   Outputs: `{BATCH}_HV23_merged.xlsx` and `{BATCH}_HV1_merged.xlsx`.

4. **Final profiles merge (Step 05, `05_06/ms_final_profiles_merge.py`)** —
   Merges the `LID` sheets of the HV23 and HV1 merged workbooks into the
   combined per-LID final MS mtDNA profile, `{BATCH}_final_profiles.xlsx`.
   No truthset is used here; concordance is performed separately in Step 06.

5. **Final profiles vs truth (Step 06, `05_06/ms_final_profiles_vs_truth.py`,
   optional)** — Compares merged profiles against a truthset by normalized LID.
   It performs full-profile, HV2-3 region, and HV1 region concordance
   (`Y` exact, `Y1` non-het IUPAC-compatible, `F` missing/extra/incompatible
   or het-vs-non-het; `Y2` is not emitted; length variants `309.xC`/`573.xC`
   are ignored in HV2-3 and full comparisons). It appends
   `Full Profile / HV2-3 / HV1 Concordance`, `Discordance reason`,
   `Discordance Category`, and truth `Range`/`Variants`/`Flags` columns, and
   exposes the `LID_no_PC_NTC`, `Summary`, and `Discordance_detail` sheets
   consumed by `merge_final_profiles_vs_truth.py`. Output:
   `{BATCH}_final_profiles_vs_truth.xlsx`.

6. **ETL + auto-pass filter (Step 08, `pipeline.process_batch` →
   `etl.process`)** — The terminal pipeline step. Reads the `LID` sheet from
   `final_profiles.xlsx`, turns each row into a `Sample`, keeps only the
   auto-pass regions, drops samples with no auto-pass regions or no variants,
   writes per-region intermediate JSONs (`HV1_{LID}.json`,
   `HV2-3_{LID}.json`) for the auto-pass regions, and writes the batch JSON
   via `Batch.write()`. This is the step that produces the `Sample` objects
   consumed by GSTmtDNAv1 downstream.

7. **AB1 sort by classification (Step 09, `sort_review_data.py`)** —
   Classifies each LID by its HV1 and HV2-3 Review Class into
   no-triggering-flag vs triggering-flag categories and copies the
   corresponding AB1 traces into categorized folders for the manual-review
   workflow (see [AB1 sort](#ab1-sort-step-09)).

A generic Excel merge utility, `01_02/merge_files.py`, concatenates the first
worksheet of multiple input files row-wise and adds `Batch` and `Source_File`
helper columns. It is used for ad-hoc workbook merging, not as a pipeline
stage.

## Module structure

```
src/tools/mutation_surveyor/
├── etl.py                              # ETL: final_profiles xlsx → per-LID Sample
├── flagging.py                         # MSFlag, FlagLevel, region status flags, flag routing
├── pipeline.py                         # process_batch() + main() CLI (Steps 08–09)
├── sort_review_data.py                 # AB1 sort by no/triggering-flag classification
├── txt_to_excel.py                     # TXT custom-report → XLSX conversion (Step 1)
├── merge_final_profiles_vs_truth.py    # Merge per-batch vs_truth workbooks
├── utils.py                            # Shared constants + helpers (IUPAC, parsing, controls)
├── 01_02/                              # Steps 01–02: trace processing + control QC
│   ├── ms_trace_process.py             # parse/clean MS report, emit polyC flags
│   ├── ms_control_QC.py                # PC/NTC profile-level control QC
│   └── merge_files.py                  # generic xlsx row-wise merge utility
├── 03_04/                              # Steps 03–04: HV2-3 + HV1 merges
│   ├── ms_hv23_merge.py                # base + adjusted HV2-3 views, HV2-3 flags
│   └── ms_hv1_merge.py                 # HV1 trace merge, HV1 flags split
└── 05_06/                              # Steps 05–06: final profiles + vs_truth
    ├── ms_final_profiles_merge.py      # combine HV23 + HV1 LID sheets
    └── ms_final_profiles_vs_truth.py   # concordance vs truthset workbook
```

The subdirectory scripts are standalone and are not imported by the top-level
pipeline; each subdirectory has its own `README_*.md` for internal logic,
column names, and algorithms.

## ETL (`etl.process`)

`etl.process()` is the heart of Step 08. It reads one LID's row from
`final_profiles.xlsx` and builds a `Sample`. The stages, in order:

1. **Load rCRS** into a `{position: base}` map (1-based) from `ref_path`
   (default `ref/rCRS.fasta`).
2. **Read the LID sheet**, strip column names, and locate the row matching
   `sample_id`. If the LID is missing, a minimal empty `Sample` is returned
   and a warning is logged.
3. **Parse variant tokens** from `Sample profile - Variants` via
   `utils.parse_variant_tokens()`, producing `(Position, allele, token)`
   triples, and map each to a `Variant` (SNP, insertion, deletion, PHP/IUPAC,
   het SNP → IUPAC, het deletion → `seq="DEL"`). Per-position quality
   (`Sample profile - Quality`) and source files (`Sample profile - Files`)
   are attached.
4. **Build intervals** from `Sample profile - Range`, falling back to the
   canonical full-region bounds from `settings.regions.REGIONS` when the text
   is empty or unparseable.
5. **Route MSFlags** into `sample_flags` and `variant_flags` — *but only when
   `combined_flags` is passed and `sample_id` is present in it*. See the
   [combined_flags guard](#the-combined_flags-guard) below.
6. **Flag variants** via core `flag_variants()` from `src/core/flagging.py`
   (polyC, het, indel, etc.), merged into `variant_flags` with
   order-preserving dedup. These are always computed from the parsed
   `Variant` list, regardless of `combined_flags`.
7. **Build region status flags** from `HV1 Review Class` and
   `HV2-3 Review Class` via `flagging.build_region_status_flags()`, appended
   to `sample_flags`.
8. **Run sample-level core flagging** via
   `SampleFlagger(variant_list, intervals).analyze()` (only when
   `variant_list` is non-empty); core sample reasons are appended to
   `sample_flags` if not already present.
9. **Deduplicate** `sample_flags` against per-variant flags via
   `core.deduplicate_sample_flags()`.
10. **Construct the `Sample`** with `source_tool=Tool.MUTATION_SURVEYOR`,
    populating `information["well"]` from the `Well` column when present.

The HV/statistic fields (`hv1`, `hv2`, `hv3`, `no_*`) are left `None` here —
they are computed by `Batch` at aggregation time.

### The `combined_flags` guard

Flag routing (`sample_flags_from()` / `variant_flags_from()`) only runs when
`combined_flags` is passed **and** `sample_id` is present in it.
`pipeline.process_batch()` calls `etl.process()` **without** `combined_flags`,
so in the standard CLI flow `sample_flags` start empty and only the core
region-status flags plus `SampleFlagger.analyze()` flags are populated.
Per-variant flags from `flag_variants()` are always computed from the parsed
`Variant` list regardless of `combined_flags`.

## Flagging architecture

The MS pipeline uses a two-layer flag system:

- **MS flags** (`src/tools/mutation_surveyor/flagging.py`) are per-LID and
  per-region. They carry a `FlagLevel` so the ETL can route them into
  `sample_flags` vs `variant_flags`. Representative names: `Autopass HV1`,
  `Review HV2 and HV3`, `No full HV1 region`, `HV1 polyC variant`, `PHP`.
- **Core flags** (`src/core/flagging.py`) are per-variant and are computed by
  `flag_variants()` during ETL. Representative names: `Insertion C at 309.1`,
  `Heteroplasmy at 16192`, `Deletion at 523`.

### `MSFlag` and `FlagLevel`

`FlagLevel` is an enum with two values: `SAMPLE` (`"sample"`) for flags about
the sample as a whole (coverage gaps, pair status, autopass/review status),
and `VARIANT` (`"variant"`) for flags about a specific variant position
(polyC, PHP, indel; het is detected via IUPAC codes in `Variant.seq`).

`MSFlag` is a Pydantic model with fields `name` (str, required), `level`
(`FlagLevel`, required), and optional `pos` (`Position | None`), `ref`
(`str | None`), `seq` (`str | None`). It exposes:

- `is_review_trigger` (property): `self.name not in _NON_REVIEW_TRIGGER_NAMES`.
- `to_text()`: returns `self.name` (the flat string used in Excel columns).

The routing helpers:

- `sample_flags_from(flags) -> list[str]` — sample-level flag names.
- `variant_flags_from(flags) -> dict[str, list[str]]` — variant-level flags
  keyed `"pos|ref|alt"`; flags without a position are omitted.
- `build_region_status_flags(hv1_review_class, hv23_review_class) -> list[MSFlag]`
  — one `SAMPLE`-level `MSFlag` per region that has a non-empty review class.
  A class containing `"Auto"` and not `"Review"` maps to the `AUTOPASS_*`
  flag; any other non-empty/non-whitespace string maps to `REVIEW_*`. Empty
  or whitespace strings produce no flag for that region.

### Region status flags

The four sample-level region status flag constants, all `FlagLevel.SAMPLE`:

| Constant | Value | Meaning |
|---|---|---|
| `AUTOPASS_HV1` | `Autopass HV1` | HV1 passed auto-review (no review-triggering flags) |
| `AUTOPASS_HV2_AND_HV3` | `Autopass HV2 and HV3` | HV2-3 passed auto-review |
| `REVIEW_HV1` | `Review HV1` | HV1 has review-triggering flags |
| `REVIEW_HV2_AND_HV3` | `Review HV2 and HV3` | HV2-3 has review-triggering flags |

### Review classification and `_NON_REVIEW_TRIGGER_NAMES`

The MS merge steps (03_04/) assign a review class per region from the flag
lists: `Auto - no flag` when there are no review-triggering flags, and
`Review` when there is at least one. `MSFlag.is_review_trigger` decides
whether a flag triggers review. Three flag names are explicitly excluded from
triggering review:

```python
_NON_REVIEW_TRIGGER_NAMES = frozenset({
    "HV1 polyC variant",
    AUTOPASS_HV1,          # "Autopass HV1"
    AUTOPASS_HV2_AND_HV3,  # "Autopass HV2 and HV3"
})
```

These are the informational HV1 polyC variant flag and the two autopass
status flags — they record information or a passing outcome, not a problem.

### Representative flag names by source step

The flag names below are emitted by the subdirectory merge scripts into the
`HV1 Flag list` / `HV2-3 Flag list` columns of the merged workbooks. ETL
consumes the resulting `HV1 Review Class` / `HV2-3 Review Class` columns, not
the raw flag lists. This list is representative; the subdirectory READMEs have
the authoritative enumeration.

| Source step | Representative flag names |
|---|---|
| 03_04 HV1 merge | `No full HV1 region`; `Missing HV1F` / `Missing HV1R` / `Missing HV1F and HV1R`; `HV1 trace mismatch`; `HV1 Trace Range Adjusted` / `HV1 Trace Variants Adjusted`; `HV1 same-position multi-call: {pos_text}`; `HV1 source het` / `HV1 source polyC het`; `HV1 indel outside polyC`; `HV1 polyC PHP` / `HV1 polyC deletion` / `HV1 polyC insertion`; `HV1 polyC variant` (non-triggering); `PHP` |
| 03_04 HV2-3 merge | `No full region`; `Missing HV2F` / `Missing HV3R` / `Missing HV2F and HV3R`; `Pair status not OK`; `Trace mismatch`; `HV2-3 same-position multi-call: {text}`; `non-C in HV2 polyC` / `non-C in HV3 polyC`; `HV3R PHP 341-437`; `Inconsistent HV2-3 truth`; `PHP` |
| 01_02 trace processing | `HV1_PolyC_flag`, `HV2-3_PolyC_flag` (C-stretch indicators `BAD`/`SUSPECT`/`LOW`, carried forward to the merges and feeding the `HV1 polyC *` / `non-C in HV* polyC` flags above) |

### Variant flagging through core

As of this version, variant-level flags for het, insertion, and deletion are
unified through `src/core/flagging.py`. The MS pipeline normalizes `_het`
tokens to IUPAC codes in `Variant.seq` (e.g. `309C_het` → `seq="Y"`), so
core's `flag_variants()` emits `Heteroplasmy at {pos}`. Het deletions
(`seq="DEL"`) and het insertions (lowercase `seq`) also produce
`Heteroplasmy at {pos}`. Insertions and deletions are likewise handled by
core, with exception rules for `315.1`, polyC C insertions, the `249`
deletion, and the `AC` `523-524` pair. Deferred decisions — whether PHP flag
routing should merge into het, be added to core, or stay MS-only, and whether
PolyC informational flags should stay MS-only or move to core — are pending
lab/clinical input.

## Auto-pass region filter (Step 08)

After ETL, `process_batch()` keeps only the auto-pass regions of each sample.
The mapping from sample-level autopass flags to regions is:

```python
AUTOPASS_FLAG_TO_REGIONS = {
    "Autopass HV1": ["HV1"],
    "Autopass HV2 and HV3": ["HV2", "HV3"],
}
```

`autopass_regions_from_flags(sample_flags)` looks up each flag in this map and
unions the mapped regions, returning a `set[str]` (e.g. `{"HV1"}` or
`{"HV2", "HV3"}`). Samples with no auto-pass regions, or no remaining
variants after filtering, are dropped. For each kept sample, per-region
intermediate JSONs are written under `{region_json_dir}/{HV1,HV2-3}/`, and
`validate_region_group_intervals()` logs any range-QC warnings. Finally
`Batch(list(results.values()), ref_path).write(output_dir, batch_id)` writes
`statistic_fullbatch.json` plus the per-sample JSON.

## AB1 sort (Step 09)

`sort_review_data.py` classifies each LID by its HV1 and HV2-3 Review Class
and copies AB1 trace files into categorized folders for manual review.
Classification uses the `HV1 Review Class` and `HV2-3 Review Class` columns: a
class starting with `"Auto"` (and not containing `"Review"`) means no
triggering flag for that region; anything else (including an empty string,
which is treated as Review with a warning) means a triggering flag.

A LID is `NO_TRIGGERING_FLAG` only when both `hv1_no_triggering` and
`hv23_no_triggering` are true; otherwise it is `TRIGGERING_FLAG`. Empty and
control LIDs (`is_control_lid`) are skipped.

The sort then:

1. Reads the `LID` sheet of `{batch_id}_final_profiles.xlsx` and calls
   `classify_lids()`.
2. Reads the `Traces` sheet from **both** `{batch_id}_HV23_merged.xlsx` and
   `{batch_id}_HV1_merged.xlsx`, concatenates them, and builds a per-LID
   per-primer filename lookup (deduplicated). Raises `FileNotFoundError` if
   neither traces file exists.
3. Creates the `no_triggering_flag/AB1/{HV1F,HV1R,HV2F,HV3R}/` and
   `triggering_flag/AB1/{...}/` directory trees.
4. For each LID, copies each region's primer traces into the matching folder
   using `REGION_TO_PRIMERS` (`"HV1" → ["HV1F","HV1R"]`,
   `"HV2-3" → ["HV2F","HV3R"]`). Missing source files are counted in
   `files_missing`.
5. Writes `no_triggering_flag/no_triggering_flag_classification.xlsx`
   (Classification + Summary sheets).

The classification xlsx has a **Classification** sheet (`LID`, `pass HV1`,
`pass HV2+HV3`, `pass both`, one row per LID, sorted by LID) and a **Summary**
sheet with aggregate counts (`Total LIDs`, `Pass HV1`, `Pass HV2+HV3`,
`Pass both`, `No triggering flag (pass both)`, `Triggering flag (fail one or
both)`).

The result is returned as a `SortResult` (Pydantic model) with
`no_triggering_flag_lids`, `triggering_flag_lids` (both sorted), `files_copied`,
`files_missing`, and `output_dir`.

## Output layout

The pipeline writes to three locations, serving different audiences.

### `results/analysis/mutation_surveyor/{BATCH}/` — reviewer-facing artifacts

All files for one batch in one place: intermediates, inputs, the control QC
report, and the AB1 review packages.

```
results/analysis/mutation_surveyor/{BATCH}/
├── input/                                        ← source files from MS software
├── {BATCH}_MS_trace_processed.xlsx               ← trace processing output
├── {BATCH}_control_QC.xlsx                       ← PC/NTC control QC (optional, Step 02)
├── {BATCH}_HV1_merged.xlsx                       ← HV1 merge output
├── {BATCH}_HV23_merged.xlsx                      ← HV2+HV3 merge output
├── {BATCH}_final_profiles.xlsx                   ← merged profiles
├── {BATCH}_final_profiles_vs_truth.xlsx          ← concordance validation
├── no_triggering_flag/
│   ├── no_triggering_flag_classification.xlsx
│   └── AB1/
│       ├── HV1F/  HV1R/  HV2F/  HV3R/            ← traces from autopassing samples
└── triggering_flag/
    └── AB1/
        ├── HV1F/  HV1R/  HV2F/  HV3R/            ← traces from review samples
```

### `analysis/json_temp/{BATCH_ID}/mutation_surveyor/` — region JSON interchange

Per-region intermediate JSONs are written to the shared `json_temp` directory
under `mutation_surveyor/`, with region-group subdirectories. Each region JSON
includes `region_sources` provenance (e.g. `{"HV1": "mutation_surveyor"}`).
These are temporary interchange files consumed by the unify step. The
`json_temp` directory also contains sibling `sequencher/` and `merge/`
subdirectories for Sequencher region JSONs and unified output respectively.

```
analysis/json_temp/{BATCH_ID}/
├── mutation_surveyor/
│   ├── HV1/      └── HV1_{LID}.json
│   └── HV2-3/    └── HV2-3_{LID}.json
├── sequencher/
│   ├── HV1/      └── HV1_{LID}.json
│   └── HV2-3/    └── HV2-3_{LID}.json
└── merge/
    ├── {LID}.json
    └── statistic_fullbatch.json
```

### `results/tools/mutation_surveyor/{BATCH}/` — downstream-facing JSON

The canonical JSON output consumed by downstream tasks.

```
results/tools/mutation_surveyor/{BATCH}/
├── {LID}/{LID}.json                ← per-sample JSON
└── statistic_fullbatch.json        ← batch-level aggregation
```

## CLI (`pipeline.py`)

```bash
python src/tools/mutation_surveyor/pipeline.py \
  -i results/analysis/mutation_surveyor/MS_300326_004/MS_300326_004_final_profiles.xlsx \
  -o results/tools/mutation_surveyor/MS_300326_004 \
  --json-dir analysis/json_temp/MS_300326_004/mutation_surveyor \
  --batch-id MS_300326_004 \
  --ref-path ref/rCRS.fasta
```

| Argument | Required | Default | Purpose |
|---|---|---|---|
| `-i` / `--input-dir` | yes | — | Path to `*_final_profiles.xlsx` (from Step 05); its parent is used as `excels_dir` for the AB1 sort |
| `-o` / `--output-dir` | yes | — | Output directory for batch JSON (per-sample + `statistic_fullbatch.json`) |
| `--json-dir` | no | same as `--output-dir` | Output directory for per-region intermediate JSON files |
| `--batch-id` | no | derived from input filename (strip `_final_profiles`) | Batch identifier |
| `--ref-path` | no | `ref/rCRS.fasta` | Path to rCRS reference FASTA |

Two paths are **not** CLI arguments and are worth calling out:

- `excels_dir` is derived as `input_dir.parent`, so the AB1 sort can locate
  the sibling `{BATCH}_HV23_merged.xlsx`, `{BATCH}_HV1_merged.xlsx`, and
  `{BATCH}_final_profiles.xlsx` workbooks.
- `raw_dir` (the AB1 source directory) is hard-coded inside `main()` as
  `/mnt/nas/bca/mtDNA/science/data/raw/{batch_id}`. The `data_dir` constant is
  `/mnt/nas/bca/mtDNA/science/data`.

Custom-report file search and TXT→XLSX conversion are performed by the shell
script (`scripts/tools/mutation_surveyor.sh`, Step 1), not by `pipeline.py`.
The Python CLI consumes the already-converted `*_final_profiles.xlsx`
directly.

### `txt_to_excel.py` CLI

```bash
python src/tools/mutation_surveyor/txt_to_excel.py \
  -i path/to/{BATCH}_df_custom_report.txt \
  -o path/to/{BATCH}_df_custom_report.xlsx
```

`-i`/`--input` (required) is the input `*_custom_report.txt`; `-o`/`--output`
(optional) is the output `.xlsx` path, derived by normalizing the stem to
`_df_custom_report.xlsx` when omitted. A valid report must contain a
trace-table header row with both `Trace #` and `Sample Name`.

### `merge_final_profiles_vs_truth.py` CLI

`merge_final_profiles_vs_truth.py` is a standalone cross-batch utility (not
part of the shell-script pipeline). It concatenates the per-batch
`LID_no_PC_NTC`, `Summary`, and `Discordance_detail` sheets across batches,
fills the `Batch` column from the batch directory name, and produces:

1. `all_batches_final_profiles_vs_truth.xlsx` — all rows with `Batch` filled.
2. `HV1_AutoNoFlag_Y1F.xlsx` — `LID_no_PC_NTC` rows where
   `HV1 Review Class == "Auto - no flag"` and
   `HV1 Concordance in ("Y1", "F")`.
3. `HV23_AutoNoFlag_Y1F.xlsx` — same filter for the HV2-3 columns.

Its argparse arguments are `--base-dir` (default
`/mnt/nas/bca/mtDNA/mtDNA_workflow_2_development/analysis/mutation_surveyor`),
`--batches` (defaulting to 7 hard-coded `MS_*` ids, `nargs="+"`), and
`--output-dir` (default `~/workspaces/mtdna_raw/temp`).

## Settings

There is **no `MutationSurveyorSettings`** class. The pipeline is CLI-driven;
all runtime inputs come from `pipeline.py` argparse arguments. The only
settings consumed at runtime are the shared `Settings` fields:

| Setting | Path in `Settings` | Default |
|---|---|---|
| `regions.REGIONS` | `GenomicRegionsSettings` | `{"HV1": [16024, 16365], "HV2": [73, 340], "HV3": [438, 576]}` |
| `regions.POLYC_REGIONS` | `GenomicRegionsSettings` | `{"HV2": [303, 315]}` |
| `regions.POLYC_QC_REGIONS` | `GenomicRegionsSettings` | `[(16183,16193),(16289,16290),(302,315),(455,463),(515,525)]` |
| `regions.SKIP_AFTER_POLYC` | `GenomicRegionsSettings` | `{"HV1": [16193, 16195]}` |
| `directories.data` | `DirectoryPathsSettings` | `Path("data")` (env prefix `MTDNA_DIRS_`) |

ETL reads `regions.REGIONS` via `_parse_ranges_from_text()` (canonical bounds
fallback). `_load_rcrs()` reads the rCRS reference from the `ref_path` CLI
argument. `pipeline.main()` reads `settings.directories.data` before
overriding `raw_dir` with the hard-coded NAS path above.

## Public API

The key entry points, with one-line signatures and purposes:

**`etl.py`**

```python
def process(
    sample_id: str,
    input_path: Path,
    ref_path: str = "ref/rCRS.fasta",
    batch_id: str | None = None,
    combined_flags: dict[str, list[MSFlag]] | None = None,
) -> Sample
```
ETL for one LID: parse variant tokens → build intervals → route MSFlags (only
when `combined_flags` is passed) → `flag_variants()` → region status flags →
`SampleFlagger.analyze()` + dedup → `Sample`.

**`pipeline.py`**

```python
def autopass_regions_from_flags(sample_flags: list[str]) -> set[str]
```
Union of regions mapped from `AUTOPASS_FLAG_TO_REGIONS` for the given flags.

```python
def process_batch(
    final_profiles: Path,
    output_dir: Path,
    ref_path: str = "ref/rCRS.fasta",
    batch_id: str | None = None,
    json_dir: Path | None = None,
) -> dict[str, Sample]
```
Run ETL + auto-pass filter over `final_profiles.xlsx`; write per-region and
batch JSON. Returns the kept `dict[LID, Sample]`.

**`flagging.py`**

```python
class FlagLevel(Enum): SAMPLE; VARIANT
class MSFlag(BaseModel): name: str; level: FlagLevel; pos; ref; seq  # + is_review_trigger, to_text()
def sample_flags_from(flags: list[MSFlag]) -> list[str]
def variant_flags_from(flags: list[MSFlag]) -> dict[str, list[str]]
def build_region_status_flags(hv1_review_class: str, hv23_review_class: str) -> list[MSFlag]
```
`MSFlag` carries the level used to route flags; `build_region_status_flags`
maps `HV1 Review Class` / `HV2-3 Review Class` to the `AUTOPASS_*` / `REVIEW_*`
sample flags. Region status constants: `AUTOPASS_HV1` (`"Autopass HV1"`),
`AUTOPASS_HV2_AND_HV3` (`"Autopass HV2 and HV3"`), `REVIEW_HV1`
(`"Review HV1"`), `REVIEW_HV2_AND_HV3` (`"Review HV2 and HV3"`).

**`sort_review_data.py`**

```python
class ReviewCategory(Enum): NO_TRIGGERING_FLAG; TRIGGERING_FLAG
class LIDClassification(BaseModel): lid: str; hv1_no_triggering: bool; hv23_no_triggering: bool; category: ReviewCategory
class SortResult(BaseModel): no_triggering_flag_lids; triggering_flag_lids; files_copied; files_missing; output_dir
def classify_lids(df: pd.DataFrame) -> dict[str, LIDClassification]
def sort_review_data(batch_id: str, results_dir, raw_dir, excels_dir, output_dir: str | Path | None = None) -> SortResult
```
Classify LIDs by HV1/HV2-3 Review Class and copy AB1 traces per
`REGION_TO_PRIMERS` (`HV1 → HV1F/HV1R`, `HV2-3 → HV2F/HV3R`) into
`no_triggering_flag/` and `triggering_flag/` folders; write the classification
xlsx.

**`utils.py`**

```python
def safe_str(x: object) -> str                 # stripped string; NaN/None -> ""
def parse_variant_tokens(text: object) -> list[tuple[Position, str, str]]
def is_control_lid(lid: object) -> bool        # NTC\d+ or PC\d+_.+ (not bare PC/NTC prefixes)
```
Shared helpers plus the `UNORDERED_MAP` IUPAC reverse map and
`EXPECTED_PRIMERS = ["HV1F", "HV1R", "HV2F", "HV3R"]`.

`parse_variant_tokens` matches each token with
`^(\d+(?:\.\d+)?)([A-Za-z]+(?:_[Hh][Ee][Tt])?(?:_[Dd][Ee][Ll])?)$`
(`_HET` normalized to `_het`), with a bare-position fallback; `Position` is
`int` for base positions and `str` for insertion positions (e.g. `"309.1"`).

## Python API example

```python
from src.tools.mutation_surveyor.pipeline import process_batch

results = process_batch(
    final_profiles="path/to/final_profiles.xlsx",
    output_dir="results/tools/mutation_surveyor/BATCH",
    ref_path="ref/rCRS.fasta",
    batch_id="BATCH_001",
    json_dir="analysis/json_temp/BATCH/mutation_surveyor",  # optional, defaults to output_dir
)
# results: dict[str, Sample] — keyed by LID (only auto-pass samples)
```

## Out-of-scope workflow context

The following steps are performed outside this Python module and are noted
for workflow context only:

- **AB1 processing** is done by the external Mutation Surveyor software, which
  loads AB1 files, runs automatic trace analysis, and exports the custom
  report xlsx that this pipeline consumes.
- **LID manifest** is resolved by the shell script from
  `/mnt/nas/bca/mtDNA/science/data/metadata/{BATCH_ID}.xlsx` (required columns
  `LID`, `MẺ CHẠY`; the batch id in the filename must match the `MẺ CHẠY`
  column values) and passed to `ms_trace_process.py --lid-manifest` for
  inventory QC.
- **Manual review**: LIDs classified as `Review` are separated for manual
  review of AB1 traces, corrected or confirmed by a reviewer, and merged back
  with auto-pass results.
- **GSTmtDNAv1 final review** happens after manual review results are merged.

## Cross-References

- [fasta.md](fasta.md) — reverse FASTA → Sample flow
- [unify.md](unify.md) — MS + Sequencher unification into unified profiles
- [../core/batch.md](../core/batch.md) — `Batch` aggregation and JSON output
- [../core/flagging.md](../core/flagging.md) — core per-variant `flag_variants()`
- [../ARCHITECTURE.md](../ARCHITECTURE.md) — pipeline overview