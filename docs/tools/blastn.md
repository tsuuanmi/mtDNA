# BLASTn Tool Module

> **Module**: `src/tools/blastn/`
> Aligns raw AB1 traces (or pre-assembled FASTA) against rCRS with BLASTN and
> calls variants into the standard `Sample` model for comparison, FASTA
> generation, and unification.

## Overview

BLASTn is the in-house variant caller for Sanger AB1 trace data. Each sample
is a set of AB1 traces (forward and reverse primer reads across the HV1, HV2,
and HV3 hypervariable regions of the mitochondrial control region). The module
converts those traces to FASTQ, trims low-quality ends, aligns them against
rCRS with NCBI BLASTN, parses the tabular alignment output, and runs a
multi-stage ETL that calls SNPs and indels, normalizes indels, handles polyC
length heteroplasmy, applies curated conversion and filter rules, and
assembles a standardized `Sample(source_tool=Tool.BLASTN)` for downstream
batch aggregation.

The AB1 pipeline has two halves sharing one ETL core: a **preprocessing** half
(`preprocessing.py`) that turns AB1 files into BLASTN tabular output via
external tools (UGENE, seqtk, BLASTN), and an **ETL** half (`etl.py`) that
turns that tabular output into a `Sample`. A **FASTA-input counterpart**
(`fasta_pipeline.py`) accepts pre-assembled consensus FASTA files instead of
AB1 traces, runs a single BLASTN alignment per file, and reuses the same ETL,
combine, region-JSON, and `Batch.write()` finalization path, so both inputs
produce an identical output layout.

```
AB1 traces ─▶ preprocessing (UGENE/seqtk/BLASTN) ─▶ BLASTN TSV ─┐
                                                                 ├─▶ ETL ─▶ Sample ─▶ combine + region JSONs + Batch.write()
FASTA files ─▶ BLASTN alignment ────────────────────▶ BLASTN TSV ┘
```

Variant calling here is alignment-based rather than peak/quality
threshold-based: `quality_threshold`, `min_peak_value`, and
`heteroplasmy_threshold` exist in the settings shape for parity with the other
tools but are not consumed by this module. Per-base Phred quality *is* used,
however, as a tie-break and quality gate during consensus and polyC handling.

## Module Structure

```
src/tools/blastn/
├── __init__.py          # empty — import directly from etl.py / pipeline.py
├── etl.py               # Pure ETL: BLASTN TSV -> Sample; public entry point process()
├── transform.py         # Variant transforms: indel right-norm, conversion rules, polyC removal, region validation, polyC warning flags
├── pipeline.py          # Batch orchestration for AB1 input: ProcessPoolExecutor, combine same-LID, region JSON + Batch.write
├── preprocessing.py     # AB1 -> BLASTN TSV subprocess pipeline (UGENE, seqtk, BLASTN)
├── fasta_pipeline.py    # FASTA-input batch pipeline (pre-assembled FASTA -> BLASTN TSV -> ETL)
└── utils.py             # BLASTN column constants, strand helpers, TSV parsing, quality extraction, subprocess runner
```

`utils.py` has no internal dependencies. `preprocessing.py` depends on
`utils`, `config`, and `loguru`. `transform.py` depends on `core/models`,
`core/variants`, and `config`, and inline-imports
`etl._load_no_filter_positions` (deferred) to avoid an `etl ↔ transform`
import cycle. `etl.py` depends on `utils`, `transform`, and the core
models/flagging/sample/variants modules. `pipeline.py` and
`fasta_pipeline.py` depend on `etl` plus the core
batch/flagging/models/paths/region/variants modules; `fasta_pipeline.py` also
reuses combine and region-JSON helpers from `pipeline.py` and uses
`Bio.SeqIO`. The module never imports the removed legacy root modules
(`blastn`, `analysis`, `postprocessing`, `write_report`); `etl.py` never
imports from `pipeline.py`, and `utils.py` never imports `Sample`/`Variant`.

## Preprocessing

`preprocessing.preprocess_sample(sample_id, list_ab1, output_dir, ref_path)`
runs a 7-step subprocess pipeline that turns one sample's AB1 traces into the
BLASTN tabular file consumed by ETL. Steps 1–4 run once per AB1 trace; steps
5–7 run once per sample after all traces are processed. Every external command
goes through `utils.run_command`, and the UGENE, BLASTN, and seqtk binaries
are validated up front with `shutil.which` (a missing binary raises
`FileNotFoundError`). All artifacts are written to
`{output_dir}/{sample_id}/`, and each BLASTN subprocess uses
`BlastnSettings.blastn_threads` threads (default 1).

1. **AB1 → FASTQ** — UGENE `convert-seq --format=fastq` produces `{name}.fastq`.
2. **Quality trim** — seqtk `trimfq -q {trim_threshold}` trims low-quality
   ends into `{name}.trim.fastq`.
3. **Extract Phred scores** — Biopython parses the trimmed FASTQ and writes
   `{name}.trim.quality` as `>id` headers followed by comma-separated Phred
   scores; this is the per-base quality used later by consensus and polyC.
4. **FASTQ → FASTA** — seqtk `seq -a` produces `{name}.trim.fasta`.
5. **Concatenate FASTA** — `_concatenate_fasta` globs `*.trim.fasta` into
   `{sample_id}.concatenate.fasta`, merging the forward/reverse primer reads
   for the sample.
6. **Concatenate quality** — `_concatenate_quality` globs `*.trim.quality`
   into `{sample_id}.concatenate.quality`, keeping per-trace quality aligned
   with the concatenated FASTA.
7. **Multi-sequence BLASTN alignment** — BLASTN aligns the concatenated FASTA
   against rCRS with `-outfmt 7` (tabular, 16 columns) into
   `{sample_id}.concatenate.blastn`. This is the primary ETL input. The former
   per-trace BLASTN text report was removed because no downstream stage
   consumed it.

`preprocess_sample` returns a `PreprocessResult` named tuple:

```python
class PreprocessResult(NamedTuple):
    blastn_tsv_path: Path
    quality_path: Path | None = None
    concatenated_fasta_path: Path | None = None
```

`blastn_tsv_path` is the concatenated BLASTN TSV; `quality_path` and
`concatenated_fasta_path` are `None` when those files were not produced. On
subprocess failure it raises `subprocess.CalledProcessError`.

### Parallelism

Sample-level parallelism is controlled by `BlastnSettings.max_workers`
(`None` → `os.cpu_count()`), driven by a `ProcessPoolExecutor` in
`pipeline.py`. Inner BLASTN threading is decoupled via
`BlastnSettings.blastn_threads` (default 1): because samples already run in
parallel, each BLASTN subprocess is single-threaded by default to avoid CPU
oversubscription (`max_workers × blastn_threads` threads). This matches the
Tracy sample-parallel pattern. Raise `MTDNA_BLASTN_BLASTN_THREADS` only when
running a single sample.

## ETL Processing

`etl.process(sample_id, input_path, ref_path, batch_id=None) -> Sample` is the
single ETL entry point. It reads a BLASTN format-7 TSV and returns a
`Sample(source_tool=Tool.BLASTN)`. The processing stages below run in order;
helper names are included so the source can be followed directly. Stages 7–10
remove variants (polyC, special positions, conversion, region validation),
while stage 12 annotates the remaining variants — removal and flagging are
distinct operations. The removal order mirrors the legacy `PostProcess.post_process()`
so the new `Sample` reaches final-output parity with the legacy pipeline.

1. **Load reference** — `ref_seq = str(SeqIO.read(ref_path, "fasta").seq)`.
2. **Parse BLASTN TSV** — `utils.parse_blastn_tsv(input_path)` returns one
   dict per alignment row, keyed by `utils.BLASTN_COLUMNS`.
3. **Normalize strand** — `utils.normalize_strand_all(rows)` augments each row
   with `forward_ref_seq`, `forward_sample_seq`, `forward_ref_start`, and
   `forward_ref_end`. Plus-strand rows copy the originals; minus-strand rows
   are reverse-complemented (`Bio.Seq.Seq.reverse_complement`) with start/end
   swapped, so all downstream logic works in forward-reference coordinates.
4. **Load quality (optional)** — looks for
   `{input}` with `.concatenate.blastn` replaced by `.concatenate.quality`;
   if it exists, `utils.extract_quality_scores()` parses it into
   `dict[str, list[int]]`, otherwise `{}`.
5. **Call variants per position** — `_call_variants(rows, quality_data, ref_seq)`:
   `_collect_variants_by_position` iterates rows and calls `_call_row_variants`
   for SNPs, insertions (`ref == "-"` at `.1` positions), and deletions
   (`seq == "-"`), with insertion-aware position arithmetic and per-variant
   polyC context via `_analyze_polyc_context`. Indels are then
   right-normalized with `transform._right_norm_indel`. Only variants inside
   `regions.REGIONS` (HV1/HV2/HV3) are kept (`get_hv_region_for_position`).
   `_compute_consensus_variants` then determines the consensus base per
   position (IUPAC codes for ties, with `-` excluded from ties), applies
   strand-aware polyC `flagged_for_filtering` consensus filtering, and removes
   heterozygous single-strand variants — matching legacy
   `Analysis.concensus_variant()`.
    BLASTn retains its independent consensus-stage polyC policy: plus-strand
    calls at `316–335` and minus-strand calls at `283–302` are directionally
    flagged. Tracy uses a separate, more conservative unbounded per-trace
    directional coverage policy before merge, so the callers do not share
    boundaries or processing stage.

6. **Compute intervals** — `_compute_intervals(rows, regions)` extracts
   `forward_ref_start`/`forward_ref_end` ranges, merges them via
   `core.sample.merge_intervals`, and intersects each region's bounds via
   `core.variants.get_overlap`, producing `Dict[str, List[List[int]]]`.
7. **Handle polyC regions (post-hoc removal)** —
   `transform._handle_polyc_removal` removes low-quality variants in
   `SKIP_AFTER_POLYC` ranges and in affected `POLYC_REGIONS`: a variant is
   kept when its max per-trace quality is ≥ 40; 309.1C/309.2C (HV2) and
   16189T>C (HV1) are always kept; HV3/other variants are removed unless
   quality ≥ 40. Kept variants are tagged `polyc_region` / `kept_by_quality`.
   This is a distinct removal pass from the consensus-level
   `flagged_for_filtering` filter in stage 5; it ports
   `PostProcess.handle_polyC_regions`.
8. **Filter positions** — `_filter_positions` loads
   `rules/filter/filter.txt` and `rules/filter/no_filter.txt`
   (`_load_filter_positions` / `_load_no_filter_positions`), removes variants
   at filter positions unless they are in the no-filter set, and applies the
   dynamic **309.x rule** (any `309.x` position is filtered unless in
   `no_filter`). This ports `PostProcess.filter_special_positions`.
9. **Apply conversion rules** — `transform._apply_conversion_rules` loads
   `*.yaml`/`*.yml` from `directories.convert_rules` and applies
   `MustHave`/`MustNotHave`-gated `Change` (`Add`/`Remove`) operations. This
   ports `PostProcess.apply_conversion_rules`.
10. **Region validation** — `transform._apply_region_validation_rules`
    rejects variants seen only in non-native HV regions (native region from
    `get_hv_region_for_position`) when max quality < 40; positions in
    `no_filter_positions` are always kept. When `no_filter_positions is None`
    it performs a **deferred inline import** of
    `etl._load_no_filter_positions` to avoid the `etl ↔ transform` import
    cycle. This ports `PostProcess.apply_region_validation_rules`.
11. **Build `Variant` objects** — `_build_variant_objects` converts the raw
    dicts to Pydantic `Variant` records (pos/ref/seq/quality/files, with
    `peaks=None`).
12. **Flag variants** — `flag_variants(variant_list)` produces `variant_flags`;
    `transform._polyc_warning_flags(filtered_variants)` merges polyC,
    after-polyC, and 16189T>C warning strings into `variant_flags` (keyed
    `pos|ref|alt`); `SampleFlagger(variant_list, intervals).analyze()` produces
    sample-level flags; `deduplicate_sample_flags` drops any sample flag that
    duplicates a per-variant flag string. This ports
    `PostProcess.warning_variants`.
13. **Construct `Sample`** — returns `Sample(source_tool=Tool.BLASTN, ...)`
    with `variants`, `intervals`, `sample_flags`, and `variant_flags`
    populated; `hv1`/`hv2`/`hv3` and statistics fields (`no_snps`, `no_ins`,
    …) are `None` (filled by `Batch.write()`), `information` is `None` (filled
    downstream, e.g. `region_sources` by `write_region_jsons`), and
    `batch_id` is passed through unchanged.

## Variant Transforms (`transform.py`)

`transform.py` holds the variant-transformation logic extracted from `etl.py`.
All functions are internal (underscore-prefixed); `__all__` re-exports the five
used by `etl.py`: `_apply_conversion_rules`, `_apply_region_validation_rules`,
`_handle_polyc_removal`, `_polyc_warning_flags`, and `_right_norm_indel`.

**Indel right-normalization.** BLASTN places indels at the leftmost position
of a tandem repeat, but the canonical mtDNA numbering (shared with the other
tool pipelines) uses the rightmost equivalent position. `_right_norm_indel`
splits insertions and deletions, sorts each by position, and runs
`_normalize_indel_blocks` per type. `_normalize_single_indel` right-shifts a
single indel along identical trailing bases (insertions land at `{pos}.1`,
deletions at the integer position); `_normalize_indel_block` normalizes a run
of consecutive indels as a unit (deletions become per-base positions from the
block start; insertions are normalized individually and re-indexed `.1`,
`.2`, …). This is used by `etl._collect_variants_by_position`.

**YAML conversion rules.** `_load_conversion_rules` loads `*.yaml`/`*.yml`
(sorted) from `directories.convert_rules` (warns if missing).
`_check_must_have` / `_check_must_not_have` evaluate rule conditions against a
position-keyed variant dict (`"*"` wildcards ref/seq), and
`_apply_rule_changes` applies `Change.Remove` (pop) and `Change.Add` (insert
with `file=["conversion_rule"]`) in place. `_parse_variant_rule` parses
strings like `'309.1TC'` or `'16193CT'` into `(position, ref, seq)` (last two
chars = ref + seq). `_apply_conversion_rules` orchestrates these and returns
the updated variant list.

**polyC handling (post-hoc removal).** `_polyc_stretch_detected` flags a
polyC stretch when HV1 has 16189T>C, or when HV2 has any of 309.1/309.2/315.1
or its reference block has ≥ 7 C bases (`_HV2_C_COUNT_THRESHOLD = 7`).
`_max_quality` / `_keep_by_quality` keep a variant when its max per-trace
quality is ≥ `_QUALITY_KEEP_THRESHOLD = 40`. `_apply_in_polyc_rule` keeps
309.1C/309.2C in HV2 and 16189T>C in HV1, and removes HV3/other variants
unless quality ≥ 40; kept variants are tagged `polyc_region` (and
`kept_by_quality`). `_handle_polyc_removal` is the post-hoc removal pass
described in ETL stage 7 — distinct from the consensus-level
`flagged_for_filtering` filter.

**Region validation.** `_regions_from_files` derives `["HV1"]`/`["HV2"]`/`["HV3"]`
labels from source file names. `_apply_region_validation_rules` rejects
variants seen only in non-native HV regions when max quality < 40, keeping
`no_filter_positions` always; it loads `no_filter_positions` via the deferred
inline import noted above.

**polyC warning flags.** `_polyc_warning_flags` produces `variant_flags`
strings keyed by `pos|ref|alt`: `"polyC region {region}: potential length
heteroplasmy"` (in `POLYC_REGIONS`), `"near polyC region {region}: poor signal
quality"` (in `SKIP_AFTER_POLYC`), and `"16189T>C creates polyC stretch"`.
These are consumed by ETL stage 12.

Module-level constants: `_QUALITY_KEEP_THRESHOLD = 40`,
`_HV1_POLYC_CREATING_POS = 16189`, `_HV2_C_COUNT_THRESHOLD = 7`.

## BLASTN TSV Format

The module expects BLASTN tabular output **format 7** (tab-separated, with
`#`-prefixed comment lines):

```
# BLASTN 2.14.0+
query1    subject1    plus    580    99.48    601    3    0    1    601    16024    16624    ...
```

The 16 columns are defined in `utils.BLASTN_COLUMNS`:

```
qseqid, sseqid, sstrand, nident, pident, length, mismatch, gapopen,
qstart, qend, sstart, send, qseq, sseq, evalue, bitscore
```

`utils.BLASTN_COMMENT_PREFIX = "#"` is used for pandas
`read_csv(..., comment=...)`. `parse_blastn_tsv` raises `FileNotFoundError` if
the path is missing and `ValueError` if the column count is not 16; it
converts the DataFrame to plain records immediately so the rest of the ETL
works with plain dicts.

`utils` also provides the strand helpers used by ETL stage 3
(`normalize_strand`, `normalize_strand_all`), `extract_quality_scores` (returns
`{}` with a warning if the file is missing; skips invalid lines with a
warning), `run_command` (wraps `subprocess.run` with `capture_output=True,
text=True, check=True`, logging stderr on failure), and `get_files_by_id`
(maps each ID from a TXT list to AB1 file paths whose filename contains the
ID, logging an error for IDs with no match).

## Output Layout

Both pipelines construct paths via `core.paths.OutputPaths(output_dir)` and
call `paths.create_dirs()`, producing three subdirectories under the base
output dir:

```
<output_dir>/
├── json/              # Batch.write() output
│   ├── statistic_fullbatch.json
│   └── <LID>/<LID>.json
├── regions/           # Per-region intermediate JSONs (write_region_jsons)
│   ├── HV1/HV1_{LID}.json
│   └── HV2-3/HV2-3_{LID}.json
└── preprocess/        # Tool-specific preprocessing artifacts
    └── <LID>/
        ├── *.fastq, *.trim.fastq, *.trim.quality, *.trim.fasta             (AB1 pipeline)
        ├── <LID>.concatenate.fasta, <LID>.concatenate.quality              (AB1 pipeline)
        ├── <LID>.concatenate.blastn                                        (AB1 + FASTA pipelines)
        └── (no .concatenate.quality for FASTA pipeline)
```

`Batch.write()` is called with `nest_batch_id=False` because the base output
dir is already batch-specific, so per-LID JSONs land directly under `json/`
as `<LID>/<LID>.json` (matching the Tracy reference pipeline).

### Region JSON pattern

Region JSON file patterns are defined in `core.region.REGION_FILE_PATTERNS`.
Each LID produces one `HV1_{LID}.json` and one `HV2-3_{LID}.json` with merged
primer data:

| Region group | Pattern | Path | Region keys (`REGION_TO_KEYS`) |
|--------------|---------|------|--------------------------------|
| `HV1` | `HV1_{lid}.json` | `regions/HV1/HV1_{LID}.json` | `["HV1"]` |
| `HV2-3` | `HV2-3_{lid}.json` | `regions/HV2-3/HV2-3_{LID}.json` | `["HV2", "HV3"]` |

`write_region_jsons()` writes the region JSONs from the **combined** samples
(after same-LID merging), using data-driven `REGION_TO_KEYS` groups. BLASTN
does not set auto-pass flags, so region groups are NOT derived from
`sample_flags` (unlike Mutation Surveyor); groups with no variants/intervals
are skipped. For each group it filters the combined `Sample` to the group's
region keys (`filter_sample_by_regions`), writes the filtered sample dict
(`sample_to_dict`), and adds `information.region_sources = {<key>: <source_tool>}`
provenance. After writing, `validate_region_group_intervals` runs range QC and
logs warnings for intervals extending beyond canonical region boundaries.

## Batch Processing

### AB1 pipeline (`pipeline.py`)

`pipeline.process_batch(input_dir, output_dir, options=None) -> dict[str, Sample]`
orchestrates the AB1 pipeline:

1. Read sample IDs from `options.samples_path` (returns `{}` if missing).
2. Map IDs → AB1 files via `utils.get_files_by_id`.
3. For each sample, in parallel via `ProcessPoolExecutor`
   (`max_workers = BlastnSettings.max_workers or os.cpu_count()`):
   preprocess with `preprocess_sample`, then ETL with `etl.process`.
4. Combine same-LID Samples (`_combine_same_tool_samples`), merging
   forward/reverse primer traces.
5. Write per-region JSONs via `write_region_jsons()` from the combined
   samples.
6. Call `Batch.write(..., nest_batch_id=False)`.
7. Return `dict[str, Sample]` mapping `sample_id` → combined `Sample`.

`BatchOptions` carries only paths and identity; tuning parameters
(`trim_threshold`, `word_size`, `max_workers`) come from `BlastnSettings`, not
this dataclass or the CLI.

```python
@dataclass
class BatchOptions:
    ref_path: str = "ref/rCRS.fasta"
    batch_id: str | None = None
    samples_path: Path | None = None
```

### FASTA pipeline (`fasta_pipeline.py`)

For pre-assembled FASTA files (one sample per file, ID = filename stem) instead
of AB1 traces, `fasta_pipeline.process_batch_fasta(input_path, output_dir,
options=None) -> dict[str, Sample]` collects FASTA files via `_get_fasta_files`
(recursive, by `FASTA_EXTENSIONS = (".fasta", ".fa", ".fas", ".fna")`,
case-insensitive via `suffix.lower()`), runs one BLASTN alignment per file
(`_blastn_fasta`, outfmt 7 — same 16 columns as the AB1 pipeline), ETLs each
via `etl.process`, and then reuses `_build_combined_samples`,
`_write_full_region_jsons`, and `Batch.write(..., nest_batch_id=False)`, so
finalization is identical to the AB1 pipeline. FASTA carries no per-base
quality, so no `.concatenate.quality` is produced and the ETL uses empty
quality scores (same as the legacy FASTA path).

```python
@dataclass
class FastaBatchOptions:
    ref_path: str = "ref/rCRS.fasta"
    batch_id: str | None = None
```

## Configuration

BLASTN tuning parameters live in `src/config.py` as
`class BlastnSettings(BaseSettings)` with env prefix `MTDNA_BLASTN_`
(`model_config = SettingsConfigDict(env_prefix="MTDNA_BLASTN_", extra="forbid")`).
External tool paths (UGENE, BLASTN, seqtk) come from `get_settings().tools`
(`ToolsSettings`, env prefix `MTDNA_TOOLS_`).

| Field | Type | Default | Env var | Description |
|-------|------|---------|---------|-------------|
| `trim_threshold` | `float` | `0.01` | `MTDNA_BLASTN_TRIM_THRESHOLD` | Quality threshold for seqtk trimming |
| `word_size` | `int` | `22` | `MTDNA_BLASTN_WORD_SIZE` | Word size for BLASTN alignment |
| `quality_threshold` | `int` | `35` | `MTDNA_BLASTN_QUALITY_THRESHOLD` | Minimum quality score (not consumed by the BLASTN tool) |
| `min_peak_value` | `int` | `125` | `MTDNA_BLASTN_MIN_PEAK_VALUE` | Minimum peak value (not consumed by the BLASTN tool) |
| `heteroplasmy_threshold` | `float` | `0.8` | `MTDNA_BLASTN_HETEROPLASMY_THRESHOLD` | Minimum heteroplasmy ratio (not consumed by the BLASTN tool) |
| `max_workers` | `int \| None` | `None` | `MTDNA_BLASTN_MAX_WORKERS` | Sample-level parallelism; `None` → `os.cpu_count()` |
| `blastn_threads` | `int` | `1` | `MTDNA_BLASTN_BLASTN_THREADS` | Threads per BLASTN subprocess (keep 1 with sample-parallel `max_workers`) |

## CLI

### AB1 input (`pipeline.py`)

```bash
python src/tools/blastn/pipeline.py -i <dir> -o <dir> -s <txt> [--batch-id ID] [--ref-path PATH]
```

| Arg | Required | Default | Purpose |
|-----|----------|---------|---------|
| `-i` / `--input-dir` | yes | — | Directory containing AB1 files |
| `-o` / `--output-dir` | yes | — | Base output dir (`json/`, `regions/`, `preprocess/` auto-created) |
| `-s` / `--samples` | yes | — | Path to TXT file with sample IDs (one per line) |
| `--batch-id` | no | `None` | Batch identifier |
| `--ref-path` | no | `ref/rCRS.fasta` | Path to rCRS reference FASTA |

Example:

```bash
python src/tools/blastn/pipeline.py \
  -i data/raw/BATCH_001 \
  -o results/tools/blastn/BATCH_001 \
  -s data/raw/BATCH_001.txt \
  --batch-id BATCH_001 \
  --ref-path ref/rCRS.fasta
```

### FASTA input (`fasta_pipeline.py`)

```bash
python -m src.tools.blastn.fasta_pipeline -i <file|dir> -o <dir> [-r PATH] [--batch-id ID]
```

| Arg | Required | Default | Purpose |
|-----|----------|---------|---------|
| `-i` / `--input` | yes | — | FASTA file or directory of FASTA files (searched recursively) |
| `-o` / `--output-dir` | yes | — | Base output dir (`json/`, `regions/`, `preprocess/` auto-created) |
| `-r` / `--ref-path` | no | `ref/rCRS.fasta` | Path to rCRS reference FASTA |
| `--batch-id` | no | `None` | Batch identifier |

Example:

```bash
python -m src.tools.blastn.fasta_pipeline \
  -i data/fasta/BATCH_001 \
  -o results/tools/blastn/BATCH_001 \
  -r ref/rCRS.fasta \
  --batch-id BATCH_001
```

Each `main()` parses argparse, builds the options dataclass, calls the batch
function, and logs the count of processed samples.

## Public API

The key entry points; anyone needing full parameter detail can read the
source.

- **`etl.process(sample_id, input_path, ref_path, batch_id=None) -> Sample`**
  — pure ETL from a BLASTN format-7 TSV to a `Sample(source_tool=Tool.BLASTN)`
  (the 13-stage pipeline above).
- **`pipeline.process_batch(input_dir, output_dir, options=None) -> dict[str, Sample]`**
  — AB1 batch: preprocess + ETL in parallel, combine same-LID, write region
  JSONs and `Batch.write()`.
- **`pipeline.BatchOptions`** — dataclass (`ref_path`, `batch_id`,
  `samples_path`) for AB1 batch options.
- **`fasta_pipeline.process_batch_fasta(input_path, output_dir, options=None) -> dict[str, Sample]`**
  — FASTA batch: one BLASTN alignment per file, then the shared ETL/combine/
  write finalization.
- **`fasta_pipeline.FastaBatchOptions`** — dataclass (`ref_path`, `batch_id`)
  for FASTA batch options.
- **`fasta_pipeline.FASTA_EXTENSIONS`** — `(".fasta", ".fa", ".fas", ".fna")`,
  the FASTA extensions recognized by `_get_fasta_files`.
- **`preprocessing.preprocess_sample(sample_id, list_ab1, output_dir, ref_path) -> PreprocessResult`**
  — the 8-step AB1 → BLASTN TSV subprocess pipeline.
- **`preprocessing.PreprocessResult`** — named tuple
  (`blastn_tsv_path`, `quality_path`, `concatenated_fasta_path`).

The `utils` helpers (`parse_blastn_tsv`, `normalize_strand` /
`normalize_strand_all`, `extract_quality_scores`, `run_command`,
`get_files_by_id`, and the `BLASTN_COLUMNS` / `BLASTN_COMMENT_PREFIX`
constants) and the `transform` helpers (indel right-norm, YAML conversion
rules, polyC removal, region validation, polyC warning flags) are internal
building blocks described in their respective sections above.

## Differences from Legacy

The legacy FASTA-to-JSON path subclassed the legacy `SampleProcessor` and used
the old root-level BLASTn/analysis/post-processing modules; those modules have
been removed and the current module is the sole implementation. Key
differences: variants are `List[Variant]` (Pydantic) instead of
`dict[float, dict]`; consensus and statistics are delegated to `Batch.write()`
rather than `Analysis`; post-processing is integrated into ETL rather than a
separate `PostProcess` class. Flagging (`flag_variants()` + `SampleFlagger`)
and per-region JSONs (`write_region_jsons()`) are new; pandas is used only in
`parse_blastn_tsv()` and converted to plain data immediately; configuration
uses `get_settings()` (Pydantic) instead of `os.getenv()`; logging uses
`loguru` instead of `print()`. Indel calling, the polyC consensus filter,
`PostProcess.handle_polyC_regions` post-hoc removal,
`apply_region_validation_rules`, and `warning_variants` were all ported into
the new ETL/transform functions with consensus-level parity. Per-variant
quality is always stored as a list and `_max_quality` handles scalar/list; the
legacy `Analysis.concensus_variant` stored scalar quality for single-alignment
variants, which made `PostProcess` crash (`max(int)`) on ~6% of real samples
during `apply_region_validation_rules`, so those samples now produce output
instead of being skipped.

## Cross-References

- [README.md](../README.md) — pipeline overview and tool list
- [ARCHITECTURE.md](../ARCHITECTURE.md) — §12 Standard Tool Module Template
- [variants.md](../variants.md) — variant model, numbering, and HV-region conventions