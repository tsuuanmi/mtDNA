# Tracy Integration

> Tracy-based variant calling from AB1 trace files for mtDNA analysis.

## Overview

`src/tools/tracy/` turns raw AB1 chromatograms into the standard `Sample`
model used by the comparison, flagging, and batch aggregation layers. Each AB1
file is first decomposed by the external **Tracy** `decompose` command, which
aligns the trace to rCRS and emits SNPs, insertions, and deletions with peak
and quality data. The module then parses Tracy's JSON output, applies a series
of position- and primer-specific transforms, filters low-quality calls, and
builds a `Sample` with `source_tool=Tool.TRACY`.

Tracy is an **alternative analysis path** to the BLASTN pipeline: both produce
per-sample JSON that feeds into the same comparison, FASTA generation, and
unification stages. Tracy additionally supports multi-sample batch processing,
where traces from several primers (forward/reverse, HV1/HV2/HV3) covering the
same LID are merged into one consensus profile across all hypervariable regions.

The module follows the [Standard Tool Module Template](README.md) and shares
its shape with Sequencher and Mutation Surveyor.

## Module Structure

```
src/tools/tracy/
├── __init__.py        # empty package marker
├── etl.py             # pure ETL: one Tracy decompose JSON -> Sample
├── pipeline.py        # batch orchestration, masking, report writing, CLI
├── preprocessing.py   # AB1 -> Tracy decompose (external subprocess)
├── quality_control.py  # trace metrics, noisy ranges, and report records
├── noise_mask.py      # remove noisy-range variants and recompute flags
├── transforms.py      # position-specific, polyC, and strand variant transforms
└── utils.py           # shared trace, coordinate, overlap, and peak helpers
```

`etl.py` is the pure core: it reads one Tracy decompose JSON, applies
transforms and filters, and returns a `Sample`, never touching files beyond its
input and never importing from `pipeline.py`. `pipeline.py` orchestrates the
batch — it reads sample IDs, decomposes and analyzes traces, runs ETL and
per-trace masking in `ProcessPoolExecutor` workers, combines same-LID reads, and
writes the `json/`, `regions/`, and `preprocess/` trees.
`preprocessing.py` is the only module that invokes the external Tracy binary,
and it is called only by `pipeline.py`. `quality_control.py` is the single
source of truth for trace metrics, 10–20-base window analysis, noisy range
detection, and QC report records. `noise_mask.py` applies those ranges to a
completed per-trace `Sample` and recomputes flags. `utils.py` owns shared peak
extraction, alignment normalization, coordinate mapping, and range-overlap
logic. Data flows `preprocessing.decompose_sample()` ->
`quality_control.analyze_trace_file()` -> `etl.process()` ->
`noise_mask.apply_noise_mask()` -> `preprocess/<sample_id>/qc_report.json` ->
trace combination -> `Batch` (`gen_seq_from_variants_regions()`,
`statistic_variants()`, `write()`) -> downstream comparison, FASTA generation,
and unification.

## Preprocessing: Tracy Decompose

`preprocessing.decompose_sample()` turns AB1 traces into Tracy JSON. For each
sample it locates every `.ab1` file in the input directory whose name contains
the sample ID, resolves the Tracy binary through
`shutil.which(get_settings().tools.tracy)`, and runs Tracy's `decompose`
subcommand once per AB1 file:

```
<tracy_path> decompose <ab1_file> \
    --genome <ref_path> \
    --callVariants \
    --trim <config.trim> \
    --pratio <config.pratio> \
    --maxindel <config.maxindel> \
    --outprefix <sample_output_dir>/<ab1_stem>
```

Tracy decomposes the chromatogram against rCRS and calls variants; with
`--outprefix` set to `preprocess/<sample>/<stem>`, it writes its artifacts
(`<stem>.align`, `<stem>.decomp`, `<stem>.bcf`, `<stem>.json`) into a per-sample
preprocess directory. The function returns only the JSON paths successfully
produced in the current run. This prevents QC and ETL from consuming stale JSON
files left in the sample directory. Missing inputs, a missing binary, or failed
subprocesses produce no returned path; failures are logged and the batch
continues.

The three decompose knobs live in `DecomposeConfig`, a `NamedTuple` of
`trim` (`int`, default `7`, bases trimmed from trace edges), `pratio` (`float`,
default `0.3`, peak ratio threshold for decomposition), and `maxindel` (`int`,
default `1000`, maximum indel length to detect). When `decompose_sample` is
called with `config=None`, it builds the `DecomposeConfig` from
`get_settings().tracy`. The Tracy executable path itself is not a
`TracySettings` field — it lives in `ToolsSettings.tracy` (env
`MTDNA_TOOLS_TRACY`, default `./tools/tracy`, validated as an executable
`FilePath`).

## ETL Processing

`etl.process(sample_id, input_path, *, ref_seq, config=None) -> Sample` is the
single entry point of the pure ETL. It takes a sample ID (LID), the path to one
Tracy decompose JSON, the rCRS reference string, and an optional
`_ProcessConfig` (when `None`, built from `get_settings().tracy` plus
`get_settings().regions.REGIONS`). It returns a `Sample` with
`source_tool=Tool.TRACY`, populated `variants`, `intervals`, `variant_flags`,
`sample_flags`, and `batch_id`; the `hv1`/`hv2`/`hv3` fields and all stats
fields (`no_snps`, `no_ins`, `no_dels`, `no_identicals`, `no_unread`) are left
`None`, and `information=None`. HV sequences and counts are owned by `Batch` —
the ETL must not compute them. The stages run in this order:

1. **Build config** — resolve `_ProcessConfig` from settings if needed;
   `variant_regions` from `config.variant_regions` or fall back to
   `get_settings().regions.REGIONS`; use `config.sample` for logging, else
   `sample_id`.
2. **Load JSON** (`_load_tracy_json`) and **normalize `ref1pos`**
   (`utils.normalize_ref_positions`) — Tracy >= 0.7.8 emits `ref1pos` as a single integer
   (the start position); expand it into a per-column array using `ref1align`
   (gap columns reuse the previous position minus one when beyond `start`);
   returns the scalar start.
3. **Detect primer type** (`detect_primer_type`) — derive `PrimerTypeInfo` from
   `input_path.stem`.
4. **Compute trim bounds** (`utils.alignment_bounds`) — measure leading/trailing
   dash counts in `alt1align`/`ref1align`; `trim_left`/`trim_right` are the max
   alt/ref dash counts at each edge and define inclusive
   `start_index`/`end_index`.
5. **Call raw variants** (`_call_variant_at_index` over `[start_index, end_index]`)
   — compare `alt1align[i]` vs `ref1align[i]`, skip matches. Peak data is read
   at `peak_index = i - deletions` (a running deletion counter offsets the
   basecall index). Each raw dict carries `pos`, `ref`, `seq`, `peaks`
   `[A,C,G,T]`, `index`, `peak_index`, `quality`; when `alt != "-"`,
   `detect_heteroplasmy` overrides `seq` with an IUPAC base when not `"N"`.
6. **Update variant positions** (`update_variant_positions`) — remap raw
   indices to genomic/decimal positions based on strand/primer type, re-detect
   heteroplasmy after reverse-strand complementation, and number consecutive
   insertions at the same base (`315.1`, `315.2`, `513.1`...`513.4`). Returns
   `file_intervals` (`{filename: [[start, end]]}`).
7. **HV2F position-73 special case** (`_maybe_add_hv2f_pos73`) — when the
   primer is HV2F and the first variant's genomic position is beyond 73, create
   an SNP at position 73 (skipped if ref==seq) and append the `[73, 75]`
   interval.
8. **Detect variant conditions** (`detect_variant_conditions`) — single pass
   over the decimal-positioned variants collecting every condition flag the
   transforms consume.
9. **Apply transforms** (`apply_all_transformations`) — position-specific
   dispatch, AC-repeat insertion remap (`513.x` -> `524.x`), HV2 polyC
   transforms (primer removal, conversion-deletion creation, `16189`/`16193`
   handling), then canonical polyC insertion creation. Mutates `raw_variants`
   in place and may append variants.
10. **Filter variants** (`_should_skip_variant`) — drop variants marked
    `remove`, with `"N"` seq, outside `variant_regions` intervals, HV1
    deletions except at `POS_16189`/`POS_16193`, low-quality variants
    (`quality < quality_threshold`, except `POS_73`), and variants failing
    `validate_peak_quality`. Strip `index`/`peak_index` from survivors.
11. **Build `Variant` objects** (`_build_variant_objects`) — peaks become
    `list[list[int|None]]`, quality `list[int]`; invalid variants are logged
    and skipped.
12. **Compute intervals** (`_compute_intervals`) — intersect `file_intervals`
    with each region's bounds and merge overlaps via `merge_intervals`.
13. **Flag** — per-variant via `flag_variants(variant_list)`, sample-level via
    `SampleFlagger(...).analyze()`, then deduplicate sample-level flags against
    per-variant flags.
14. **Construct and return the `Sample`.**

The critical transform order is `update_variant_positions()` ->
`_maybe_add_hv2f_pos73()` -> `detect_variant_conditions()` ->
`apply_all_transformations()`. Position update runs before condition detection
so insertion positions are decimal strings, which the AC-repeat remap inside
`apply_all_transformations` relies on.

`config` is the `_ProcessConfig` `NamedTuple`, carrying the filtering parameters
`process()` does not pull from settings directly: `quality_threshold` (`int`),
`min_peak_value` (`int`), `pratio` (`float`), `heteroplasmy_threshold`
(`float`), `variant_regions` (`dict[str, list[int]] | None`, `None` ->
`get_settings().regions.REGIONS`), `batch_id` (`str | None`, propagated to
`Sample.batch_id`), and `sample` (`str | None`, a free-text label for log
messages, `None` -> `sample_id`).

### Shared directional polyC suppression

Tracy calls candidate variants before applying directional polyC policy; it
does not trim or mask a sequence flank before calling. `src.core.polyc` owns
the directional classifications and reasons, while Tracy maps its known primer
provenance to a direction and applies the result before same-LID merge.

| Policy | Forward suppression | Reverse suppression | Affected primers |
|---|---|---|---|
| HV2 directional coverage | `>=304` | `<=315` | HV2F/HV3F, HV2R/HV3R |
| HV1 16189-created polyC | strictly above `16189` | strictly below `16189` | HV1F, HV1R |

For normal trace coverage, HV2F supplies raw candidates only through `73–303`,
and HV3R supplies raw candidates from `316–576`; canonical HV2 polyC insertions
remain separate. The shared directional filters stay unbounded so expanded traces
receive the same protection.

The HV1 policy is triggered by a called `C` anchored at 16189; the trigger
position itself remains usable. Its separate 16189-deletion conversion and
synthetic 16193 deletion remain Tracy-specific nomenclature transforms. The
HV2 internal polyC region (`303–315`) likewise retains dedicated canonical
insertion normalization; directional coverage governs its trusted raw calls.

A candidate in an applicable direction is marked removed with the shared
reason, logged with its trace, allele, and quality, and excluded from the
trusted per-trace variant list before same-LID merge. Here, removed means
*detected but excluded from consensus*; it does not mean that the input
sequence was discarded before calling. Unknown primer provenance is left
unchanged. These policies are deliberately conservative: future chromatogram
quality, peak separation, phasing, noise, or opposite-trace evidence may allow
reliable flagged calls to be recovered without changing the raw calling stage.

## Variant Transforms

`transforms.py` holds the pure functions that rewrite variant dicts in place.
It is the mechanism behind the canonical mtDNA numbering — the rightmost
equivalent position within a tandem-repeat run, matching the regenerate,
Sequencher, and FASTA pipelines so that variants round-trip the consensus.

`apply_all_transformations(variants, conditions, primers, params) -> None` is the single entry point. For each variant it runs, in order: position-specific dispatch (`_apply_position_specific_transforms`, driven by the `_POSITION_TRANSFORMS` table mapping integer positions to handlers — e.g. `248`/`250` -> 249 ref=`A`, `456` -> 459, `16189` deletion seq -> `C`), then AC-repeat insertion remap (`_apply_insertion_remap`, sliding `513.x` insertions to `524.x`), then polyC transforms (`_apply_polyc_transforms`: insertion counting, shared-policy primer filtering, `16189`/`16193` handling, conversion-deletion creation). After the loop it appends any additional variants created by transforms and then creates canonical HV2 polyC insertions from the `HV2_polyC` count via `_create_polyc_insertions` — patterns keyed by count (`1` -> `["315.1"]`, `2` -> `["309.1","315.1"]`, `3` -> `["309.1","309.2","315.1"]`), with the count increased by one when both `has_309_C_T` and `has_310_T_C` hold.

`update_variant_positions(variants, position_params, start_index, end_index)` does the strand-aware position remap and returns `file_intervals`. Forward and HV2F traces map `pos = normalize_position(ref1pos + pos_base(pos) - insertions)`; reverse traces go through `_transform_reverse_strand`, which computes `pos = normalize_position(ref1pos + ref_align_len - ref_end_dash - 1 - index - insertions)`, complements `ref`/`seq` via `Seq.reverse_complement()`, and swaps the peak channels A<->T and G<->C. After the transform, `_redetect_heteroplasmy` re-runs `detect_heteroplasmy` on SNPs (deletions/insertions skipped) and marks the variant `remove` when `ref == seq`. Insertions (`ref == "-"`) are numbered with a shared `insertion_count_at_pos` counter as `f"{base_pos - 1}.{count}"`, producing the decimal-string positions the AC-repeat remap depends on; reverse-strand variants are processed in descending `index` order. When `start_index`/`end_index` is `None` (no variants in the trimmed region) it returns an empty dict.

## Utilities

`utils.py` is a pure helper library — no file I/O, no side effects, no
`Sample`/`Variant` dependency. It carries the mtDNA position constants and the
primer/peak/condition logic the ETL relies on.

The position constants encode the special sites referenced across filters and transforms: `POS_73` (HV2F position-73 special case — a pass-through in `validate_peak_quality`, a quality-filter exemption in `_should_skip_variant`, and the HV2F pos-73 insertion via `_maybe_add_hv2f_pos73`), `POS_16000` (the HV1-region guard `pos > 16000` for deletion skipping), `POS_16189` and `POS_16193` (HV1 polyC handling — the 16189 deletion/`T>C` conditions and the deletion skip-exemption, with 16193 created alongside the 16189 deletion). Additional constants cover the other transform/condition sites (`POS_248`, `POS_250`, `POS_300`, `POS_309`, `POS_310`, `POS_315`, `POS_316`, `POS_320`, `POS_456`, `POS_459`, `POS_513`-`POS_515`, `POS_523`-`POS_525`), and `NUM_PEAK_CHANNELS = 4` (A, C, G, T).

`PrimerTypeInfo` is a `NamedTuple` of boolean primer flags derived by substring checks on the trace filename: `is_forward` (`HV1F` or `HV3F`), `is_reverse` (`HV1R`, `HV2R`, or `HV3R`), `is_hv2f` (`HV2F`), `is_hv2f_hv3f` (`HV2F` or `HV3F`), `is_hv2r_hv3r` (`HV2R` or `HV3R`), `is_hv1f` (`HV1F`), `is_hv1r` (`HV1R`); `detect_primer_type(filename) -> PrimerTypeInfo` builds it.

`detect_heteroplasmy(peaks, position, heteroplasmy_threshold=None, sample=None)` inspects the four peak heights and returns an IUPAC code: `"N"` when all peaks are zero or no base reaches `max_peak * threshold`, the single passing base when exactly one passes, or `nucleotide_to_iupac(detected_bases)` when several pass (a heteroplasmic call, logged at INFO). A `None` threshold falls back to `get_settings().tracy.heteroplasmy_threshold`.

`validate_peak_quality(variant, min_peak_value=None, pratio=None, sample=None)` returns `True` when the variant passes. In order: position 73 with `ref="A"`, `seq="G"` always passes; deletions always pass; an all-`None` peak list fails; `max_peak < min_peak_value` fails (logged); heteroplasmic variants (IUPAC ambiguity codes) auto-pass; for homoplasmic variants `highest_peak_ratio = max_peak / sum(peaks)` must exceed `pratio` or the variant fails (logged). Missing `min_peak_value`/`pratio` fall back to settings.

`detect_variant_conditions(variants) -> dict` is the single pass that collects the condition flags the transforms consume — booleans such as `has_248_deletion`, `has_309_C_T`, `has_310_T_C`, `has_456_deletion`, `has_513_*_insertion`, `has_16189_T_C`, `has_16189_deletion`, and the integer counter `HV2_polyC` (incremented for each insertion with `ref == "-"` and `309 <= pos_int <= 316`). Insertion positions are matched as exact decimal strings (e.g. `"513.1"`), so it must run after `update_variant_positions`.

## Merging Variants Across Reads

`etl.merge_variants(variants: list[Variant]) -> list[Variant]` collapses the
variants from multiple Tracy reads — forward and reverse primers, or several HV
region groups — for one LID into one IUPAC consensus per `(pos, ref)`, sorted by
`pos_sort_key`. It groups variants by `(pos_sort_key(pos), ref)`. A
single-element group is kept as-is. For a multi-element group, if the group
mixes deletions (`"-"`) with substitutions the variants are kept separate
(IUPAC does not cover deletions). Otherwise the group collapses into one
variant: `seq = "-"` if any member is a deletion, else
`nucleotide_to_iupac(seqs)` — which also unifies identical calls
(`T + T -> "T"`), so no separate deduplication pass is needed. The `files`,
`quality`, and `peaks` of the merged variant are concatenated from every source
in first-appearance order. This is the function behind the combine-same-LID
flow in the batch pipeline.

## Batch Processing

`pipeline.process_batch(input_dir, output_dir, options=None) -> dict[str, Sample]`
runs the full batch. It defaults `options` to `BatchOptions()`, coerces the
directories to `Path`, builds `OutputPaths(output_dir)` and calls
`create_dirs()` (creating `json/`, `regions/`, `preprocess/`), and loads the
shared reference once via `SeqIO.read(options.ref_path, "fasta")`.

Sample IDs come exclusively from a TXT file (`options.samples_path`, one
non-empty stripped line per ID, read by `_read_sample_ids`) — there is no AB1
glob fallback. If the path is `None` or the file is empty, the batch logs and
returns `{}`. Parallelism is controlled by
`max_workers = get_settings().tracy.max_workers or max(1, os.cpu_count() or 1)`.

Processing has two parallel phases. First, `_prepare_sample` runs decompose,
keeps each newly produced `*.json` whose name contains `"HV"`, and analyzes
those exact files with `quality_control.analyze_trace_file()`. Second,
`_process_prepared_sample` runs `etl.process()` and, when
`noise_mask_enabled=true`, applies only that trace's noisy ranges before any
same-sample merge. The parent adds exclusions to each trace result and writes
`preprocess/<sample_id>/qc_report.json`, then combines masked samples. Workers
never write reports. QC records malformed or unreadable trace data as `error` or
`unavailable`; ETL failures in the handled `ValueError`, `KeyError`, `OSError`,
and `RuntimeError` classes are logged while the batch continues.

Same-LID samples are then combined by `_build_combined_samples` ->
`_combine_same_tool_samples`: results are grouped by `sample_id`, variants are
merged via `etl.merge_variants`, intervals via `merge_intervals` per region,
flags are merged, and `SampleFlagger.analyze()` plus `flag_variants()` are
re-run on the merged variant set so that cross-region sample-level flags and
IUPAC-consensus flags are detected. This yields one combined `Sample` per LID.

For each combined sample, `_write_full_region_jsons` emits the region groups
`{"HV1": REGION_TO_KEYS["HV1"], "HV2-3": REGION_TO_KEYS["HV2-3"]}` via
`write_region_jsons(...)` into `regions/` (groups without data are skipped),
then runs `validate_region_group_intervals` per group and logs any QC warnings.
Finally, `Batch(list(combined.values())).write(str(paths.json_dir), options.batch_id or "unknown", nest_batch_id=False)` writes `statistic_fullbatch.json` plus per-LID JSONs flat under `json/`. The base output dir is already batch-specific, so `nest_batch_id=False` avoids a duplicated `<batch_id>` path segment. Completion counts are logged, and a warning is emitted if fewer samples were processed than loaded.

## Output Layout

```
<output_dir>/                         # base output dir (typically results/tools/tracy/<batch_id>)
├── json/                             # Batch.write(..., nest_batch_id=False)
│   ├── statistic_fullbatch.json      # full-batch statistic
│   └── <LID>/<LID>.json              # per-LID JSON (flat under json_dir)
├── regions/                          # write_region_jsons() output
│   ├── HV1/
│   │   └── HV1_<LID>.json            # per-LID HV1 region JSON (skipped if no data)
│   └── HV2-3/
│       └── HV2-3_<LID>.json          # per-LID HV2+HV3 region JSON (skipped if no data)
└── preprocess/                       # per-sample Tracy artifacts and QC reports
    └── <sample_id>/                  # one subdir per sample
        ├── qc_report.json            # noise QC for this sample's traces
        ├── <stem>.align
        ├── <stem>.decomp
        ├── <stem>.bcf
        └── <stem>.json               # decompose JSON consumed by etl.process()
```

Each `preprocess/<sample_id>/qc_report.json` describes that sample's Tracy
traces, including clean, suspicious, likely-noisy, unavailable, and error
results. It contains mask settings, window thresholds, per-trace metrics,
merged reference/alignment ranges, and variants excluded from each trace. A
mask range requires at least two overlapping or directly adjacent likely-noisy
windows; an isolated window is reported as suspicious and cannot exclude
variants.
Likely-noisy variants are absent from sample, region, and batch JSON output by
default, and the matching reference spans are removed from the originating
trace's intervals. A clean trace for the same LID can still provide coverage for
that span after merging. Set `MTDNA_TRACY_NOISE_MASK_ENABLED=false` to keep
variants and intervals while still producing QC reports. The pipeline removes
the obsolete batch-level `preprocess/qc_report.json` when reusing an output
directory.

`OutputPaths` exposes `base_dir`, `json_dir` (`base_dir / "json"`), `regions_dir` (`base_dir / "regions"`), and `preprocess_dir` (`base_dir / "preprocess"`); `create_dirs()` creates all three with `parents=True, exist_ok=True`. Region JSON filenames come from `REGION_FILE_PATTERNS` in `src/core/region.py` (`{"HV1": "HV1_{lid}.json", "HV2-3": "HV2-3_{lid}.json"}`), formatted with the sample's `sample_id`; region groups are mapped via `REGION_TO_KEYS` (`{"HV1": ["HV1"], "HV2-3": ["HV2", "HV3"]}`). Each region JSON is written into a per-group subdirectory and contains only that region's filtered data plus an `information.region_sources` provenance map; groups without intervals/variants are skipped.

## Configuration

`TracySettings` is a Pydantic `BaseSettings` subclass in `src/config.py` with
`env_prefix="MTDNA_TRACY_"` and `extra="forbid"`. Every field maps 1:1 to an
environment variable by upper-cased name under the prefix. The settings
instance is accessed via `get_settings().tracy`, and `TracySettings` is
attached to the root settings as `tracy: TracySettings` (default factory
`TracySettings`).

| Field | Type | Default | Env Var | Description |
|-------|------|---------|---------|-------------|
| `trim` | `int` | `7` | `MTDNA_TRACY_TRIM` | Tracy decompose `--trim`; bases trimmed from trace edges |
| `pratio` | `float` | `0.3` | `MTDNA_TRACY_PRATIO` | Tracy decompose `--pratio`; peak ratio threshold for decomposition |
| `maxindel` | `int` | `1000` | `MTDNA_TRACY_MAXINDEL` | Tracy decompose `--maxindel`; maximum indel length to detect |
| `quality_threshold` | `int` | `35` | `MTDNA_TRACY_QUALITY_THRESHOLD` | Minimum Phred-like quality score for variant filtering |
| `min_peak_value` | `int` | `125` | `MTDNA_TRACY_MIN_PEAK_VALUE` | Minimum chromatogram peak height for variant filtering |
| `heteroplasmy_threshold` | `float` | `0.8` | `MTDNA_TRACY_HETEROPLASMY_THRESHOLD` | Minimum peak ratio (relative to max peak) for IUPAC/heteroplasmic calls |
| `max_workers` | `int \| None` | `None` | `MTDNA_TRACY_MAX_WORKERS` | Max parallel `ProcessPoolExecutor` workers; `None` -> `os.cpu_count() or 1` |
| `noise_mask_enabled` | `bool` | `true` | `MTDNA_TRACY_NOISE_MASK_ENABLED` | Mask variants and subtract coverage that overlap a trace's likely-noisy ranges; QC reporting always remains enabled |
| `noise_window_size` | `int` | `15` | `MTDNA_TRACY_NOISE_WINDOW_SIZE` | Bases per QC window; valid range 10–20 |
| `noise_window_step` | `int` | `5` | `MTDNA_TRACY_NOISE_WINDOW_STEP` | Bases between successive windows; cannot exceed `noise_window_size` |
| `noise_min_valid_bases` | `int` | `10` | `MTDNA_TRACY_NOISE_MIN_VALID_BASES` | Valid peak measurements required to emit a window; cannot exceed `noise_window_size` |
| `noise_min_supporting_windows` | `int` | `2` | `MTDNA_TRACY_NOISE_MIN_SUPPORTING_WINDOWS` | Likely-noisy windows required to emit a mask range; minimum 2 |
| `noise_snr_threshold` | `float` | `3.0` | `MTDNA_TRACY_NOISE_SNR_THRESHOLD` | Threshold for median and per-base peak SNR checks |
| `noise_low_snr_threshold` | `float` | `2.0` | `MTDNA_TRACY_NOISE_LOW_SNR_THRESHOLD` | Low-end SNR required with sustained low-SNR evidence |
| `noise_purity_threshold` | `float` | `0.60` | `MTDNA_TRACY_NOISE_PURITY_THRESHOLD` | Minimum acceptable dominant-peak purity |
| `noise_background_threshold` | `float` | `0.35` | `MTDNA_TRACY_NOISE_BACKGROUND_THRESHOLD` | Maximum acceptable median-background-to-signal ratio |
| `noise_second_peak_threshold` | `float` | `0.50` | `MTDNA_TRACY_NOISE_SECOND_PEAK_THRESHOLD` | Maximum acceptable second-to-first peak ratio |
| `noise_quality_threshold` | `float` | `20.0` | `MTDNA_TRACY_NOISE_QUALITY_THRESHOLD` | Minimum acceptable `basecallQual` value when available |
| `noise_signal_threshold` | `float` | `75.0` | `MTDNA_TRACY_NOISE_SIGNAL_THRESHOLD` | Minimum acceptable dominant-peak signal |
| `noise_bad_base_fraction` | `float` | `0.40` | `MTDNA_TRACY_NOISE_BAD_BASE_FRACTION` | Fraction of bases required for sustained likely-noisy evidence |
| `noise_suspicious_base_fraction` | `float` | `0.20` | `MTDNA_TRACY_NOISE_SUSPICIOUS_BASE_FRACTION` | Fraction of bases required for suspicious evidence; cannot exceed `noise_bad_base_fraction` |

QC analyzes only valid base measurements. Missing individual peak measurements are skipped; a trace is `unavailable` if required trace fields are absent, no valid bases remain, or fewer than `noise_min_valid_bases` are available. A likely-noisy range requires the configured number of overlapping or directly adjacent likely-noisy windows. With masking enabled (the default), both variants and coverage intervals from that trace are subtracted at the resulting canonical reference range before forward/reverse traces are combined. See [Tracy Noise QC](tracy-noise-qc.md) for the metrics, classification, report schema, and masking rules.

The Tracy executable path is configured separately in `ToolsSettings.tracy`
(env `MTDNA_TOOLS_TRACY`, default `./tools/tracy`, validated as an executable
`FilePath`), accessed via `get_settings().tools.tracy`.

## CLI

`pipeline.main()` builds an `argparse` parser, constructs a `BatchOptions`,
calls `process_batch(...)`, and logs the processed-sample count. Tracy-specific
parameters (`trim`, `pratio`, `maxindel`, variant-quality thresholds,
`max_workers`, and all `noise_*` settings) are **not** CLI flags — they come
from `TracySettings` / the `MTDNA_TRACY_*` env vars.
`BatchOptions` is a dataclass with `ref_path: str = "ref/rCRS.fasta"`,
`batch_id: str | None = None`, and `samples_path: Path | None = None` (the TXT
path, required at runtime).

| Arg | Required | Default | Purpose |
|-----|----------|---------|---------|
| `-i` / `--input-dir` | yes | — | Directory containing AB1 trace files |
| `-o` / `--output-dir` | yes | — | Base output dir (`json/`, `regions/`, `preprocess/` auto-created) |
| `-s` / `--samples` | yes | — | Path to TXT file with sample IDs (one per line) |
| `--batch-id` | no | `None` | Batch identifier (passed to `Batch.write`; falls back to `"unknown"`) |
| `--ref-path` | no | `ref/rCRS.fasta` | Path to rCRS reference FASTA |

```bash
python src/tools/tracy/pipeline.py \
    -i data/raw/20241218_mtDNA_11 \
    -o results/automate_pipeline \
    -s data/raw/20241218_mtDNA_11/sample_ids.txt \
    --batch-id 20241218_mtDNA_11 \
    --ref-path ref/rCRS.fasta
```

## Public API

Entry points and key types, with one-line signatures and purposes. The pure
ETL and helper modules expose only what is listed here; everything else is an
internal helper.

- **`etl.process(sample_id, input_path, *, ref_seq, config=None) -> Sample`** — parse one Tracy decompose JSON into a `Sample` (pure ETL entry point).
- **`etl.merge_variants(variants: list[Variant]) -> list[Variant]`** — IUPAC consensus per `(pos, ref)` across multiple reads for one LID.
- **`etl._ProcessConfig`** — `NamedTuple` of `quality_threshold`, `min_peak_value`, `pratio`, `heteroplasmy_threshold`, `variant_regions`, `batch_id`, `sample`; the `config` argument type of `process()`.
- **`pipeline.process_batch(input_dir, output_dir, options=None) -> dict[str, Sample]`** — batch orchestration: decompose, QC, ETL, per-trace masking, trace merging, and output writing.
- **`pipeline.BatchOptions`** — dataclass (`ref_path`, `batch_id`, `samples_path`) for `process_batch()`.
- **`pipeline.main() -> None`** — CLI entry point.
- **`preprocessing.decompose_sample(sample, input_dir, output_dir, ref_path, *, config=None) -> tuple[Path, ...]`** — run the external `tracy decompose` on a sample's AB1 files; returns JSON files produced successfully in the current run.
- **`preprocessing.DecomposeConfig`** — `NamedTuple` of `trim`, `pratio`, `maxindel`; the `config` argument type of `decompose_sample()`.
- **`quality_control.NoiseConfig.from_settings() -> NoiseConfig`** — read all trace-QC and masking settings from `TracySettings`.
- **`quality_control.analyze_trace_data(data, *, sample_id, filename, config=None) -> TraceQCResult`** — calculate per-trace peak metrics, windows, and merged noisy ranges from parsed Tracy JSON.
- **`quality_control.analyze_trace_file(path, *, sample_id, config=None) -> TraceQCResult`** — load a trace JSON and calculate its QC result.
- **`quality_control.build_qc_report(results, config=None) -> dict[str, Any]`** — serialize deterministic per-sample trace-QC report data.
- **`noise_mask.apply_noise_mask(sample, ranges) -> NoiseMaskResult`** — exclude noisy-range variants, subtract the same spans from that trace's intervals, and recompute flags.
- **`transforms.apply_all_transformations(variants, conditions, primers, params) -> None`** — position-specific dispatch, AC-repeat insertion remap, HV2 polyC transforms, canonical polyC insertion creation (mutates in place).
- **`transforms.update_variant_positions(variants, position_params, start_index, end_index) -> dict[str, list[list[int]]]`** — strand-aware position remap and insertion decimal numbering; returns `file_intervals`.
- **`utils.PrimerTypeInfo`** — `NamedTuple` of primer boolean flags.
- **`utils.detect_primer_type(filename) -> PrimerTypeInfo`** — primer flags from the trace filename.
- **`utils.detect_heteroplasmy(peaks, position, heteroplasmy_threshold=None, sample=None) -> str`** — IUPAC code from peak ratios.
- **`utils.validate_peak_quality(variant, min_peak_value=None, pratio=None, sample=None) -> bool`** — peak-height and dominant-ratio quality check.
- **`utils.detect_variant_conditions(variants) -> dict`** — single-pass condition flags for the transforms.
- **`utils.variant_overlaps_range(position, ref, range_start, range_end) -> bool`** — canonical SNP/insertion/deletion overlap check.

Module communication is one-directional: `pipeline.py` coordinates
`preprocessing.py`, `quality_control.py`, `etl.py`, and `noise_mask.py`;
`noise_mask.py` consumes QC ranges and core models; `etl.py` imports from
`transforms.py`, `utils.py`, and `src/core`; QC never imports ETL or core models;
`utils.py` never depends on `Sample`/`Variant`; and the CLI wrapper holds no
business logic.

## Cross-References

- [tracy-noise-qc.md](tracy-noise-qc.md) — Trace-QC metrics, report schema, and masking rules

- [README.md](README.md) — Standard Tool Module Template
- [ARCHITECTURE.md](../ARCHITECTURE.md) — System architecture overview
- [variants.md](../variants.md) — Position normalization and interval utilities
- [config.md](../config.md) — Tracy settings and genomic regions