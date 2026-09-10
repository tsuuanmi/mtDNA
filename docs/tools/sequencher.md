# Sequencher Tool

> Sequencher TXT → Sample ETL pipeline for the mtDNA raw data processing system.

## Overview

The Sequencher module imports variant tables exported by the Sequencher
software (plain-text `*.TXT` format) and converts each export into a
standardized `Sample` object compatible with the downstream comparison,
flagging, and batch-aggregation layers. A batch run discovers every `*.TXT`
in an input directory, runs a pre-ETL quality-control gate on each file,
parses the survivors in parallel, and writes both per-region intermediate
JSONs and a flat batch-level `json/` tree plus a QC report.

Following the standard tool-pipeline contract from
[`ARCHITECTURE.md`](../ARCHITECTURE.md), the Sequencher ETL is **pure**: it
populates only `variants`, `intervals`, and flags. The HV sequences
(`hv1`/`hv2`/`hv3`) and per-sample statistics (`no_snps`, `no_ins`,
`no_dels`, `no_identicals`, `no_unread`) are left `None` on the `Sample` and
are computed later by `Batch` at aggregation time. Variants are emitted with
`source_tool = Tool.SEQUENCHER`.

## Module structure

```
src/tools/sequencher/
├── __init__.py          # empty (no re-exports; public symbols imported from their modules)
├── etl.py               # TXT -> Sample ETL (pure, no HV/stats computation)
├── pipeline.py          # Batch orchestration, pre-ETL QC gate, parallel processing, CLI
├── quality_control.py   # Filename parsing, range validation, and the pre-ETL QC gate
└── utils.py             # Region/interval helpers and region-group detection
```

`etl.py` owns the single-file ETL. `pipeline.py` owns batch orchestration and
the CLI. `quality_control.py` is the single source of truth for filename
parsing, analysis-range validation, and the pre-ETL pass/fail gate. `utils.py`
holds region/interval helpers plus region-group detection (`NUM_REGIONS = 3`).
Private helpers (leading `_`) are used internally and are not re-exported.

## Quality control and the filename standard

Quality control is the defining feature of this tool. `quality_control.run_qc`
is a **pre-ETL hard gate**: a file that fails it never enters `etl.process`,
never produces a region JSON, and never appears in
`statistic_fullbatch.json`. A failure does not abort the batch — the remaining
files still process normally. The gate runs in `process_batch` on every
discovered TXT *before* any ETL work, and a full report is written to
`preprocess/qc_report.json` immediately after discovery (before processing) so
it always reflects every discovered file, pass or fail.

### The filename standard

There is exactly one correct standard format for a Sequencher TXT filename.
Every other format fails QC and produces no JSON.

- **Full coverage** — `LID.TXT` (e.g. `2513827.TXT`). The file is treated as
  covering all three HV regions.
- **Explicit partial coverage** —
  `LID-a1-b1 a2-b2 ... an-bn.TXT` (e.g.
  `2513827-73-340 438-576 16024-16365.TXT`). Each `start-end` pair names a
  sequenced interval; pairs are separated by single spaces and the LID is
  joined to the first pair with a single hyphen `-`.
- **Optional region-group prefix** — `HV1_` or `HV2-3_` may precede the LID
  (e.g. `HV1_2513827-16024-16365.TXT`, `HV2-3_2513827-73-340 438-576.TXT`).
  The prefix declares which region group the file covers and is enforced by
  rule 6 below.

The standard format is matched by `VALID_INTERVAL_PATTERN`
(`^\d+-\d+( \d+-\d+)*$`) for the range part, with the LID-to-range separator
fixed to `-`.

### What makes a filename non-standard (hard QC fail)

Any deviation from the above is unsafe and hard-fails the gate:

- A **non-standard LID-to-range separator** (space, dot, underscore, comma,
  ...) — e.g. `2510768-73-340_438-573_16024-16365.TXT` or
  `...16024.16365.TXT`. A malformed separator can silently drop a region's
  range, so it is never emitted. The parser still extracts whatever ranges it
  can find (so the flag is informative), but `format_error` is set to
  `SAMPLE_FLAG_RANGE_POSSIBLY_WRONG` and QC fails.
- A **non-space separator between range pairs**.
- **No LID extractable** from the filename (the reason string is
  `"No sample ID (LID) could be extracted from filename"`).
- A **missing Sequencher variant-table header** in the TXT content — the file
  must contain the column header `Pos`, `Seq`, `Con`, and `Required Edit`
  (matched by `SEQUENCHER_TABLE_HEADER_PATTERN`). A file without it fails with
  reason `SEQUENCHER_TABLE_HEADER_ERROR` =
  `"Invalid Sequencher TXT format: missing variant table header"`. Without
  this header the row parser would treat prose lines as data and emit a
  misleading variant set.

### Range validation rules

When ranges are present, `validate_analysis_ranges` checks the parsed pairs
through six rules. The rule numbers below match the source docstring and the
`_check_*` helper names in `quality_control.py`; rules 1-4 are structural and
checked first (call order: direction → position bounds → sorted order →
region consistency), then rule 6 (prefix), then rule 5 (completeness) last.

| # | Rule | Condition checked | Failure flag | Gates JSON? |
|---|------|-------------------|--------------|-------------|
| 1 | Sorted order | `a1 < b1 < a2 < b2 < ... < an < bn` (i.e. `ranges[i][1] < ranges[i+1][0]`) | `Range - Possibly wrong` | Yes |
| 2 | Position bounds | For every pair, at least one of `start`/`end` falls inside HV2 `[73,340]`, HV3 `[438,576]`, or HV1 `[16024,16365]` | `Range - Possibly wrong` | Yes |
| 3 | Region consistency | Each pair's `start` and `end` map to the **same** HV region | `Range - Possibly wrong` | Yes |
| 4 | Range direction | For every pair `a-b`, `a < b` (`start >= end` fails) | `Range - Possibly wrong` | Yes |
| 5 | Completeness | With a prefix: all regions expected by that prefix are present. Without a prefix: all of HV1/HV2/HV3 are present | `Range - Missing` | **No** (soft) |
| 6 | Prefix consistency | When `region_prefix` is set, every found region ⊆ the prefix's expected regions (`HV1`→{HV1}, `HV2-3`→{HV2,HV3}) | `Range - Possibly wrong` | Yes (only when a prefix is set) |

The distinction between hard and soft is the key behavior. **Hard** rules
(1-4 and 6) and any non-standard format return
`(should_display=False, SAMPLE_FLAG_RANGE_POSSIBLY_WRONG)` and block JSON
output entirely. The **soft** rule 5 (completeness) is the only one that
returns `(True, SAMPLE_FLAG_RANGE_MISSING)`: a correctly-formatted filename
that is simply missing a region still passes QC and produces JSON, with
`Range - Missing` carried onto the `Sample` as a sample-level flag by
`etl.process`. `"FULL REGION"` and empty ranges short-circuit to
`(True, None)`.

### Flag constants

The three quality-control constants, all defined in `quality_control.py`:

- `SAMPLE_FLAG_RANGE_POSSIBLY_WRONG` = `"Range - Possibly wrong"` — the
  hard-fail flag/reason for rules 1-4 and 6 and for non-standard separators.
- `SAMPLE_FLAG_RANGE_MISSING` = `"Range - Missing"` — the soft sample flag
  for incomplete HV coverage (rule 5); does not gate JSON.
- `SEQUENCHER_TABLE_HEADER_ERROR` =
  `"Invalid Sequencher TXT format: missing variant table header"` — the
  hard-fail reason when the TXT lacks the variant-table header.

Internal regexes: `REGION_PREFIX_PATTERN = ^(HV1|HV2-3)_`,
`SEQUENCHER_TABLE_HEADER_PATTERN =
\bPos\b\s+\bSeq\b\s+\bCon\b\s+Required\s+Edit\b`, and `LID_RANGE_PAIR_RE`
(a lookbehind-anchored overlapping scan that robustly separates the LID from
the ranges even when the separator is non-standard, without mistaking
trailing LID digits for a range start).

## ETL processing

`etl.process(sample_id, input_path, batch_id=None) -> Sample` is the single-file
ETL entry point. It runs only for files that have already passed the QC gate.
The stages, in order:

1. **Parse filename ranges** — `quality_control.parse_filename_ranges(input_path.name)`
   → `FilenameParseResult` (`sample_id`, `ranges`, `format_error`,
   `region_prefix`).
2. **Resolve sequencing regions** — `_resolve_seq_ranges(ranges, filename, region_prefix)`
   → `(seq_regions, range_flags)`. `"FULL REGION"` or an empty list selects
   `_default_seq_regions()` (full HV1/HV2/HV3 from settings, no flags).
   Otherwise `validate_analysis_ranges` is called: if
   `should_display=False`, a warning is logged, the flag is appended, and the
   defaults are used as a fallback; if `should_display=True`, the parsed
   ranges are used directly and any returned flag (e.g. `Range - Missing`) is
   appended.
3. **Carry format errors** — if `parse_result.has_format_error` and the
   `format_error` string is not already in `range_flags`, append it. (QC gates
   format errors before ETL, so this is defensive.)
4. **Parse TXT variants** — `_parse_sequencher_txt(input_path)` locates the
   variant-table header line (the first line containing `Pos`, `Seq`, **and**
   `Con`), treats data as starting at `i + 2` (skipping the header and a blank
   line), and splits each non-empty line on whitespace. Lines with fewer than
   three parts are skipped. `pos_str = parts[0].replace(",", "")` is normalized
   via `src.core.variants.normalize_position` (a `None` result skips the line);
   `ref_base`/`alt_base` are `parts[1]`/`parts[2]` with `":"` mapped to `"-"`.
   Each entry is `{pos, ref, seq, file: [], quality: []}`.
5. **Build intervals** — `utils.get_intervals_from_manual(seq_regions)` maps
   each HV region name (`HV1`, `HV2`, `HV3`) to its overlap with the manual
   ranges (via `src.core.variants.get_overlap` against
   `settings.regions.REGIONS`); each span is converted to `list[int]`.
6. **Build variant list** — each parsed dict becomes a `Variant` with `pos`,
   `ref`, `seq`, `quality` (coerced to `list[int]`), `files` (default `[]`),
   `peaks` (default `None`). **All** parsed variants are kept, including
   positions outside the analyzed intervals, so the comparison layer sees the
   full picture.
7. **Variant flags** — `src.core.flagging.flag_variants(variant_list)` → the
   per-variant flag dict.
8. **Sample-level analyze flags** — the sequenced regions are filtered to
   `_expected_regions_for_prefix(region_prefix or "")` and
   `SampleFlagger(variant_list, intervals=flagging_intervals or None).analyze()`
   is run; the analyze reasons (e.g. `"No 315.1 variant"`,
   `"Consecutive indels (5+)"`) are appended to the sample flag list.
9. **Deduplicate** — `deduplicate_sample_flags(sample_flags_list, variant_flags_dict)`
   drops any sample-level flag that duplicates a per-variant flag string.
10. **Construct `Sample`** — pure ETL output: `variants`, `intervals`,
    `sample_flags`, `variant_flags`, `source_tool=Tool.SEQUENCHER`, `batch_id`,
    and `hv1`/`hv2`/`hv3`/`no_snps`/`no_ins`/`no_dels`/`no_identicals`/`no_unread`
    and `information` all `None`.

## Batch processing

`pipeline.process_batch(input_dir, output_dir, options=None) -> dict[str, Sample]`
orchestrates a full batch. It returns a dict mapping LID → combined `Sample`
(empty if no tasks passed QC or none processed successfully). The flow:

1. **Defaults & paths** — `options = options or BatchOptions()`;
   `paths = OutputPaths(output_dir)`; `json/` and `preprocess/` are created
   explicitly, and `regions/` is created **only when** `options.json_dir is
   None` (the override case leaves it uncreated).
2. **Workers** — `max_workers = settings.sequencher.max_workers or
   max(1, os.cpu_count() or 1)`.
3. **Discover & QC** — `txt_files = sorted(input_dir.glob("*.TXT"))`. Each
   file is run through `run_qc(name, path)` → `QCResult`; a `qc_results`
   entry `{sample_id, filename, status("pass"|"fail"), reason}` is recorded
   and failing files are skipped (logged as a warning) and never reach
   `process()`.
4. **Write QC report** — immediately after discovery (before processing),
   `preprocess/qc_report.json` is written with
   `{tool:"sequencher", total, passed, failed, samples:[...]}` so it always
   reflects every discovered TXT.
5. **Empty guard** — if no tasks passed QC, logs a warning and returns `{}`.
6. **Parallel ETL** — a `ProcessPoolExecutor(max_workers)` submits
   `process(sample_id, input_path, options.batch_id)` per passing task;
   results are collected via `as_completed`. Exceptions
   `(ValueError, KeyError, OSError, RuntimeError)` are logged and the sample
   is skipped (does not abort the batch).
7. **Per-region JSONs** — the output directory is `options.json_dir` if set,
   otherwise `paths.regions_dir`. For each processed item,
   `_write_region_jsons_for_sample` dispatches by
   `utils.detect_region_group_from_filename(filename)`: a detected group
   (`HV1` or `HV2-3`) writes a single region group; an indeterminate (full
   region) filename writes both `HV1` and `HV2-3`. Both paths call
   `src.core.region.write_region_jsons` and run
   `validate_region_group_intervals` (warnings logged, not gating).
8. **Combine same-LID** — `_build_combined_samples` groups processed items by
   `sample_id`; singletons pass through; multiple region Samples for one LID
   are merged by `_combine_same_tool_samples` (concat variants, merge
   intervals, union flags, dedup variants by position, re-run
   `SampleFlagger.analyze` + `flag_variants` + `deduplicate_sample_flags` on
   the merged set).
9. **Batch write** —
   `Batch(list(combined.values())).write(str(paths.json_dir), options.batch_id or "unknown", nest_batch_id=False)`,
   producing flat output (the base dir is already batch-specific).
10. **Logging** — logs the unique-LID count vs TXT count and warns if any TXT
    failed processing.

`pipeline.main()` is the CLI entry point (see [CLI](#cli)).

## Output layout

`process_batch` builds `OutputPaths(output_dir)` from `src/core/paths.py`:

```
{output_dir}/
├── json/                 # Batch.write() output (nest_batch_id=False)
│   ├── statistic_fullbatch.json
│   └── {LID}/{LID}.json
├── regions/              # per-region intermediate JSONs (default; --json-dir overrides)
│   ├── HV1/HV1_{LID}.json
│   └── HV2-3/HV2-3_{LID}.json
└── preprocess/
    └── qc_report.json    # pre-ETL QC report (pass + fail)
```

`OutputPaths` exposes `json_dir = base/"json"`,
`regions_dir = base/"regions"`, `preprocess_dir = base/"preprocess"`; its
`create_dirs()` creates all three, but Sequencher creates `json/` and
`preprocess/` explicitly and `regions/` only when `--json-dir` is not set.

### Region JSON pattern

Per-region intermediate JSONs follow the patterns in `src/core/region.py`:

- `REGION_FILE_PATTERNS = {"HV1": "HV1_{lid}.json", "HV2-3": "HV2-3_{lid}.json"}`
- `REGION_TO_KEYS = {"HV1": ["HV1"], "HV2-3": ["HV2", "HV3"]}`

`write_region_jsons(sample, output_dir, ref_path, region_groups)` writes, for
each group with variants/intervals,
`output_dir/{group}/{group}_{LID}.json` (filtering the sample to the group's
region keys via `filter_sample_by_regions`) and injects
`information.region_sources = dict.fromkeys(region_keys, sample.source_tool.value)`
for provenance. It returns `dict[group_name, Path]` of written files. Each
region JSON therefore carries `region_sources` provenance in `information`.

### QC report

`qc_report.json` records every discovered TXT (pass and fail) with a reason,
regardless of any later per-sample processing failure:

```json
{
  "tool": "sequencher",
  "total": 3,
  "passed": 1,
  "failed": 2,
  "samples": [
    {"sample_id": "2513827", "filename": "2513827.TXT", "status": "pass", "reason": null},
    {"sample_id": "", "filename": "73-340.TXT", "status": "fail", "reason": "No sample ID (LID) could be extracted from filename"}
  ]
}
```

## Configuration

`SequencherSettings` is defined in `src/config.py` as a Pydantic `BaseSettings`
with `env_prefix="MTDNA_SEQUENCHER_"` and `extra="forbid"`, accessed via
`get_settings().sequencher`. Because the only configurable field is an
environment-variable-backed worker count, settings are presented as a single
table:

| Field | Type | Default | Env var | Description |
|-------|------|---------|---------|-------------|
| `max_workers` | `int \| None` | `None` | `MTDNA_SEQUENCHER_MAX_WORKERS` | Max parallel `ProcessPoolExecutor` workers. `None` → `max(1, os.cpu_count() or 1)`. |

The HV region bounds used by validation and the interval helpers come from
`get_settings().regions.REGIONS` (`GenomicRegionsSettings`, env_prefix
`MTDNA_REGIONS_`): `HV1=[16024, 16365]`, `HV2=[73, 340]`, `HV3=[438, 576]`.

## CLI

`src/tools/sequencher/pipeline.py` exposes `main()`:

```bash
python src/tools/sequencher/pipeline.py \
  --input-dir <input_dir> \
  --output-dir <output_dir> \
  [--json-dir <dir>] \
  [--ref-path ref/rCRS.fasta] \
  [--batch-id <id>]
```

| Arg | Type | Required | Default | Purpose |
|-----|------|----------|---------|---------|
| `--input-dir` | `Path` | yes | — | Directory containing Sequencher `*.TXT` files. |
| `--output-dir` | `Path` | yes | — | Base directory for batch-level JSON output (`json/`, `regions/`, `preprocess/` created under it). |
| `--ref-path` | `str` | no | `ref/rCRS.fasta` | Path to rCRS reference FASTA (passed to `write_region_jsons`). |
| `--batch-id` | `str` | no | `None` | Batch identifier (falls back to `"unknown"` in `Batch.write`). |
| `--json-dir` | `Path` | no | `None` | Override directory for per-region intermediate JSONs (otherwise `{output-dir}/regions`). |

Concrete example:

```bash
python src/tools/sequencher/pipeline.py \
  --input-dir data/sequencher/20260615_batch \
  --output-dir results/tools/sequencher/20260615_batch \
  --ref-path ref/rCRS.fasta \
  --batch-id 20260615_batch
```

## Public API

Public symbols are imported directly from their modules (`process` from
`etl.py`, `process_batch` from `pipeline.py`); `__init__.py` is empty. The
key entry points:

- `etl.process(sample_id, input_path, batch_id=None) -> Sample` — single-file
  ETL: parse filename ranges, validate, parse TXT, build intervals/variants,
  flag, construct a pure `Sample` (HV/stats `None`).
- `pipeline.process_batch(input_dir, output_dir, options=None) -> dict[str, Sample]`
  — batch orchestration: discover, QC-gate, parallel ETL, per-region JSONs,
  combine same-LID, batch write.
- `pipeline.BatchOptions` — dataclass
  (`ref_path: str = "ref/rCRS.fasta"`, `batch_id: str | None = None`,
  `json_dir: Path | None = None`) controlling batch behavior.
- `pipeline.main() -> None` — CLI entry point (builds argparse, constructs
  `BatchOptions`, calls `process_batch`).
- `quality_control.run_qc(filename, file_path=None) -> QCResult` — the
  pre-ETL hard gate; returns `QCResult(sample_id, passed, reason)`.
- `quality_control.QCResult` — `NamedTuple` (`sample_id: str`, `passed: bool`,
  `reason: str | None`); `sample_id` is `""` when no LID could be extracted.
- `quality_control.validate_analysis_ranges(ranges, region_prefix=None, filename=None) -> tuple[bool, str | None]`
  — the 6-rule range validator; returns `(should_display, flag)`.
- `quality_control.parse_filename_ranges(filename) -> FilenameParseResult` —
  extract LID, ranges, format error, and region prefix from a filename.
- `quality_control.FilenameParseResult` — `NamedTuple`
  (`sample_id`, `ranges: str | list[list[int]]`, `format_error: str | None`,
  `region_prefix: str | None = None`) with `is_full_region` and
  `has_format_error` properties.
- `quality_control.extract_region_prefix(filename) -> tuple[str | None, str]`
  — strip an `HV1_` / `HV2-3_` prefix; returns `(prefix, remainder)`.
- `quality_control.has_sequencher_variant_table_header(file_path) -> bool` —
  detect the required `Pos`/`Seq`/`Con`/`Required Edit` header (UTF-8 with
  latin-1 fallback; returns `False` on `OSError`).
- `utils.get_intervals_from_manual(seq_regions) -> dict` — map manual ranges
  to per-HV-region intervals via `get_overlap` against `REGIONS`.
- `utils.is_full_region(intervals) -> str` — return `"FULL REGION"` when all
  three HV regions match the canonical boundaries exactly (sorted), else a
  space-separated `start-end` string (`"None"` for empty input).
- `utils.parse_intervals_string(intervals_string) -> list[list[int]]` — parse
  comma- or space-separated interval strings (old and new formats); empty or
  `"FULL REGION"` → `[]`; unparseable parts skipped with a warning.
- `utils.detect_region_group_from_filename(filename) -> str | None` — return
  `"HV1"`, `"HV2-3"`, or `None` (full region / indeterminate). Decision order:
  region prefix → `None` when full-region or format-error → otherwise inspect
  range boundaries against canonical HV bounds.

`_parse_sequencher_txt(file_path) -> dict[str, dict[str, Any]]` and
`_default_seq_regions() -> list[list[int]]` are private `etl.py` helpers
described in the [ETL processing](#etl-processing) stages above.

## Cross-References

- [../ARCHITECTURE.md](../ARCHITECTURE.md) — tool pipeline contract and the pure-ETL rule
- [../core/batch.md](../core/batch.md) — `Batch` aggregation, HV/stats computation, and `write()`
- [../core/region.md](../core/region.md) — per-region JSON I/O and naming conventions
- [../core/paths.md](../core/paths.md) — `OutputPaths` layout
- [unify.md](unify.md) — consumes Sequencher per-region JSONs for batch unification