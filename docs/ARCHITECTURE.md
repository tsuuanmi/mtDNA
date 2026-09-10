# mtDNA Pipeline — System Architecture

> **Scope**: Current architecture of the mtDNA raw data processing pipeline.
> This document replaces the former ROADMAP.md and reflects the implemented codebase as of 2026-07-02.

## 1. Architecture Overview

```
  Tracy (.ab1)   BLASTn (.ab1)   Sequencher (.TXT)   Mutation Surveyor (.xlsx/.tsv)
       │               │                │                      │
       └───────────────┴────────────────┴──────────────────────┘
                               │
               Tool Entry Point (src/tools/<tool>/pipeline.py main())
                               │  delegates to ETL module
                               ▼
  ┌──────────────────────────────────────────────────────────────┐
  │  ETL Module (src/tools/<tool>/etl.py)                       │
  │                                                              │
  │  1. Parse Raw Variants     ← tool-specific parsing           │
  │  2. Validate & Transform  ← ranges, positions, conditions   │
  │  3. Variant Flagging      ← flag_variants()                  │
  │                                                              │
  │  → produces Sample (variants + intervals, no HV sequences)   │
  └──────────────────────────────┬───────────────────────────────┘
                                 │
                                 ▼
                      ┌─────────────────────┐
                      │  Sample (canonical) │
                      │  (src/core/models)  │
                      │  hv1/hv2/hv3: None  │  ← populated by Batch
                      └──────────┬──────────┘
                                 │
                                 ▼
  ┌──────────────────────────────────────────────────────────────┐
  │  Batch (src/core/batch)                                      │
  │                                                              │
  │  1. generate_sequence() via sample_to_dict() → hv1, hv2, hv3 │
  │  2. statistic_variants()          → no_snps, no_ins, ...     │
  │  3. Write HV sequences back to Sample                        │
  │  4. Serialize Sample → JSON                                  │
  └──────────────────────────────┬───────────────────────────────┘
                                 │
                 ┌───────────────┼───────────────┐
                 │               │               │
                 ▼               ▼               ▼
           ┌──────────┐   ┌──────────┐   ┌──────────┐
           │Comparison│   │FASTA Gen │   │  Other   │
           │(pairwise)│   │          │   │ (future) │
           └──────────┘   └──────────┘   └──────────┘

  All downstream consumers accept List[Sample] as input.
  `generate_sequence()` runs inside `sample_to_dict()` (called by Batch).
```

### 1.1 Standard Tool Pipeline — Named Steps

ETL modules handle parsing and validation only. Steps 5–6 (sequence generation
and variant statistics) are shared computations owned by `Batch`. ETL produces
a `Sample` with variants and intervals; Batch computes HV sequences, writes them
back to `Sample.hv1`/`hv2`/`hv3`, then serializes.

| Step | Function | Input | Output | Module |
|------|----------|-------|--------|--------|
| 1. Parse Raw Variants | Tool-specific parsing of raw trace or text data | Raw file (.ab1, .TXT, .xlsx) | Variant dicts | `src/tools/<tool>/etl.py` |
| 2. Validate & Transform | Range validation, position normalization, condition detection | Variant dicts + intervals | Cleaned variants + sample_flags | `src/tools/<tool>/etl.py` |
| 3. Variant Flagging | `flag_variants()` | Variant objects | Per-variant flag reasons (`pos\|ref\|alt` → reasons) | `src/core/flagging.py` |
| 4. Build Sample | Construct `Sample()` model | Variants + intervals + flags | `Sample` (hv1/hv2/hv3 = None) | `src/core/models.py` |
| 5. Sequence Generation | `generate_sequence()` | Sample.variants + Sample.intervals | hv1/hv2/hv3 + hv_seq_refs | `src/core/sample.py` (called inside `sample_to_dict()` by Batch) |
| 6. Variant Statistics | `statistic_variants()` | Variants + HV seqs + refs | SNP/ins/del/identical/unread counts | `src/core/variants.py` (called inside `sample_to_dict()` by Batch) |
| 7. Batch Serialization | `sample_to_dict()` — writes hv1/hv2/hv3 back to Sample | `Sample` | batch JSON entry | `src/core/sample.py` |

**Tool-specific ETL details** (parsing, strand handling, polyC logic, etc.)
are documented in `docs/tools/`.

### 1.2 Downstream Consumption

Downstream tasks accept `List[Sample]` as input, loaded from `statistic_fullbatch.json`
via `load_sample_batch()`. They read `Sample.hv1`/`hv2`/`hv3` directly — they do
NOT call `generate_sequence()` directly (it runs inside `sample_to_dict()` invoked by Batch). They do not import tool-specific code.

#### Comparison (`src/core/comparison.py`)

```
batch_a.json + batch_b.json
       │              │
       ▼              ▼
  load_sample_batch()  load_sample_batch()
       │              │
       ▼              ▼
  Sample A            Sample B
       │              │
       └──────┬───────┘
              ▼
     SampleComparator.compare()
              │
              ▼
     ComparisonResult
     ├── concordance (Y/N)
     ├── unique_variants_a
     ├── unique_variants_b
     └── per-variant flags
              │
              ▼
     JSON + TSV + Excel output
```

#### FASTA Generation (`src/generate_fasta.py`)

```
statistic_fullbatch.json
              │
              ▼
  load_sample_batch() → List[Sample]
              │
              ▼
  FastaGenerator.generate_from_samples()
     │  For each sample:
     │    Read sample.hv1/hv2/hv3 (computed by Batch)
     │    Write .fasta file (HV1/HV2/HV3)
              │
              ▼
  FastaValidator.validate_all()
     │  Re-runs canonical code path and diffs
     │  alignment report
              ▼
     .fasta files + validation report
```

### Key Principles

1. **Tool modules live in `src/tools/<tool>/`** — Each tool has a `pipeline.py` with a `main()` CLI entry point that serves as the callable entry point. Shell scripts invoke it directly: `python -m src.tools.<tool>.pipeline`.
2. **ETL modules live in `src/tools/<tool_name>/`** — The ETL component that parses raw input, validates, and converts to a `Sample` object lives inside the tool directory. The pipeline module delegates to it.
3. **`Sample` is the canonical model** — Every ETL module produces a `Sample` object (defined in `src/core/models.py`). This is the single interchange format between ETL layers, comparison, and flagging.
4. **Downstream tasks are tool-agnostic** — Comparison, FASTA generation, reports, flagging, statistics, and concordance consume `Sample` objects or batch JSON produced by `Batch`. They do not import tool-specific code and live in `src/core/`.
5. **Modules are callable** — Each tool module exposes a **Python API** (importable functions) plus a **thin CLI wrapper**. Pipeline scripts call the Python API directly; the CLI provides shell convenience.
6. **Two-level flag system** — Sample-level flags (`Sample.sample_flags`) are tool-specific and set by ETL. Variant-level flags (`Sample.variant_flags`) are general and computed by `src/core/flagging.py` via `flag_variants()`.
7. **ETL contract** — Every tool ETL must call `flag_variants()` before constructing a `Sample`. `generate_sequence()` and `statistic_variants()` are called inside `sample_to_dict()` (invoked by Batch) only — ETL modules must NOT call them directly. `sample_to_dict()` writes HV sequences back to `Sample.hv1`/`hv2`/`hv3`.

---

## 2. Core Layer (`src/core/`)

The core layer contains tool-agnostic logic shared by all pipelines. Downstream modules import from `src/core/` and never from tool modules.

### 2.1 Models (`src/core/models.py`)

Pydantic models that define the canonical data contract.

| Model | Purpose |
|-------|---------|
| `Variant` | Single variant entry — atomic unit of comparison. Fields: `pos` (int/str), `ref`, `seq`, `quality` (List[float]), `files` (List[str]), `peaks` |
| `Sample` | Standardized result that every tool must produce. The single interchange format. |
| `ToolResult` | Per-tool data extracted from comparison A/B pairs. |
| `ComparisonResult` | Result of comparing two Sample objects. |
| `Tool` | Enum: `sequencher`, `mutation_surveyor`, `tracy`, `blastn`, `fasta`, `fis`, `sanger`, `unified`, `unknown` |

**Sample key fields:**

| Field | Type | Description |
|-------|------|-------------|
| `sample_id` | `str` | Sample identifier (LID from lab) |
| `variants` | `List[Variant]` | Detected variants |
| `source_tool` | `Tool` | Which tool produced this sample; `unified` for merged results |
| `sample_flags` | `List[str]` | Sample-level flags (range problems, QC) |
| `variant_flags` | `Dict[str, List[str]]` | Per-variant flags keyed by "pos\|ref\|alt" |
| `intervals` | `Optional[Dict]` | Per-HV-region intervals |
| `hv1/hv2/hv3` | `Optional[str]` | HV consensus sequences |
| `no_snps/no_ins/no_dels/no_identicals/no_unread` | `Optional[int]` | Variant statistics |
| `peaks` | `Optional[Dict]` | Per-position peak data from traces |
| `information` | `Optional[Dict]` | Sample metadata. Always includes `source_tool` and `batch_id` at serialization. `region_sources` (per-region tool provenance) included for unified/merged results. Tool-specific fields (well, type_sample, etc.) populated by ETL. |
| `batch_id` | `Optional[str]` | Batch identifier |

**Constants:**

| Constant | Value | Purpose |
|----------|-------|---------|
| `PEAK_BASE_ORDER` | `["A", "C", "G", "T"]` | Peak data index convention |
| `SAMPLE_FLAG_RANGE_POSSIBLY_WRONG` | `"Range - Possibly wrong"` | Sequencher range validation failed |
| `SAMPLE_FLAG_RANGE_MISSING` | `"Range - Missing"` | Valid ranges but missing HV region(s) |

### 2.2 Comparison (`src/core/comparison.py`)

Tool-agnostic pairwise comparison of `Sample` objects. Produces `ComparisonResult` with binary concordance (Y/N/N/A) and per-variant flagging.

**Key class:** `SampleComparator`

| Method | Purpose |
|--------|---------|
| `compare(result_a, result_b)` | Full pairwise comparison → ComparisonResult |
| `is_concordant(result_a, result_b)` | Binary concordance check |
| `find_unique_variants(a, b, ignore_special)` | Variants in a not present in b |
| `format_intervals(intervals)` | Format intervals dict to display string |
| `compare_batches(samples_a, samples_b)` | Batch-level comparison with summary/variant rows |

**Output formats:** JSON, summary TSV, variant-level TSV via `write_comparison_json()`, `write_comparison_tsv()`, `write_variant_level_tsv()`.

**CLI:** `python -m src.core.comparison -a batch_a.json -b batch_b.json -o output_base [--batch-id ID]`

### 2.3 Flagging (`src/core/flagging.py`)

Two-level flag system:

- **`flag_variants(variants: List[Variant]) → Dict[str, List[str]]`** — Per-variant flags keyed by "pos|ref|alt". General to all tools.
- **`flag_sample(sample_flags: List[str]) → List[str]`** — Passthrough for sample-level flags (currently; will merge tool-specific + general in future).

**VariantFlag** dataclass: `position`, `ref`, `alt`, `flags: List[str]`, `key() → "pos|ref|alt"`.

**VariantAnalyzer** class detects: special positions, consecutive indels, lowercase variants, heteroplasmy outside polyC, combined conditions (459 del + 16192 issue).

**SampleFlagger** assigns severity levels (Critical/High/Medium/Low) and pattern-based flags (309.1 insertion, 310 variant, 460 variant, polyC region, nomenclature issues, etc.).

### 2.3.1 Directional polyC Policy (`src/core/polyc.py`)

`polyc.py` centralizes tool-agnostic directional policy: HV2 forward exclusion
at `>=304` and reverse exclusion at `<=315`, plus 16189-C-triggered HV1
forward/reverse suppression. Its helpers own classification, exclusion
decisions, and reasons; tools only map native primer metadata to a direction
and apply results to already-called candidates.

### 2.4 Sample I/O (`src/core/sample.py`) and Batch (`src/core/batch.py`)

Sample serialization/deserialization lives in `src/core/sample.py`:

| Function | Purpose |
|----------|---------|
| `sample_to_dict(sample, ref_path)` | Compute HV + stats, write back to Sample, serialize to canonical dict |
| `load_sample_batch(path)` | Load a batch JSON file into a list of Sample objects |
| `_normalize_keyed_sample(raw)` | Convert grouped variant format to flat list for Sample construction |
| `_parse_variant_dict(v, sample_id)` | Normalize a single variant dict for Sample model validation |

Batch (`src/core/batch.py`) delegates per-sample serialization to `sample_to_dict()` and owns batch-level concerns only: multi-sample aggregation and file writing.

**Key class:** `Batch`

| Method | Purpose |
|--------|---------|
| `to_json()` | Convert all samples to batch dict (delegates to `sample_to_dict`) |
| `write(path, batch_id)` | Write batch JSON + per-sample JSONs to file |

**JSON schema:**

```json
{
  "SAMPLE_ID": {
    "no_snps": 0, "no_ins": 0, "no_dels": 0,
    "no_identicals": 0, "no_unread": 0,
    "HV1": "<consensus>", "HV2": "<consensus>", "HV3": "<consensus>",
    "variants": {
      "snps": [{"pos": ..., "ref": ..., "seq": ..., "file": [...], "quality": [...]}],
      "insertions": [...], "deletions": [...]
    },
    "intervals": {"HV1": [[16024, 16365]], "HV2": [[73, 340]], "HV3": [[438, 576]]},
    "sample_flags": ["<flag reason>", ...],
    "variant_flags": {"pos|ref|alt": ["<flag reason>", ...]},
    "peaks": {"pos": [[A, C, G, T], ...]},
    "information": {"source_tool": "...", "batch_id": "...", ...}
  }
}
```

### 2.5 Merger (`src/core/mtdna_merger.py`)

Single-pass mtDNA profile merger. `mtDNAMerger` merges multiple per-sample JSON files with strict overlap consistency checks. Uses Sample objects internally and Batch.sample_to_dict() for output serialization. Enforces: variants inside declared intervals, ref allele matches rCRS, overlap regions have identical alleles. Returns MergeResult(sample=Sample, warnings=list) instead of raw dict. Also provides merge_samples() for merging Sample objects directly.

See [`docs/core/mtdna_merger.md`](core/mtdna_merger.md) for full documentation of the merger implementation, including validation rules, merge flow, error handling, and CLI usage.

### 2.6 Uploader (`src/modules/sheets/uploader.py`)

Uploads Sequencher comparison results (Excel/TSV) to Google Sheets via `gspread`. Not part of the core analysis pipeline — used for reporting/sharing results.

---

## 3. Tool Modules (`src/tools/<tool>/`)

### 3.1 Sequencher (`src/tools/sequencher/pipeline.py`)

Entry point for Sequencher TXT processing. Delegates to `src/tools/sequencher/etl.py`.

**API:**
- `process(sample_id, input_path, batch_id=None)` → `Sample`
- `process_batch(input_dir, output_dir, ..., json_dir=None)` → `Dict[str, Sample]`

**CLI:** `python -m src.tools.sequencher.pipeline --input-dir input_dir --output-dir output_dir [--ref-path] [--batch-id]`

### 3.2 Tracy (`src/tools/tracy/`)

Direct variant calling from AB1 files using the Tracy tool. Decomposes AB1 files, calls variants with quality filtering and heteroplasmy detection. Follows the [Standard Tool Module Template](tools/README.md).

**Module Structure:**

| File | Purpose |
|------|---------|
| `src/tools/tracy/etl.py` | Pure ETL: Tracy decompose JSON → Sample |
| `src/tools/tracy/pipeline.py` | Batch orchestration, parallel processing via ProcessPoolExecutor |
| `src/tools/tracy/preprocessing.py` | AB1 → Tracy decompose (subprocess) |
| `src/tools/tracy/quality_control.py` | Per-trace peak-noise metrics, window classification, and QC report records |
| `src/tools/tracy/noise_mask.py` | Per-trace likely-noisy-range masking and flag recomputation |
| `src/tools/tracy/transforms.py` | Position-specific, polyC, and strand variant transformations |
| `src/tools/tracy/utils.py` | Alignment, coordinate, primer, and peak helpers |
| `src/tools/tracy/__init__.py` | Empty; import `process` from etl.py, `process_batch` from pipeline.py |

**Public API:**

- `process(sample_id, input_path, *, ref_seq, ...) → Sample` — Pure ETL for a single Tracy decompose JSON file
- `process_batch(input_dir, output_dir, options) → dict[str, Sample]` — Batch orchestration with parallel processing

**Configuration:** Tracy parameters loaded from `get_settings().tracy` (Pydantic settings with `MTDNA_TRACY_` env prefix). Per-trace noise QC evaluates overlapping peak-metric windows after decomposition. By default, variants and the same coverage spans are removed from the originating trace before forward/reverse traces are combined; a clean mate trace can retain coverage for the final LID profile. `MTDNA_TRACY_NOISE_MASK_ENABLED=false` retains calls and intervals while reports remain enabled. See [Tracy Noise QC](tools/tracy-noise-qc.md).

### 3.3 BLASTn (`src/tools/blastn/pipeline.py`)

BLASTN-based pipeline: AB1 → FASTQ → FASTA → BLASTN alignment → variant calling → consensus.

**Key entry points:** `process_batch(...)` (batch orchestration with `ProcessPoolExecutor`), `main()` CLI, configured via `BatchOptions`.

Pipeline: `preprocessing.py` (AB1 → BLASTN TSV) → `etl.py` (TSV → Sample) → `Batch.write()` + per-region JSONs. A parallel `fasta_pipeline.py` handles pre-assembled FASTA input. See [tools/blastn.md](tools/blastn.md).

---

## 4. ETL & Feature Modules

### 4.1 Sequencher ETL (`src/tools/sequencher/`)

| File | Purpose |
|------|---------|
| `etl.py` | Main ETL: TXT → Sample. Pure ETL: parses filename ranges, validates, flags variants. HV sequences and statistics computed by Batch. |
| `utils.py` | Shared utilities: `parse_filename_ranges()`, `validate_analysis_ranges()`, `parse_manual_txt()`, `safe_extract_sample_id()`, range/interval helpers. |

**ETL flow:** Parse TXT → Extract variants → Validate ranges → Compute flags → Build `Sample` (HV sequences and statistics are computed by Batch).

### 4.2 LAB (`src/modules/LAB/`)

| File | Purpose |
|------|---------|
| `manual.py` | Manual pipeline processing for LAB data |
| `generate_mtdna_html_report.py` | HTML report generation |

### 4.3 NGS (`src/modules/NGS/`)

| File | Purpose |
|------|---------|
| `src/tools/fis/` | Reusable FIS workbook ETL and canonical Sample JSON output; see [tools/fis.md](tools/fis.md) |
| `prepare_fis_comparison.py` | Map and normalize FIS samples for core comparison |
| `prepare_sanger_comparison.py` | Select and normalize mapped Sanger samples |
| `enrich_fis_sanger_comparison.py` | Add FIS metadata to shared comparison results |
| `compare_fis_fasta_direct.py` | Compare FIS FASTA results directly |
| `compare_fis_sanger_direct.py` | Legacy standalone FIS vs Sanger comparison |
| `FIS_merge_fasta.py` | Merge FIS FASTA files |
| `FIS_merge_variants.py` | Merge FIS variant files |
| `merge_comparison.py` | Merge comparison results |
| `merge_json_files.py` | Merge JSON files |
| `calculate_variant_freq.py` | Calculate variant frequencies |
| `convert_variants_to_tsv.py` | Convert variants to TSV format |
| `variant_frequency_analyzer.py` | Analyze variant frequencies |

### 4.4 TNLS (`src/modules/TNLS/`)

| File | Purpose |
|------|---------|
| `extract.py` | Data extraction |
| `merge_family_info.py` | Merge family information |
| `merge_str_mtdna_ids.py` | Merge STR and mtDNA IDs |
| `map_cccd.py` | Map CCCD identifiers |
| `fix_barcode_merge.py` | Fix barcode merge issues |
| `filter_family_variants.py` | Filter family variants |
| `kinship_analysis.py` | Kinship analysis |
| `upload_mapping_to_gsheets.py` | Upload mapping data to Google Sheets |
| `verification/` | TNLS verification subpackage |

### 4.5 Statistics (`src/modules/statistics/`)

| File | Purpose |
|------|---------|
| `query_multi_filters.py` | Multi-filter statistical queries |
| `merge.py` | Merge data and per-batch results |
| `blind_copy.py` | Blind-copy AB1/FASTA files and generate blinded metrics JSON |
| `add_variants_to_samples.py` | Add variant data and analyzed intervals to the sample mapping TSV |
| `add_variants_to_metadata.py` | Add variant data to metadata |
| `data/modules/statistics/sample_mapping.tsv` | Local ignored blind-id to sample-id mapping input |

---

## 5. Pipeline Orchestrators

The former standalone Python orchestrators for automate, manual, and comparison stages have been removed. End-to-end orchestration is now driven by `scripts/pipeline.sh`, which calls the tool modules and core steps directly. See [pipeline.md](pipeline.md) for the full step-by-step breakdown.

### 5.1 Shell Pipeline (`scripts/pipeline.sh`)

6-step pipeline script:

| Step | Name | Module | Description |
|------|------|--------|-------------|
| 1 | validate | `src/modules/quality_control/pc_ntc.py` | Validate Positive Control / NTC samples |
| 2 | prepare | (shell + pandas) | Organize AB1 files, extract sample list from metadata Excel |
| 3 | automate | `src/tools/blastn/pipeline.py` or `src/tools/tracy/pipeline.py` | Run analysis engine (`-p`/`--pipeline` mode) |
| 4 | manual | `src/tools/sequencher/pipeline.py` + `src/core/comparison.py` | Sequencher ETL and pairwise comparison vs automated results |
| 5 | fasta | `src/generate_fasta.py` | Generate FASTA files |
| 6 | report | `src/generate_reports.py` | Generate PDF reports |

Supports two modes: **tracy** (default) and **blastn**.

---

## 6. Supporting Modules

### 6.4 Variant Utilities (`src/core/variants.py`)

Centralized position and variant domain functions. All variant/position operations should import from here.

**Key functions:**

| Function | Purpose |
|----------|---------|
| `normalize_position(pos)` | Normalize to int (bases) or str (insertions like "309.1") |
| `pos_base(pos)` | Extract integer base position |
| `pos_sort_key(pos)` | Sort key for position ordering |
| `statistic_variants(...)` | Compute variant counts (SNPs, ins, del, identical, unread) |
| `split_variant_types(variants)` | Categorize into {snps, insertions, deletions} |
| `variant_to_dict(v)` | Convert Variant to canonical {pos, ref, seq, file, quality} dict |
| `get_hv_region_for_position(pos)` | Determine which HV region a position belongs to |
| `is_position_in_intervals(pos, intervals)` | Check if position falls within given intervals |
| `get_overlap(a, b)` | Compute overlap between interval lists |

> Sequence generation (`generate_sequence()`) lives in `src/core/sample.py` and is called inside `sample_to_dict()` by Batch. The older `gen_seq_from_variants_regions()` has been removed.

**Position contract:** `int` for base positions, `str` for insertions. `normalize_position()` enforces this.

### 6.5 Configuration (`src/config.py`)

Pydantic Settings-based configuration with `.env` support.

| Settings class | Purpose |
|---------------|---------|
| `Settings` | Main settings (project root, nested settings) |
| `ToolsSettings` | External tool paths (tracy, blastn, seqtk, ugene) |
| `DirectoryPathsSettings` | Project directory paths |
| `GenomicRegionsSettings` | HV1/HV2/HV3 regions, polyC regions, QC regions |
| `ParametersSettings` | Default thresholds and file extensions |

**Singleton access:** `get_settings()` (cached).


### 6.6 Report Generation (`src/generate_reports.py`)

`ReportGenerator` — Extracts sample info from Excel/TSV, manages auth tokens, generates PDF reports via API.

### 6.7 FASTA Generation (`src/generate_fasta.py`)

`FastaGenerator` — Generates FASTA files from per-sample JSON (loaded as `Sample`) and regenerates HV consensus sequences via `generate_sequence()`. Supports `--force` for clean-slate regeneration. Validates output via `FastaValidator`.

---

## 7. Validation (`src/validation/`)

### 7.1 FASTA Validation (`src/validation/validation_fasta.py`)

`FastaValidator` — Re-runs the canonical `generate_sequence()` code path (via the standard `Sample`) and compares output against generated FASTA. Produces position-level alignment diffs.

### 7.2 Regenerate Results Validation (`src/validation/validation_regenerate_results.py`)

`RegenerateResultsValidator` — Validates JSON files from the manual regeneration step. Checks: batch JSON structure, per-sample JSON existence, HV sequence integrity, variant consistency, interval correctness. Outputs comprehensive JSON report.

---


## 8. Data Flow

The standard data flow is documented in §1 (Architecture Overview) and §1.1 (Standard Tool Pipeline — Named Steps). Tool-specific data flows are documented in `docs/tools/`.

### 8.1 Comparison Pipeline (Core)

```
batch_a.json (Sample[]) + batch_b.json (Sample[])
    ↓ load_sample_batch()
    ↓ SampleComparator.compare_batches()
    ↓
ComparisonResult → JSON + summary TSV + variant-level TSV + Excel
```

### 8.2 Batch Serialization

```
List[Sample]
    ↓ Batch(samples).to_json()
    ↓
statistic_fullbatch.json (v4 schema)
```

See §1.2 for detailed consumption diagrams.

---

## 9. Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Canonical model | `Sample` (Pydantic) | Single interchange format; validated, serializable, typed |
| Position format | `int` for bases, `str` for insertions | Matches `normalize_position()` convention; preserves precision |
| Two-level flags | `sample_flags` (tool-specific) + `variant_flags` (general) | Clean separation of concerns; ETL owns sample-level, core owns variant-level |
| Variant flags key | `"pos\|ref\|alt"` string | Stable identity for variant matching across tools |
| Batch JSON | Grouped variants + flags + peaks + information | Historical format with modern extensions |
| Native format handling | Each module parses its own format | JSON, TSV, Excel, TXT — no universal parser |
| Callable interface | Python API + thin CLI wrapper | Testability (API), composability (import), convenience (CLI) |
| ETL granularity | One `etl.py` per module | Start simple; split Extract/Transform/Load only when complexity demands |
| Downstream imports | Import from `src/core/`, never from tool modules | Enforces tool-agnostic downstream |
| Tool-specific logic | Stays in module directory | Each module owns its ETL and flagging; comparison lives in `src/core/` |
| Config | Pydantic Settings + `.env` | Type-safe, validated, environment-configurable |
| Logging | Loguru throughout | Consistent, structured logging |

---

## 10. Adding a New Tool

Adding a new tool requires:

1. Add `main()` CLI entry point to `src/tools/<new_tool>/pipeline.py`
2. Create `src/tools/<new_tool>/` with `etl.py` and `__init__.py`
3. Ensure ETL follows the standard pipeline contract (see §1.1):
   - Parse raw variants
   - Validate & transform
   - Call `flag_variants()` for per-variant flags
   - Construct a `Sample` object with variants, intervals, and flags
   - **Do NOT** call `generate_sequence()` or `statistic_variants()` — Batch owns those (via `sample_to_dict()`)
4. Ensure ETL output is a `Sample` object (with `source_tool`, `batch_id`, `intervals`, `variants`)
5. Ensure batch aggregation via `sample_to_dict()` produces correct output
6. Add to `Tool` enum in `src/core/models.py`
7. No downstream code changes required — they consume `Sample` objects or batch JSON

---


## 10. Source Module Index

The docs are organized to match the current high-level `src/` tree:

| Source Module | Doc | Purpose |
|---|---|---|
| `src/config.py` | [config.md](config.md) | Centralized Pydantic Settings configuration |
| `src/fasta_to_sample.py` | [tools/fasta.md](tools/fasta.md) | FASTA-to-Sample conversion helper |
| `src/generate_fasta.py` | [generation.md](generation.md), [tools/fasta.md](tools/fasta.md) | FASTA file generation |
| `src/generate_reports.py` | [generation.md](generation.md) | Report generation and PDF download |
| `src/core/batch.py` | [core/batch.md](core/batch.md) | Batch processing and HV sequence generation |
| `src/core/comparison.py` | [core/comparison.md](core/comparison.md) | Pairwise sample comparison |
| `src/core/flagging.py` | [core/flagging.md](core/flagging.md) | Sample and variant flagging system |
| `src/core/models.py` | [core/models.md](core/models.md) | Pydantic data models (Position, Variant, Sample, etc.) |
| `src/core/mtdna_merger.py` | [core/mtdna_merger.md](core/mtdna_merger.md) | Multi-tool mtDNA profile merger |
| `src/core/paths.py` | [core/paths.md](core/paths.md) | Centralized output path helpers |
| `src/core/polyc.py` | [core/polyc.md](core/polyc.md) | Shared directional HV1/HV2 polyC policy |
| `src/core/regenerate.py` | [core/regenerate.md](core/regenerate.md) | Regenerate gate and final JSON filtering |
| `src/core/region.py` | [core/region.md](core/region.md) | Per-region JSON I/O and range/interval QC |
| `src/core/sample.py` | [core/sample.md](core/sample.md) | Sample serialization and `generate_sequence()` |
| `src/core/unify.py` | [core/unify.md](core/unify.md) | Batch unification helpers |
| `src/core/variants.py` | [core/variants.md](core/variants.md) | Variant and position domain utilities |
| `src/tools/blastn/` | [tools/blastn.md](tools/blastn.md) | BLASTN pipeline (AB1 → BLASTN TSV → Sample) |
| `src/tools/tracy/` | [tools/tracy.md](tools/tracy.md) | Tracy variant caller (standard module) |
| `src/tools/sequencher/` | [tools/sequencher.md](tools/sequencher.md) | Sequencher ETL pipeline |
| `src/tools/mutation_surveyor/` | [tools/mutation_surveyor.md](tools/mutation_surveyor.md) | Mutation Surveyor pipeline and step scripts |
| `src/modules/LAB/` | [modules/LAB/generate_mtdna_html_report.md](modules/LAB/generate_mtdna_html_report.md) | LAB report/manual workflow helpers |
| `src/modules/NGS/` | [modules/NGS/README.md](modules/NGS/README.md) | NGS comparison, merging, and variant-frequency utilities |
| `src/modules/TNLS/` | [modules/TNLS/README.md](modules/TNLS/README.md) | TNLS metadata, kinship, and verification utilities |
| `src/modules/quality_control/pc_ntc.py` | — | Positive Control / NTC validation |
| `src/modules/sheets/merger.py` | [modules/sheets/merger.md](modules/sheets/merger.md) | Excel file merger for comparison results |
| `src/modules/sheets/uploader.py` | [modules/sheets/uploader.md](modules/sheets/uploader.md) | Google Sheets upload |
| `src/modules/statistics/` | [modules/statistics/README.md](modules/statistics/README.md) | Multi-filter statistics queries |
| `src/validation/validation_fasta.py` | [validation/validation_fasta.md](validation/validation_fasta.md) | FASTA validation |
| `src/validation/validation_regenerate_results.py` | [validation/validation_regenerate_results.md](validation/validation_regenerate_results.md) | Result regeneration validation |
| `scripts/pipeline.sh` | [pipeline.md](pipeline.md) | End-to-end pipeline orchestration |

## 11. Future Work

### Regenerate Gate (Partially Implemented)

> **See also**: [`docs/tools/mutation_surveyor.md`](tools/mutation_surveyor.md) — MS AutoRun + MS1 target workflow extending the Regenerate Gate with region-based merge, tool priority rules, auto-pass criteria, review packages, and TXT-to-Excel preprocessing.

`mtDNAMerger` implements Step 1 (strict validation) and produces unified output with `source_tool = "unified"` and `region_sources` provenance tracking. Step 2 implements a **Sequencher-wins priority rule**: when both Sequencher and Mutation Surveyor cover the same region, Sequencher's variants are used (manual review is more reliable than auto-pass). Configurable YAML-based priority (ROADMAP §5) remains deferred.

```
┌─────────────────────────────────────────────────────────────────┐
│  REGENERATE GATE                                                │
│                                                                 │
│  Step 1: Strict validation (mtDNAMerger) ✅ IMPLEMENTED          │
│          - All variants inside declared intervals               │
│          - Ref allele matches rCRS                              │
│          - Overlapping regions have identical alleles           │
│          - source_tool = "unified" for merged results           │
│          - region_sources = {"HV1": "...", "HV2": "...", "HV3": "..."} in information │
│                                                                 │
│  Step 2: Sequencher-wins priority ✅ IMPLEMENTED                │
│          - When both tools cover same region, Sequencher wins   │
│          - MS variants for that region are discarded            │
│          - region_sources shows only the winning tool           │
│          - Tool priority: sequencher > mutation_surveyor >     │
│            tracy > blastn                                       │
│                                                                 │
│  Step 3: Configurable YAML priority ❌ DEFERRED                 │
│          (deferred — planned YAML config schema)                │
│                                                                 │
│  → Final Unified JSON (same schema, richer data)                │
│    source_tool = "unified", region_sources = {"HV1": "...", "HV2": "...", "HV3": "..."} in information  │
└─────────────────────────────────────────────────────────────────┘
```

The Sequencher-wins priority rule applies because MS auto-pass is not 100% reliable (can miss variants), while Sequencher manual review is more trustworthy. When both tools cover the same region, Sequencher's variants replace MS's entirely (no union or supplementing). Allele conflicts at the same position from different tools still raise `MergeValidationError`.

### Module Protocol Extraction

When 3+ tool modules share a clear pattern, extract a `ToolModule` Protocol from the convention. Tracy, BLASTn, Sequencher, and Mutation Surveyor now all follow the standard tool module template (see §12), so this extraction is feasible when a new tool is added.

### FASTA Generation Alignment

`generate_fasta.py` loads each per-sample JSON as a `Sample` (via `src.core.sample.load_sample`) and regenerates HV consensus sequences via `generate_sequence()`. Full alignment to the `Sample` model is complete.

### Sequence Generation Ownership

`gen_seq_from_variants_regions()` (the older tuple-returning generator) has been removed. HV consensus generation is now owned by `generate_sequence()` in `src/core/sample.py`, called inside `sample_to_dict()` by Batch.

| Module | Status | Notes |
|--------|--------|-------|
| `src/core/sample.py` | ✅ Canonical | `generate_sequence()` — the structured `Sequence` generator |
| `src/core/batch.py` | ✅ Canonical | Delegates to `sample_to_dict()` for serialization; HV computation is inside `sample_to_dict()` |
| `src/core/mtdna_merger.py` | ✅ Done | Uses `Sample` objects + `sample_to_dict()` for output |
| `src/generate_fasta.py` | ✅ Done | Reads `Sample.hv1`/`hv2`/`hv3` via `generate_sequence()` |
| `src/validation/validation_fasta.py` | ✅ Done | Replays `generate_sequence()` for validation

## 12. Standard Tool Module Template

> **Full specification**: [tools/README.md](tools/README.md)

All tool modules follow the standard template defined in docs/tools/README.md.
Key points:

| Aspect | Standard |
|--------|----------|
| **Module location** | src/tools/<tool>/ (implementation + CLI entry point via pipeline.py main()) |
| **Mandatory files** | __init__.py, etl.py, pipeline.py |
| **Optional files** | utils.py, preprocessing.py, flagging.py, step subdirectories |
| **CLI convention** | Thin wrapper only: argparse to pipeline.process_batch() |
| **ETL contract** | Pure function, returns Sample with hv*=None, stats=None |
| **Pipeline contract** | Orchestrator, delegates to etl, uses ProcessPoolExecutor, calls Batch.write() |
| **Preprocessing** | preprocessing.py called by pipeline.py, never by etl.py |
| **Documentation** | docs/tools/<tool>.md with mandatory sections (Purpose, API, Structure, ETL Contract, Pipeline Contract, CLI, Pipeline Position) |

Current conformance:

| Tool | Status | Notes |
|------|--------|-------|
| Sequencher | Conforms | Canonical example of the standard template |
| Mutation Surveyor | Partially conforms | Has `__init__.py`, full module suite; `txt_to_excel.py` as preprocessing; CLI has file discovery logic |
| Tracy | Conforms | Follows standard template with etl.py, pipeline.py, preprocessing.py, utils.py |
| BLASTn | Conforms | Standard module at src/tools/blastn/ (etl, pipeline, preprocessing, utils, transform, fasta_pipeline) |

