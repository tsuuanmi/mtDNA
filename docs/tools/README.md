# Standard Tool Module Template

> **Status**: Active standard for all mtDNA pipeline tool modules.
> **Applies to**: FIS, Sequencher, Mutation Surveyor, Tracy, BLASTn, and any future tool modules.
> **See also**: [ARCHITECTURE.md](../ARCHITECTURE.md) — §12 Standard Tool Module Template

## 1. Purpose

This document defines the standard folder structure, naming conventions, and
contracts for tool modules in the mtDNA raw data processing pipeline. Every
tool module must follow this template to ensure consistency, testability, and
maintainability.

## 2. Module Location

Each tool module lives in a single location:

    src/tools/<tool>/            # Implementation module (ETL, pipeline, utils, etc.)

The CLI entry point is the ``main()`` function in ``pipeline.py``, invoked
directly by shell scripts (e.g. ``python -m src.tools.<tool>.pipeline``).
All business logic lives inside ``src/tools/<tool>/``. The ``pipeline.py``
owns both orchestration and the CLI interface.

## 3. Mandatory Files

Every tool module **must** contain these files:

    src/tools/<tool>/
    +-- __init__.py       # empty (no re-exports; import directly from etl.py/pipeline.py)
    +-- etl.py            # Pure ETL: raw input -> Sample
    +-- pipeline.py       # Batch orchestration and parallel processing

### 3.1 __init__.py -- Package marker

- **Empty** in the current tools; it only marks the directory as a Python package.
- **Import public symbols directly from their modules**: `process` from `etl.py`,
  `process_batch` from `pipeline.py` (e.g. `from src.tools.<tool>.etl import process`).
- **Must not import** private (_-prefixed) functions.

### 3.2 etl.py -- Pure ETL

The ETL module is **pure**: it parses raw input and produces a Sample object.
It does not compute HV sequences or statistics -- those are the responsibility of
Batch.

| Contract | Detail |
|----------|--------|
| **Pure function** | No file I/O (beyond reading the input file), no network, no mutation of external state |
| **Single entry point** | process(sample_id, input_path, **kwargs) -> Sample |
| **Returns** | Sample with source_tool=<Tool.TOOL_NAME>, hv1=None, hv2=None, hv3=None, all stats fields None |
| **Variant flagging** | Must call src/core/flagging.flag_variants() for per-variant flags |
| **Sample flags** | May call src/core/flagging.SampleFlagger for sample-level flags |
| **Data model** | Must use src/core/models.Sample and src/core/models.Variant -- no tool-specific model classes |
| **Intervals** | Must compute and set Sample.intervals from tool-specific range logic |
| **No HV/stats** | Batch owns gen_seq_from_variants_regions() and statistic_variants() -- ETL must not compute these |
| **Error handling** | Raise on invalid input; log warnings for recoverable issues. Do not silently drop data. |

### 3.3 pipeline.py -- Batch Orchestration

The pipeline module orchestrates batch processing. It imports and calls other
modules (etl, preprocessing, utils) as needed.

| Contract | Detail |
|----------|--------|
| **Delegates** | Must call etl.process() for single-sample processing, not reimplement ETL |
| **Parallel** | Must use ProcessPoolExecutor for batch processing (or equivalent) |
| **Output** | Must call Batch(samples).write() to produce statistic_fullbatch.json |
| **Region JSONs** | Must write per-region intermediate JSONs via src/core/region.write_region_jsons() |
| **Entry point** | process_batch(input_dir, output_dir, **kwargs) -> Dict[str, Sample] |
| **Error resilience** | Log per-sample failures, continue batch. Do not abort entire batch on one failure. |
| **Sample combination** | Must handle same-LID multi-region Samples via _combine_same_tool_samples() or equivalent |
| **Orchestrator** | Imports and calls etl, preprocessing, utils -- never the reverse |

## 4. Optional Files

Optional files are added when a tool needs them. They follow standardized
naming conventions.

    src/tools/<tool>/
    +-- ... (mandatory files)
    +-- utils.py            # Tool-specific utilities
    +-- preprocessing.py    # Input format conversion
    +-- flagging.py         # Tool-specific quality control flags
    +-- <step_dirs>/        # Step subdirectories (tool-specific)

### 4.1 utils.py -- Tool-Specific Utilities

- **Created when needed** -- tools with no shared helper functions may omit this file.
- **Pure functions** -- no side effects, testable in isolation.
- **No model dependency** -- should not import Sample/Variant directly. Stays
  at the utility level (parsing, validation, format conversion).
- **Canonical source** -- each tool utils.py is the single source of truth
  for that tool filename parsing, range validation, and format-specific logic.

Examples: parse_filename_ranges() (Sequencher), safe_str() (MS),
position normalization helpers.

### 4.2 preprocessing.py -- Input Format Conversion

- **Purpose**: Converts tool-specific raw input into the standardized format
  that etl.py consumes.
- **Called by**: pipeline.py (before ETL), **never** by etl.py directly.
- **Idempotent**: Running twice on the same input must produce the same output.
- **Standardized name**: All input conversion modules use the name
  preprocessing.py, regardless of the specific conversion. The current
  Mutation Surveyor helper is still `txt_to_excel.py`; treat that as the
  existing preprocessing implementation until it is renamed at refactor time.
- **File discovery**: If the tool needs to locate input files by pattern, that
  logic belongs in pipeline.py or preprocessing.py, **not** in the CLI wrapper.

Examples:
- txt_to_excel.py -> preprocessing.py (MS: TXT custom report -> XLSX)
- retired root-level preprocessing helper -> src/tools/blastn/preprocessing.py (BLASTn: AB1->FASTQ->FASTA)
- Tracy decompose_sample() -> preprocessing.py (AB1->Tracy decompose output)

### 4.3 flagging.py -- Tool-Specific Quality Control Flags

- **Purpose**: QC flags specific to this tool that go beyond src/core/flagging.py.
- **Optional**: src/core/flagging.py provides shared flagging (heteroplasmy,
  insertion/deletion rules, etc.). Tool-specific flagging.py is only needed
  when a tool has flags not covered by core.
- **Compatibility**: Must produce flags compatible with Sample.sample_flags
  (list[str]) and Sample.variant_flags (dict[str, list[str]]).
- **No duplication**: Must not duplicate flags already defined in src/core/flagging.py.
- **Documentation**: Each flag must be documented with its meaning and trigger condition.

Example: Mutation Surveyor flagging.py defines AUTOPASS_HV1,
AUTOPASS_HV2_AND_HV3, REVIEW_HV1, REVIEW_HV2_AND_HV3, and region status
classification.

### 4.4 Step Subdirectories

Some tools have multi-step processing orchestrated by external scripts. These
step subdirectories are **tool-specific and not part of the standard template**.

Example: Mutation Surveyor 01_02/, 03_04/, 05_06/ directories contain
standalone scripts called by scripts/tools/mutation_surveyor.sh.

When a tool has step subdirectories, they live inside src/tools/<tool>/
and each step has its own README documenting its logic.

## 5. CLI Entry Point -- pipeline.py main()

The CLI entry point is the ``main()`` function in ``pipeline.py``. It handles
argument parsing and delegation to ``process_batch()``. Shell scripts invoke
it directly: ``python -m src.tools.<tool>.pipeline ...``.

| Contract | Detail |
|----------|--------|
| **main() in pipeline.py** | argparse -> process_batch(). No business logic beyond arg parsing. |
| **No file discovery** | File discovery belongs in preprocessing.py or the shell script. |
| **No preprocessing triggers** | Preprocessing is called by pipeline.py, not the CLI. |
| **No data processing** | All processing happens inside the module. |

## 6. Documentation -- docs/tools/<tool>.md

Every tool must have a documentation file in docs/tools/. The documentation
must include these **mandatory sections**:

1. **Purpose** -- What the tool processes, input/output formats, position in the pipeline.
2. **Public API** -- Function signatures with types and return values.
3. **Module Structure** -- File-by-file description of each module file.
4. **ETL Contract** -- What etl.py guarantees (pure, no HV/stats, Sample model).
5. **Pipeline Contract** -- What pipeline.py guarantees (parallel processing, Batch.write).
6. **CLI Entry Point** -- Command-line interface and arguments.
7. **Pipeline Position** -- Diagram showing upstream inputs and downstream consumers.

**Optional sections** (add when applicable):

8. **Preprocessing** -- Input format conversion details.
9. **Quality Control** -- Tool-specific flags and their meanings.
10. **Tool-Specific Details** -- Filename formats, range validation rules, configuration.

## 7. Pipeline Position

Every tool module sits in the same position in the overall pipeline:

    Tool-specific raw input (.TXT, .xlsx, .ab1, etc.)
           |
           v
    src/tools/<tool>/pipeline.py main() (CLI entry point)
           |
           v
    src/tools/<tool>/pipeline.py -> process_batch()
           |  may call preprocessing.py for format conversion
           |  calls etl.process() for each sample
           |  uses utils.py for parsing/validation
           |  uses core/flagging.py (and optionally <tool>/flagging.py)
           |
           v
    src/tools/<tool>/etl.py -> process()
           |  produces Sample (variants + intervals + flags)
           |  hv1/hv2/hv3 = None, stats = None
           |
           v
    Sample (src/core/models.py)
           |
           v
    src/core/batch.py -> Batch
      - gen_seq_from_variants_regions() -> hv1, hv2, hv3
      - statistic_variants() -> no_snps, no_ins, ...
      - write() -> batch JSON
           |
           v
    Downstream: comparison, FASTA gen, unification, etc.

## 8. Reference Implementations

### Sequencher (canonical example)

    src/tools/sequencher/
    +-- __init__.py     # empty (no re-exports)
    +-- etl.py          # TXT -> Sample (pure ETL)
    +-- pipeline.py     # Batch orchestration, ProcessPoolExecutor + CLI entry point
    +-- quality_control.py  # Pre-ETL QC gate: filename/range validation, variant-table-header check
    +-- utils.py        # Filename parsing, range validation, interval utilities
    docs/tools/sequencher.md  # Full documentation

### Mutation Surveyor (complex example)

    src/tools/mutation_surveyor/
    +-- __init__.py             # empty (no re-exports)
    +-- etl.py                  # XLSX -> Sample (pure ETL)
    +-- pipeline.py             # Batch orchestration + auto-pass filter + AB1 sort + CLI entry point
    +-- utils.py                # Token parsing, safe_str helpers
    +-- flagging.py             # MS-specific flags (AUTOPASS_*, REVIEW_*, region status)
    +-- txt_to_excel.py         # TXT -> XLSX conversion (current MS preprocessing helper)
    +-- sort_review_data.py     # AB1 file sorting by classification
    +-- merge_final_profiles_vs_truth.py  # Merge per-batch vs_truth workbooks into combined + filtered outputs
    +-- 01_02/                  # Step 1-2: trace processing + control QC
    +-- 03_04/                  # Step 3-4: HV region merging
    +-- 05_06/                  # Step 5-6: final profiles merge + truth comparison
    docs/tools/mutation_surveyor.md # Full documentation

### Tracy (conforms)

    src/tools/tracy/
    +-- __init__.py              # empty (no re-exports)
    +-- etl.py                   # Pure ETL: Tracy decompose JSON → Sample
    +-- pipeline.py              # Batch orchestration, parallel processing + CLI entry point
    +-- preprocessing.py         # AB1 → Tracy decompose (subprocess)
    +-- quality_control.py       # Per-trace noise metrics, ranges, and report records
    +-- noise_mask.py            # Per-trace noisy-range masking and flag recomputation
    +-- transforms.py            # Position-specific, polyC, and strand transforms
    +-- utils.py                 # Alignment, coordinate, primer, and peak helpers
    docs/tools/tracy.md          # Pipeline and API documentation
    docs/tools/tracy-noise-qc.md # QC metrics, reports, and masking behavior

### BLASTn (conforms)

    src/tools/blastn/
    +-- __init__.py              # empty (no re-exports)
    +-- etl.py                   # BLASTN output -> Sample (pure ETL)
    +-- pipeline.py              # Batch orchestration, parallel processing + CLI entry point
    +-- preprocessing.py         # AB1 -> FASTQ -> FASTA -> BLASTN alignment
    +-- utils.py                 # BLASTN column constants, strand helpers, TSV parsing
    +-- transform.py             # Indel right-norm, conversion rules, polyC handling
    +-- fasta_pipeline.py        # FASTA-input batch pipeline (pre-assembled FASTA -> BLASTN TSV)
    docs/tools/blastn.md         # Full documentation

## 9. Module Communication Rules

    # Direction of imports -- pipeline.py is the orchestrator:

    # pipeline.py imports from:
    from src.tools.<tool>.etl import process
    from src.tools.<tool>.utils import ...
    # pipeline.py may import from:
    from src.tools.<tool>.preprocessing import ...  # if preprocessing exists
    from src.tools.<tool>.flagging import ...        # if flagging exists

    # etl.py imports from:
    from src.core.flagging import flag_variants, SampleFlagger
    from src.core.models import Sample, Tool, Variant
    from src.tools.<tool>.utils import ...           # if utils exists

    # NEVER: etl.py imports from pipeline.py
    # NEVER: utils.py imports from Sample/Variant
    # NEVER: pipeline.py main() contains business logic beyond argparse + delegation

## 10. Checklist for New Tool Modules

When adding a new tool module, verify:

- [ ] src/tools/<tool>/__init__.py exists (empty; import process/process_batch directly from etl.py/pipeline.py)
- [ ] src/tools/<tool>/etl.py implements process(sample_id, input_path, ...) -> Sample
- [ ] src/tools/<tool>/pipeline.py implements process_batch(...) -> Dict[str, Sample]
- [ ] ETL returns Sample with hv1=None, hv2=None, hv3=None and all stats None
- [ ] ETL uses src/core/flagging.flag_variants() for per-variant flags
- [ ] ETL uses src/core/models.Sample and Variant -- no tool-specific model classes
- [ ] Pipeline calls Batch(samples).write() for output
- [ ] Pipeline writes per-region JSONs via write_region_jsons()
- [ ] Pipeline uses ProcessPoolExecutor for parallel processing
- [ ] src/tools/<tool>/pipeline.py has main() CLI entry point (argparse -> process_batch)
- [ ] docs/tools/<tool>.md follows the mandatory sections template
- [ ] Tool-specific preprocessing.py is called by pipeline.py, not by etl.py
- [ ] Optional files (utils.py, flagging.py) follow naming conventions
