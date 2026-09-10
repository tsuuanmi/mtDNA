# Core Data Models

> Canonical Pydantic models for the mtDNA pipeline: Position, Variant, Sample, ToolResult, ComparisonResult, Sequence, and Tool.

## Purpose

`src/core/models.py` defines the single source of truth for all data structures shared across tools. Every tool (automate, Sequencher, Mutation Surveyor, Tracy, BLASTN) must convert its output to these models so the comparison and flagging layers are unified.

## Pipeline Position

Models are the interchange format at every pipeline boundary:

```
ETL Module (tool-specific) → Sample → Batch → Comparison / FASTA / Reports
                                    ↑
                              Tool enum
```

## Public API

### `Position`

```python
Position = int | str
```

Canonical genomic position type:
- Base positions: `int` (e.g., `73`, `309`, `16569`)
- Insertion positions: `str` (e.g., `"309.1"`, `"217.10"`)

Every module must import `Position` from `src.core.models` rather than re-defining the alias locally. Use `validate_position()` at external boundaries to normalize and enforce the type.

### `validate_position(pos: float | str) → Position`

Normalise and validate a position value into canonical Position type. Enforces the Position contract at external boundaries (Excel reads, API inputs, etc.). Rejects `float` — callers must convert explicitly (`int(x)` for base positions, or `str` for insertions) so that precision loss is a conscious decision.

| Input | Result |
|-------|--------|
| `int` | Returned as-is |
| `str` like `"309"` | Converted to `int(309)` |
| `str` like `"309.1"` | Returned as `"309.1"` |
| `float` | **Rejected** — raises `TypeError` |

### `Tool`

```python
class Tool(StrEnum):
    SEQUENCHER = "sequencher"
    MUTATION_SURVEYOR = "mutation_surveyor"
    TRACY = "tracy"
    BLASTN = "blastn"
    FASTA = "fasta"
    FIS = "fis"
    SANGER = "sanger"
    UNIFIED = "unified"
    UNKNOWN = "unknown"
```

Enumeration of supported analysis tools. `FASTA` marks FASTA-derived samples. `FIS` and `SANGER` identify the two canonical batches prepared by the staged NGS comparison workflow. `UNIFIED` marks merged results from `MtDnaMerger`/`merge_by_regions()` that combine source regions; the contributing tool for each region is tracked in `information.region_sources`. `UNKNOWN` is a sentinel for a batch/sample that carries no tool attribution (e.g. a batch JSON without `information.source_tool`); loaders (`load_sample_batch`, `read_region_json`) set it instead of silently defaulting to a concrete tool, and comparison output displays it as `Unknown`.

### `Variant`

```python
class Variant(BaseModel):
    pos: int | str        # 1-based genomic position (int for bases, str for insertions)
    ref: str                    # Reference base ("-" for insertions)
    seq: str                    # Alternate sequence base ("-" for deletions)
    quality: list[int] = []     # Per-trace Phred-like quality scores
    files: list[str] = []       # Source trace filenames
    peaks: list[list[int | None]] | None = None  # Per-variant peak data
```

Single variant entry — the atomic unit of comparison. All fields use snake_case. `pos` must be pre-normalized before construction.

### `Sample`

```python
class Sample(BaseModel):
    sample_id: str                            # Sample identifier (LID from lab)
    variants: list[Variant] = []              # List of detected variants
    source_tool: Tool                   # Primary analysis tool
    sample_flags: list[str] = []              # Sample-level flag reasons
    hv1: str | None = None                # HV1 consensus sequence (16024-16365)
    hv2: str | None = None                # HV2 consensus sequence (73-340)
    hv3: str | None = None                # HV3 consensus sequence (438-576)
    intervals: dict[str, list[list[int]]] | None = None  # Per-HV-region intervals
    no_snps: int | None = None            # Number of SNP variants
    no_ins: int | None = None              # Number of insertion variants
    no_dels: int | None = None            # Number of deletion variants
    no_identicals: int | None = None      # Number of positions identical to reference
    no_unread: int | None = None          # Number of unread/N positions
    batch_id: str | None = None           # Batch identifier
    variant_flags: dict[str, list[str]] = {} # Per-variant flag reasons keyed by "pos|ref|alt"
    information: dict[str, Any] | None = None  # Sample metadata dict
```

Standardized result that every tool must produce. This is the single interchange format between ETL layers, comparison, and flagging.

### `ToolResult`

```python
class ToolResult(BaseModel):
    tool: Tool                           # Source tool
    analyzed_range: str = ""                   # Analyzed range string
    variants: list[dict[str, Any]] = []       # Variant dicts in canonical format
    flagged: bool = False                      # Whether this tool's result is flagged
    flag_reasons: list[str] = []              # Flag reasons for this tool
    variant_flags: dict[str, list[str]] = {}  # Per-variant flag reasons
```

Per-tool data extracted from ComparisonResult A/B pairs.

### `ComparisonResult`

```python
class ComparisonResult(BaseModel):
    sample_id: str                              # Sample identifier
    tool_a: ToolResult                         # Tool A result data
    tool_b: ToolResult                         # Tool B result data
    unique_variants_a: list[dict[str, Any]] = []  # Variants unique to tool A
    unique_variants_b: list[dict[str, Any]] = []  # Variants unique to tool B
    concordance: str = "N"                     # Y/N concordance between the two results
```

Result of comparing two Sample objects.

### `Sequence`

```python
class Sequence(BaseModel):
    reference_seq: str = ""                    # Full reference sequence string
    consensus_seq: str = ""                    # Full consensus sequence string
    hv_seq_refs: dict[str, Any] = {}           # Per-region reference sequences
    hv_seqs: dict[str, Any] = {}               # Per-region consensus sequences
    ref_dict: dict[int | str, str] = {}  # Position-keyed reference dict
    con_dict: dict[int | str, str] = {}  # Position-keyed consensus dict
```

Structured result from `generate_sequence()`. Replaces the old tuple return, making each field self-documenting.

## Key Concepts

### Position Contract

The mtDNA pipeline uses a strict position type system:
1. All external inputs must be normalized via `validate_position()` or `normalize_position()`
2. `float` positions are **rejected** at validation boundaries to prevent precision loss
3. `pos_base()` extracts the integer base for comparisons
4. `pos_sort_key()` provides deterministic ordering

### Tool Source Requirement

Every `Sample` must set `source_tool` to a valid `Tool` enum value. For merged results from `MtDnaMerger`/`merge_by_regions()`, `source_tool` is `Tool.UNIFIED` and the contributing source for each region is listed in `information.region_sources` (e.g., `{"HV1": "mutation_surveyor", "HV2": "sequencher", "HV3": "sequencher"}`). When a batch/region JSON carries no `source_tool`, loaders set `Tool.UNKNOWN` (displayed as `Unknown`) rather than defaulting to a concrete tool. This enables:
- Tool-specific comparison logic
- Traceability in comparison results
- Correct flagging behavior

### Variant Flags Key Format

`variant_flags` uses the key format `"pos|ref|alt"` (e.g., `"73|A|G"`) to uniquely identify each variant. Only variants that have at least one flag are included.

## Cross-References

- [core/flagging.md](flagging.md) — Flag building and severity classification
- [core/batch.md](batch.md) — Batch serialization uses `Sample` → JSON
- [core/sample.md](sample.md) — `sample_to_dict()` and `load_sample_batch()`
- [variants.md](variants.md) — position helpers and `validate_position()` re-export
