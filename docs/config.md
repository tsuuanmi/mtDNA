# Configuration

> Centralized pipeline configuration using Pydantic Settings with environment variable support.

## Purpose

`src/config.py` provides a single source of truth for all pipeline configuration — external tool paths, directory locations, genomic region definitions, and analysis parameters. It uses Pydantic v2 `BaseSettings` for type-safe, validated configuration that can be overridden via `.env` files or environment variables.

## Pipeline Position

Configuration is the **first module loaded** by every other module in the pipeline. All paths, regions, and parameters flow from `get_settings()`.

```
.env files / environment variables
        │
        ▼
   Settings (Pydantic)
        │
        ▼
   All pipeline modules
```

## Public API

### `get_settings() → Settings`

Cached settings instance (recommended access pattern). Returns a fully-validated `Settings` object. Never create `Settings()` directly — always use `get_settings()`.

### `Settings`

Main configuration class composing all sub-settings.

| Field | Type | Description |
|-------|------|-------------|
| `project_root` | `Path \| None` | Auto-detected project root (defaults to `src/` parent) |
| `tools` | `ToolsSettings` | External tool executable paths |
| `directories` | `DirectoryPathsSettings` | Input/output directory paths |
| `regions` | `GenomicRegionsSettings` | Genomic region definitions |
| `parameters` | `ParametersSettings` | Analysis parameters |

### `ToolsSettings`

Validated paths to external tools. All paths are `FilePath` — validation fails if the file doesn't exist.

| Field | Default | Env Prefix |
|-------|---------|------------|
| `ugene` | `./tools/ugene-50.0/ugene` | `MTDNA_TOOLS_UGENE` |
| `blastn` | `./tools/ugene-50.0/tools/blast/blastn` | `MTDNA_TOOLS_BLASTN` |
| `seqtk` | `./tools/seqtk/seqtk` | `MTDNA_TOOLS_SEQTK` |
| `tracy` | `./tools/tracy` | `MTDNA_TOOLS_TRACY` |

### `DirectoryPathsSettings`

Project directory paths. Input directories use `DirectoryPath` (must exist). Output directories are auto-created.

| Field | Default | Type | Description |
|-------|---------|------|-------------|
| `filter_rules` | `rules/filter` | `DirectoryPath` | Filter rule files |
| `convert_rules` | `rules/convert` | `DirectoryPath` | Convert rule files |
| `data` | `data` | `Path` | Input data directory |
| `ref` | `ref` | `DirectoryPath` | Reference sequences |
| `tools` | `tools` | `DirectoryPath` | Bundled tools directory |
| `results` | `results` | `Path` | Analysis results (auto-created) |
| `logs` | `logs` | `Path` | Log files (auto-created) |
| `credentials` | `credentials` | `Path` | Credentials (auto-created) |

### `GenomicRegionsSettings`

Hypervariable region coordinates and polyC handling.

| Field | Description |
|-------|-------------|
| `REGIONS` | `{HV1: [16024, 16365], HV2: [73, 340], HV3: [438, 576]}` |
| `POLYC_REGIONS` | `{HV2: [303, 315]}` — known polyC stretch positions |
| `POLYC_QC_REGIONS` | List of tuples for polyC QC checking |
| `SKIP_AFTER_POLYC` | `{HV1: [16193, 16195]}` — positions skipped after polyC |

### `ParametersSettings`

Default analysis parameters and file extensions.

| Field | Default | Description |
|-------|---------|-------------|
| `DEFAULT_TRIM_THRESHOLD` | `0.01` | Quality trim threshold for seqtk |
| `DEFAULT_WORD_SIZE` | `22` | BLASTN word size |
| `DEFAULT_NUM_THREADS` | `os.cpu_count()` | Thread count |
| `AB1_EXTENSION` | `.ab1` | AB1 file extension |
| `FASTQ_EXTENSION` | `.fastq` | FASTQ file extension |
| `FASTA_EXTENSION` | `.fasta` | FASTA file extension |
| `BLASTN_EXTENSION` | `.blastn` | BLASTN output extension |
| `REPORT_BASE_URL` | `https://genestory.ai/...` | PDF report API endpoint |
| `REPORT_RETRIES` | `5` | API retry count |
| `REPORT_RETRY_DELAY` | `5` | Seconds between retries |



### `TracySettings`

Tracy tool runtime parameters with environment variable override support.

| Field | Default | Env Prefix | Description |
|-------|---------|------------|-------------|
| `trim` | `7` | `MTDNA_TRACY_` | Tracy decompose trim parameter |
| `pratio` | `0.3` | `MTDNA_TRACY_` | Tracy decompose pratio parameter |
| `maxindel` | `1000` | `MTDNA_TRACY_` | Tracy decompose maxindel parameter |
| `quality_threshold` | `35` | `MTDNA_TRACY_` | Minimum quality score for variant calling |
| `min_peak_value` | `125` | `MTDNA_TRACY_` | Minimum peak height for variant calling |
| `heteroplasmy_threshold` | `0.8` | `MTDNA_TRACY_` | Minimum heteroplasmy ratio for detection |

### `Settings` (updated)

The `Settings` class now includes a `tracy` field:

| Field | Type | Description |
|-------|------|-------------|
| `project_root` | `Path | None` | Auto-detected project root |
| `tools` | `ToolsSettings` | External tool executable paths |
| `directories` | `DirectoryPathsSettings` | Input/output directory paths |
| `regions` | `GenomicRegionsSettings` | Genomic region definitions |
| `parameters` | `ParametersSettings` | Analysis parameters |
| `tracy` | `TracySettings` | Tracy tool runtime parameters |

Environment variable override for Tracy settings uses the `MTDNA_TRACY_` prefix (single underscore, matching the `env_prefix` in `TracySettings`):

```bash
MTDNA_TRACY_TRIM=7
MTDNA_TRACY_PRATIO=0.3
MTDNA_TRACY_MAXINDEL=1000
MTDNA_TRACY_QUALITY_THRESHOLD=35
MTDNA_TRACY_MIN_PEAK_VALUE=125
MTDNA_TRACY_HETEROPLASMY_THRESHOLD=0.8
```

## Key Concepts

### Path Resolution

Relative paths are resolved against `project_root` (auto-detected as `src/` parent). This ensures consistent behavior regardless of working directory.

### Auto-Created Directories

The `create_output_directories` validator automatically creates `results/`, `logs/`, `credentials/`, and `data/` directories on settings instantiation.

### Environment Variable Override

All settings can be overridden via environment variables. Each nested settings class has its own `env_prefix`:

```bash
MTDNA_TOOLS_BLASTN=/usr/local/bin/blastn
MTDNA_DIRS_RESULTS=/tmp/mtdna_results
MTDNA_DEFAULT_TRIM_THRESHOLD=0.05
MTDNA_TRACY_TRIM=7
MTDNA_TRACY_QUALITY_THRESHOLD=35
```

## Configuration

- `.env` — Primary environment file (checked first)
- `.env.local` — Local overrides (checked second, takes precedence)
- Environment variables — Highest precedence

## Cross-References

- [ARCHITECTURE.md](ARCHITECTURE.md) — System overview and configuration design decisions
- [variants.md](variants.md) — Uses `get_settings().regions` for region lookups
- [tools/tracy.md](tools/tracy.md) — Uses Tracy executable path from Settings
