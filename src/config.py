"""
Configuration module for mtDNA analysis pipeline (2026 best practices — enhanced).

- Strong FilePath/DirectoryPath validation
- Automatic creation of output directories
- Multiple .env file support
- Self-documenting fields
- Compatible with Pydantic v2 (no more deprecation warnings)
"""

import os
from functools import lru_cache
from pathlib import Path

from loguru import logger
from pydantic import Field, field_validator, model_validator
from pydantic.types import DirectoryPath, FilePath
from pydantic_settings import BaseSettings, SettingsConfigDict


class TracySettings(BaseSettings):
    """Tracy tool runtime parameters."""

    trim: int = Field(default=7, description="Tracy decompose trim parameter")
    pratio: float = Field(default=0.3, description="Tracy decompose pratio parameter")
    maxindel: int = Field(default=1000, description="Tracy decompose maxindel parameter")
    quality_threshold: int = Field(default=35, description="Minimum quality score for variant calling")
    min_peak_value: int = Field(default=125, description="Minimum peak height for variant calling")
    heteroplasmy_threshold: float = Field(default=0.8, description="Minimum heteroplasmy ratio for detection")
    max_workers: int | None = Field(default=None, description="Max parallel workers (default: system CPU count)")
    qc_enabled: bool = Field(default=True, description="Enable optional per-trace evidence QC")
    signal_noise_enabled: bool = Field(default=True, description="Evaluate signal-noise evidence")
    polyc_enabled: bool = Field(default=False, description="Evaluate directional polyC and homopolymer evidence")
    read_edge_enabled: bool = Field(default=True, description="Report exploratory read-edge uncertainty")
    noise_mask_enabled: bool = Field(default=True, description="Exclude variants and coverage in likely-noisy ranges")
    noise_window_size: int = Field(default=15, ge=10, le=20, description="Bases per trace-noise window")
    noise_window_step: int = Field(default=5, ge=1, le=20, description="Bases between trace-noise windows")
    noise_min_valid_bases: int = Field(default=10, ge=1, le=20, description="Valid bases required per window")
    noise_min_supporting_windows: int = Field(
        default=2,
        ge=2,
        description="Likely-noisy windows required to create a mask range",
    )
    noise_snr_threshold: float = Field(default=3.0, gt=0, description="Minimum acceptable peak SNR")
    noise_low_snr_threshold: float = Field(default=2.0, gt=0, description="Strong low-end SNR threshold")
    noise_purity_threshold: float = Field(default=0.60, ge=0, le=1, description="Minimum acceptable signal purity")
    noise_background_threshold: float = Field(default=0.35, ge=0, description="Maximum background-to-signal ratio")
    noise_second_peak_threshold: float = Field(default=0.50, ge=0, description="Maximum second-to-first peak ratio")
    noise_quality_threshold: float = Field(default=20.0, ge=0, description="Minimum acceptable basecall quality")
    noise_signal_threshold: float = Field(default=75.0, ge=0, description="Minimum acceptable peak signal")
    noise_bad_base_fraction: float = Field(default=0.40, ge=0, le=1, description="Fraction required for noisy status")
    noise_suspicious_base_fraction: float = Field(
        default=0.20,
        ge=0,
        le=1,
        description="Fraction required for suspicious status",
    )

    model_config = SettingsConfigDict(env_prefix="MTDNA_TRACY_", extra="forbid")

    @model_validator(mode="after")
    def validate_noise_windows(self) -> "TracySettings":
        """Ensure window settings form usable overlapping windows."""
        if self.noise_window_step > self.noise_window_size:
            msg = "noise_window_step cannot exceed noise_window_size"
            raise ValueError(msg)
        if self.noise_min_valid_bases > self.noise_window_size:
            msg = "noise_min_valid_bases cannot exceed noise_window_size"
            raise ValueError(msg)
        if self.noise_suspicious_base_fraction > self.noise_bad_base_fraction:
            msg = "noise_suspicious_base_fraction cannot exceed noise_bad_base_fraction"
            raise ValueError(msg)
        return self


class BlastnSettings(BaseSettings):
    """BLASTN tool runtime parameters."""

    trim_threshold: float = Field(default=0.01, description="Quality threshold for seqtk trimming")
    word_size: int = Field(default=22, description="Word size for BLASTN alignment")
    quality_threshold: int = Field(default=35, description="Minimum quality score for variant calling")
    min_peak_value: int = Field(default=125, description="Minimum peak value for variant calling")
    heteroplasmy_threshold: float = Field(default=0.8, description="Minimum heteroplasmy ratio for detection")
    max_workers: int | None = Field(default=None, description="Max parallel workers (default: system CPU count)")
    blastn_threads: int = Field(
        default=1,
        description="Threads per BLASTN subprocess. Use 1 with sample-parallel max_workers",
    )

    model_config = SettingsConfigDict(env_prefix="MTDNA_BLASTN_", extra="forbid")


class SequencherSettings(BaseSettings):
    """Sequencher tool runtime parameters."""

    max_workers: int | None = Field(default=None, description="Max parallel workers (default: system CPU count)")

    model_config = SettingsConfigDict(env_prefix="MTDNA_SEQUENCHER_", extra="forbid")


class ToolsSettings(BaseSettings):
    """External tool paths (validated as executable files)."""

    ugene: FilePath = Field(
        default=Path("./tools/ugene-50.0/ugene"),
        description="UGENE executable",
    )
    blastn: FilePath = Field(
        default=Path("./tools/ugene-50.0/tools/blast/blastn"),
        description="BLASTN executable",
    )
    seqtk: FilePath = Field(
        default=Path("./tools/seqtk/seqtk"),
        description="seqtk executable",
    )
    tracy: FilePath = Field(
        default=Path("./tools/tracy"),
        description="Tracy executable",
    )

    model_config = SettingsConfigDict(env_prefix="MTDNA_TOOLS_", extra="forbid")


class DirectoryPathsSettings(BaseSettings):
    """Project directory paths."""

    # Input/reference directories - must exist
    filter_rules: DirectoryPath = Field(default=Path("rules/filter"), description="Filter rules directory")
    convert_rules: DirectoryPath = Field(default=Path("rules/convert"), description="Convert rules directory")
    data: Path = Field(default=Path("data"), description="Input data directory")
    ref: DirectoryPath = Field(default=Path("ref"), description="Reference sequences directory")
    tools: DirectoryPath = Field(default=Path("tools"), description="Bundled tools directory")

    # Output/working directories - will be auto-created
    results: Path = Field(default=Path("results"), description="Analysis results directory")
    logs: Path = Field(default=Path("logs"), description="Log files directory")
    credentials: Path = Field(default=Path("credentials"), description="Credentials directory")

    model_config = SettingsConfigDict(env_prefix="MTDNA_DIRS_", extra="forbid")


class GenomicRegionsSettings(BaseSettings):
    """Genomic regions of interest for mtDNA analysis."""

    REGIONS: dict[str, list[int]] = Field(
        default_factory=lambda: {
            "HV1": [16024, 16365],
            "HV2": [73, 340],
            "HV3": [438, 576],
        },
    )
    POLYC_REGIONS: dict[str, list[int]] = Field(default_factory=lambda: {"HV2": [303, 315]})
    POLYC_QC_REGIONS: list[tuple[int, int]] = Field(
        default_factory=lambda: [
            (16183, 16193),
            (16289, 16290),
            (302, 315),
            (455, 463),
            (515, 525),
        ],
    )
    SKIP_AFTER_POLYC: dict[str, list[int]] = Field(default_factory=lambda: {"HV1": [16193, 16195]})

    model_config = SettingsConfigDict(extra="forbid")


class ParametersSettings(BaseSettings):
    """Default parameters for analysis."""

    DEFAULT_TRIM_THRESHOLD: float = Field(default=0.01, description="Default quality trim threshold")
    DEFAULT_WORD_SIZE: int = 22
    DEFAULT_NUM_THREADS: int = Field(default_factory=lambda: os.cpu_count() or 1)

    AB1_EXTENSION: str = ".ab1"
    FASTQ_EXTENSION: str = ".fastq"
    FASTA_EXTENSION: str = ".fasta"
    BLASTN_EXTENSION: str = ".blastn"

    REPORT_BASE_URL: str = "https://genestory.ai/pdfexporter/pdf/generate/mt-adn"
    REPORT_RETRIES: int = 5
    REPORT_RETRY_DELAY: int = 5

    model_config = SettingsConfigDict(env_prefix="MTDNA_", extra="forbid")


class Settings(BaseSettings):
    """Main settings class."""

    project_root: Path | None = None

    tools: ToolsSettings = Field(default_factory=ToolsSettings)
    directories: DirectoryPathsSettings = Field(default_factory=DirectoryPathsSettings)
    regions: GenomicRegionsSettings = Field(default_factory=GenomicRegionsSettings)
    parameters: ParametersSettings = Field(default_factory=ParametersSettings)
    tracy: TracySettings = Field(default_factory=TracySettings)
    blastn: BlastnSettings = Field(default_factory=BlastnSettings)
    sequencher: SequencherSettings = Field(default_factory=SequencherSettings)

    model_config = SettingsConfigDict(
        env_file=[".env", ".env.local"],
        env_file_encoding="utf-8",
        env_nested_delimiter="__",
        extra="ignore",
    )

    def model_post_init(self, __context: object, /) -> None:
        """Convert relative paths to absolute based on project root."""
        if self.project_root is None:
            self.project_root = Path(__file__).resolve().parent.parent

        # Resolve tool paths
        for field_name in ToolsSettings.model_fields:
            tool_path: Path = getattr(self.tools, field_name)
            if not tool_path.is_absolute():
                setattr(self.tools, field_name, (self.project_root / tool_path).resolve())

        # Resolve directory paths
        for field_name in DirectoryPathsSettings.model_fields:
            dir_path: Path = getattr(self.directories, field_name)
            if not dir_path.is_absolute():
                setattr(self.directories, field_name, (self.project_root / dir_path).resolve())

    @field_validator("directories", mode="after")
    @classmethod
    def create_output_directories(cls, dirs: DirectoryPathsSettings) -> DirectoryPathsSettings:
        """Automatically create output/working directories."""
        for name in ("results", "logs", "credentials", "data"):
            path: Path = getattr(dirs, name)
            path.mkdir(parents=True, exist_ok=True)
        return dirs

    def show(self) -> None:
        """Pretty-print current configuration (handy for debugging)."""
        logger.info("=== mtDNA Pipeline Configuration ===")
        logger.info("Project root : {}", self.project_root)
        logger.info("Tools        : {}", self.tools.model_dump())
        logger.info("Directories  : {}", self.directories.model_dump())
        logger.info("====================================")


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Cached settings (recommended access pattern)."""
    return Settings()
