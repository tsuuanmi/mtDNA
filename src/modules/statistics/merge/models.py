"""Typed data contracts for merge paths, batch summaries, and source signatures."""

from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict


@dataclass(frozen=True)
class MergePaths:
    """Input and output path contract for one merge run."""

    data_dir: Path
    results_dir: Path
    output_dir: Path

    @property
    def metadata_dir(self) -> Path:
        """Per-batch metadata Excel workbooks and TSVs (sample_id, barcode)."""
        return self.data_dir / "metadata"

    @property
    def regenerate_dir(self) -> Path:
        """Per-batch regenerate results, keyed by statistic_fullbatch.json."""
        return self.results_dir / "manual_pipeline" / "regenerate"

    @property
    def fasta_dir(self) -> Path:
        """Per-batch FASTA outputs."""
        return self.results_dir / "fasta"

    @property
    def archive_dir(self) -> Path:
        """Archived per-batch FASTA outputs."""
        return self.results_dir / "archive"

    @property
    def merged_metadata_tsv(self) -> Path:
        """Merged batch/sample/barcode TSV output."""
        return self.output_dir / "merged_metadata.tsv"

    @property
    def merged_statistics_json(self) -> Path:
        """Sample-ID-keyed merged statistics JSON output."""
        return self.output_dir / "merged_statistics.json"

    @property
    def merged_statistics_manifest(self) -> Path:
        """Freshness manifest for the merged statistics JSON."""
        return self.output_dir / "merged_statistics.manifest.json"

    @property
    def batch_statistics_tsv(self) -> Path:
        """Per-batch Successful/Metadata/Regenerate/Fasta/Status table."""
        return self.output_dir / "batch_statistics.tsv"

    @property
    def sample_discrepancies_tsv(self) -> Path:
        """Per-sample TSV-vs-JSON discrepancy table."""
        return self.output_dir / "sample_discrepancies.tsv"

    @property
    def scan_cache(self) -> Path:
        """Shared scan and metadata source cache."""
        return self.output_dir / "scan_cache.json"


@dataclass
class BatchSummary:
    """Summary of samples in a batch across different processing stages."""

    batch_name: str
    successful_count: int = 0
    metadata_count: int = 0
    regenerate_count: int = 0
    fasta_count: int = 0

    @property
    def status(self) -> str:
        """Determine processing status based on counts."""
        if self.metadata_count == 0:
            return "No Meta"
        if self.fasta_count == 0:
            return "No FASTA"
        if self.successful_count == 0:
            return "No Pass"
        if self.metadata_count == self.successful_count:
            return "OK"
        return "Mismatch"


@dataclass
class DirectoryScanResult:
    """Result of scanning directories for samples."""

    regenerate_samples: dict[str, set[str]]
    fasta_samples: dict[str, set[str]]


class SourceSignatureEntry(TypedDict):
    """File metadata used to validate merged-statistics cache freshness."""

    path: str
    mtime_ns: int
    ctime_ns: int
    size: int
