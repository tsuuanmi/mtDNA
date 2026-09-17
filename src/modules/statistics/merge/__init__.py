"""Merge per-batch metadata and results into consolidated TSV/JSON exports."""

from src.modules.statistics.merge.config import DEFAULT_DATA_DIR, DEFAULT_EXCLUDED_BATCHES, DEFAULT_RESULTS_DIR
from src.modules.statistics.merge.metadata import convert_excel_to_tsv, process_excel_files, read_metadata_files
from src.modules.statistics.merge.models import (
    BatchSummary,
    DirectoryScanResult,
    MergePaths,
    SourceSignatureEntry,
)
from src.modules.statistics.merge.reporting import analyze_sample_discrepancies, batch_statistics, create_merged_tsv
from src.modules.statistics.merge.scanner import DirectoryScanner, get_successful_samples
from src.modules.statistics.merge.statistics_json import merge_statistics_json

__all__ = [
    "DEFAULT_DATA_DIR",
    "DEFAULT_EXCLUDED_BATCHES",
    "DEFAULT_RESULTS_DIR",
    "BatchSummary",
    "DirectoryScanResult",
    "DirectoryScanner",
    "MergePaths",
    "SourceSignatureEntry",
    "analyze_sample_discrepancies",
    "batch_statistics",
    "convert_excel_to_tsv",
    "create_merged_tsv",
    "get_successful_samples",
    "merge_statistics_json",
    "process_excel_files",
    "read_metadata_files",
]
