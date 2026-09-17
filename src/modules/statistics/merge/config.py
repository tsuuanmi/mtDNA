"""Default merge input/output locations and shared parsing constants."""

from pathlib import Path

# Minimum tab-separated columns required to identify a sample (sample_id, barcode)
MIN_TSV_PARTS = 2

# Production data/results roots used as the CLI defaults, like the population module
DEFAULT_DATA_DIR = Path("/mnt/nas/bca/mtDNA/science/data")
DEFAULT_RESULTS_DIR = Path("/mnt/nas/bca/mtDNA/science/results")

# Default list of batches to exclude
DEFAULT_EXCLUDED_BATCHES = [
    "20250429_mtDNA_100",
]
