"""Standard output directory structure for tool pipelines.

Provides ``OutputPaths`` for constructing the canonical subdirectory layout
used by BLASTn and Tracy (and eventually all tools):

    base_dir/
    ├── json/              # Batch.write() output
    ├── regions/           # Per-region intermediate JSONs
    └── preprocess/         # Tool-specific preprocessing artifacts

Pipeline code constructs these paths from the top-level ``--output-dir``
argument, then passes the appropriate subdirectory to each write function.
This keeps ``Batch.write()`` and ``write_region_jsons()`` signatures unchanged
while ensuring clean separation of output categories.

MS and Sequencher are not yet migrated to this pattern (deferred).
"""

from pathlib import Path


class OutputPaths:
    """Standard output directory structure for tool pipelines.

    Creates and exposes three subdirectories under a base directory:

    - ``json_dir`` — for ``Batch.write()`` output (per-sample JSONs and
      ``statistic_fullbatch.json``)
    - ``regions_dir`` — for ``write_region_jsons()`` output (per-region
      intermediate JSONs like ``HV1_{LID}.json`` and ``HV2-3_{LID}.json``)
    - ``preprocess_dir`` — for tool-specific preprocessing artifacts
      (``.blastn``, ``.fasta``, ``.fastq``, ``.quality`` for BLASTn;
      (per-sample ``.align``, ``.decomp``, ``.bcf``, and ``qc_report.json`` for Tracy)

    Usage::

        paths = OutputPaths(Path("results/tools/blastn/20260615_batch"))
        paths.create_dirs()
        # paths.json_dir       -> results/tools/blastn/20260615_batch/json
        # paths.regions_dir    -> results/tools/blastn/20260615_batch/regions
        # paths.preprocess_dir -> results/tools/blastn/20260615_batch/preprocess
    """

    def __init__(self, base_dir: Path) -> None:
        self.base_dir = Path(base_dir)
        self.json_dir = self.base_dir / "json"
        self.regions_dir = self.base_dir / "regions"
        self.preprocess_dir = self.base_dir / "preprocess"

    def create_dirs(self) -> None:
        """Create all three subdirectories if they don't exist."""
        self.json_dir.mkdir(parents=True, exist_ok=True)
        self.regions_dir.mkdir(parents=True, exist_ok=True)
        self.preprocess_dir.mkdir(parents=True, exist_ok=True)
