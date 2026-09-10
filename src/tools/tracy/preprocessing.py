"""Tracy preprocessing — AB1 → Tracy decompose output.

Runs the Tracy decompose CLI tool on AB1 trace files to produce JSON
output files for ETL consumption. Called by pipeline.py only; never by
etl.py (per TOOLS_STANDARD.md §4.2).

This module is idempotent: running decompose twice on the same input
produces the same output (assuming Tracy itself is deterministic).
"""

import shutil
import subprocess
from pathlib import Path
from typing import NamedTuple

from loguru import logger

from src.config import get_settings


class DecomposeConfig(NamedTuple):
    """Configuration parameters for Tracy decompose."""

    trim: int
    pratio: float
    maxindel: int


def decompose_sample(
    sample: str,
    input_dir: str | Path,
    output_dir: str | Path,
    ref_path: str,
    *,
    config: DecomposeConfig | None = None,
) -> tuple[Path, ...]:
    """Run Tracy decompose and return JSON files produced in this run.

    Validates the binary path with ``shutil.which()`` before execution. Returning
    explicit outputs prevents callers from consuming stale files in the sample
    directory.
    """
    if config is None:
        config = DecomposeConfig(
            trim=get_settings().tracy.trim,
            pratio=get_settings().tracy.pratio,
            maxindel=get_settings().tracy.maxindel,
        )

    input_path = Path(input_dir)
    sample_output_dir = Path(output_dir) / sample
    sample_output_dir.mkdir(parents=True, exist_ok=True)

    # Find AB1 files for this sample
    ab1_files = sorted(str(path) for path in input_path.iterdir() if path.suffix == ".ab1" and sample in path.name)

    if not ab1_files:
        logger.warning("No AB1 files found for sample {}", sample)
        return ()

    tracy_binary = get_settings().tools.tracy

    # Validate binary path
    tracy_path = shutil.which(str(tracy_binary))
    if tracy_path is None:
        logger.error("Tracy binary not found: {}", tracy_binary)
        return ()

    json_outputs: list[Path] = []
    # Process each AB1 file
    for ab1_file in ab1_files:
        basename = Path(ab1_file).stem
        output_prefix = str(sample_output_dir / basename)

        command: list[str] = [
            tracy_path,
            "decompose",
            ab1_file,
            "--genome",
            ref_path,
            "--callVariants",
            "--trim",
            str(config.trim),
            "--pratio",
            str(config.pratio),
            "--maxindel",
            str(config.maxindel),
            "--outprefix",
            output_prefix,
        ]

        result = subprocess.run(command, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            logger.error("Tracy decompose failed for {}", basename)
            if result.stderr:
                logger.error("Error: {}", result.stderr)
            continue

        json_output = Path(f"{output_prefix}.json")
        if json_output.is_file():
            json_outputs.append(json_output)
        else:
            logger.error("Tracy decompose did not produce {}", json_output)

    return tuple(json_outputs)
