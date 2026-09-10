"""Batch aggregator: List[Sample] → statistic_fullbatch.json.

Delegates Sample serialization to :func:`src.core.sample.sample_to_dict`.
This module owns batch-level concerns only: aggregating multiple samples,
assembling the top-level batch dict, and writing JSON files.
"""

import json
from pathlib import Path
from typing import Any

from loguru import logger

from src.core.models import Sample
from src.core.sample import sample_to_dict


class Batch:
    """Aggregates a list of Samples into the canonical statistic_fullbatch.json format.

    Delegates per-sample serialization to :func:`sample_to_dict` from
    :mod:`src.core.sample`. Batch owns batch-level concerns: multi-sample
    aggregation and file writing.
    """

    def __init__(self, samples: list[Sample], ref_path: str = "ref/rCRS.fasta") -> None:
        self.samples = samples
        self.ref_path = ref_path

    def to_json(self) -> dict[str, Any]:
        """Convert all samples to the statistic_fullbatch.json dict.

        Returns:
            Dict mapping sample_id → canonical sample dict with HV sequences,
            counts, and intervals.

        """
        output: dict[str, Any] = {}
        for sample in sorted(self.samples, key=lambda s: s.sample_id):
            output[sample.sample_id] = sample_to_dict(sample, self.ref_path)
        return output

    def write(
        self,
        output_dir: str | Path,
        batch_id: str,
        *,
        nest_batch_id: bool = True,
    ) -> Path:
        """Write the batch JSON to a file in a batch-specific subdirectory.

        Args:
            output_dir: Base output directory.
            batch_id: Batch identifier used as the subdirectory name when
                ``nest_batch_id`` is True; otherwise unused for path construction.
            nest_batch_id: When True (default), write under
                ``{output_dir}/{batch_id}/`` so multiple batches can share one
                ``output_dir``. When False, write directly into ``output_dir`` —
                use this when ``output_dir`` is already batch-specific (for
                example ``OutputPaths.json_dir`` under a per-batch base such as
                ``results/tools/<tool>/<batch_id>/json``), to avoid a duplicated
                ``<batch_id>`` path segment.

        Returns:
            The resolved path that was written.

        """
        data = self.to_json()
        base = Path(output_dir)
        batch_dir = base / batch_id if nest_batch_id else base
        out_path = (batch_dir / "statistic_fullbatch.json").resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2)
        logger.success("Wrote batch JSON ({} samples) to {}", len(data), out_path)

        # Write per-sample JSON files
        for sample_id, sample_dict in data.items():
            sample_dir = batch_dir / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            sample_path = sample_dir / f"{sample_id}.json"
            with sample_path.open("w", encoding="utf-8") as fh:
                json.dump(sample_dict, fh, indent=2)
        logger.success("Wrote per-sample JSONs ({} samples) to {}", len(data), batch_dir)

        return out_path
