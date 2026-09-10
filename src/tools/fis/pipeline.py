"""FIS batch orchestration and canonical JSON writing."""

import argparse
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from loguru import logger

from src.core.batch import Batch
from src.core.models import Sample
from src.core.paths import OutputPaths
from src.core.region import REGION_TO_KEYS, write_region_jsons
from src.tools.fis.etl import FisRun, build_samples, load_run, raw_records

MAX_FASTA_N_FRACTION = 0.9


def default_region_workers() -> int:
    """Return the bounded default for independent region JSON writes."""
    return max(1, math.ceil((os.cpu_count() or 1) * 0.75))


def _write_sample_regions(sample: Sample, regions_dir: Path, reference: str) -> None:
    """Write all configured region JSONs for one FIS sample."""
    write_region_jsons(sample, regions_dir, reference, REGION_TO_KEYS)


def write_canonical_samples(samples: dict[str, Sample], output_dir: Path, reference: Path, batch_id: str) -> Path:
    """Write standard FIS batch and concurrently generated region JSON outputs."""
    paths = OutputPaths(output_dir)
    paths.create_dirs()
    workers = default_region_workers()
    logger.info("Writing FIS region JSONs for {} samples with {} workers", len(samples), workers)
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="fis-region") as executor:
        futures = [
            executor.submit(_write_sample_regions, sample, paths.regions_dir, str(reference))
            for sample in samples.values()
        ]
        for future in futures:
            future.result()
    return Batch(list(samples.values()), ref_path=str(reference)).write(paths.json_dir, batch_id, nest_batch_id=False)


def _write_ngs_artifacts(run: FisRun, samples: dict, output_dir: Path) -> None:
    """Write established NGS FIS artifacts from the same parsed ETL run."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "FIS_QCFail_samples.txt").write_text(
        "".join(f"{sample}\n" for sample in run.qc_fail_samples), encoding="utf-8"
    )
    filtered_t2 = run.t2.loc[~run.t2["SampleCode"].astype(str).isin(run.qc_fail_samples)].copy()
    filtered_t2.rename(columns={"SampleCode": "sample_id_base"}).to_csv(output_dir / "FIS.tsv", sep="\t", index=False)
    transformed = [((sample.information or {}).get("fis_metadata", {})) for sample in samples.values()]
    pd.DataFrame(transformed).to_csv(output_dir / "FIS_transformed.tsv", sep="\t", index=False)
    (output_dir / "FIS.json").write_text(
        json.dumps(raw_records(run, samples), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if run.t3 is None:
        msg = "R2 workbook has no T3.Seq sheet required for NGS artifacts"
        raise ValueError(msg)
    sequences = run.t3.loc[~run.t3["SampleCode"].astype(str).isin(run.qc_fail_samples)].copy()
    sequences.to_csv(output_dir / "FIS.seq.tsv", sep="\t", index=False)
    fasta_dir = output_dir / "FASTA"
    fasta_dir.mkdir(exist_ok=True)
    region_order = ("HV1", "HV2", "HV3")
    for sample_id, rows in sequences.groupby("SampleCode", sort=False):
        regions: dict[str, str] = {}
        for _, row in rows.iterrows():
            marker = str(row.get("Marker", ""))
            consensus = row.get("Consensus")
            if pd.isna(consensus):
                continue
            region = (
                "HV1"
                if "HVS-I[" in marker
                else "HV2"
                if "HVS-II[" in marker
                else "HV3"
                if "HVS-III[" in marker
                else None
            )
            if region is None:
                logger.warning("Skipping unrecognized FIS sequence marker for {}: {}", sample_id, marker)
                continue
            start = marker.split(":", maxsplit=1)[1].split("-", maxsplit=1)[0] if ":" in marker else ""
            sequence = f"GGGGT{consensus}" if region == "HV2" and start == "73" else str(consensus)
            if region in regions:
                logger.warning("Ignoring duplicate FIS {} sequence for {}", region, sample_id)
                continue
            regions[region] = sequence
        missing = [region for region in region_order if region not in regions]
        if missing:
            logger.warning("{} missing FIS sequence regions: {}", sample_id, ", ".join(missing))
        if not regions or any(
            not sequence or sequence.upper().count("N") / len(sequence) > MAX_FASTA_N_FRACTION
            for sequence in regions.values()
        ):
            logger.warning(
                "Excluded {} from FIS FASTA output due to missing sequence or excessive N content", sample_id
            )
            continue
        content = "".join(f">{region}\n{regions[region]}\n" for region in region_order if region in regions)
        (fasta_dir / f"{sample_id}.fasta").write_text(content, encoding="utf-8")


def process_batch(
    r1_file: Path,
    r2_file: Path,
    output_dir: Path,
    reference: Path,
    batch_id: str,
    *,
    write_ngs_artifacts: bool = False,
) -> dict:
    """Process FIS workbooks into canonical Samples and optional NGS artifacts."""
    run = load_run(r1_file, r2_file)
    samples = build_samples(run, reference, batch_id)
    write_canonical_samples(samples, output_dir, reference, batch_id)
    if write_ngs_artifacts:
        _write_ngs_artifacts(run, samples, output_dir)
    logger.success("Processed FIS batch: {} samples", len(samples))
    return samples


def main() -> None:
    """Run the FIS tool pipeline."""
    parser = argparse.ArgumentParser(description="Convert FIS Excel workbooks to canonical Sample JSON")
    parser.add_argument("--r1-file", type=Path, required=True, help="R1 DataProduction workbook")
    parser.add_argument("--r2-file", type=Path, required=True, help="R2 MT workbook")
    parser.add_argument("--output-dir", type=Path, required=True, help="Batch output directory")
    parser.add_argument("--reference", type=Path, default=Path("ref/rCRS.fasta"), help="rCRS FASTA")
    parser.add_argument("--batch-id", required=True, help="FIS batch identifier")
    parser.add_argument(
        "--write-ngs-artifacts", action="store_true", help="Write established root-level NGS FIS artifacts"
    )
    args = parser.parse_args()
    process_batch(
        args.r1_file,
        args.r2_file,
        args.output_dir,
        args.reference,
        args.batch_id,
        write_ngs_artifacts=args.write_ngs_artifacts,
    )


if __name__ == "__main__":
    main()
