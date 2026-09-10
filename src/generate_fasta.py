"""
Generate FASTA - Position-Based Dictionary Approach

This module generates FASTA files from variant data using a position-based dictionary.
Variant processing is delegated to the canonical implementation in src.core.variants.
Post-generation validation is performed via src.validation.validation_fasta.

Key approach:
1. Read sample JSON files with variant data
2. Regenerate consensus HV sequences via generate_sequence (canonical flow)
3. Write FASTA output files per sample
4. Validate generated FASTA files and write alignment report

Default behavior: skip samples that already have a .fasta file in the output
directory. Existing files are NEVER modified or overwritten. This ensures
idempotent re-runs without losing previously generated results.

With --force: wipe the output directory clean first, then regenerate all
samples. This gives a clean-slate result identical to the old unconditional
wipe behavior.
"""

import shutil
from argparse import ArgumentParser, Namespace
from os import environ
from pathlib import Path

import dotenv
from loguru import logger

from src.core.sample import generate_sequence, load_sample
from src.validation.validation_fasta import FastaValidator

dotenv.load_dotenv()


class FastaGenerator:
    """Generate FASTA files from variant data using centralized variant processing."""

    def __init__(self, input_dir: str, output_dir: str, ref_path: str, *, force: bool = False) -> None:
        """Initialize FastaGenerator with input and output directories and reference path.

        Args:
            input_dir: Path to input directory containing per-sample JSON files.
            output_dir: Path to output directory for FASTA files.
            ref_path: Path to the reference FASTA file.
            force: If True, wipe output directory before writing (clean-slate mode).
                   If False (default), preserve existing .fasta files.
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.ref_path = ref_path
        self.force = force

        if force:
            # --force: wipe-and-recreate (same as old unconditional behavior)
            if self.output_dir.exists():
                shutil.rmtree(self.output_dir)
            self.output_dir.mkdir(parents=True)
        else:
            # Default: ensure output directory exists, preserve existing files
            self.output_dir.mkdir(parents=True, exist_ok=True)

    def write_fasta_file(self, output_path: str, region_sequences: dict[str, str]) -> None:
        """Write sequences to FASTA format file."""
        with Path(output_path).open("w") as f:
            f.writelines(f">{region}\n{seq}\n" for region, seq in region_sequences.items())

    def generate_reviewed_fasta(self) -> tuple[int, int, int]:
        """Generate FASTA files for all samples in the input directory.

        Returns:
            Tuple of (success_count, failure_count, skip_count).
            - success_count: samples newly generated successfully
            - failure_count: samples that failed processing
            - skip_count: samples skipped because .fasta already exists
        """
        success_count = 0
        failure_count = 0
        skip_count = 0

        for entry in self.input_dir.iterdir():
            sample_name = entry.name
            sample_path = self.input_dir / sample_name
            if not sample_path.is_dir():
                continue

            output_filename = f"{sample_name.split('.')[0]}.fasta"
            output_path = self.output_dir / output_filename

            # Default: skip samples that already have a .fasta file
            if not self.force and output_path.exists():
                logger.info(f"Skipping existing: {output_filename}")
                skip_count += 1
                continue

            try:
                json_path = sample_path / f"{sample_name}.json"
                if not json_path.exists():
                    logger.warning(f"JSON file not found for sample: {sample_name}")
                    failure_count += 1
                    continue

                # Load the standard Sample from the per-sample JSON, then
                # regenerate consensus HV sequences via generate_sequence
                # (canonical flow). Stored sequences on the Sample are not
                # used; the canonical generator is always invoked.
                sample_obj = load_sample(json_path)

                result = generate_sequence(sample_obj, str(self.ref_path))
                hv1 = "".join(result.hv_seqs.get("HV1", {}).get("seq", []))
                hv2 = "".join(result.hv_seqs.get("HV2", {}).get("seq", []))
                hv3 = "".join(result.hv_seqs.get("HV3", {}).get("seq", []))

                # Prepare sequences for saving
                save_regions = {
                    "HV1": hv1,
                    "HV2": hv2,
                    "HV3": hv3,
                }

                # Write to FASTA file
                self.write_fasta_file(str(output_path), save_regions)
                success_count += 1

            except (ValueError, KeyError, OSError) as e:
                logger.error(f"Processing sample {sample_name}: {e}")
                failure_count += 1

        return success_count, failure_count, skip_count


def parse_args() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Generate FASTA files from variant data (V2)")
    parser.add_argument("-i", "--input", required=True, help="Path to input directory")
    parser.add_argument("-o", "--output", required=True, help="Path to output directory")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA file")
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Wipe output directory before writing (clean-slate mode). "
        "Without this flag, existing .fasta files are preserved and skipped.",
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help="Results root for the validation report (overrides RESULTS_DIR). "
        "The report is written to <results-dir>/validation/validation_fasta_<batch>.txt.",
    )
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()

    results_dir = Path(args.results_dir) if args.results_dir else Path(environ.get("RESULTS_DIR", "results"))

    # Create generator with required reference path
    generator = FastaGenerator(args.input, args.output, args.reference, force=args.force)
    success, failure, skipped = generator.generate_reviewed_fasta()

    logger.info(f"FASTA generation complete. Successful: {success}, Failed: {failure}, Skipped: {skipped}")

    # Post-step: validate generated FASTA files
    batch_id = Path(args.output).name
    logger.info("Running post-generation FASTA validation …")
    validator = FastaValidator(
        fasta_dir=args.output,
        json_dir=args.input,
        ref_path=args.reference,
        batch_id=batch_id,
    )
    all_pass = validator.validate_all()

    validator.write_report(results_dir)

    if not all_pass:
        logger.warning("FASTA validation detected mismatches — review the report")


if __name__ == "__main__":
    main()
