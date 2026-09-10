"""Positive/Negative Control (PC/NTC) quality-control module.

Runs the modular BLASTN pipeline on PC/NTC control samples and, when an
expected-profile directory is provided, compares each PC sample against its
matched Sequencher expected profile to decide pass/fail. NTC samples pass
when they produce no variants. Emits a JSON QC report.

This replaces the legacy ``src/process_pc_samples.py`` Tracy cross-comparison
path: heteroplasmy removal is intentionally dropped (PC output is validation
evidence only; no downstream step consumes it).

Control sample tokens recognized in filenames/sample IDs: ``PC1``, ``PC2``,
``NTC``. Add more tokens here if additional controls are introduced.

Output is written under ``<results_root>/validation/<batch_id>/`` (default
``<cwd>/results/validation/``), matching the results-layout convention used by
``src/validation/validation_fasta.py``.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import dotenv
from loguru import logger

from src.core.comparison import SampleComparator
from src.core.models import Sample
from src.core.variants import variant_to_dict
from src.tools.blastn.pipeline import BatchOptions as BlastnBatchOptions
from src.tools.blastn.pipeline import process_batch as blastn_process_batch
from src.tools.sequencher.pipeline import BatchOptions as SequencherBatchOptions
from src.tools.sequencher.pipeline import process_batch as sequencher_process_batch

dotenv.load_dotenv()

# Tokens that identify control samples. The order matters for classification:
# NTC is checked before PC so an NTC id is never misclassified as PC.
CONTROL_TOKENS: tuple[str, ...] = ("PC1", "PC2", "NTC")
PC_TOKENS: tuple[str, ...] = ("PC1", "PC2")
NTC_LABEL: str = "NTC"

# Default reference path (matches the rest of the project).
DEFAULT_REF_PATH = "ref/rCRS.fasta"

# Minimum "_"-separated parts expected in an AB1 filename to extract a sample id.
MIN_FILENAME_PARTS = 3


def _classify(sample_id: str) -> str:
    """Return the control type for a sample id: 'PC', 'NTC', or ''."""
    upper = sample_id.upper()
    if NTC_LABEL in upper:
        return "NTC"
    for token in PC_TOKENS:
        if token.upper() in upper:
            return "PC"
    return ""


def _pc_token(sample_id: str) -> str | None:
    """Return the PC token ('PC1' or 'PC2') carried by a PC sample id, or None."""
    upper = sample_id.upper()
    for token in PC_TOKENS:
        if token.upper() in upper:
            return token
    return None


def _ids_for_token(sample_ids: list[str], token: str) -> list[str]:
    """Sample ids that carry a control token (substring match, case-insensitive)."""
    up = token.upper()
    return [sid for sid in sample_ids if up in sid.upper()]


def _variant_dicts(sample: Sample) -> list[dict[str, Any]]:
    """Serialize a Sample's variants to canonical dicts, dropping None entries."""
    return [d for d in (variant_to_dict(v) for v in sample.variants) if d is not None]


class PCNTCProcessor:
    """Process PC/NTC control samples through the BLASTN-only QC flow."""

    def __init__(
        self,
        batch_id: str,
        data_dir: str | None = None,
        results_dir: str | None = None,
        ref_path: str = DEFAULT_REF_PATH,
        pc_reference_dir: str | None = None,
    ) -> None:
        self.batch_id = batch_id
        self.data_dir = Path(data_dir) if data_dir else Path(os.environ.get("DATA_DIR", ""))
        # Results root: explicit --results-dir, else RESULTS_DIR env, else <cwd>/results.
        # Validation output is written under <results_root>/validation/<batch_id>/ to match
        # the results-layout convention used by src/validation/validation_fasta.py.
        results_root = (
            Path(results_dir) if results_dir else Path(os.environ.get("RESULTS_DIR", str(Path.cwd() / "results")))
        )
        self.results_root = results_root
        self.ref_path = ref_path
        self.pc_reference_dir = Path(pc_reference_dir) if pc_reference_dir else None

        # Resolve the reference FASTA: absolute paths are used as-is; relative
        # paths are resolved against the current working directory (project root
        # when invoked via `python -m src.modules.quality_control.pc_ntc`).
        self.reference = Path(self.ref_path) if Path(self.ref_path).is_absolute() else Path.cwd() / self.ref_path

        # PC/NTC validation directories (results/validation/<batch_id>/...). The
        # copied-control working area (pc_raw_dir) lives under results so it never
        # writes into the AB1 source (the --data-dir batch directory).
        self.validation_base = results_root / "validation"
        self.pc_results_base = self.validation_base / self.batch_id
        self.pc_raw_dir = self.pc_results_base / "raw"
        self.pc_blastn_dir = self.pc_results_base / "blastn"
        self.pc_expected_dir = self.pc_results_base / "expected"
        self.report_path = self.pc_results_base / "pc_ntc_qc_report.json"
        self.sample_list_file = Path(f"{self.pc_raw_dir}.txt")

    # -- Step 1: prepare -----------------------------------------------------
    def prepare_pc_data(self) -> list[str]:
        """Copy PC/NTC AB1 files from the data dir and write the sample-list TXT.

        Returns the list of control sample ids. Raises on missing source data.
        The AB1 source is the ``--data-dir`` value: the batch directory that
        contains the HV subdirectories holding the .ab1 files. For a normal NAS
        run, ``scripts/pipeline.sh`` passes ``LAB_DATA_DIR/<batch_id>`` here; for
        an isolated/local run it passes the user-provided ``--data-dir``.
        """
        logger.info("[Step 1] Organizing PC/NTC samples")

        source_dir = self.data_dir
        if not source_dir.exists():
            logger.error("Source directory not found: {}", source_dir)
            raise FileNotFoundError(source_dir)

        if self.pc_raw_dir.exists():
            shutil.rmtree(self.pc_raw_dir)
        self.pc_raw_dir.mkdir(parents=True, exist_ok=True)

        # Copy all .ab1 files (same approach as scripts/pipeline.sh). pc_raw_dir
        # lives under results (outside the source in the normal layout); the
        # exclusion is a safety net for configs where results nest under the
        # data dir, so the copy never reads from its own destination.
        find_cmd = [
            "find",
            str(source_dir),
            "-type",
            "f",
            "-name",
            "*.ab1",
            "!",
            "-path",
            f"{self.pc_raw_dir}{os.sep}*",
            "-exec",
            "cp",
            "{}",
            str(self.pc_raw_dir),
            ";",
        ]
        subprocess.run(find_cmd, capture_output=True, text=True, check=True)

        # Filter: keep only control files (by token substring)
        kept = 0
        for ab1 in list(self.pc_raw_dir.glob("*.ab1")):
            name = ab1.name.upper()
            if any(tok in name for tok in CONTROL_TOKENS):
                kept += 1
            else:
                ab1.unlink(missing_ok=True)

        if kept == 0:
            logger.warning("No PC/NTC files found after filtering")
            return []

        # Derive sample ids from filenames: prefix_region_sampleid_etc -> parts[2]
        # (same convention as legacy step1) e.g. H03_20260410_PC1_mtDNA_247547_HV2F_24
        sample_ids: set[str] = set()
        for ab1 in self.pc_raw_dir.glob("*.ab1"):
            parts = ab1.stem.split("_")
            if len(parts) >= MIN_FILENAME_PARTS and _classify(parts[2]):
                sample_ids.add(parts[2])
            elif _classify(ab1.stem):
                sample_ids.add(ab1.stem)

        ordered = sorted(sample_ids)
        with self.sample_list_file.open("w") as f:
            for sid in ordered:
                f.write(f"{sid}\n")

        logger.info("Found {} control samples: {}", len(ordered), ordered)
        return ordered

    # -- Step 2: BLASTN -------------------------------------------------------
    def run_blastn(self) -> dict[str, Sample]:
        """Run the modular BLASTN pipeline on the control samples."""
        logger.info("[Step 2] Running BLASTN on control samples")

        if not self.sample_list_file.exists():
            raise FileNotFoundError(self.sample_list_file)

        if self.pc_blastn_dir.exists():
            shutil.rmtree(self.pc_blastn_dir)
        self.pc_blastn_dir.mkdir(parents=True, exist_ok=True)

        options = BlastnBatchOptions(
            ref_path=str(self.reference),
            batch_id=self.batch_id,
            samples_path=self.sample_list_file,
        )
        samples = blastn_process_batch(
            input_dir=self.pc_raw_dir,
            output_dir=self.pc_blastn_dir,
            options=options,
        )
        logger.info("BLASTN produced {} control samples", len(samples))
        return samples

    # -- Step 3 (optional): expected profiles --------------------------------
    def load_expected(self) -> dict[str, Sample] | None:
        """Convert the expected Sequencher TXT directory to Samples, if provided."""
        if self.pc_reference_dir is None:
            logger.info("[Step 3] No expected-profile directory provided; skipping comparison")
            return None

        logger.info("[Step 3] Loading expected profiles from {}", self.pc_reference_dir)
        if not self.pc_reference_dir.exists():
            logger.warning("Expected-profile directory not found: {}", self.pc_reference_dir)
            return None

        self.pc_expected_dir.mkdir(parents=True, exist_ok=True)
        options = SequencherBatchOptions(
            ref_path=str(self.reference),
            batch_id=self.batch_id,
            json_dir=self.pc_expected_dir,
        )
        expected = sequencher_process_batch(
            input_dir=self.pc_reference_dir,
            output_dir=self.pc_expected_dir,
            options=options,
        )
        logger.info("Loaded {} expected profiles", len(expected))
        return expected

    # -- Step 4: evaluate + report -------------------------------------------
    def evaluate(
        self,
        blastn_samples: dict[str, Sample],
        expected_samples: dict[str, Sample] | None,
        control_ids: list[str],
    ) -> dict[str, Any]:
        """Evaluate PC/NTC pass/fail and build the report payload."""
        entries: list[dict[str, Any]] = []

        ntc_ids = _ids_for_token(control_ids, NTC_LABEL)
        pc_ids = [sid for sid in control_ids if _classify(sid) == "PC"]

        # NTC: pass = no variants (absent or empty); fail = >=1 variant
        for sid in ntc_ids:
            sample = blastn_samples.get(sid)
            if sample is None:
                entries.append(self._ntc_entry(sid, verdict="Y", observed=0, note="absent_from_blastn_output"))
            elif len(sample.variants) == 0:
                entries.append(self._ntc_entry(sid, verdict="Y", observed=0, note=None))
            else:
                entries.append(
                    self._ntc_entry(sid, verdict="N", observed=len(sample.variants), note="ntc_has_variants")
                )

        # PC: priority pass rule (see PRD revision)
        entries.extend(self._pc_entry(sid, blastn_samples, expected_samples) for sid in pc_ids)

        return {
            "batch_id": self.batch_id,
            "pc_reference_dir": str(self.pc_reference_dir) if self.pc_reference_dir else None,
            "ref_path": str(self.reference),
            "samples": sorted(entries, key=lambda e: e["sample_id"]),
        }

    def _ntc_entry(self, sid: str, verdict: str, observed: int, note: str | None) -> dict[str, Any]:
        return {
            "sample_id": sid,
            "control_type": "NTC",
            "pass": verdict,
            "concordance": "N/A",
            "observed_variant_count": observed,
            "expected_variant_count": 0,
            "unique_observed": [],
            "unique_expected": [],
            "note": note,
        }

    def _pc_entry(
        self,
        sid: str,
        blastn_samples: dict[str, Sample],
        expected_samples: dict[str, Sample] | None,
    ) -> dict[str, Any]:
        sample = blastn_samples.get(sid)

        # Rule 1: empty/absent PC -> FAIL (assay produced no data)
        if sample is None or len(sample.variants) == 0:
            return {
                "sample_id": sid,
                "control_type": "PC",
                "pass": "N",
                "concordance": "N/A",
                "observed_variant_count": 0 if sample is None else len(sample.variants),
                "expected_variant_count": 0,
                "unique_observed": [],
                "unique_expected": [],
                "note": "pc_no_variants",
            }

        observed_dicts = _variant_dicts(sample)

        # Rule 4: no expected dir provided -> N/A
        if expected_samples is None:
            return {
                "sample_id": sid,
                "control_type": "PC",
                "pass": "N/A",
                "concordance": "N/A",
                "observed_variant_count": len(observed_dicts),
                "expected_variant_count": 0,
                "unique_observed": observed_dicts,
                "unique_expected": [],
                "note": "no_expected_profile",
            }

        # Match expected by PC token substring
        token = _pc_token(sid)
        expected_sid = None
        if token is not None:
            for key in expected_samples:
                if token.upper() in key.upper():
                    expected_sid = key
                    break

        # Rule 2: expected dir provided but no match -> N/A
        if expected_sid is None:
            return {
                "sample_id": sid,
                "control_type": "PC",
                "pass": "N/A",
                "concordance": "N/A",
                "observed_variant_count": len(observed_dicts),
                "expected_variant_count": 0,
                "unique_observed": observed_dicts,
                "unique_expected": [],
                "note": "no_expected_profile",
            }

        # Rule 3: matched -> is_concordant decides pass
        expected_sample = expected_samples[expected_sid]
        expected_dicts = _variant_dicts(expected_sample)
        concordance = SampleComparator.is_concordant(sample, expected_sample)
        comparison = SampleComparator.compare(sample, expected_sample)
        return {
            "sample_id": sid,
            "control_type": "PC",
            "pass": "Y" if concordance == "Y" else "N",
            "concordance": concordance,
            "observed_variant_count": len(observed_dicts),
            "expected_variant_count": len(expected_dicts),
            "unique_observed": comparison.unique_variants_a,
            "unique_expected": comparison.unique_variants_b,
            "note": None,
        }

    def write_report(self, report: dict[str, Any]) -> Path:
        self.pc_results_base.mkdir(parents=True, exist_ok=True)
        with self.report_path.open("w") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        logger.info("QC report written to {}", self.report_path)
        return self.report_path

    def run(self) -> dict[str, Any]:
        control_ids = self.prepare_pc_data()
        if not control_ids:
            logger.warning("No control samples to process; writing empty report")
            empty_report = {
                "batch_id": self.batch_id,
                "pc_reference_dir": str(self.pc_reference_dir) if self.pc_reference_dir else None,
                "samples": [],
            }
            self.write_report(empty_report)
            return empty_report
        blastn_samples = self.run_blastn()
        expected_samples = self.load_expected()
        report = self.evaluate(blastn_samples, expected_samples, control_ids)
        self.write_report(report)
        return report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PC/NTC quality control: BLASTN-only validation with optional expected-profile comparison.",
    )
    parser.add_argument("-b", "--batch-id", required=True, help="Batch identifier (e.g. 20250220_mtDNA_19)")
    parser.add_argument("--data-dir", default=None, help="Data directory (overrides DATA_DIR)")
    parser.add_argument(
        "--results-dir",
        default=None,
        help="Results root (overrides RESULTS_DIR); validation output goes under <results-root>/validation/<batch_id>/",
    )
    parser.add_argument("--ref-path", default=DEFAULT_REF_PATH, help="Path to rCRS reference FASTA")
    parser.add_argument(
        "--pc-reference-dir", default=None, help="Optional directory of Sequencher TXT expected profiles"
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    processor = PCNTCProcessor(
        batch_id=args.batch_id,
        data_dir=args.data_dir,
        results_dir=args.results_dir,
        ref_path=args.ref_path,
        pc_reference_dir=args.pc_reference_dir,
    )
    try:
        processor.run()
    except Exception as e:  # noqa: BLE001 - top-level CLI guard
        logger.error("PC/NTC QC failed: {}", e)
        # Non-blocking in the pipeline: exit non-zero only on crash so the
        # caller can decide (pipeline.sh treats failure as a warning).
        sys.exit(1)


if __name__ == "__main__":
    main()
