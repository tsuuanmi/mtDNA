"""
Validate Regenerate Results - Simplified Version

Validates JSON files from manual_pipeline.py regenerate_results() and outputs comprehensive JSON report.
This module is integrated into manual_pipeline.py as the final validation step.
"""

import json
from pathlib import Path

from loguru import logger

from src.config import get_settings
from src.core.variants import (
    is_position_in_intervals,
    normalize_position,
    pos_base,
)
from src.tools.sequencher.utils import get_intervals_from_manual, parse_intervals_string


class RegenerateResultsValidator:
    """Validator for regenerate_results output - validates structure and outputs comprehensive JSON report."""

    def __init__(
        self,
        regenerate_dir: Path,
        comparison_dir: Path | None = None,
        strict: bool = False,
    ):
        self.regenerate_dir = Path(regenerate_dir)
        self.comparison_dir = Path(comparison_dir) if comparison_dir else None
        self.strict = strict
        self.batch_id = self.regenerate_dir.name

        # Results storage
        self.batch_data = {}
        self.expected_samples = set()
        self.valid_samples: list[str] = []
        self.invalid_samples: list[dict] = []
        self.missing_jsons: set[str] = set()
        self.sample_validations: dict[str, dict] = {}
        self.errors: list[str] = []
        self.compatibility_issues: list[str] = []

    def validate_batch_json(self) -> bool:
        """Load statistic_fullbatch.json as source of truth for expected samples."""
        batch_path = self.regenerate_dir / "statistic_fullbatch.json"

        if not batch_path.exists():
            self.errors.append("Missing statistic_fullbatch.json")
            logger.error("Batch JSON not found")
            return False

        try:
            with open(batch_path) as f:
                self.batch_data = json.load(f)

            self.expected_samples = set(self.batch_data.keys())
            logger.success(f"Batch JSON loaded: {len(self.expected_samples)} samples")
            return True

        except json.JSONDecodeError as e:
            self.errors.append(f"Invalid JSON: {e!s}")
            return False
        except Exception as e:
            self.errors.append(f"Error loading batch JSON: {e!s}")
            return False

    def validate_sample_json(self, sample_id: str) -> dict:
        """Check if sample JSON file exists and is parseable."""
        result = {"sample_id": sample_id, "valid": False, "errors": []}

        json_path = self.regenerate_dir / sample_id / f"{sample_id}.json"

        if not json_path.exists():
            result["errors"].append("JSON not found")
            self.missing_jsons.add(sample_id)
            return result

        try:
            with open(json_path) as f:
                json.load(f)  # Just check if parseable
            result["valid"] = True

        except json.JSONDecodeError:
            result["errors"].append("Invalid JSON format")
        except Exception as e:
            result["errors"].append(f"Error: {e!s}")

        return result

    def validate_all_samples(self) -> bool:
        """Validate all sample JSONs."""
        for sample_id in sorted(self.expected_samples):
            validation = self.validate_sample_json(sample_id)
            self.sample_validations[sample_id] = validation

            if validation["valid"]:
                self.valid_samples.append(sample_id)
            else:
                self.invalid_samples.append(validation)
                for error in validation["errors"][:2]:
                    logger.error(f"{sample_id}: {error}")

        return len(self.invalid_samples) == 0

    def _validate_intervals_match(self, sample_id: str, manual, regenerated) -> list[str]:
        """Check if manual intervals match regenerated intervals.

        Args:
            sample_id: Sample identifier
            manual: Manual intervals dict with HV1/HV2/HV3 keys (or "FULL REGION" string)
            regenerated: Regenerated intervals dict with HV1/HV2/HV3 keys (or "FULL REGION" string)

        Returns:
            List of issues found (empty if all match)
        """
        issues = []

        # Handle "FULL REGION" case - this means use default full regions
        # If either is "FULL REGION", consider it acceptable (regeneration may expand it)
        if (isinstance(manual, str) and manual == "FULL REGION") or (
            isinstance(regenerated, str) and regenerated == "FULL REGION"
        ):
            return issues

        # Both are dicts - compare region by region
        for region in ["HV1", "HV2", "HV3"]:
            manual_interval = manual.get(region, [])
            regenerated_interval = regenerated.get(region, [])

            # Empty intervals are valid (sequencing failed)
            if manual_interval == regenerated_interval:
                continue

            # Convert to comparable format for better checking
            manual_str = ",".join(manual_interval) if manual_interval else "empty"
            regen_str = ",".join(regenerated_interval) if regenerated_interval else "empty"

            issues.append(
                f"{sample_id}: {region} interval mismatch - manual: [{manual_str}], regenerated: [{regen_str}]",
            )

        return issues

    def validate_manual_compatibility(self) -> bool:
        """Validate compatibility with manual pipeline (if comparison dir exists)."""
        if not self.comparison_dir or not self.comparison_dir.exists():
            return True

        intervals_file = self.comparison_dir / "all_intervals.json"
        if not intervals_file.exists():
            self.compatibility_issues.append("Missing all_intervals.json")
            return False

        try:
            with open(intervals_file) as f:
                manual_intervals = json.load(f)
        except Exception as e:
            self.compatibility_issues.append(f"Failed to load intervals: {e!s}")
            return False

        # Validate intervals consistency
        issues_before = len(self.compatibility_issues)

        for sample_id in self.valid_samples:
            if sample_id not in manual_intervals:
                continue

            sample_path = self.regenerate_dir / sample_id / f"{sample_id}.json"
            try:
                with open(sample_path) as f:
                    data = json.load(f)

                regenerated_intervals = data.get("intervals", {})
                manual_interval_data = manual_intervals[sample_id]

                # Convert manual intervals string to dict format if needed
                if isinstance(manual_interval_data, str) and manual_interval_data != "FULL REGION":
                    # Parse intervals string and organize by region using existing utility function
                    try:
                        manual_interval_data = get_intervals_from_manual(
                            parse_intervals_string(manual_interval_data),
                        )
                    except Exception as e:
                        logger.warning(f"{sample_id}: Could not parse intervals string: {e}")
                        continue

                issues = self._validate_intervals_match(sample_id, manual_interval_data, regenerated_intervals)

                if issues:
                    for issue in issues:
                        logger.warning(issue)
                        self.compatibility_issues.append(issue)

            except Exception as e:
                issue = f"{sample_id}: Failed to check intervals - {e!s}"
                logger.error(issue)
                self.compatibility_issues.append(issue)

        intervals_errors = len(self.compatibility_issues) - issues_before
        if intervals_errors > 0:
            logger.warning(f"{intervals_errors} interval consistency issues found")
        else:
            logger.success("All intervals match manual data")

        # Load diff file if exists
        diff_file = self.comparison_dir / "diff_formatted.tsv"
        manual_edits = {}
        if diff_file.exists():
            try:
                with open(diff_file) as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        headers = lines[0].strip().split("\t")
                        for line in lines[1:]:
                            parts = line.strip().split("\t")
                            if len(parts) >= len(headers):
                                row = dict(zip(headers, parts))
                                sample_id = row.get("Sample", "")
                                if sample_id not in manual_edits:
                                    manual_edits[sample_id] = []
                                manual_edits[sample_id].append(
                                    {
                                        "pos": row.get("POS_REF", ""),
                                        "type": row.get("edited_type", ""),
                                        "seq_manual": row.get("SEQ_manual", ""),
                                    },
                                )
            except Exception as e:
                logger.warning(f"Could not load diff file: {e}")

        # Validate edits were applied
        issues_before = len(self.compatibility_issues)

        for sample_id in self.valid_samples:  # Check all samples
            if sample_id not in manual_edits:
                continue

            sample_path = self.regenerate_dir / sample_id / f"{sample_id}.json"
            try:
                with open(sample_path) as f:
                    data = json.load(f)

                variants = data.get("variants", {})
                all_vars = variants.get("HV1", []) + variants.get("HV2", []) + variants.get("HV3", [])
                var_by_pos = {str(normalize_position(v.get("pos", 0))): v for v in all_vars}

                for edit in manual_edits[sample_id]:
                    pos = str(normalize_position(edit["pos"]))
                    edit_type = edit["type"]
                    seq_manual = edit.get("seq_manual", "")

                    # Skip edits outside genomic regions
                    try:
                        if not is_position_in_intervals(pos_base(pos), list(get_settings().regions.REGIONS.values())):
                            continue
                    except (ValueError, TypeError):
                        continue

                    if edit_type == "ADD" and pos not in var_by_pos:
                        self.compatibility_issues.append(f"{sample_id}: ADD at {pos} not applied")
                    elif edit_type == "REMOVE" and pos in var_by_pos:
                        self.compatibility_issues.append(f"{sample_id}: REMOVE at {pos} not applied")
                    elif edit_type == "EDIT":
                        if pos not in var_by_pos:
                            self.compatibility_issues.append(
                                f"{sample_id}: EDIT at {pos} not applied (variant missing)",
                            )
                        elif seq_manual and var_by_pos[pos].get("seq") != seq_manual:
                            self.compatibility_issues.append(
                                f"{sample_id}: EDIT at {pos} has wrong sequence - "
                                f"expected '{seq_manual}', found '{var_by_pos[pos].get('seq')}'",
                            )

            except Exception as e:
                self.compatibility_issues.append(f"{sample_id}: Failed to validate - {e!s}")
                logger.error(f"{sample_id}: Validation error - {e!s}")

        errors = len(self.compatibility_issues) - issues_before
        if errors > 0:
            logger.warning(f"{errors} compatibility issues")
            return False

        logger.success("Manual compatibility OK")
        return True

    def generate_report(self) -> dict:
        """Generate comprehensive JSON validation report."""
        batch_valid = len(self.batch_data) > 0 and len(self.errors) == 0
        samples_valid = len(self.invalid_samples) == 0
        compatibility_valid = len(self.compatibility_issues) == 0

        overall_pass = batch_valid and samples_valid and compatibility_valid

        # Include valid sample IDs (without detailed flags)
        valid_samples_details = sorted(self.valid_samples)

        report = {
            "batch_id": self.batch_id,
            "regenerate_dir": str(self.regenerate_dir),
            "overall_status": "PASS" if overall_pass else "FAIL",
            "validation_status": {
                "batch_json_valid": batch_valid,
                "all_samples_valid": samples_valid,
                "manual_compatibility_valid": compatibility_valid,
            },
            "statistics": {
                "total_samples": len(self.expected_samples),
                "valid_samples": len(self.valid_samples),
                "invalid_samples": len(self.invalid_samples),
                "missing_json_files": len(self.missing_jsons),
            },
            "samples": {
                "valid": valid_samples_details,
                "invalid": self.invalid_samples,
                "missing_json_files": sorted(list(self.missing_jsons)),
            },
            "issues": {
                "validation_errors": self.errors,
                "compatibility_issues": self.compatibility_issues,
            },
        }

        # Save single validation report directly in regenerate_dir
        json_path = self.regenerate_dir / "validation.json"
        with open(json_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.success(f"Validation report: {json_path}")

        return report

    def run_validation(self) -> bool:
        """Run complete validation workflow."""

        # Step 1: Load batch JSON (source of truth)
        if not self.validate_batch_json():
            logger.error("Batch JSON validation failed")
            self.generate_report()
            return False

        # Step 2: Check sample files exist
        self.validate_all_samples()

        # Step 3: Validate manual compatibility
        self.validate_manual_compatibility()

        # Generate report
        report = self.generate_report()

        # Summary
        if report["overall_status"] == "PASS":
            logger.success(
                f"ALL VALIDATIONS PASSED - "
                f"{report['statistics']['valid_samples']}/{report['statistics']['total_samples']} samples valid",
            )
        else:
            logger.error("VALIDATION FAILED")
            if report["issues"]["validation_errors"]:
                logger.error(f"  Errors: {len(report['issues']['validation_errors'])}")
            if report["issues"]["compatibility_issues"]:
                logger.error(f"  Compatibility issues: {len(report['issues']['compatibility_issues'])}")

        return report["overall_status"] == "PASS"
