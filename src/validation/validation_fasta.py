"""
FASTA Validation — single-path validation with alignment report.

Validates generated FASTA files by re-running the canonical
generate_sequence code path (via the standard Sample) and comparing output.
Produces a text report with position-level alignment diffs.
"""

import sys
from argparse import ArgumentParser
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from loguru import logger

from src.config import get_settings
from src.core.sample import generate_sequence, load_sample
from src.core.variants import pos_base, pos_sort_key

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class AlignmentData:
    """Position-level alignment data for a single region."""

    ref_aligned: str
    alt_aligned: str
    pos_markers: list[int]
    diff_line: str
    variants_applied: list[dict]


@dataclass
class RegionResult:
    """Validation result for a single HV region within a sample."""

    region: str
    start: int
    end: int
    passed: bool
    expected_seq: str
    actual_seq: str
    ref_aligned: str = ""
    alt_aligned: str = ""
    diff_line: str = ""
    pos_markers: list[int] = field(default_factory=list)
    variants_applied: list[dict] = field(default_factory=list)
    error: str = ""


@dataclass
class SampleResult:
    """Validation result for a single sample."""

    sample_id: str
    regions: dict[str, RegionResult] = field(default_factory=dict)
    passed: bool = True
    error: str = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _merge_con_dict(con_dict: dict) -> dict:
    """Merge con_dict integer and string keys into a unified dict.

    con_dict has two kinds of keys:
      - int keys from the reference genome (original values)
      - str keys from variant application (variant values)

    For SNP and deletion positions, both ``16129`` (int, ref base) and
    ``"16129"`` (str, variant base) coexist.  The variant value must win.

    Insertion keys like ``"217.1"`` are kept as string keys; they are
    handled separately in ``_build_alignment``.

    Returns a new dict with int keys for base positions and string keys
    only for insertion positions.
    """
    merged: dict = {}
    for key, value in con_dict.items():
        if isinstance(key, str) and "." in key:
            # Insertion position — keep as string key
            merged[key] = value
        elif isinstance(key, int):
            # Integer key — store as int
            merged[key] = value
        elif isinstance(key, str):
            # String key that is NOT an insertion (e.g. "16129")
            # Override the integer key for this base position.
            merged[int(key)] = value
        else:
            merged[key] = value
    return merged


# ---------------------------------------------------------------------------
# FastaValidator
# ---------------------------------------------------------------------------


class FastaValidator:
    """Validate generated FASTA files against the canonical code path."""

    def __init__(self, fasta_dir: str | Path, json_dir: str | Path, ref_path: str, batch_id: str) -> None:
        self.fasta_dir = Path(fasta_dir)
        self.json_dir = Path(json_dir)
        self.ref_path = ref_path
        self.batch_id = batch_id
        self.settings = get_settings()
        self.results: list[SampleResult] = []

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def validate_all(self) -> bool:
        """Validate every sample in json_dir. Return True if all pass."""
        if not self.json_dir.exists():
            logger.error(f"JSON directory not found: {self.json_dir}")
            return False

        sample_ids = sorted(d.name for d in self.json_dir.iterdir() if d.is_dir())
        if not sample_ids:
            logger.warning(f"No sample directories found in {self.json_dir}")
            return True

        for sample_id in sample_ids:
            result = self.validate_sample(sample_id)
            self.results.append(result)

        all_pass = all(r.passed for r in self.results)
        logger.info(
            f"Validation {'PASSED' if all_pass else 'FAILED'}: "
            f"{sum(r.passed for r in self.results)}/{len(self.results)} samples pass",
        )
        return all_pass

    def validate_sample(self, sample_id: str) -> SampleResult:
        """Validate a single sample by replaying the canonical FASTA generation path."""
        result = SampleResult(sample_id=sample_id)

        try:
            # 1. Load JSON (same as generate_fasta.py)
            json_path = self.json_dir / sample_id / f"{sample_id}.json"
            if not json_path.exists():
                result.passed = False
                result.error = f"JSON file not found: {json_path}"
                return result

            sample_obj = load_sample(json_path)

            # 2. Replay the canonical sequence generation (generate_sequence)
            seq_result = generate_sequence(sample_obj, str(self.ref_path))
            hv_seqs = seq_result.hv_seqs
            ref_dict = seq_result.ref_dict
            con_dict = seq_result.con_dict

            # 3. Build expected FASTA content from the canonical HV sequences
            save_regions = {
                "HV1": "".join(hv_seqs.get("HV1", {}).get("seq", [])),
                "HV2": "".join(hv_seqs.get("HV2", {}).get("seq", [])),
                "HV3": "".join(hv_seqs.get("HV3", {}).get("seq", [])),
            }

            # 4. Parse actual FASTA
            fasta_path = self.fasta_dir / f"{sample_id.split('.', maxsplit=1)[0]}.fasta"
            if not fasta_path.exists():
                result.passed = False
                result.error = f"FASTA file not found: {fasta_path}"
                return result

            actual_regions = self._parse_fasta(fasta_path)

            # 5. Merge con_dict so that string keys (variants) override int keys (ref)
            merged_con = _merge_con_dict(con_dict)

            # 6. Compare and build alignment per region
            regions = self.settings.regions.REGIONS
            for region_name, bounds in regions.items():
                start, end = bounds
                expected = save_regions.get(region_name, "")
                actual = actual_regions.get(region_name, "")

                region_result = RegionResult(
                    region=region_name,
                    start=start,
                    end=end,
                    passed=(expected == actual),
                    expected_seq=expected,
                    actual_seq=actual,
                )

                if expected != actual:
                    logger.warning(f"Sample {sample_id} {region_name}: FASTA mismatch")

                # Build alignment report regardless of pass/fail
                alignment = self._build_alignment(
                    ref_dict,
                    merged_con,
                    hv_seqs[region_name]["ranges"],
                )
                region_result.ref_aligned = alignment.ref_aligned
                region_result.alt_aligned = alignment.alt_aligned
                region_result.pos_markers = alignment.pos_markers
                region_result.variants_applied = alignment.variants_applied
                region_result.diff_line = alignment.diff_line

                result.regions[region_name] = region_result

            result.passed = all(r.passed for r in result.regions.values())

        except (ValueError, KeyError, OSError) as exc:
            logger.error(f"Validation error for {sample_id}: {exc}")
            result.passed = False
            result.error = str(exc)

        return result

    # ------------------------------------------------------------------
    # FASTA parsing
    # ------------------------------------------------------------------

    def _parse_fasta(self, fasta_path: Path) -> dict[str, str]:
        """Parse a FASTA file. Region headers are plain >REGION_NAME format."""
        regions: dict[str, str] = {}
        current_region = None
        current_seq: list[str] = []

        with fasta_path.open() as f:
            for raw_line in f:
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if stripped.startswith(">"):
                    # Save previous region
                    if current_region is not None:
                        regions[current_region] = "".join(current_seq)
                    # Parse header — plain region name (e.g., >HV1, >HV2, >HV3)
                    header = stripped[1:]  # strip '>'
                    # Handle legacy >HV2|68-340 format for backward compat
                    current_region = header.split("|")[0]
                    current_seq = []
                else:
                    current_seq.append(stripped)

        # Save last region
        if current_region is not None:
            regions[current_region] = "".join(current_seq)

        return regions

    # ------------------------------------------------------------------
    # Alignment construction (single-path from canonical dicts)
    # ------------------------------------------------------------------

    @staticmethod
    def _classify_position(pos: int, ref_base: str, con_base: str) -> dict | None:
        """Classify a single position as deletion, unread, SNP, or None."""
        if con_base == "-" and ref_base != "-":
            return {"pos": pos, "ref": ref_base, "alt": con_base, "type": "deletion"}
        if con_base == "N" and ref_base != "N":
            return {"pos": pos, "ref": ref_base, "alt": con_base, "type": "unread"}
        if con_base != ref_base and con_base not in ("-", "N") and ref_base != "-":
            return {"pos": pos, "ref": ref_base, "alt": con_base, "type": "SNP"}
        return None

    @staticmethod
    def _emit_insertions(
        insertion_keys: list[str],
        merged_con: dict,
        ref_aligned: list[str],
        alt_aligned: list[str],
        pos_markers: list[int],
        variants_applied: list[dict],
    ) -> None:
        """Append insertion positions to the alignment lists."""
        inserted_bases: list[str] = []
        for ins_key in insertion_keys:
            ins_base = merged_con[ins_key]
            inserted_bases.append(ins_base)
            ref_aligned.append("-")
            alt_aligned.append(ins_base)
            pos_markers.append(0)  # insertion gap — no rCRS position

        variants_applied.append(
            {
                "pos_start": insertion_keys[0],
                "pos_end": insertion_keys[-1],
                "ref": "-",
                "alt": "".join(inserted_bases),
                "type": "insertion",
                "bases": len(inserted_bases),
            },
        )

    def _build_alignment(
        self,
        ref_dict: dict,
        merged_con: dict,
        extraction_ranges: list[list[int]],
    ) -> AlignmentData:
        """Build position-level alignment from the canonical ref_dict and merged con_dict.

        *merged_con* is the con_dict after ``_merge_con_dict`` — integer keys
        for base positions (with variant values winning over reference), and
        string keys only for insertion positions like ``"217.1"``.

        *extraction_ranges* is a list of [start, end] intervals (from
        hv_seqs[region]["ranges"]) defining the actual positions included
        in the FASTA content.

        Returns AlignmentData with ref_aligned, alt_aligned, pos_markers,
        diff_line, and variants_applied.
        """
        # Sort extraction ranges by start position
        ranges = sorted(extraction_ranges, key=lambda r: r[0])

        ref_aligned: list[str] = []
        alt_aligned: list[str] = []
        pos_markers: list[int] = []  # rCRS position at each column (0 = insertion gap)
        variants_applied: list[dict] = []

        for range_start, range_end in ranges:
            for pos in range(range_start, range_end + 1):
                # Collect insertion keys after this base position
                insertion_keys = sorted(
                    [k for k in merged_con if isinstance(k, str) and "." in k and pos_base(k) == pos],
                    key=pos_sort_key,
                )

                # Emit the base position
                ref_base = ref_dict.get(pos, "N")
                con_base = merged_con.get(pos, "N")

                ref_aligned.append(ref_base)
                alt_aligned.append(con_base)
                pos_markers.append(pos)

                # Classify this position
                variant_entry = self._classify_position(pos, ref_base, con_base)
                if variant_entry is not None:
                    variants_applied.append(variant_entry)

                # Emit insertion positions after this base
                if insertion_keys:
                    self._emit_insertions(
                        insertion_keys,
                        merged_con,
                        ref_aligned,
                        alt_aligned,
                        pos_markers,
                        variants_applied,
                    )

        # Build diff line
        diff_chars: list[str] = []
        for ref_c, alt_c in zip(ref_aligned, alt_aligned, strict=False):
            if ref_c == alt_c:
                diff_chars.append(".")
            elif ref_c == "-":
                diff_chars.append("^")  # insertion
            elif alt_c == "-":
                diff_chars.append("D")  # deletion
            elif alt_c == "N":
                diff_chars.append("N")  # unread
            else:
                diff_chars.append(ref_c)  # SNP — show actual ref base

        return AlignmentData(
            ref_aligned="".join(ref_aligned),
            alt_aligned="".join(alt_aligned),
            pos_markers=pos_markers,
            diff_line="".join(diff_chars),
            variants_applied=variants_applied,
        )

    # ------------------------------------------------------------------
    # Report formatting
    # ------------------------------------------------------------------

    def _format_report(self) -> str:
        """Format the full text validation report."""
        lines: list[str] = []
        sep = "=" * 60
        now = datetime.now(tz=UTC).strftime("%Y-%m-%d %H:%M:%S")

        lines.append(sep)
        lines.append(f"FASTA Validation Report — Batch: {self.batch_id}")
        lines.append(f"Reference: {self.ref_path}")
        lines.append(f"Date: {now}")
        lines.append(sep)

        for sample_result in self.results:
            lines.append("")
            lines.append(f"--- Sample: {sample_result.sample_id} ---")

            if sample_result.error and not sample_result.regions:
                lines.append(f"  ERROR: {sample_result.error}")
                continue

            for region_name in ["HV1", "HV2", "HV3"]:
                rr = sample_result.regions.get(region_name)
                if rr is None:
                    continue
                self._format_region_report(lines, rr)

        # Summary
        pass_count = sum(1 for r in self.results if r.passed)
        total_count = len(self.results)
        lines.append("")
        lines.append(sep)
        lines.append(f"Summary: {pass_count}/{total_count} samples PASS")
        lines.append(sep)

        return "\n".join(lines)

    @staticmethod
    def _format_region_report(lines: list[str], rr: RegionResult) -> None:
        """Append a single region's report section to *lines*."""
        status = "PASS" if rr.passed else "FAIL"
        lines.append("")
        lines.append(f"## {rr.region} ({rr.start}-{rr.end}) — {status}")

        # Variant summary
        snp_count = sum(1 for v in rr.variants_applied if v["type"] == "SNP")
        ins_list = [v for v in rr.variants_applied if v["type"] == "insertion"]
        ins_count = len(ins_list)
        ins_bases = sum(v.get("bases", 0) for v in ins_list)
        del_count = sum(1 for v in rr.variants_applied if v["type"] == "deletion")
        unread_count = sum(1 for v in rr.variants_applied if v["type"] == "unread")

        lines.append(
            f"  Variants: {snp_count} SNP, {ins_count} insertion "
            f"({ins_bases} bases), {del_count} deletion, {unread_count} unread",
        )
        lines.append("  List:")

        for v in rr.variants_applied:
            if v["type"] == "SNP":
                lines.append(f"    {v['pos']}  {v['ref']}>{v['alt']}   SNP")
            elif v["type"] == "insertion":
                lines.append(
                    f"    {v['pos_start']}-{v['pos_end']}  ->  {v['alt']}   insertion ({v['bases']} bases)",
                )
            elif v["type"] == "deletion":
                lines.append(f"    {v['pos']}  {v['ref']}>-   deletion")
            elif v["type"] == "unread":
                lines.append(f"    {v['pos']}  {v['ref']}>N   unread")

        # Alignment block with line wrapping at ~70 chars
        if rr.ref_aligned:
            FastaValidator._format_alignment_block(lines, rr)

    @staticmethod
    def _format_alignment_block(lines: list[str], rr: RegionResult) -> None:
        """Append wrapped alignment lines (pos/ref/alt/dif) to *lines*."""
        width = 70
        total_len = len(rr.ref_aligned)

        for offset in range(0, total_len, width):
            chunk_end = min(offset + width, total_len)
            chunk_len = chunk_end - offset

            # Build position ruler: show rCRS position numbers at regular intervals
            # Place the number so its LAST digit aligns with the column
            ruler = [" "] * chunk_len
            last_placed = -20  # ensure first label can be placed
            min_label_spacing = 8
            for i in range(chunk_len):
                pm = rr.pos_markers[offset + i]
                if pm != 0 and (i - last_placed) >= min_label_spacing:
                    label = str(pm)
                    # Right-align: last digit at column i
                    for j, ch in enumerate(label):
                        ci = i - len(label) + 1 + j
                        if 0 <= ci < chunk_len and ruler[ci] == " ":
                            ruler[ci] = ch
                    last_placed = i

            ref_chunk = rr.ref_aligned[offset:chunk_end]
            alt_chunk = rr.alt_aligned[offset:chunk_end]
            dif_chunk = rr.diff_line[offset:chunk_end]

            lines.append("")
            lines.append(f"  pos   {''.join(ruler)}")
            lines.append(f"  ref   {ref_chunk}")
            lines.append(f"  alt   {alt_chunk}")
            lines.append(f"  dif   {dif_chunk}")

    # ------------------------------------------------------------------
    # Report writing
    # ------------------------------------------------------------------

    def write_report(self, output_dir: str | Path) -> Path:
        """Write report to results/validation/validation_fasta_{batch_id}.txt."""
        output_dir = Path(output_dir)
        report_dir = output_dir / "validation"
        report_dir.mkdir(parents=True, exist_ok=True)
        report_path = report_dir / f"validation_fasta_{self.batch_id}.txt"

        report_text = self._format_report()
        report_path.write_text(report_text, encoding="utf-8")
        logger.info(f"Report written to {report_path}")
        return report_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point for FASTA validation."""
    parser = ArgumentParser(description="Validate generated FASTA files")
    parser.add_argument("-f", "--fasta_dir", required=True, help="FASTA output directory")
    parser.add_argument("-j", "--json_dir", required=True, help="Source JSON directory")
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA path")
    parser.add_argument("-b", "--batch_id", required=True, help="Batch identifier")
    args = parser.parse_args()

    validator = FastaValidator(
        fasta_dir=args.fasta_dir,
        json_dir=args.json_dir,
        ref_path=args.reference,
        batch_id=args.batch_id,
    )
    all_pass = validator.validate_all()
    validator.write_report(Path("results"))
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
