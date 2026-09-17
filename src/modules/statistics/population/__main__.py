"""CLI for current approved-cohort population reporting."""

import argparse
from datetime import UTC, datetime
from pathlib import Path

from loguru import logger

from src.config import get_settings
from src.modules.statistics.population import AnalysisOptions, analyze_population, load_cohort, write_reports

DEFAULT_INPUT = Path("/mnt/nas/bca/mtDNA/science/results/merged/merged_statistics.json")


def _default_output_dir() -> Path:
    results = get_settings().directories.results
    return results / "modules" / "statistics" / "population" / f"current-cohort-{datetime.now(UTC):%Y%m%d-%H%M%S}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Current upstream-approved sample-ID-keyed export (default: %(default)s)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="New or empty analysis directory "
        "(default: <results>/modules/statistics/population/current-cohort-<timestamp>)",
    )
    parser.add_argument("--confidence-level", type=float, default=0.95)
    parser.add_argument("--no-confidence-intervals", action="store_true")
    parser.add_argument("--precision", type=int, default=6, help="Executive text decimal places only")
    args = parser.parse_args()
    output_dir = args.output_dir or _default_output_dir()
    logger.info("Population analysis: input {}, output {}", args.input, output_dir)
    try:
        options = AnalysisOptions(
            confidence_level=None if args.no_confidence_intervals else args.confidence_level,
            presentation_precision=args.precision,
        )
        result = analyze_population(load_cohort(args.input), options=options)
        write_reports(result, output_dir)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.exit(1, f"Population analysis failed: {error}\n")
    haplotypes = result.summary.haplotypes
    logger.success(
        "Population analysis complete: {} samples, {} complete / {} incomplete combined profiles, "
        "{} distinct variants, {} distinct haplotypes, RMP {}, DP {}",
        result.summary.n_samples_analyzed,
        result.summary.n_samples_complete_haplotype,
        result.summary.n_samples_incomplete_haplotype,
        result.summary.variants["all"].n_distinct_variants,
        haplotypes.n_distinct_haplotypes,
        haplotypes.empirical_random_match_probability,
        haplotypes.discrimination_power,
    )
    logger.success("Executive summary: {}", (output_dir / "mtDNA_population_executive_summary.txt").resolve())


if __name__ == "__main__":
    main()
