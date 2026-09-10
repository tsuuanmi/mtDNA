"""Split a FIS-to-Sanger comparison report into the three defined sample sets."""

from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd
from loguru import logger

SET_COLUMNS = ["Sample_ID", "Sanger_Batch", "Set"]
SET_VALUES = ("1", "2", "3")


def arg_parser() -> Namespace:
    """Parse FIS-to-Sanger set-splitting arguments."""
    parser = ArgumentParser(description="Split a FIS-to-Sanger comparison TSV into three sample-set reports")
    parser.add_argument("-i", "--input-file", type=Path, required=True, help="FIS-to-Sanger comparison TSV")
    parser.add_argument("-s", "--sets-file", type=Path, required=True, help="Headerless Sample_ID/batch/set TSV")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Directory for set reports; defaults to the input TSV directory",
    )
    return parser.parse_args()


def load_sets(sets_file: Path) -> dict[str, str]:
    """Load the headerless Sample_ID, Sanger_Batch, Set definitions."""
    sets = pd.read_csv(sets_file, sep="\t", header=None, names=SET_COLUMNS, dtype=str, keep_default_na=False)
    if sets.empty:
        msg = f"Sets file is empty: {sets_file}"
        raise ValueError(msg)
    invalid_sets = sorted(set(sets["Set"]) - set(SET_VALUES))
    if invalid_sets:
        msg = f"Sets file has unsupported set values: {', '.join(invalid_sets)}"
        raise ValueError(msg)
    if (sets["Sample_ID"] == "").any():
        msg = f"Sets file has empty sample IDs: {sets_file}"
        raise ValueError(msg)
    duplicates = sets.loc[sets["Sample_ID"].duplicated(), "Sample_ID"].tolist()
    if duplicates:
        msg = f"Sets file has duplicate sample IDs: {', '.join(duplicates[:5])}"
        raise ValueError(msg)
    return dict(zip(sets["Sample_ID"], sets["Set"], strict=True))


def resolve_set_id(fis_sample: str, set_ids: set[str]) -> str | None:
    """Resolve the longest set ID matching a complete FIS sample-ID prefix."""
    matches = [set_id for set_id in set_ids if fis_sample == set_id or fis_sample.startswith(f"{set_id}_")]
    return max(matches, key=len) if matches else None


def split_fis_sanger_comparison_sets(input_file: Path, sets_file: Path, output_dir: Path | None = None) -> None:
    """Write TSV and Excel reports for sets 1, 2, and 3 without altering the source report."""
    report = pd.read_csv(input_file, sep="\t", dtype=str, keep_default_na=False)
    if "FIS_Sample" not in report.columns:
        msg = f"Comparison file {input_file} is missing FIS_Sample"
        raise ValueError(msg)

    sample_to_set = load_sets(sets_file)
    set_ids = set(sample_to_set)
    fis_samples: list[str] = report["FIS_Sample"].astype(str).tolist()
    resolved_ids = [resolve_set_id(sample, set_ids) for sample in fis_samples]
    report_sets = pd.Series(
        [sample_to_set.get(sample_id) if sample_id else None for sample_id in resolved_ids],
        index=report.index,
        dtype="string",
    )
    output_directory = output_dir or input_file.parent
    output_directory.mkdir(parents=True, exist_ok=True)

    for set_value in SET_VALUES:
        subset = report.loc[report_sets == set_value]
        output_file = output_directory / f"{input_file.stem}_set_{set_value}.tsv"
        subset.to_csv(output_file, sep="\t", index=False)
        subset.to_excel(output_file.with_suffix(".xlsx"), index=False)
        logger.success("Saved set {} comparison: {} rows to {}", set_value, len(subset), output_file)

    unresolved = len(report_sets[report_sets.isna()])
    if unresolved:
        logger.warning("Excluded {} comparison rows not assigned by {}", unresolved, sets_file)


def main() -> None:
    """Run FIS-to-Sanger set splitting."""
    args = arg_parser()
    split_fis_sanger_comparison_sets(args.input_file, args.sets_file, args.output_dir)


if __name__ == "__main__":
    main()
