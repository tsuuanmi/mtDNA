#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Any

from loguru import logger


def load_json_data(file_path: Path) -> dict[str, Any]:
    """
    Load data from a JSON file.

    Args:
        file_path: Path to the JSON file

    Returns:
        Dictionary containing the loaded JSON data

    Raises:
        FileNotFoundError: If the file doesn't exist
        json.JSONDecodeError: If the file contains invalid JSON
    """
    try:
        with file_path.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in {file_path}: {e}")
        raise


def merge_json_files(json_files: list[Path], output_file: Path) -> dict[str, Any]:
    """
    Merge multiple JSON files with the same structure.

    Args:
        json_files: List of paths to JSON files to merge
        output_file: Path where to save the merged JSON

    Returns:
        The merged dictionary
    """
    merged_data: dict[str, Any] = {}

    for file_path in json_files:
        logger.info(f"Processing {file_path}")
        try:
            data = load_json_data(file_path)
            # Merge the data
            merged_data.update(data)
        except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
            logger.error(f"Error processing {file_path}: {e}")
            continue

    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Save merged data
    try:
        with output_file.open("w") as f:
            json.dump(merged_data, f, indent=2)
        logger.success(f"Successfully merged {len(json_files)} files to {output_file}")
    except OSError as e:
        logger.error(f"Error writing merged file: {e}")
        raise

    return merged_data


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge multiple JSON files with the same structure")
    parser.add_argument("input_files", nargs="+", help="Input JSON files to merge")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file")

    args = parser.parse_args()

    # Convert input files to Path objects
    input_paths = [Path(f) for f in args.input_files]
    output_path = Path(args.output)

    merge_json_files(input_paths, output_path)


if __name__ == "__main__":
    main()
