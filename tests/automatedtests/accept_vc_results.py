#!/usr/bin/env python3
"""Accept the latest VC test results as the new targets."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile


HEADER = ("test", "max_vc", "max_gradient")


def read_results(path: Path) -> tuple[list[str], dict[str, tuple[str, str]]]:
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != HEADER:
            raise ValueError(f"{path} must have the header {','.join(HEADER)}")

        order: list[str] = []
        results: dict[str, tuple[str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            test_name = row["test"]
            if not test_name or test_name in results:
                raise ValueError(f"Invalid or duplicate test name at {path}:{line_number}")

            max_vc_text = row["max_vc"]
            max_gradient_text = row["max_gradient"]
            try:
                max_vc = float(max_vc_text)
                max_gradient = float(max_gradient_text)
            except ValueError as error:
                raise ValueError(f"Invalid numeric value at {path}:{line_number}") from error
            if not math.isfinite(max_vc) or max_vc <= 0:
                raise ValueError(f"max_vc must be finite and positive at {path}:{line_number}")
            if not math.isfinite(max_gradient) or max_gradient < 0:
                raise ValueError(f"max_gradient must be finite and non-negative at {path}:{line_number}")

            order.append(test_name)
            results[test_name] = (max_vc_text, max_gradient_text)

    return order, results


def accept_results(target_path: Path, results_path: Path) -> None:
    target_order, targets = read_results(target_path)
    _, results = read_results(results_path)

    missing = targets.keys() - results.keys()
    unexpected = results.keys() - targets.keys()
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing results: " + ", ".join(sorted(missing)))
        if unexpected:
            details.append("unexpected results: " + ", ".join(sorted(unexpected)))
        raise ValueError("Cannot accept results; " + "; ".join(details))

    target_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            encoding="utf-8",
            dir=target_path.parent,
            prefix=target_path.name + ".",
            suffix=".tmp",
            delete=False,
        ) as destination:
            temporary_name = destination.name
            writer = csv.writer(destination, lineterminator="\n")
            writer.writerow(HEADER)
            for test_name in target_order:
                writer.writerow((test_name, *results[test_name]))
        os.replace(temporary_name, target_path)
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def main() -> None:
    automated_tests_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--targets",
        type=Path,
        default=automated_tests_dir / "vc_targets.csv",
        help="target CSV to overwrite",
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=automated_tests_dir / "vc_results.csv",
        help="results CSV produced by limatest",
    )
    args = parser.parse_args()

    accept_results(args.targets.resolve(), args.results.resolve())
    print(f"Accepted {args.results} into {args.targets}")


if __name__ == "__main__":
    main()
