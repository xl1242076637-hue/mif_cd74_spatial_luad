#!/usr/bin/env python
"""Add harmonized stage and inclusion labels to parsed GEO sample metadata."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.metadata import annotate_sample_record  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata.csv",
        help="Parsed sample metadata CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
        help="Annotated metadata CSV.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.input.open("r", encoding="utf-8-sig", newline="") as handle:
        records = [annotate_sample_record(row) for row in csv.DictReader(handle)]

    fieldnames = sorted({key for record in records for key in record})
    preferred = [
        "series_accession",
        "sample_accession",
        "title",
        "histological_type",
        "group",
        "interpreted_stage",
        "interpreted_condition",
        "include_in_luad_progression",
        "radiological_type",
        "gender",
        "source_name",
        "tissue",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    print(f"Annotated {len(records)} samples.")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

