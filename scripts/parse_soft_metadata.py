#!/usr/bin/env python
"""Parse downloaded GEO family SOFT files into one sample metadata table."""

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

from luad_niche.soft import collect_sample_metadata  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Project root containing data/raw and results/tables.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata.csv",
        help="Output CSV path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    soft_paths = sorted(args.project_root.glob("data/raw/GSE*/GSE*_family.soft.gz"))
    records = collect_sample_metadata(soft_paths)
    if not records:
        raise SystemExit("No sample metadata records parsed.")

    fieldnames = sorted({key for record in records for key in record})
    preferred = [
        "series_accession",
        "sample_accession",
        "title",
        "source_name",
        "histological_type",
        "radiological_type",
        "group",
        "gender",
        "tissue",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    print(f"Parsed {len(records)} samples from {len(soft_paths)} SOFT files.")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

