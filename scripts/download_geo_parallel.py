#!/usr/bin/env python
"""Download GEO supplementary files with parallel HTTP byte-range requests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.download import (  # noqa: E402
    annotate_download_records,
    build_dataset_targets,
    build_download_plan,
    download_file_by_parallel_ranges,
    fetch_geo_directory_listing,
    filter_records,
    write_csv,
    write_jsonl,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=PROJECT_ROOT / "config" / "datasets.yaml")
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--dataset", action="append", required=True)
    parser.add_argument("--include-regex", required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--range-chunk-mb", type=int, default=8)
    parser.add_argument("--max-attempts", type=int, default=20)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_parallel_download_plan.jsonl",
    )
    parser.add_argument(
        "--csv-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_parallel_file_inventory.csv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    records: list[dict] = []
    for dataset in args.dataset:
        targets = build_dataset_targets(args.manifest, args.project_root, dataset=dataset)
        for target in targets:
            print(f"Listing {target['accession']} from {target['supplementary_url']}")
            remote_files = fetch_geo_directory_listing(target["supplementary_url"])
            dataset_records = filter_records(build_download_plan(target, remote_files), args.include_regex)
            records.extend(dataset_records)
            print(f"  {target['accession']}: {len(dataset_records)} matching files")

    records = annotate_download_records(records, fetch_sizes=True, timeout=args.timeout)
    write_jsonl(records, args.output)
    write_csv(records, args.csv_output)
    print(f"Wrote JSONL plan: {args.output}")
    print(f"Wrote CSV inventory: {args.csv_output}")

    for record in records:
        local_path = Path(record["local_path"])
        print(f"Parallel download: {record['accession']} {record['filename']}")
        size = download_file_by_parallel_ranges(
            record["url"],
            local_path,
            timeout=args.timeout,
            workers=args.workers,
            range_chunk_size=args.range_chunk_mb * 1024 * 1024,
            max_attempts=args.max_attempts,
        )
        print(f"Downloaded {local_path} ({size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
