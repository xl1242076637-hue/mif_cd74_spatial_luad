#!/usr/bin/env python
"""List or download GEO supplementary files from the project manifest."""

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
    download_file_resumable,
    fetch_geo_directory_listing,
    filter_records,
    write_csv,
    write_jsonl,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=PROJECT_ROOT / "config" / "datasets.yaml",
        help="Path to datasets.yaml.",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Project root containing data/, docs/, and results/.",
    )
    parser.add_argument(
        "--dataset",
        action="append",
        help="GSE accession to process. May be supplied more than once. Defaults to all manifest datasets.",
    )
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="Only list remote files and write a plan. This is the default when --download is absent.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download files in addition to writing the plan.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_download_plan.jsonl",
        help="JSONL output path for the download plan.",
    )
    parser.add_argument(
        "--csv-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_file_inventory.csv",
        help="CSV output path for the download plan.",
    )
    parser.add_argument(
        "--include-regex",
        help="Only keep files whose filename matches this regular expression.",
    )
    parser.add_argument(
        "--fetch-sizes",
        action="store_true",
        help="Fetch remote Content-Length values and add local download status columns to the CSV/JSONL plan.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Per-request timeout in seconds.",
    )
    parser.add_argument(
        "--range-chunk-mb",
        type=int,
        default=64,
        help="Byte-range chunk size in MB for resumable large-file downloads.",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=100,
        help="Maximum retry attempts across byte-range chunks.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = args.dataset or [None]
    all_records: list[dict] = []

    for dataset in selected:
        targets = build_dataset_targets(args.manifest, args.project_root, dataset=dataset)
        for target in targets:
            print(f"Listing {target['accession']} from {target['supplementary_url']}")
            remote_files = fetch_geo_directory_listing(target["supplementary_url"])
            records = build_download_plan(target, remote_files)
            records = filter_records(records, args.include_regex)
            all_records.extend(records)
            print(f"  {target['accession']}: {len(records)} supplementary files")

    if args.fetch_sizes or args.download:
        print("Annotating remote sizes and local download status.")
        all_records = annotate_download_records(all_records, fetch_sizes=True, timeout=args.timeout)

    write_jsonl(all_records, args.output)
    write_csv(all_records, args.csv_output)
    print(f"Wrote JSONL plan: {args.output}")
    print(f"Wrote CSV inventory: {args.csv_output}")

    if args.download:
        for record in all_records:
            local_path = Path(record["local_path"])
            remote_bytes = record.get("remote_bytes")
            local_bytes = local_path.stat().st_size if local_path.exists() else 0
            if remote_bytes is not None and local_bytes == int(remote_bytes):
                print(f"Skipping existing file: {local_path}")
                continue
            size = download_file_resumable(
                record["url"],
                local_path,
                timeout=args.timeout,
                range_chunk_size=args.range_chunk_mb * 1024 * 1024,
                max_attempts=args.max_attempts,
            )
            print(f"Downloaded {local_path} ({size} bytes)")
    else:
        print("Download not requested; rerun with --download to fetch files.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
