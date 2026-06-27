#!/usr/bin/env python
"""Download GEO family SOFT metadata files for datasets in the project manifest."""

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
    add_local_download_status,
    download_file_resumable,
    fetch_remote_content_length,
    load_manifest,
    write_csv,
    write_jsonl,
)
from luad_niche.geo import geo_family_soft_url  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=PROJECT_ROOT / "config" / "datasets.yaml",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
    )
    parser.add_argument(
        "--dataset",
        action="append",
        help="GSE accession to process. May be supplied more than once. Defaults to all manifest datasets.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download SOFT files after writing the status table.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_soft_download_plan.jsonl",
    )
    parser.add_argument(
        "--csv-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_soft_file_inventory.csv",
    )
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--range-chunk-mb", type=int, default=64)
    parser.add_argument("--max-attempts", type=int, default=100)
    return parser.parse_args()


def build_records(manifest_path: Path, project_root: Path, datasets: list[str] | None, timeout: int) -> list[dict]:
    manifest = load_manifest(manifest_path)
    wanted = {dataset.upper() for dataset in datasets} if datasets else None
    records: list[dict] = []
    for entry in manifest["datasets"]:
        accession = str(entry["accession"]).upper()
        if wanted and accession not in wanted:
            continue
        url = geo_family_soft_url(accession)
        local_path = project_root / "data" / "raw" / accession / f"{accession}_family.soft.gz"
        record = {
            "accession": accession,
            "filename": f"{accession}_family.soft.gz",
            "url": url,
            "local_path": str(local_path),
            "remote_bytes": fetch_remote_content_length(url, timeout=timeout),
        }
        records.append(add_local_download_status(record))
    if wanted:
        found = {record["accession"] for record in records}
        missing = sorted(wanted - found)
        if missing:
            raise ValueError(f"Dataset(s) not found in manifest: {', '.join(missing)}")
    return records


def main() -> int:
    args = parse_args()
    records = build_records(args.manifest, args.project_root, args.dataset, args.timeout)
    write_jsonl(records, args.output)
    write_csv(records, args.csv_output)
    print(f"Wrote JSONL plan: {args.output}")
    print(f"Wrote CSV inventory: {args.csv_output}")

    if args.download:
        for record in records:
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
        print("Download not requested; rerun with --download to fetch SOFT files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
