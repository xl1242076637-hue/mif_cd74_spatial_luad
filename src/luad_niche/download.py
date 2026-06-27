"""Download planning helpers for GEO supplementary files."""

from __future__ import annotations

import csv
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote, urljoin

import requests
from requests import RequestException
import yaml


def load_manifest(manifest_path: Path) -> dict:
    """Load the dataset manifest from YAML."""
    with Path(manifest_path).open("r", encoding="utf-8") as handle:
        manifest = yaml.safe_load(handle)
    if not isinstance(manifest, dict) or "datasets" not in manifest:
        raise ValueError(f"Manifest {manifest_path} must contain a 'datasets' list.")
    return manifest


def build_dataset_targets(
    manifest_path: Path,
    project_root: Path,
    dataset: str | None = None,
) -> list[dict]:
    """Build deterministic local raw-data targets from the manifest."""
    manifest = load_manifest(Path(manifest_path))
    wanted = dataset.upper() if dataset else None
    targets: list[dict] = []

    for entry in manifest["datasets"]:
        accession = str(entry["accession"]).upper()
        if wanted and accession != wanted:
            continue
        targets.append(
            {
                "accession": accession,
                "supplementary_url": entry["supplementary_url"],
                "raw_dir": Path(project_root) / "data" / "raw" / accession,
                "priority": entry.get("priority"),
                "modality": entry.get("modality"),
                "role": entry.get("role"),
            }
        )

    if wanted and not targets:
        raise ValueError(f"Dataset {wanted} not found in {manifest_path}.")
    return targets


def parse_geo_directory_listing(html: str, base_url: str) -> list[dict]:
    """Extract downloadable file links from an NCBI GEO FTP HTML listing."""
    files: list[dict] = []
    seen: set[str] = set()

    for href in re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.IGNORECASE):
        if href.startswith("?") or href.startswith("../") or href in {"/", ""}:
            continue
        if href.endswith("/"):
            continue
        resolved_url = urljoin(base_url, href)
        if not resolved_url.startswith(base_url):
            continue
        filename = unquote(href.split("/")[-1])
        if not filename or filename in seen:
            continue
        seen.add(filename)
        files.append({"filename": filename, "url": resolved_url})

    return files


def fetch_geo_directory_listing(url: str, timeout: int = 60) -> list[dict]:
    """Fetch and parse an NCBI GEO supplementary directory listing."""
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return parse_geo_directory_listing(response.text, url)


def build_download_plan(target: dict, remote_files: Iterable[dict]) -> list[dict]:
    """Map remote files to deterministic local paths."""
    plan: list[dict] = []
    raw_dir = Path(target["raw_dir"])
    accession = str(target["accession"]).upper()

    for remote in remote_files:
        filename = remote["filename"]
        plan.append(
            {
                "accession": accession,
                "filename": filename,
                "url": remote["url"],
                "local_path": str(raw_dir / filename),
            }
        )
    return plan


def filter_records(records: Iterable[dict], include_regex: str | None = None) -> list[dict]:
    """Filter records by filename using an optional regular expression."""
    records = list(records)
    if not include_regex:
        return records
    pattern = re.compile(include_regex)
    return [record for record in records if pattern.search(record["filename"])]


def write_jsonl(records: Iterable[dict], output_path: Path) -> None:
    """Write records as UTF-8 JSON Lines."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")


def write_csv(records: Iterable[dict], output_path: Path) -> None:
    """Write records as a small UTF-8 CSV table."""
    records = list(records)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    preferred = [
        "accession",
        "filename",
        "url",
        "local_path",
        "remote_bytes",
        "local_bytes",
        "download_complete",
        "download_status",
    ]
    extra = sorted({key for record in records for key in record if key not in preferred})
    fieldnames = [field for field in preferred if any(field in record for record in records)]
    fieldnames.extend(extra)
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def fetch_remote_content_length(url: str, timeout: int = 60) -> int | None:
    """Return the remote Content-Length, or None when the server omits it."""
    response = requests.head(url, allow_redirects=True, timeout=timeout)
    response.raise_for_status()
    content_length = response.headers.get("Content-Length")
    if content_length is None:
        return None
    return int(content_length)


def add_local_download_status(record: dict) -> dict:
    """Add local byte counts and a coarse completion status to one record."""
    annotated = dict(record)
    local_path = Path(str(record["local_path"]))
    local_bytes = local_path.stat().st_size if local_path.exists() else 0
    remote_bytes = annotated.get("remote_bytes")
    complete = bool(remote_bytes is not None and local_bytes == int(remote_bytes))
    if complete:
        status = "complete"
    elif local_bytes > 0:
        status = "partial"
    else:
        status = "missing"

    annotated.update(
        {
            "local_bytes": local_bytes,
            "download_complete": complete,
            "download_status": status,
        }
    )
    return annotated


def annotate_download_records(records: Iterable[dict], fetch_sizes: bool = False, timeout: int = 60) -> list[dict]:
    """Add optional remote sizes and local download status to records."""
    annotated: list[dict] = []
    for record in records:
        enriched = dict(record)
        if fetch_sizes and "remote_bytes" not in enriched:
            enriched["remote_bytes"] = fetch_remote_content_length(enriched["url"], timeout=timeout)
        annotated.append(add_local_download_status(enriched))
    return annotated


def build_byte_ranges(total_bytes: int, range_chunk_size: int) -> list[tuple[int, int, int]]:
    """Build inclusive byte ranges as (index, start, end)."""
    if total_bytes <= 0:
        raise ValueError("total_bytes must be positive.")
    if range_chunk_size <= 0:
        raise ValueError("range_chunk_size must be positive.")
    ranges: list[tuple[int, int, int]] = []
    start = 0
    index = 0
    while start < total_bytes:
        end = min(start + range_chunk_size - 1, total_bytes - 1)
        ranges.append((index, start, end))
        start = end + 1
        index += 1
    return ranges


def download_file(url: str, local_path: Path, timeout: int = 120) -> int:
    """Download one file and return its final byte size."""
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with local_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return local_path.stat().st_size


def download_byte_range(
    url: str,
    part_path: Path,
    start: int,
    end: int,
    timeout: int = 120,
    chunk_size: int = 1024 * 1024,
    max_attempts: int = 20,
) -> int:
    """Download one inclusive byte range to a complete part file."""
    part_path = Path(part_path)
    part_path.parent.mkdir(parents=True, exist_ok=True)
    expected_bytes = end - start + 1
    if part_path.exists() and part_path.stat().st_size == expected_bytes:
        return expected_bytes

    tmp_path = part_path.with_suffix(part_path.suffix + ".tmp")
    last_error: Exception | None = None
    for _ in range(max_attempts):
        if tmp_path.exists():
            tmp_path.unlink()
        try:
            headers = {"Range": f"bytes={start}-{end}"}
            with requests.get(url, stream=True, timeout=timeout, headers=headers) as response:
                response.raise_for_status()
                if response.status_code != 206:
                    raise RuntimeError(f"Server did not honor byte range for {url}: HTTP {response.status_code}")
                with tmp_path.open("wb") as handle:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            handle.write(chunk)
            downloaded_bytes = tmp_path.stat().st_size
            if downloaded_bytes != expected_bytes:
                raise RuntimeError(
                    f"Downloaded range size mismatch for {part_path}: "
                    f"{downloaded_bytes} != {expected_bytes}"
                )
            tmp_path.replace(part_path)
            return expected_bytes
        except (OSError, RequestException, RuntimeError) as error:
            last_error = error
            continue
    if last_error is not None:
        raise last_error
    raise RuntimeError(f"Could not download byte range {start}-{end} for {url}")


def download_file_by_parallel_ranges(
    url: str,
    local_path: Path,
    timeout: int = 120,
    workers: int = 8,
    range_chunk_size: int = 8 * 1024 * 1024,
    max_attempts: int = 20,
) -> int:
    """Download a file through parallel HTTP Range requests and assemble it."""
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)
    remote_bytes = fetch_remote_content_length(url, timeout=timeout)
    if remote_bytes is None:
        raise ValueError(f"Remote Content-Length is required for parallel range download: {url}")
    if local_path.exists() and local_path.stat().st_size == remote_bytes:
        return remote_bytes

    part_dir = local_path.parent / f"{local_path.name}.parts"
    ranges = build_byte_ranges(remote_bytes, range_chunk_size)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                download_byte_range,
                url,
                part_dir / f"{index:05d}.part",
                start,
                end,
                timeout=timeout,
                max_attempts=max_attempts,
            ): (index, start, end)
            for index, start, end in ranges
        }
        for future in as_completed(futures):
            index, start, end = futures[future]
            bytes_downloaded = future.result()
            print(f"Downloaded part {index + 1}/{len(ranges)} ({start}-{end}, {bytes_downloaded} bytes)")

    tmp_output = local_path.with_suffix(local_path.suffix + ".assembling")
    if tmp_output.exists():
        tmp_output.unlink()
    with tmp_output.open("wb") as output_handle:
        for index, _start, _end in ranges:
            part_path = part_dir / f"{index:05d}.part"
            with part_path.open("rb") as part_handle:
                for chunk in iter(lambda: part_handle.read(1024 * 1024), b""):
                    output_handle.write(chunk)
    if tmp_output.stat().st_size != remote_bytes:
        raise RuntimeError(f"Assembled file size mismatch: {tmp_output.stat().st_size} != {remote_bytes}")
    tmp_output.replace(local_path)
    return local_path.stat().st_size


def download_file_resumable(
    url: str,
    local_path: Path,
    timeout: int = 120,
    chunk_size: int = 1024 * 1024,
    range_chunk_size: int = 64 * 1024 * 1024,
    max_attempts: int = 100,
) -> int:
    """Download a file, resuming from the current local byte count when possible."""
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)
    remote_bytes = fetch_remote_content_length(url, timeout=timeout)

    if remote_bytes is None:
        last_error: Exception | None = None
        for _ in range(max_attempts):
            try:
                return download_file(url, local_path, timeout=timeout)
            except RequestException as error:
                last_error = error
        if last_error is not None:
            raise last_error
        return local_path.stat().st_size

    attempts = 0
    while True:
        existing_bytes = local_path.stat().st_size if local_path.exists() else 0
        if existing_bytes == remote_bytes:
            return existing_bytes
        if existing_bytes > remote_bytes:
            raise ValueError(f"Local file is larger than remote file: {local_path}")
        if attempts >= max_attempts:
            raise RuntimeError(
                f"Could not finish {url} after {max_attempts} attempts "
                f"({existing_bytes}/{remote_bytes} bytes)."
            )

        end_byte = min(existing_bytes + range_chunk_size - 1, remote_bytes - 1)
        headers = {"Range": f"bytes={existing_bytes}-{end_byte}"}
        try:
            with requests.get(url, stream=True, timeout=timeout, headers=headers) as response:
                response.raise_for_status()
                if response.status_code != 206:
                    raise RuntimeError(f"Server did not honor byte range for {url}: HTTP {response.status_code}")
                with local_path.open("ab") as handle:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        if chunk:
                            handle.write(chunk)
                if local_path.stat().st_size == existing_bytes:
                    attempts += 1
        except RequestException:
            attempts += 1
            last_size = local_path.stat().st_size if local_path.exists() else 0
            if last_size == existing_bytes:
                continue
            continue
