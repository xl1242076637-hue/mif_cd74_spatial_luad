"""Small parser for GEO family SOFT sample metadata."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterable


def normalize_characteristic_key(key: str) -> str:
    """Normalize a GEO characteristic key to a stable snake_case column name."""
    normalized = key.strip().lower()
    if normalized == "histolgical type":
        normalized = "histological type"
    normalized = re.sub(r"[^a-z0-9]+", "_", normalized).strip("_")
    return normalized


def parse_soft_samples(text: str, series_accession: str) -> list[dict]:
    """Parse sample-level metadata from GEO SOFT text."""
    samples: list[dict] = []
    current: dict | None = None

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.startswith("^SAMPLE = "):
            if current:
                samples.append(current)
            sample_accession = line.split(" = ", 1)[1].strip()
            current = {
                "series_accession": series_accession,
                "sample_accession": sample_accession,
            }
            continue
        if current is None:
            continue
        if line.startswith("!Sample_title = "):
            current["title"] = line.split(" = ", 1)[1].strip()
        elif line.startswith("!Sample_geo_accession = "):
            current["sample_accession"] = line.split(" = ", 1)[1].strip()
        elif line.startswith("!Sample_source_name_ch1 = "):
            current["source_name"] = line.split(" = ", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split(" = ", 1)[1].strip()
            if ":" not in value:
                continue
            key, characteristic_value = value.split(":", 1)
            current[normalize_characteristic_key(key)] = characteristic_value.strip()

    if current:
        samples.append(current)
    return samples


def read_soft_samples(path: Path, series_accession: str | None = None) -> list[dict]:
    """Read a gzipped or plain SOFT file and parse sample metadata."""
    path = Path(path)
    inferred_series = series_accession or path.name.split("_family", 1)[0]
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    else:
        text = path.read_text(encoding="utf-8", errors="replace")
    return parse_soft_samples(text, inferred_series)


def collect_sample_metadata(paths: Iterable[Path]) -> list[dict]:
    """Parse and concatenate metadata from multiple SOFT files."""
    records: list[dict] = []
    for path in paths:
        records.extend(read_soft_samples(Path(path)))
    return records

