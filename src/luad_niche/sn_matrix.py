"""Readers for tabular single-nucleus count matrices."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


GSE308103_FILENAME_RE = re.compile(r"^(GSM\d+)_(.+)\.raw_counts\.mtx\.txt\.gz$")


def parse_gse308103_filename(path: Path) -> dict[str, str]:
    """Parse sample accession and sample name from a GSE308103 raw-count filename."""
    match = GSE308103_FILENAME_RE.match(Path(path).name)
    if not match:
        raise ValueError(f"Unexpected GSE308103 filename: {path}")
    return {"sample_accession": match.group(1), "sample_name": match.group(2)}


def read_tabular_count_barcodes(path: Path) -> list[str]:
    """Read the first-line cell barcodes from a gzipped tabular count matrix."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
    if not header:
        raise ValueError(f"Empty count matrix: {path}")
    return header.split("\t")


def read_tabular_selected_genes(path: Path, genes: list[str]) -> pd.DataFrame:
    """Read selected gene counts without computing full per-cell totals.

    Use this when library sizes are already available from a previous pass.
    Unrequested rows are skipped before numeric parsing, which is much faster
    for targeted follow-up panels.
    """
    requested = list(dict.fromkeys(genes))
    wanted = set(requested)
    selected: dict[str, np.ndarray] = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
        if not header:
            raise ValueError(f"Empty count matrix: {path}")
        barcodes = header.split("\t")

        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            gene, sep, values = line.partition("\t")
            if not sep or gene not in wanted:
                continue
            counts = np.fromstring(values, sep="\t", dtype=float)
            if counts.size != len(barcodes):
                raise ValueError(
                    f"Gene {gene!r} in {path} has {counts.size} counts; "
                    f"expected {len(barcodes)}."
                )
            if gene in selected:
                selected[gene] = selected[gene] + counts
            else:
                selected[gene] = counts

    return pd.DataFrame(
        {gene: selected[gene] for gene in requested if gene in selected},
        index=barcodes,
    )


def read_tabular_selected_genes_with_totals(path: Path, genes: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    """Read selected gene counts and full per-cell totals from a gzipped matrix.

    The expected format is:

    - first row: tab-separated cell barcodes
    - subsequent rows: gene symbol followed by one raw count per cell
    """
    requested = list(dict.fromkeys(genes))
    wanted = set(requested)
    selected: dict[str, np.ndarray] = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
        if not header:
            raise ValueError(f"Empty count matrix: {path}")
        barcodes = header.split("\t")
        totals = np.zeros(len(barcodes), dtype=float)

        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            gene, sep, values = line.partition("\t")
            if not sep:
                continue
            counts = np.fromstring(values, sep="\t", dtype=float)
            if counts.size != len(barcodes):
                raise ValueError(
                    f"Gene {gene!r} in {path} has {counts.size} counts; "
                    f"expected {len(barcodes)}."
                )
            totals += counts
            if gene in wanted:
                if gene in selected:
                    selected[gene] = selected[gene] + counts
                else:
                    selected[gene] = counts

    matrix = pd.DataFrame(
        {gene: selected[gene] for gene in requested if gene in selected},
        index=barcodes,
    )
    return matrix, pd.Series(totals, index=barcodes, name="total_counts")
