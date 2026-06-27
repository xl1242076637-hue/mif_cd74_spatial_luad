"""Streaming readers for wide gene-by-cell text matrices."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def read_wide_matrix_header(path: Path) -> list[str]:
    """Return cell IDs from a tab-delimited gene-by-cell matrix header."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
    if not header:
        raise ValueError(f"Empty matrix file: {path}")
    columns = header.split("\t")
    if len(columns) < 2:
        raise ValueError(f"Matrix header must include a gene column and at least one cell: {path}")
    return columns[1:]


def read_wide_selected_genes(path: Path, genes: list[str]) -> pd.DataFrame:
    """Read selected genes from a wide gene-by-cell matrix.

    The source matrix is expected to have one header row with a first gene-name
    column followed by cell IDs. Returned data are cells-by-genes. Duplicate gene
    rows are summed, matching the behavior used for 10x feature duplicates.
    """
    requested = list(dict.fromkeys(genes))
    requested_set = set(requested)
    values_by_gene: dict[str, np.ndarray] = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
        if not header:
            raise ValueError(f"Empty matrix file: {path}")
        cell_ids = header.split("\t")[1:]
        n_cells = len(cell_ids)

        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            first_tab = line.find("\t")
            if first_tab < 0:
                continue
            gene = line[:first_tab]
            if gene not in requested_set:
                continue
            values = np.fromstring(line[first_tab + 1 :], sep="\t", dtype=np.float32)
            if len(values) != n_cells:
                raise ValueError(f"Gene {gene} has {len(values)} values but header has {n_cells} cells.")
            if gene in values_by_gene:
                values_by_gene[gene] = values_by_gene[gene] + values
            else:
                values_by_gene[gene] = values

    present = [gene for gene in requested if gene in values_by_gene]
    missing = [gene for gene in requested if gene not in values_by_gene]
    expr = pd.DataFrame({gene: values_by_gene[gene] for gene in present}, index=cell_ids)
    expr.attrs["present_genes"] = present
    expr.attrs["missing_genes"] = missing
    return expr
