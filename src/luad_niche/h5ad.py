"""Lightweight h5ad readers for selected genes and spatial coordinates."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def decode_value(value):
    """Decode h5py bytes into Python strings."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def decode_array(values) -> list:
    """Decode an HDF5 dataset or array into a Python list."""
    return [decode_value(value) for value in values]


def read_csr_group(group: h5py.Group) -> csr_matrix:
    """Read an AnnData CSR matrix group."""
    shape = tuple(int(x) for x in group.attrs["shape"])
    return csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)


def read_h5ad_var_names(path: Path) -> list[str]:
    """Read var index gene names from an h5ad file."""
    with h5py.File(path, "r") as handle:
        index_key = decode_value(handle["var"].attrs.get("_index", "_index"))
        return [str(value) for value in decode_array(handle["var"][index_key][:])]


def read_h5ad_obs(path: Path) -> pd.DataFrame:
    """Read obs metadata and spatial coordinates from an h5ad file."""
    with h5py.File(path, "r") as handle:
        obs_group = handle["obs"]
        index_key = decode_value(obs_group.attrs.get("_index", "spots"))
        obs = pd.DataFrame({"spot": decode_array(obs_group[index_key][:])})

        for key, obj in obs_group.items():
            if key == index_key:
                continue
            if isinstance(obj, h5py.Group) and {"categories", "codes"}.issubset(obj.keys()):
                categories = decode_array(obj["categories"][:])
                codes = obj["codes"][:]
                obs[key] = [categories[code] if code >= 0 else None for code in codes]
            elif isinstance(obj, h5py.Dataset) and len(obj.shape) == 1 and obj.shape[0] == len(obs):
                obs[key] = decode_array(obj[:])

        if "obsm" in handle and "spatial" in handle["obsm"]:
            spatial = np.asarray(handle["obsm/spatial"][:])
            obs["x"] = spatial[:, 0]
            obs["y"] = spatial[:, 1]
        return obs


def read_h5ad_selected_genes(
    path: Path,
    genes: list[str],
    matrix_key: str = "X",
) -> pd.DataFrame:
    """Read selected genes into a spots-by-genes dataframe."""
    requested = list(dict.fromkeys(genes))
    with h5py.File(path, "r") as handle:
        var_group = handle["var"]
        index_key = decode_value(var_group.attrs.get("_index", "_index"))
        var_names = [str(value) for value in decode_array(var_group[index_key][:])]
        obs_group = handle["obs"]
        obs_index_key = decode_value(obs_group.attrs.get("_index", "spots"))
        spots = [str(value) for value in decode_array(obs_group[obs_index_key][:])]
        gene_to_index = {gene: idx for idx, gene in enumerate(var_names)}
        present_genes = [gene for gene in requested if gene in gene_to_index]
        if not present_genes:
            return pd.DataFrame(index=spots)
        matrix = read_csr_group(handle[matrix_key])
        selected = matrix[:, [gene_to_index[gene] for gene in present_genes]].toarray()
    return pd.DataFrame(selected, index=spots, columns=present_genes)


def compute_panel_scores(
    expr: pd.DataFrame,
    panels: Mapping[str, list[str]],
) -> pd.DataFrame:
    """Compute per-spot mean expression scores for gene panels."""
    scores = pd.DataFrame(index=expr.index)
    genes_used: dict[str, list[str]] = {}
    for panel_name, genes in panels.items():
        present = [gene for gene in genes if gene in expr.columns]
        genes_used[panel_name] = present
        column = f"{panel_name}_score"
        scores[column] = expr[present].mean(axis=1) if present else np.nan
    scores.attrs["panel_genes_used"] = genes_used
    return scores

