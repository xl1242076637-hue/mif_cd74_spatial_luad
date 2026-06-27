from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from luad_niche.h5ad import (
    compute_panel_scores,
    read_h5ad_obs,
    read_h5ad_selected_genes,
)


def write_tiny_h5ad(path: Path) -> None:
    matrix = csr_matrix(np.array([[1.0, 0.0, 3.0], [0.0, 2.0, 4.0]], dtype=np.float32))
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        x = handle.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["shape"] = matrix.shape
        x.create_dataset("data", data=matrix.data)
        x.create_dataset("indices", data=matrix.indices)
        x.create_dataset("indptr", data=matrix.indptr)

        obs = handle.create_group("obs")
        obs.attrs["_index"] = "spots"
        obs.create_dataset("spots", data=np.array(["spot1", "spot2"], dtype=object), dtype=string_dtype)
        cat = obs.create_group("cell_type")
        cat.create_dataset("categories", data=np.array(["Epithelial", "Macrophage"], dtype=object), dtype=string_dtype)
        cat.create_dataset("codes", data=np.array([0, 1], dtype=np.int8))

        var = handle.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.array(["EPCAM", "IL1B", "SPP1"], dtype=object), dtype=string_dtype)

        obsm = handle.create_group("obsm")
        obsm.create_dataset("spatial", data=np.array([[10.0, 20.0], [30.0, 40.0]]))


def test_read_h5ad_obs_decodes_categorical_and_spatial_coordinates(tmp_path):
    path = tmp_path / "tiny.h5ad"
    write_tiny_h5ad(path)

    obs = read_h5ad_obs(path)

    assert obs[["spot", "cell_type", "x", "y"]].to_dict("records") == [
        {"spot": "spot1", "cell_type": "Epithelial", "x": 10.0, "y": 20.0},
        {"spot": "spot2", "cell_type": "Macrophage", "x": 30.0, "y": 40.0},
    ]


def test_read_h5ad_selected_genes_returns_spots_by_genes(tmp_path):
    path = tmp_path / "tiny.h5ad"
    write_tiny_h5ad(path)

    expr = read_h5ad_selected_genes(path, ["SPP1", "MISSING", "EPCAM"])

    assert list(expr.index) == ["spot1", "spot2"]
    assert list(expr.columns) == ["SPP1", "EPCAM"]
    assert expr.loc["spot1", "SPP1"] == 3.0
    assert expr.loc["spot2", "EPCAM"] == 0.0


def test_compute_panel_scores_averages_present_genes():
    expr = pd.DataFrame(
        {"EPCAM": [2.0, 4.0], "IL1B": [0.0, 6.0]},
        index=["spot1", "spot2"],
    )

    scores = compute_panel_scores(
        expr,
        {"epithelial": ["EPCAM", "KRT8"], "macrophage": ["IL1B"]},
    )

    assert scores.loc["spot1", "epithelial_score"] == 2.0
    assert scores.loc["spot2", "macrophage_score"] == 6.0
    assert scores.attrs["panel_genes_used"] == {
        "epithelial": ["EPCAM"],
        "macrophage": ["IL1B"],
    }

