import gzip

import numpy as np
import pandas as pd

from luad_niche.wide_matrix import read_wide_matrix_header, read_wide_selected_genes


def write_tiny_wide_matrix(path):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("Index\tcell1\tcell2\tcell3\n")
        handle.write("EPCAM\t1\t0\t2\n")
        handle.write("SPP1\t0\t3\t1\n")
        handle.write("EPCAM\t4\t1\t0\n")
        handle.write("LYZ\t2\t2\t2\n")


def test_read_wide_matrix_header_returns_cell_ids(tmp_path):
    path = tmp_path / "matrix.txt.gz"
    write_tiny_wide_matrix(path)

    assert read_wide_matrix_header(path) == ["cell1", "cell2", "cell3"]


def test_read_wide_selected_genes_preserves_requested_order_and_sums_duplicates(tmp_path):
    path = tmp_path / "matrix.txt.gz"
    write_tiny_wide_matrix(path)

    expr = read_wide_selected_genes(path, ["SPP1", "EPCAM", "MISSING"])

    expected = pd.DataFrame(
        {
            "SPP1": np.array([0.0, 3.0, 1.0], dtype=np.float32),
            "EPCAM": np.array([5.0, 1.0, 2.0], dtype=np.float32),
        },
        index=["cell1", "cell2", "cell3"],
    )
    pd.testing.assert_frame_equal(expr, expected)
    assert expr.attrs["present_genes"] == ["SPP1", "EPCAM"]
    assert expr.attrs["missing_genes"] == ["MISSING"]
