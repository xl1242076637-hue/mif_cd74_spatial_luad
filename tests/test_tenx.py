import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

from luad_niche.tenx import (
    discover_10x_samples,
    read_10x_obs,
    read_10x_selected_genes,
    read_10x_selected_genes_with_totals,
)


def write_gzip_text(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def write_tiny_10x_sample(root: Path, sample: str = "GSMTEST_TD1") -> None:
    write_gzip_text(root / f"{sample}_barcodes.tsv.gz", "spot1\nspot2\n")
    write_gzip_text(
        root / f"{sample}_features.tsv.gz",
        "\n".join(
            [
                "ENSG000001\tEPCAM\tGene Expression",
                "ENSG000002\tSPP1\tGene Expression",
                "ENSG000003\tEPCAM\tGene Expression",
                "",
            ]
        ),
    )
    write_gzip_text(
        root / f"{sample}_tissue_positions_list.csv.gz",
        "\n".join(
            [
                "spot1,1,0,0,20,10",
                "spot2,1,0,1,40,30",
                "not_in_matrix,0,1,0,60,50",
                "",
            ]
        ),
    )
    matrix = coo_matrix(np.array([[1, 0], [2, 3], [4, 5]], dtype=np.int32))
    matrix_path = root / f"{sample}_matrix.mtx"
    mmwrite(matrix_path, matrix)
    with matrix_path.open("rb") as source, gzip.open(root / f"{sample}_matrix.mtx.gz", "wb") as target:
        target.write(source.read())
    matrix_path.unlink()


def write_tiny_10x_expression_sample(root: Path, sample: str = "GSMTEST_TD1") -> None:
    write_gzip_text(root / f"{sample}_barcodes.tsv.gz", "cell1\ncell2\n")
    write_gzip_text(
        root / f"{sample}_features.tsv.gz",
        "ENSG000001\tEPCAM\tGene Expression\nENSG000002\tSPP1\tGene Expression\n",
    )
    matrix = coo_matrix(np.array([[1, 0], [2, 3]], dtype=np.int32))
    matrix_path = root / f"{sample}_matrix.mtx"
    mmwrite(matrix_path, matrix)
    with matrix_path.open("rb") as source, gzip.open(root / f"{sample}_matrix.mtx.gz", "wb") as target:
        target.write(source.read())
    matrix_path.unlink()


def write_tiny_dot_10x_expression_sample(root: Path, sample: str = "GSMTEST_DOT") -> None:
    write_gzip_text(root / f"{sample}.barcodes.tsv.gz", "cell1\ncell2\n")
    write_gzip_text(
        root / f"{sample}.genes.tsv.gz",
        "ENSG000001\tEPCAM\nENSG000002\tSPP1\n",
    )
    matrix = coo_matrix(np.array([[1, 0], [2, 3]], dtype=np.int32))
    matrix_path = root / f"{sample}.matrix.mtx"
    mmwrite(matrix_path, matrix)
    with matrix_path.open("rb") as source, gzip.open(root / f"{sample}.matrix.mtx.gz", "wb") as target:
        target.write(source.read())
    matrix_path.unlink()


def write_tiny_spaceranger_sample(root: Path, sample: str = "P1_AAH") -> None:
    matrix_dir = root / sample / "filtered_feature_bc_matrix"
    spatial_dir = root / sample / "spatial"
    matrix_dir.mkdir(parents=True)
    spatial_dir.mkdir(parents=True)
    write_gzip_text(matrix_dir / "barcodes.tsv.gz", "spot1\nspot2\n")
    write_gzip_text(
        matrix_dir / "features.tsv.gz",
        "ENSG000001\tEPCAM\tGene Expression\nENSG000002\tSPP1\tGene Expression\n",
    )
    (spatial_dir / "tissue_positions.csv").write_text(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
        "spot1,1,0,0,20,10\n"
        "spot2,0,0,1,40,30\n",
        encoding="utf-8",
    )
    matrix = coo_matrix(np.array([[1, 0], [2, 3]], dtype=np.int32))
    matrix_path = matrix_dir / "matrix.mtx"
    mmwrite(matrix_path, matrix)
    with matrix_path.open("rb") as source, gzip.open(matrix_dir / "matrix.mtx.gz", "wb") as target:
        target.write(source.read())
    matrix_path.unlink()


def test_discover_10x_samples_groups_required_files(tmp_path):
    write_tiny_10x_sample(tmp_path)

    samples = discover_10x_samples(tmp_path)

    assert list(samples) == ["GSMTEST_TD1"]
    assert samples["GSMTEST_TD1"].sample_accession == "GSMTEST"
    assert samples["GSMTEST_TD1"].matrix.name == "GSMTEST_TD1_matrix.mtx.gz"


def test_discover_10x_samples_allows_expression_only_samples_when_positions_not_required(tmp_path):
    write_tiny_10x_expression_sample(tmp_path)

    samples = discover_10x_samples(tmp_path, require_tissue_positions=False)

    assert list(samples) == ["GSMTEST_TD1"]
    assert samples["GSMTEST_TD1"].tissue_positions is None


def test_discover_10x_samples_supports_dot_named_genes_files(tmp_path):
    write_tiny_dot_10x_expression_sample(tmp_path)

    samples = discover_10x_samples(tmp_path, require_tissue_positions=False)

    assert list(samples) == ["GSMTEST_DOT"]
    assert samples["GSMTEST_DOT"].features.name == "GSMTEST_DOT.genes.tsv.gz"


def test_discover_10x_samples_supports_spaceranger_directory_layout(tmp_path):
    write_tiny_spaceranger_sample(tmp_path)

    samples = discover_10x_samples(tmp_path)

    assert list(samples) == ["P1_AAH"]
    assert samples["P1_AAH"].sample_accession == "P1"
    assert samples["P1_AAH"].matrix.name == "matrix.mtx.gz"
    assert samples["P1_AAH"].tissue_positions.name == "tissue_positions.csv"


def test_read_10x_obs_returns_barcodes_without_spatial_columns_for_expression_only_sample(tmp_path):
    sample_name = "GSMTEST_TD1"
    write_tiny_10x_expression_sample(tmp_path, sample_name)
    sample = discover_10x_samples(tmp_path, require_tissue_positions=False)[sample_name]

    obs = read_10x_obs(sample)

    assert obs.to_dict("records") == [{"spot": "cell1"}, {"spot": "cell2"}]


def test_read_10x_obs_merges_barcodes_with_spatial_positions(tmp_path):
    sample_name = "GSMTEST_TD1"
    write_tiny_10x_sample(tmp_path, sample_name)
    sample = discover_10x_samples(tmp_path)[sample_name]

    obs = read_10x_obs(sample)

    assert obs[["spot", "in_tissue", "array_row", "array_col", "x", "y"]].to_dict("records") == [
        {"spot": "spot1", "in_tissue": 1, "array_row": 0, "array_col": 0, "x": 10, "y": 20},
        {"spot": "spot2", "in_tissue": 1, "array_row": 0, "array_col": 1, "x": 30, "y": 40},
    ]


def test_read_10x_obs_handles_spaceranger_headered_tissue_positions(tmp_path):
    write_tiny_spaceranger_sample(tmp_path)
    sample = discover_10x_samples(tmp_path)["P1_AAH"]

    obs = read_10x_obs(sample)

    assert obs[["spot", "in_tissue", "array_row", "array_col", "x", "y"]].to_dict("records") == [
        {"spot": "spot1", "in_tissue": 1, "array_row": 0, "array_col": 0, "x": 10, "y": 20},
        {"spot": "spot2", "in_tissue": 0, "array_row": 0, "array_col": 1, "x": 30, "y": 40},
    ]


def test_read_10x_selected_genes_returns_spots_by_gene_and_sums_duplicate_symbols(tmp_path):
    sample_name = "GSMTEST_TD1"
    write_tiny_10x_sample(tmp_path, sample_name)
    sample = discover_10x_samples(tmp_path)[sample_name]

    expr = read_10x_selected_genes(sample, ["SPP1", "MISSING", "EPCAM"])

    expected = pd.DataFrame(
        {"SPP1": [2, 3], "EPCAM": [5, 5]},
        index=["spot1", "spot2"],
    )
    pd.testing.assert_frame_equal(expr, expected)


def test_read_10x_selected_genes_with_totals_returns_total_counts_per_spot(tmp_path):
    sample_name = "GSMTEST_TD1"
    write_tiny_10x_sample(tmp_path, sample_name)
    sample = discover_10x_samples(tmp_path)[sample_name]

    expr, totals = read_10x_selected_genes_with_totals(sample, ["EPCAM"])

    assert expr["EPCAM"].tolist() == [5, 5]
    assert totals.to_dict() == {"spot1": 7, "spot2": 8}


def test_read_10x_selected_genes_supports_two_column_genes_file(tmp_path):
    sample_name = "GSMTEST_DOT"
    write_tiny_dot_10x_expression_sample(tmp_path, sample_name)
    sample = discover_10x_samples(tmp_path, require_tissue_positions=False)[sample_name]

    expr, totals = read_10x_selected_genes_with_totals(sample, ["EPCAM", "SPP1"])

    assert expr["EPCAM"].tolist() == [1, 0]
    assert expr["SPP1"].tolist() == [2, 3]
    assert totals.to_dict() == {"cell1": 3, "cell2": 3}
