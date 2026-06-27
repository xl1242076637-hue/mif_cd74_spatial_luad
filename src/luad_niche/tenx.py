"""Lightweight readers for 10x Visium MatrixMarket spatial outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


@dataclass(frozen=True)
class TenXSampleFiles:
    """Paths belonging to one extracted 10x/Visium sample."""

    sample_name: str
    sample_accession: str
    barcodes: Path
    features: Path
    matrix: Path
    tissue_positions: Path | None = None


def discover_10x_samples(
    root: Path,
    require_tissue_positions: bool = True,
) -> dict[str, TenXSampleFiles]:
    """Discover complete 10x sample file groups under an extracted GEO directory."""
    root = Path(root)
    samples: dict[str, TenXSampleFiles] = {}

    flat_layouts = [
        {
            "barcode_glob": "*_barcodes.tsv.gz",
            "barcode_suffix": "_barcodes.tsv.gz",
            "feature_suffixes": ["_features.tsv.gz", "_genes.tsv.gz"],
            "matrix_suffix": "_matrix.mtx.gz",
            "position_suffixes": ["_tissue_positions_list.csv.gz"],
        },
        {
            "barcode_glob": "*.barcodes.tsv.gz",
            "barcode_suffix": ".barcodes.tsv.gz",
            "feature_suffixes": [".features.tsv.gz", ".genes.tsv.gz"],
            "matrix_suffix": ".matrix.mtx.gz",
            "position_suffixes": [".tissue_positions_list.csv.gz"],
        },
    ]
    for layout in flat_layouts:
        for barcode_path in sorted(root.glob(layout["barcode_glob"])):
            sample_name = barcode_path.name.removesuffix(layout["barcode_suffix"])
            feature_path = next(
                (root / f"{sample_name}{suffix}" for suffix in layout["feature_suffixes"] if (root / f"{sample_name}{suffix}").exists()),
                root / f"{sample_name}{layout['feature_suffixes'][0]}",
            )
            tissue_positions = next(
                (
                    root / f"{sample_name}{suffix}"
                    for suffix in layout["position_suffixes"]
                    if (root / f"{sample_name}{suffix}").exists()
                ),
                None,
            )
            files = {
                "barcodes": barcode_path,
                "features": feature_path,
                "matrix": root / f"{sample_name}{layout['matrix_suffix']}",
                "tissue_positions": tissue_positions,
            }
            required_paths = [files["barcodes"], files["features"], files["matrix"]]
            if require_tissue_positions:
                required_paths.append(files["tissue_positions"] or root / f"{sample_name}{layout['position_suffixes'][0]}")
            missing = [str(path) for path in required_paths if not Path(path).exists()]
            if missing:
                raise FileNotFoundError(f"Incomplete 10x sample {sample_name}: missing {missing}")
            samples[sample_name] = TenXSampleFiles(
                sample_name=sample_name,
                sample_accession=sample_name.split("_", maxsplit=1)[0],
                barcodes=files["barcodes"],
                features=files["features"],
                matrix=files["matrix"],
                tissue_positions=tissue_positions,
            )

    for barcode_path in sorted(root.glob("*_barcodes.tsv.gz")):
        sample_name = barcode_path.name.removesuffix("_barcodes.tsv.gz")
        if sample_name in samples:
            continue
        files = {
            "barcodes": barcode_path,
            "features": root / f"{sample_name}_features.tsv.gz",
            "matrix": root / f"{sample_name}_matrix.mtx.gz",
            "tissue_positions": root / f"{sample_name}_tissue_positions_list.csv.gz",
        }
        required_paths = [files["barcodes"], files["features"], files["matrix"]]
        if require_tissue_positions:
            required_paths.append(files["tissue_positions"])
        missing = [str(path) for path in required_paths if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Incomplete 10x sample {sample_name}: missing {missing}")
        tissue_positions = files["tissue_positions"] if files["tissue_positions"].exists() else None
        samples[sample_name] = TenXSampleFiles(
            sample_name=sample_name,
            sample_accession=sample_name.split("_", maxsplit=1)[0],
            barcodes=files["barcodes"],
            features=files["features"],
            matrix=files["matrix"],
            tissue_positions=tissue_positions,
        )

    for barcode_path in sorted(root.glob("*/filtered_feature_bc_matrix/barcodes.tsv.gz")):
        sample_dir = barcode_path.parents[1]
        sample_name = sample_dir.name
        matrix_dir = sample_dir / "filtered_feature_bc_matrix"
        spatial_dir = sample_dir / "spatial"
        candidate_positions = [
            spatial_dir / "tissue_positions.csv",
            spatial_dir / "tissue_positions_list.csv",
            spatial_dir / "tissue_positions.csv.gz",
            spatial_dir / "tissue_positions_list.csv.gz",
        ]
        tissue_positions = next((path for path in candidate_positions if path.exists()), None)
        files = {
            "barcodes": barcode_path,
            "features": matrix_dir / "features.tsv.gz",
            "matrix": matrix_dir / "matrix.mtx.gz",
        }
        required_paths = [files["barcodes"], files["features"], files["matrix"]]
        if require_tissue_positions:
            required_paths.append(tissue_positions or spatial_dir / "tissue_positions.csv")
        missing = [str(path) for path in required_paths if not Path(path).exists()]
        if missing:
            raise FileNotFoundError(f"Incomplete 10x sample {sample_name}: missing {missing}")
        samples[sample_name] = TenXSampleFiles(
            sample_name=sample_name,
            sample_accession=sample_name.split("_", maxsplit=1)[0],
            barcodes=files["barcodes"],
            features=files["features"],
            matrix=files["matrix"],
            tissue_positions=tissue_positions,
        )
    return samples


def read_10x_barcodes(sample: TenXSampleFiles) -> list[str]:
    """Read 10x spot barcodes in matrix column order."""
    return pd.read_csv(sample.barcodes, header=None, compression="gzip")[0].astype(str).tolist()


def read_10x_features(sample: TenXSampleFiles) -> pd.DataFrame:
    """Read 10x feature metadata in matrix row order."""
    return pd.read_csv(
        sample.features,
        sep="\t",
        header=None,
        names=["feature_id", "gene", "feature_type"],
        compression="gzip",
    )


def read_10x_tissue_positions(sample: TenXSampleFiles) -> pd.DataFrame:
    """Read Visium tissue positions and expose pixel coordinates as x/y."""
    if sample.tissue_positions is None:
        raise ValueError(f"Sample {sample.sample_name} has no tissue positions file.")
    positions = pd.read_csv(
        sample.tissue_positions,
        header=None,
        compression="infer",
        names=[
            "spot",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ],
    )
    if not positions.empty and str(positions.iloc[0]["spot"]).lower() == "barcode":
        positions = positions.iloc[1:].reset_index(drop=True)
    for column in ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]:
        positions[column] = pd.to_numeric(positions[column], errors="coerce")
    positions["x"] = positions["pxl_col_in_fullres"]
    positions["y"] = positions["pxl_row_in_fullres"]
    return positions


def read_10x_obs(sample: TenXSampleFiles) -> pd.DataFrame:
    """Return matrix barcodes merged with Visium spatial coordinates."""
    obs = pd.DataFrame({"spot": read_10x_barcodes(sample)})
    if sample.tissue_positions is None:
        return obs
    positions = read_10x_tissue_positions(sample)
    return obs.merge(
        positions[["spot", "in_tissue", "array_row", "array_col", "x", "y"]],
        on="spot",
        how="left",
    )


def read_10x_selected_genes(sample: TenXSampleFiles, genes: list[str]) -> pd.DataFrame:
    """Read selected gene symbols into a spots-by-genes dataframe."""
    expr, _ = read_10x_selected_genes_with_totals(sample, genes)
    return expr


def read_10x_selected_genes_with_totals(
    sample: TenXSampleFiles,
    genes: list[str],
) -> tuple[pd.DataFrame, pd.Series]:
    """Read selected gene symbols into a spots-by-genes dataframe.

    The 10x MatrixMarket matrix is features-by-barcodes. If the same gene symbol
    appears more than once, rows are summed so panel scoring uses the full signal.
    The returned totals are all-feature library sizes per spot, useful for
    normalizing a selected-gene subset without using subset totals as the
    denominator.
    """
    requested = list(dict.fromkeys(genes))
    barcodes = read_10x_barcodes(sample)
    features = read_10x_features(sample)
    gene_to_rows: dict[str, list[int]] = {}
    for row_index, gene in enumerate(features["gene"].astype(str)):
        if gene in requested:
            gene_to_rows.setdefault(gene, []).append(row_index)

    matrix = mmread(sample.matrix).tocsr()
    total_counts = pd.Series(np.asarray(matrix.sum(axis=0)).ravel(), index=barcodes, name="total_counts")
    present_genes = [gene for gene in requested if gene in gene_to_rows]
    if not present_genes:
        return pd.DataFrame(index=barcodes), total_counts

    columns = {}
    for gene in present_genes:
        summed = matrix[gene_to_rows[gene], :].sum(axis=0)
        columns[gene] = np.asarray(summed).ravel()
    return pd.DataFrame(columns, index=barcodes), total_counts
