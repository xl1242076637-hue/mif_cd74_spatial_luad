#!/usr/bin/env python
"""Score selected LUAD niche signatures in the GSE131907 lung cancer atlas."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.wide_matrix import read_wide_selected_genes  # noqa: E402


REFINED_PANEL_NAMES = {
    "epithelial_progenitor_like_vs_other_epithelial": "refined_epithelial_progenitor_like",
    "proliferating_epithelial_vs_other_epithelial": "refined_proliferating_epithelial",
    "spp1_macrophage_vs_other_macrophage": "refined_spp1_macrophage",
    "c1q_macrophage_vs_other_macrophage": "refined_c1q_macrophage",
    "inflammatory_macrophage_vs_other_macrophage": "refined_inflammatory_macrophage",
}
SELECTED_FIGURE_SCORES = [
    "epithelial_progenitor_like",
    "refined_epithelial_progenitor_like",
    "spp1_macrophage",
    "refined_spp1_macrophage",
    "inflammatory_macrophage",
    "refined_inflammatory_macrophage",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        type=Path,
        default=PROJECT_ROOT
        / "data"
        / "raw"
        / "GSE131907"
        / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=PROJECT_ROOT
        / "data"
        / "raw"
        / "GSE131907"
        / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    )
    parser.add_argument(
        "--markers",
        type=Path,
        default=PROJECT_ROOT / "config" / "cell_state_markers.yaml",
    )
    parser.add_argument(
        "--refined-signatures",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_state_signature_genes.json",
    )
    parser.add_argument(
        "--candidate-mechanisms",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_mechanisms.yaml",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    parser.add_argument("--min-cells", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--limit-cells", type=int, default=None)
    return parser.parse_args()


def load_candidate_mechanism_genes(mechanism_path: Path | None) -> list[str]:
    if mechanism_path is None or not Path(mechanism_path).exists():
        return []
    with Path(mechanism_path).open("r", encoding="utf-8") as handle:
        mechanism_config = yaml.safe_load(handle) or {}
    genes: list[str] = []
    for axis in mechanism_config.get("axes", []):
        for field in ("source_genes", "target_genes", "bulk_genes", "perturbation_genes"):
            for gene in axis.get(field, []) or []:
                if gene not in genes:
                    genes.append(gene)
    return genes


def load_panels(marker_path: Path, refined_path: Path, mechanism_path: Path | None = None) -> dict[str, list[str]]:
    with Path(marker_path).open("r", encoding="utf-8") as handle:
        marker_config = yaml.safe_load(handle)
    panels = {
        **{name: list(genes) for name, genes in marker_config["broad_classes"].items()},
        **{name: list(genes) for name, genes in marker_config["state_panels"].items()},
    }
    refined = json.loads(Path(refined_path).read_text(encoding="utf-8"))
    panels.update(
        {
            output_name: list(refined[contrast])
            for contrast, output_name in REFINED_PANEL_NAMES.items()
            if contrast in refined
        }
    )
    mechanism_genes = load_candidate_mechanism_genes(mechanism_path)
    if mechanism_genes:
        panels["candidate_mechanism_genes"] = mechanism_genes
    return panels


def flatten_genes(panels: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panel_genes in panels.values():
        genes.extend(gene for gene in panel_genes if gene not in genes)
    return genes


def read_annotation(path: Path) -> pd.DataFrame:
    annotation = pd.read_csv(path, sep="\t", compression="gzip")
    annotation = annotation.rename(columns={annotation.columns[0]: "cell_id"})
    annotation["cell_id"] = annotation["cell_id"].astype(str)
    return annotation.set_index("cell_id", drop=False)


def summarize_scores(
    table: pd.DataFrame,
    group_columns: list[str],
    score_columns: list[str],
    min_cells: int = 0,
) -> pd.DataFrame:
    records = []
    for score_column in score_columns:
        grouped = (
            table.groupby(group_columns, dropna=False)[score_column]
            .agg(n_cells="size", mean_score="mean", median_score="median", sd_score="std")
            .reset_index()
        )
        grouped["score"] = score_column.removesuffix("_score")
        records.append(grouped)
    summary = pd.concat(records, ignore_index=True)
    if min_cells:
        summary = summary[summary["n_cells"] >= min_cells].copy()
    return summary


def summarize_genes(
    expr: pd.DataFrame,
    annotation: pd.DataFrame,
    group_columns: list[str],
    min_cells: int = 0,
) -> pd.DataFrame:
    records = []
    working = pd.concat([annotation[group_columns].reset_index(drop=True), expr.reset_index(drop=True)], axis=1)
    for gene in expr.columns:
        grouped = (
            working.groupby(group_columns, dropna=False)[gene]
            .agg(n_cells="size", mean_expression="mean", median_expression="median")
            .reset_index()
        )
        grouped["gene"] = gene
        records.append(grouped)
    summary = pd.concat(records, ignore_index=True)
    if min_cells:
        summary = summary[summary["n_cells"] >= min_cells].copy()
    return summary


def top_celltype_by_gene(gene_summary: pd.DataFrame) -> pd.DataFrame:
    ordered = gene_summary.sort_values(["gene", "mean_expression"], ascending=[True, False])
    return ordered.groupby("gene", as_index=False).head(1).reset_index(drop=True)


def top_groups_by_score(summary: pd.DataFrame, score_column: str, group_columns: list[str], top_n: int) -> pd.DataFrame:
    selected = summary[summary["score"] == score_column].copy()
    selected = selected.sort_values("mean_score", ascending=False)
    return selected[group_columns + ["score", "n_cells", "mean_score", "median_score", "sd_score"]].head(top_n)


def plot_refined_score_heatmap(celltype_summary: pd.DataFrame, output: Path) -> None:
    subset = celltype_summary[celltype_summary["score"].isin(SELECTED_FIGURE_SCORES)].copy()
    pivot = subset.pivot_table(index="Cell_type.refined", columns="score", values="mean_score")
    present_columns = [score for score in SELECTED_FIGURE_SCORES if score in pivot.columns]
    pivot = pivot[present_columns].fillna(0.0)
    scaled = pivot.copy()
    for column in scaled.columns:
        values = scaled[column]
        value_range = values.max() - values.min()
        scaled[column] = 0.0 if value_range == 0 else (values - values.min()) / value_range

    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.35 * len(scaled))))
    image = ax.imshow(scaled.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(scaled.columns)))
    ax.set_xticklabels(scaled.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(scaled.index)))
    ax.set_yticklabels(scaled.index, fontsize=8)
    ax.set_title("GSE131907 selected signature specificity by refined cell type")
    fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="Column-scaled mean score")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    panels = load_panels(args.markers, args.refined_signatures, args.candidate_mechanisms)
    all_genes = flatten_genes(panels)
    annotation = read_annotation(args.annotation)

    expr = read_wide_selected_genes(args.matrix, all_genes)
    matrix_present_genes = expr.attrs.get("present_genes", list(expr.columns))
    matrix_missing_genes = expr.attrs.get("missing_genes", [])
    if args.limit_cells:
        expr = expr.iloc[: args.limit_cells].copy()
    expr = np.log1p(expr)

    common_cells = expr.index.intersection(annotation.index)
    if len(common_cells) == 0:
        raise SystemExit("No matrix cells matched annotation cell IDs.")
    expr = expr.loc[common_cells]
    annotation = annotation.loc[common_cells]
    scores = compute_panel_scores(expr, panels)

    output_columns = ["cell_id", "Sample", "Sample_Origin", "Cell_type", "Cell_type.refined", "Cell_subtype"]
    compact = pd.concat([annotation[output_columns].reset_index(drop=True), scores.reset_index(drop=True)], axis=1)
    score_columns = [column for column in compact.columns if column.endswith("_score")]

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    compact.to_csv(
        args.table_dir / "gse131907_selected_signature_cell_scores.csv",
        index=False,
        encoding="utf-8-sig",
    )
    origin_summary = summarize_scores(compact, ["Sample_Origin"], score_columns)
    celltype_summary = summarize_scores(compact, ["Cell_type.refined"], score_columns, min_cells=args.min_cells)
    subtype_summary = summarize_scores(
        compact,
        ["Sample_Origin", "Cell_type.refined", "Cell_subtype"],
        score_columns,
        min_cells=args.min_cells,
    )
    origin_celltype_summary = summarize_scores(
        compact,
        ["Sample_Origin", "Cell_type.refined"],
        score_columns,
        min_cells=args.min_cells,
    )
    gene_celltype_summary = summarize_genes(
        expr,
        annotation,
        ["Cell_type.refined"],
        min_cells=args.min_cells,
    )
    gene_origin_celltype_summary = summarize_genes(
        expr,
        annotation,
        ["Sample_Origin", "Cell_type.refined"],
        min_cells=args.min_cells,
    )

    origin_summary.to_csv(
        args.table_dir / "gse131907_selected_signature_origin_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    celltype_summary.to_csv(
        args.table_dir / "gse131907_selected_signature_celltype_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    origin_celltype_summary.to_csv(
        args.table_dir / "gse131907_selected_signature_origin_celltype_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    subtype_summary.to_csv(
        args.table_dir / "gse131907_selected_signature_subtype_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    gene_celltype_summary.to_csv(
        args.table_dir / "gse131907_selected_gene_celltype_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    gene_origin_celltype_summary.to_csv(
        args.table_dir / "gse131907_selected_gene_origin_celltype_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top_celltype_by_gene(gene_celltype_summary).to_csv(
        args.table_dir / "gse131907_selected_gene_top_celltype.csv",
        index=False,
        encoding="utf-8-sig",
    )

    top_records = []
    for score in SELECTED_FIGURE_SCORES:
        top_records.append(
            top_groups_by_score(
                subtype_summary,
                score,
                ["Sample_Origin", "Cell_type.refined", "Cell_subtype"],
                top_n=args.top_n,
            )
        )
    pd.concat(top_records, ignore_index=True).to_csv(
        args.table_dir / "gse131907_selected_signature_top_subtypes.csv",
        index=False,
        encoding="utf-8-sig",
    )

    genes_used = {
        "matrix_present_genes": matrix_present_genes,
        "matrix_missing_genes": matrix_missing_genes,
        "panel_genes_used": scores.attrs["panel_genes_used"],
        "transform": "log1p_raw_UMI",
        "n_cells_scored": int(len(compact)),
    }
    (args.table_dir / "gse131907_selected_signature_genes_used.json").write_text(
        json.dumps(genes_used, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    plot_refined_score_heatmap(
        celltype_summary,
        args.figure_dir / "gse131907_selected_signature_celltype_heatmap.png",
    )
    print(f"Scored {len(compact)} cells with {len(expr.columns)} selected genes.")
    print(f"Wrote GSE131907 selected-signature summaries to {args.table_dir}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
