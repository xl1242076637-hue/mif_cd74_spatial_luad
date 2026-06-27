#!/usr/bin/env python
"""Plot focused GSE308103 snRNA and GSE282617 bulk orthogonal-validation panels."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.orthogonal_validation import (  # noqa: E402
    BULK_STAGE_ORDER,
    FOCUSED_GENE_CONTEXTS,
    SNRNA_STAGE_ORDER,
    prepare_bulk_focused_stage_source,
    prepare_snrna_focused_stage_source,
)


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


GENE_COLORS = {
    "MIF": "#B64342",
    "CD74": "#C95A49",
    "CD44": "#E28E2C",
    "CXCR4": "#7C6CCF",
    "SPP1": "#42949E",
    "TREM2": "#7BAA5B",
    "PLA2G7": "#5B8FD6",
    "IL1B": "#8F8F8F",
    "TNF": "#A8A8A8",
    "CXCL8": "#C9CCD1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snrna-stage-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_snrna_candidate_gene_stage_summary.csv",
    )
    parser.add_argument(
        "--bulk-stage-means",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse282617_candidate_marker_stage_means.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    return parser.parse_args()


def save_figure(fig: plt.Figure, base_path: Path) -> list[Path]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    saved = []
    for suffix, kwargs in {
        ".svg": {},
        ".pdf": {},
        ".tiff": {"dpi": 600},
        ".png": {"dpi": 300},
    }.items():
        output = base_path.with_suffix(suffix)
        fig.savefig(output, bbox_inches="tight", **kwargs)
        saved.append(output)
    plt.close(fig)
    return saved


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, fontsize=8, fontweight="bold")


def _heatmap_matrix(source: pd.DataFrame, stages: tuple[str, ...]) -> pd.DataFrame:
    matrix = source.pivot_table(index="gene", columns="stage", values="row_zscore", aggfunc="mean")
    matrix = matrix.reindex(index=list(FOCUSED_GENE_CONTEXTS), columns=list(stages))
    return matrix


def plot_heatmap(
    ax: plt.Axes,
    source: pd.DataFrame,
    stages: tuple[str, ...],
    title: str,
    *,
    context_labels: bool,
) -> None:
    matrix = _heatmap_matrix(source, stages)
    values = matrix.to_numpy(dtype=float)
    limit = max(float(np.nanmax(np.abs(values))), 1.0)
    image = ax.imshow(
        values,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
    )
    labels = []
    for gene in matrix.index:
        context = FOCUSED_GENE_CONTEXTS[gene]
        short_context = "epi" if context == "epithelial" else "macro"
        labels.append(f"{gene} ({short_context})" if context_labels else gene)
    ax.set_xticks(range(len(matrix.columns)), matrix.columns)
    ax.set_yticks(range(len(matrix.index)), labels)
    ax.set_title(title, loc="left", fontsize=7.4, fontweight="bold")
    ax.tick_params(length=0)
    colorbar = ax.figure.colorbar(image, ax=ax, shrink=0.72, pad=0.02)
    colorbar.set_label("row z-score", fontsize=6.4)
    colorbar.ax.tick_params(labelsize=5.8)


def plot_bulk_delta(ax: plt.Axes, bulk: pd.DataFrame) -> None:
    delta = bulk[["gene", "delta_late_vs_normal"]].drop_duplicates("gene").copy()
    delta = delta.set_index("gene").reindex(list(FOCUSED_GENE_CONTEXTS)).reset_index()
    y = np.arange(len(delta))
    values = delta["delta_late_vs_normal"].to_numpy(dtype=float)
    colors = [GENE_COLORS[gene] for gene in delta["gene"]]
    ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="#4D4D4D", linewidth=0.7)
    ax.set_yticks(y, delta["gene"])
    ax.invert_yaxis()
    ax.set_xlabel("IAC minus Normal mean expression")
    ax.set_title("c  Bulk progression deltas", loc="left", fontsize=7.4, fontweight="bold")


def plot_snrna_trends(ax: plt.Axes, snrna: pd.DataFrame) -> None:
    focus_genes = ("MIF", "CD74", "SPP1", "IL1B")
    x = np.arange(len(SNRNA_STAGE_ORDER))
    for gene in focus_genes:
        subset = snrna[snrna["gene"].eq(gene)].set_index("stage").reindex(SNRNA_STAGE_ORDER)
        ax.plot(
            x,
            subset["row_zscore"],
            marker="o",
            linewidth=1.2,
            markersize=3.2,
            label=gene,
            color=GENE_COLORS[gene],
        )
    ax.axhline(0, color="#A8A8A8", linewidth=0.6)
    ax.set_xticks(x, SNRNA_STAGE_ORDER)
    ax.set_ylabel("Compartment-matched row z-score")
    ax.set_title("d  Selected snRNA trends", loc="left", fontsize=7.4, fontweight="bold")
    ax.legend(ncol=2, fontsize=6.2, loc="upper left")


def main() -> int:
    args = parse_args()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    snrna_stage = pd.read_csv(args.snrna_stage_summary)
    bulk_stage = pd.read_csv(args.bulk_stage_means)
    snrna = prepare_snrna_focused_stage_source(snrna_stage)
    bulk = prepare_bulk_focused_stage_source(bulk_stage)
    snrna_output = args.table_dir / "supplementary_figure_focused_orthogonal_validation_snrna_source.csv"
    bulk_output = args.table_dir / "supplementary_figure_focused_orthogonal_validation_bulk_source.csv"
    snrna.to_csv(snrna_output, index=False, encoding="utf-8-sig")
    bulk.to_csv(bulk_output, index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 6.7), constrained_layout=True)
    grid = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.18, 0.82], width_ratios=[1.08, 0.92])
    ax_snrna = fig.add_subplot(grid[0, 0])
    ax_bulk = fig.add_subplot(grid[0, 1])
    ax_bulk_delta = fig.add_subplot(grid[1, 0])
    ax_snrna_trend = fig.add_subplot(grid[1, 1])
    plot_heatmap(
        ax_snrna,
        snrna,
        SNRNA_STAGE_ORDER,
        "a  GSE308103 snRNA compartment-matched expression",
        context_labels=True,
    )
    plot_heatmap(
        ax_bulk,
        bulk,
        BULK_STAGE_ORDER,
        "b  GSE282617 bulk expression",
        context_labels=False,
    )
    plot_bulk_delta(ax_bulk_delta, bulk)
    plot_snrna_trends(ax_snrna_trend, snrna)
    fig.suptitle(
        "Focused non-spatial cohorts provide orthogonal trend support",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    saved = save_figure(fig, args.figure_dir / "supplementary_figure_focused_orthogonal_validation")
    print(f"Wrote source data: {snrna_output}")
    print(f"Wrote source data: {bulk_output}")
    print("Wrote figure exports:")
    for path in saved:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
