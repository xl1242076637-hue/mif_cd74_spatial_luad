#!/usr/bin/env python
"""Plot the GSE308103 GRN-level virtual perturbation supplementary figure."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.figure_data import (  # noqa: E402
    prepare_grn_target_ranking_source,
    prepare_grn_top_gene_source,
)
from luad_niche.nature_figure_style import (  # noqa: E402
    NATURE_PALETTE,
    add_panel_label,
    apply_nature_style,
    save_publication_figure,
)


SIGNATURE_COLORS = {
    "C1Q macrophage": NATURE_PALETTE["green_support"],
    "SPP1 macrophage": NATURE_PALETTE["macrophage_axis"],
    "Proliferating epi": NATURE_PALETTE["mif_axis"],
    "Epi progenitor": "#D66A58",
    "Inflammatory macrophage": NATURE_PALETTE["gold_support"],
    "Resident macrophage": NATURE_PALETTE["immune_blue"],
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target-ranking",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_grn_virtual_perturbation_target_ranking.csv",
    )
    parser.add_argument(
        "--network-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_grn_virtual_perturbation_network_summary.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    return parser.parse_args()


def _panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.04) -> None:
    add_panel_label(ax, label, x=x, y=y)


def _signature_color(label: str) -> str:
    return SIGNATURE_COLORS.get(str(label), "#8F8F8F")


def plot_target_ranking(ax: plt.Axes, source: pd.DataFrame) -> None:
    table = source.sort_values("display_rank", ascending=True).copy()
    y = np.arange(len(table))
    colors = [_signature_color(label) for label in table["signature_label"]]
    values = table["top_signature_mean_impact"].to_numpy(dtype=float)
    ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.55)
    ax.set_yticks(y, table["target_gene"])
    ax.invert_yaxis()
    ax.set_xlabel("mean propagated impact score")
    ax.set_title("GRN target ranking", loc="left", fontsize=7.6, fontweight="bold")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    xmax = max(values.max() * 1.22, 0.001)
    ax.set_xlim(0, xmax)
    for idx, row in table.iterrows():
        value = float(row["top_signature_mean_impact"])
        ax.text(
            value + xmax * 0.015,
            int(row["display_rank"]) - 1,
            row["signature_label"],
            va="center",
            ha="left",
            fontsize=6.0,
            color=NATURE_PALETTE["neutral_dark"],
        )
    _panel_label(ax, "a")


def plot_top_gene_tiles(ax: plt.Axes, target_source: pd.DataFrame, gene_source: pd.DataFrame) -> None:
    targets = target_source.sort_values("display_rank")["target_gene"].tolist()
    target_y = {target: i for i, target in enumerate(targets)}
    max_position = int(gene_source["gene_position"].max()) if not gene_source.empty else 0
    ax.set_xlim(0.5, max_position + 0.5)
    ax.set_ylim(len(targets) - 0.5, -0.5)
    ax.set_xticks(range(1, max_position + 1), [f"Top {i}" for i in range(1, max_position + 1)])
    ax.set_yticks(range(len(targets)), targets)
    ax.set_title("Top propagated gene neighborhoods", loc="left", fontsize=7.6, fontweight="bold")
    for row in gene_source.itertuples(index=False):
        y = target_y[row.target_gene]
        x = int(row.gene_position)
        color = _signature_color(row.signature_label)
        rect = Rectangle((x - 0.48, y - 0.38), 0.96, 0.76, facecolor=color, alpha=0.16, edgecolor=color, linewidth=0.55)
        ax.add_patch(rect)
        ax.text(x, y, row.impacted_gene, ha="center", va="center", fontsize=5.8, color=NATURE_PALETTE["neutral_black"])
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    _panel_label(ax, "b")


def plot_network_summary(ax: plt.Axes, network_summary: pd.DataFrame) -> None:
    ax.axis("off")
    lines = ["GSE308103 MIA/LUAD GRNs"]
    for row in network_summary.sort_values("broad_class").itertuples(index=False):
        label = str(row.broad_class).capitalize()
        lines.append(f"{label}: {row.n_cells:,} cells, {row.n_network_genes} genes, {row.n_edges:,} edges")
    ax.text(
        0.0,
        0.98,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=6.4,
        color=NATURE_PALETTE["neutral_dark"],
        linespacing=1.35,
    )
    _panel_label(ax, "c", x=-0.02, y=1.08)


def main() -> int:
    apply_nature_style(font_size=7.0)
    args = parse_args()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    target_ranking = pd.read_csv(args.target_ranking)
    network_summary = pd.read_csv(args.network_summary)
    target_source = prepare_grn_target_ranking_source(target_ranking)
    gene_source = prepare_grn_top_gene_source(target_source, top_n=5)

    target_output = args.table_dir / "supplementary_figure_grn_virtual_perturbation_target_source.csv"
    gene_output = args.table_dir / "supplementary_figure_grn_virtual_perturbation_top_gene_source.csv"
    target_source.to_csv(target_output, index=False, encoding="utf-8-sig")
    gene_source.to_csv(gene_output, index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 3.85), constrained_layout=True)
    grid = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.18], width_ratios=[0.95, 1.25])
    ax_rank = fig.add_subplot(grid[0, 0])
    ax_tiles = fig.add_subplot(grid[0, 1])
    ax_summary = fig.add_subplot(grid[1, :])
    plot_target_ranking(ax_rank, target_source)
    plot_top_gene_tiles(ax_tiles, target_source, gene_source)
    plot_network_summary(ax_summary, network_summary)
    fig.suptitle(
        "GRN-level virtual perturbation prioritizes target-linked expression neighborhoods",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )

    saved = save_publication_figure(fig, args.figure_dir / "supplementary_figure_grn_virtual_perturbation")
    print(f"Wrote source data: {target_output}")
    print(f"Wrote source data: {gene_output}")
    print("Wrote figure exports:")
    for path in saved:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
