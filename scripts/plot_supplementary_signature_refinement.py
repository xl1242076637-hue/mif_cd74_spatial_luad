#!/usr/bin/env python
"""Plot the SPP1 macrophage signature specificity-refinement supplementary figure."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.figure_data import (  # noqa: E402
    prepare_signature_refinement_source,
    prepare_signature_refinement_status_summary,
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


STATUS_COLORS = {
    "expected": "#42949E",
    "off_target": "#C95A49",
    "missing": "#C9CCD1",
}

RETAINED_COLORS = {
    "Retained": "#42949E",
    "Removed": "#BDBDBD",
}

KEY_GENES = {"SPP1", "TREM2", "PLA2G7", "HBB", "FABP3"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--audit",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_signature_gse131907_specificity_audit.csv",
    )
    parser.add_argument(
        "--filtered-genes",
        type=Path,
        default=PROJECT_ROOT
        / "results"
        / "tables"
        / "gse189357_refined_signature_genes_gse131907_specificity_filtered.json",
    )
    parser.add_argument("--signature", default="spp1_macrophage_vs_other_macrophage")
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    return parser.parse_args()


def save_figure(fig: plt.Figure, base_path: Path) -> list[Path]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
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


def _load_filtered_genes(path: Path, signature: str) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    genes = data.get(signature, [])
    if not isinstance(genes, list):
        raise ValueError(f"Filtered gene entry for {signature!r} is not a list")
    return [str(gene) for gene in genes]


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.04, label, transform=ax.transAxes, fontsize=8, fontweight="bold")


def plot_status_bars(ax: plt.Axes, summary: pd.DataFrame) -> None:
    rows = [
        ("Original refined\nsignature", "original_n"),
        ("After specificity\naudit", "retained_n"),
    ]
    status_order = ["expected", "off_target", "missing"]
    for y, (label, value_column) in enumerate(rows):
        left = 0
        total = int(summary[value_column].sum())
        for status in status_order:
            value = int(summary.loc[summary["specificity_status"].eq(status), value_column].sum())
            if value == 0:
                continue
            ax.barh(
                y,
                value,
                left=left,
                color=STATUS_COLORS[status],
                edgecolor="white",
                linewidth=0.4,
                label=status.replace("_", " ") if y == 0 else None,
            )
            ax.text(left + value / 2, y, str(value), ha="center", va="center", fontsize=6.2, color="white")
            left += value
        ax.text(left + 0.6, y, f"n={total}", va="center", fontsize=6.4, color="#4D4D4D")
    ax.set_yticks(range(len(rows)), [item[0] for item in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Number of genes")
    ax.set_title("Specificity audit removes off-target genes", loc="left", fontsize=7.4, fontweight="bold")
    ax.legend(loc="lower right", fontsize=6.0)
    _panel_label(ax, "a")


def plot_gene_rank_strip(ax: plt.Axes, source: pd.DataFrame) -> None:
    y_map = {"expected": 0, "off_target": 1, "missing": 2}
    y_values = source["specificity_status"].map(y_map).fillna(3).astype(float)
    colors = [STATUS_COLORS.get(status, "#8F8F8F") for status in source["specificity_status"]]
    sizes = np.where(source["retained_after_audit"], 50, 30)
    edgecolors = np.where(source["retained_after_audit"], "#1F1F1F", "#FFFFFF")
    ax.scatter(
        source["display_rank"],
        y_values,
        s=sizes,
        c=colors,
        edgecolors=edgecolors,
        linewidths=0.7,
        zorder=3,
    )
    for _, row in source[source["gene"].isin(KEY_GENES)].iterrows():
        ax.text(
            row["display_rank"],
            y_map.get(row["specificity_status"], 3) + 0.18,
            row["gene"],
            ha="center",
            va="bottom",
            fontsize=5.7,
            rotation=45,
        )
    ax.set_yticks([0, 1, 2], ["expected\nmyeloid", "off-target", "missing"])
    ax.set_ylim(2.55, -0.45)
    ax.set_xlim(0.3, max(source["display_rank"]) + 0.7)
    ax.set_xlabel("Original marker rank")
    ax.set_title("Gene-level audit of the original SPP1 macrophage signature", loc="left", fontsize=7.4, fontweight="bold")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    _panel_label(ax, "b")


def plot_top_celltype_summary(ax: plt.Axes, celltype_source: pd.DataFrame) -> None:
    totals = celltype_source.groupby("top_celltype")["n_genes"].sum().sort_values(ascending=False)
    ordered = totals.index.tolist()
    y = np.arange(len(ordered))
    left = np.zeros(len(ordered))
    for label in ("Retained", "Removed"):
        values = (
            celltype_source[celltype_source["retained_label"].eq(label)]
            .set_index("top_celltype")
            .reindex(ordered)["n_genes"]
            .fillna(0)
            .to_numpy(dtype=float)
        )
        ax.barh(
            y,
            values,
            left=left,
            color=RETAINED_COLORS[label],
            edgecolor="white",
            linewidth=0.4,
            label=label,
        )
        left += values
    ax.set_yticks(y, ordered)
    ax.invert_yaxis()
    ax.set_xlabel("Number of genes")
    ax.set_title("Removed genes expose mixed lineage signal", loc="left", fontsize=7.4, fontweight="bold")
    ax.legend(loc="lower right", fontsize=6.0)
    _panel_label(ax, "c")


def main() -> int:
    args = parse_args()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    audit = pd.read_csv(args.audit)
    filtered_genes = _load_filtered_genes(args.filtered_genes, args.signature)
    source = prepare_signature_refinement_source(audit, filtered_genes, signature=args.signature)
    status = prepare_signature_refinement_status_summary(source)
    celltype_source = (
        source.groupby(["top_celltype", "retained_label"], dropna=False)
        .size()
        .rename("n_genes")
        .reset_index()
    )

    gene_output = args.table_dir / "supplementary_figure_spp1_signature_refinement_gene_source.csv"
    status_output = args.table_dir / "supplementary_figure_spp1_signature_refinement_status_source.csv"
    celltype_output = args.table_dir / "supplementary_figure_spp1_signature_refinement_celltype_source.csv"
    source.to_csv(gene_output, index=False, encoding="utf-8-sig")
    status.to_csv(status_output, index=False, encoding="utf-8-sig")
    celltype_source.to_csv(celltype_output, index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    grid = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[0.82, 1.18], width_ratios=[0.95, 1.05])
    ax_status = fig.add_subplot(grid[0, 0])
    ax_celltype = fig.add_subplot(grid[0, 1])
    ax_gene = fig.add_subplot(grid[1, :])

    plot_status_bars(ax_status, status)
    plot_top_celltype_summary(ax_celltype, celltype_source)
    plot_gene_rank_strip(ax_gene, source)
    fig.suptitle(
        "Specificity audit reframes the original SPP1 macrophage signature",
        x=0.02,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )

    saved = save_figure(fig, args.figure_dir / "supplementary_figure_spp1_signature_refinement")
    print(f"Wrote source data: {gene_output}")
    print(f"Wrote source data: {status_output}")
    print(f"Wrote source data: {celltype_output}")
    print("Wrote figure exports:")
    for path in saved:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
