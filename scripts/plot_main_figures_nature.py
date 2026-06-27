#!/usr/bin/env python
"""Generate Nature-style redesigned main figures for the early LUAD niche project."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.figure_data import (  # noqa: E402
    CELLTYPE_ORDER,
    SPATIAL_STAGE_ORDER,
    prepare_axis_priority,
    prepare_axis_stage_heatmap,
    prepare_candidate_gene_specificity,
    prepare_evidence_heatmap,
    prepare_perturbation_dose_response,
    prepare_perturbation_method_concordance,
    prepare_perturbation_stage_loss,
    prepare_signature_celltype_heatmap,
    prepare_specificity_status_summary,
    prepare_top_perturbation_effects,
    stage_columns,
    summarize_dataset_composition,
)
from luad_niche.nature_figure_style import (  # noqa: E402
    NATURE_PALETTE,
    add_panel_label,
    apply_nature_style,
    save_publication_figure,
    semantic_axis_color,
)


STAGE_COLORS = {
    "Normal": "#CACDD2",
    "AAH": "#BBD7E8",
    "AIS": "#6ABDB9",
    "MIA": NATURE_PALETTE["gold_support"],
    "LUAD": NATURE_PALETTE["mif_axis"],
    "IAC": "#8E3B36",
    "Adjacent": "#93B57C",
    "Tumor": NATURE_PALETTE["violet_support"],
    "Primary tumor": "#9A6FA9",
    "Metastasis/effusion": "#767676",
    "Normal lymph node": "#D8D8D8",
    "LUSC": "#4D4D4D",
    "Unknown": "#EFEFEF",
}

GENE_COLORS = {
    "MIF": NATURE_PALETTE["mif_axis"],
    "CD74": "#D66A58",
    "CD44": NATURE_PALETTE["gold_support"],
    "CXCR4": NATURE_PALETTE["violet_support"],
    "SPP1": NATURE_PALETTE["macrophage_axis"],
    "TREM2": NATURE_PALETTE["green_support"],
    "PLA2G7": NATURE_PALETTE["immune_blue"],
    "IL1B": "#8F8F8F",
    "TNF": "#A8A8A8",
    "CXCL8": "#C9CCD1",
}

CELLTYPE_COLORS = {
    "Epithelial cells": NATURE_PALETTE["mif_axis"],
    "Myeloid cells": NATURE_PALETTE["macrophage_axis"],
    "Fibroblasts": NATURE_PALETTE["green_support"],
    "Endothelial cells": NATURE_PALETTE["immune_blue"],
    "MAST cells": NATURE_PALETTE["gold_support"],
    "B lymphocytes": NATURE_PALETTE["violet_support"],
    "T/NK cells": "#8F8F8F",
    "Unassigned": "#C9CCD1",
}

ROLE_DISPLAY = {
    "GSE307534": ("Visium", "spatial signal + score priority"),
    "GSE308103": ("snRNA-seq", "progression + GRN priority"),
    "GSE131907": ("scRNA ref.", "specificity audit"),
    "GSE282617": ("bulk RNA-seq", "orthogonal expression trend"),
    "GSE189357": ("scRNA-seq", "early-state context"),
    "GSE189487": ("Visium", "pilot spatial context"),
    "GSE164789": ("scRNA-seq", "tumor-adjacent context"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    table_dir = PROJECT_ROOT / "results" / "tables"
    parser.add_argument("--metadata", type=Path, default=table_dir / "geo_sample_metadata_annotated.csv")
    parser.add_argument("--evidence-matrix", type=Path, default=table_dir / "main_axis_evidence_matrix.csv")
    parser.add_argument(
        "--perturbation-ranking",
        type=Path,
        default=table_dir / "gse307534_continuous_perturbation_mia_luad_ranking.csv",
    )
    parser.add_argument(
        "--continuous-perturbation-effects",
        type=Path,
        default=table_dir / "gse307534_continuous_perturbation_effects.csv",
    )
    parser.add_argument(
        "--virtual-perturbation-ranking",
        type=Path,
        default=table_dir / "gse307534_virtual_perturbation_mia_luad_ranking.csv",
    )
    parser.add_argument(
        "--signature-celltype-summary",
        type=Path,
        default=table_dir / "gse131907_selected_signature_celltype_summary.csv",
    )
    parser.add_argument(
        "--specificity-summary",
        type=Path,
        default=table_dir / "gse189357_refined_signature_gse131907_specificity_summary.csv",
    )
    parser.add_argument("--gene-top-celltype", type=Path, default=table_dir / "gse131907_selected_gene_top_celltype.csv")
    parser.add_argument(
        "--axis-spatial-by-stage",
        type=Path,
        default=table_dir / "gse307534_candidate_axis_spatial_by_stage.csv",
    )
    parser.add_argument(
        "--paired-stats",
        type=Path,
        default=table_dir / "gse307534_candidate_axis_paired_patient_stats.csv",
    )
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures" / "nature_redesign")
    parser.add_argument("--table-dir", type=Path, default=table_dir)
    return parser.parse_args()


def add_figure_title(fig: plt.Figure, title: str, subtitle: str | None = None) -> None:
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=10.2, fontweight="bold")
    if subtitle:
        fig.text(0.012, 0.948, subtitle, ha="left", va="top", fontsize=6.2, color=NATURE_PALETTE["neutral_dark"])


def hide_all_spines(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)


def rounded_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    label: str,
    subtitle: str,
    facecolor: str,
    edgecolor: str,
    *,
    fontsize: float = 6.8,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=0.75,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height * 0.62,
        label,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight="bold",
        color=NATURE_PALETTE["neutral_black"],
    )
    ax.text(
        x + width / 2,
        y + height * 0.34,
        subtitle,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize - 1.0,
        color=NATURE_PALETTE["neutral_dark"],
    )


def arrow_between(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=0.75,
            color=NATURE_PALETTE["neutral_dark"],
        )
    )


def draw_spatial_schematic(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Specificity-audited spatial niche mining", loc="left", fontsize=7.4, fontweight="bold", pad=2)
    ax.add_patch(Rectangle((0.0, 0.03), 1.0, 0.78, transform=ax.transAxes, facecolor="#171717", edgecolor="none"))
    rng = np.random.default_rng(5)
    xs = np.linspace(0.08, 0.92, 9)
    ys = np.linspace(0.13, 0.72, 6)
    for y in ys:
        for x in xs:
            jitter_x, jitter_y = rng.normal(0, 0.008, 2)
            dist_epi = ((x - 0.40) ** 2 + (y - 0.45) ** 2) ** 0.5
            dist_myeloid = ((x - 0.60) ** 2 + (y - 0.43) ** 2) ** 0.5
            color = "#575757"
            radius = 0.012
            if dist_epi < 0.18:
                color = NATURE_PALETTE["mif_axis"]
                radius = 0.017
            elif dist_myeloid < 0.21:
                color = NATURE_PALETTE["macrophage_axis"]
                radius = 0.016
            elif 0.20 <= dist_epi < 0.31:
                color = "#E7B35A"
                radius = 0.013
            ax.add_patch(Circle((x + jitter_x, y + jitter_y), radius, transform=ax.transAxes, facecolor=color, edgecolor="white", linewidth=0.15))
    ax.add_patch(
        FancyBboxPatch(
            (0.30, 0.25),
            0.44,
            0.38,
            boxstyle="round,pad=0.01,rounding_size=0.02",
            transform=ax.transAxes,
            facecolor="none",
            edgecolor="#F6CFCB",
            linewidth=1.0,
            linestyle="--",
        )
    )
    ax.text(0.315, 0.66, "MIF-high epithelial neighborhood", transform=ax.transAxes, color="white", fontsize=5.8, ha="left")
    ax.text(0.055, 0.86, "Epithelial score", transform=ax.transAxes, fontsize=5.6, color=NATURE_PALETTE["mif_axis"])
    ax.text(0.335, 0.86, "Myeloid target score", transform=ax.transAxes, fontsize=5.6, color=NATURE_PALETTE["macrophage_axis"])
    ax.text(0.655, 0.86, "support spots", transform=ax.transAxes, fontsize=5.6, color="#E7B35A")
    add_panel_label(ax, "a", x=-0.04, y=1.02)


def plot_cohort_stage_composition(ax: plt.Axes, composition: pd.DataFrame) -> None:
    stages = stage_columns(composition)
    y = np.arange(len(composition))
    left = np.zeros(len(composition))
    for stage in stages:
        values = composition[stage].to_numpy(dtype=float)
        if values.sum() == 0:
            continue
        ax.barh(y, values, left=left, height=0.62, color=STAGE_COLORS.get(stage, "#D8D8D8"), edgecolor="white", linewidth=0.35, label=stage)
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(composition["series_accession"], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel("metadata samples", fontsize=6)
    ax.set_title("Cohort composition", loc="left", fontsize=7.2, fontweight="bold")
    max_count = float(composition["n_samples"].max())
    ax.set_xlim(0, max_count * 1.18)
    for yi, total, extent in zip(y, composition["n_samples"], composition["data_extent"]):
        ax.text(float(total) + max_count * 0.012, yi, str(extent).replace("expression matrices", "matrices"), va="center", fontsize=5.1, color=NATURE_PALETTE["neutral_dark"])
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.48), ncol=5, fontsize=5.0, handlelength=0.9, columnspacing=0.8)
    add_panel_label(ax, "c", x=-0.12)


def plot_role_map(ax: plt.Axes, composition: pd.DataFrame) -> None:
    ax.set_axis_off()
    ax.set_title("Evidence layers", loc="left", fontsize=7.2, fontweight="bold", pad=2)
    y = 0.86
    for i, row in composition.iterrows():
        dataset = row["series_accession"]
        modality, use = ROLE_DISPLAY.get(dataset, (row["modality"], row["analysis_use"]))
        if dataset == "GSE307534":
            face = NATURE_PALETTE["mif_pale"]
            edge = NATURE_PALETTE["mif_axis"]
        elif dataset == "GSE308103":
            face = NATURE_PALETTE["macrophage_pale"]
            edge = NATURE_PALETTE["macrophage_axis"]
        elif dataset == "GSE131907":
            face = NATURE_PALETTE["immune_soft"]
            edge = NATURE_PALETTE["immune_blue"]
        else:
            face = NATURE_PALETTE["neutral_pale"]
            edge = "#D0D0D0"
        ax.add_patch(FancyBboxPatch((0.00, y - 0.065), 1.0, 0.082, transform=ax.transAxes, boxstyle="round,pad=0.004,rounding_size=0.012", facecolor=face, edgecolor=edge, linewidth=0.5))
        ax.text(0.03, y, dataset, transform=ax.transAxes, ha="left", va="center", fontsize=5.8, fontweight="bold")
        ax.text(0.31, y, modality, transform=ax.transAxes, ha="left", va="center", fontsize=5.2, color=NATURE_PALETTE["neutral_dark"])
        ax.text(0.59, y, use, transform=ax.transAxes, ha="left", va="center", fontsize=5.1, color=NATURE_PALETTE["neutral_dark"])
        y -= 0.115
    add_panel_label(ax, "d", x=-0.08, y=1.02)


def draw_workflow_ladder(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Public-data workflow", loc="left", fontsize=7.2, fontweight="bold", pad=2)
    steps = [
        ("Cohorts", "7 datasets", NATURE_PALETTE["neutral_pale"], "#BDBDBD"),
        ("Audit", "GSE131907 ref.", NATURE_PALETTE["immune_soft"], NATURE_PALETTE["immune_blue"]),
        ("Spatial", "GSE307534", NATURE_PALETTE["mif_pale"], NATURE_PALETTE["mif_axis"]),
        ("Score", "score loss", NATURE_PALETTE["gold_soft"], NATURE_PALETTE["gold_support"]),
        ("GRN", "GSE308103", NATURE_PALETTE["macrophage_pale"], NATURE_PALETTE["macrophage_axis"]),
        ("Target", "MIF-CD74", NATURE_PALETTE["mif_pale"], NATURE_PALETTE["mif_axis"]),
    ]
    xs = np.linspace(0.01, 0.84, len(steps))
    box_width = 0.125
    for i, (title, subtitle, facecolor, edgecolor) in enumerate(steps):
        rounded_box(
            ax,
            (xs[i], 0.35),
            box_width,
            0.38,
            title,
            subtitle,
            facecolor,
            edgecolor,
            fontsize=5.15,
        )
        if i < len(steps) - 1:
            arrow_between(ax, (xs[i] + box_width, 0.54), (xs[i + 1] - 0.01, 0.54))
    ax.text(
        0.02,
        0.12,
        "Boundary: public-data prioritization; no causal perturbation claim.",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=5.3,
        color=NATURE_PALETTE["neutral_dark"],
    )
    add_panel_label(ax, "b", x=-0.05, y=1.02)


def plot_figure1(metadata: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    composition = summarize_dataset_composition(metadata)
    fig = plt.figure(figsize=(7.25, 5.6), constrained_layout=False)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.25], width_ratios=[1.0, 1.0], hspace=0.40, wspace=0.30)
    ax_spatial = fig.add_subplot(gs[0, 0])
    ax_workflow = fig.add_subplot(gs[0, 1])
    ax_stage = fig.add_subplot(gs[1, 0])
    ax_roles = fig.add_subplot(gs[1, 1])
    draw_spatial_schematic(ax_spatial)
    draw_workflow_ladder(ax_workflow)
    plot_cohort_stage_composition(ax_stage, composition)
    plot_role_map(ax_roles, composition)
    add_figure_title(
        fig,
        "Public multi-cohort framework for early LUAD niche mining",
        "Specificity audit, patient-aware spatial scoring, score-level target prioritization and GRN-level target-neighborhood prioritization.",
    )
    fig.subplots_adjust(left=0.08, right=0.985, top=0.88, bottom=0.12)
    return save_publication_figure(fig, args.figure_dir / "nature_figure1_workflow_dataset_composition")


def plot_priority_bars(ax: plt.Axes, priority: pd.DataFrame) -> None:
    display = priority.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [semantic_axis_color(axis_id) for axis_id in display["axis_id"]]
    ax.barh(y, display["priority_score"], height=0.62, color=colors, edgecolor="black", linewidth=0.25)
    ax.set_yticks(y)
    ax.set_yticklabels(display["axis_short_label"], fontsize=6)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("integrated priority score", fontsize=6)
    ax.set_title("Axis priority", loc="left", fontsize=7.2, fontweight="bold")
    for yi, value in zip(y, display["priority_score"]):
        ax.text(value + 0.012, yi, f"{value:.3f}", va="center", fontsize=5.4)
    add_panel_label(ax, "b", x=-0.13)


def plot_evidence_landscape(ax: plt.Axes, heatmap_table: pd.DataFrame) -> None:
    columns = ["Spatial niche", "Specificity audit", "Bulk trend", "snRNA program", "Tumor-adjacent scRNA", "Target prioritization"]
    matrix = heatmap_table[columns].to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("nature_evidence", ["#F7F7F7", "#E4CCD8", "#B64342"])
    image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns, rotation=32, ha="right", fontsize=5.4)
    ax.set_yticks(np.arange(len(heatmap_table)))
    ax.set_yticklabels(heatmap_table["axis_short_label"], fontsize=6)
    ax.set_title("Integrated evidence landscape", loc="left", fontsize=7.4, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.0, color="white" if value > 0.62 else NATURE_PALETTE["neutral_black"])
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("0-1 display score", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.2, length=2)
    hide_all_spines(ax)
    ax.tick_params(length=0)
    add_panel_label(ax, "a", x=-0.06)


def plot_perturbation_callout(ax: plt.Axes, perturbation: pd.DataFrame) -> None:
    display = perturbation.head(6).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [GENE_COLORS.get(str(g).split(",")[0], semantic_axis_color(axis_id)) for g, axis_id in zip(display["perturbed_genes"], display["axis_id"])]
    ax.barh(y, display["coupling_loss"], height=0.60, color=colors, edgecolor="black", linewidth=0.25)
    ax.set_yticks(y)
    ax.set_yticklabels(display["perturbation_label"], fontsize=5.7)
    ax.set_xlim(0, 1.06)
    ax.set_xlabel("relative coupling loss", fontsize=6)
    ax.set_title("Score-level target-prioritization callout", loc="left", fontsize=7.2, fontweight="bold")
    for yi, value in zip(y, display["coupling_loss"]):
        ax.text(value + 0.015, yi, f"{value:.2f}", va="center", fontsize=5.3)
    ax.text(0.98, -0.30, "not causal proof", transform=ax.transAxes, ha="right", va="top", fontsize=5.4, color=NATURE_PALETTE["neutral_dark"])
    add_panel_label(ax, "c", x=-0.09)


def plot_figure2(evidence_matrix: pd.DataFrame, perturbation_ranking: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    priority = prepare_axis_priority(evidence_matrix)
    heatmap = prepare_evidence_heatmap(evidence_matrix)
    perturbation = prepare_top_perturbation_effects(perturbation_ranking, top_n=8)
    fig = plt.figure(figsize=(7.25, 5.55), constrained_layout=False)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.05, 0.95], width_ratios=[1.35, 0.90], hspace=0.55, wspace=0.45)
    ax_heat = fig.add_subplot(gs[0, :])
    ax_priority = fig.add_subplot(gs[1, 0])
    ax_perturb = fig.add_subplot(gs[1, 1])
    plot_evidence_landscape(ax_heat, heatmap)
    plot_priority_bars(ax_priority, priority)
    plot_perturbation_callout(ax_perturb, perturbation)
    add_figure_title(fig, "Integrated evidence prioritizes the MIF-CD74 axis", "Multi-cohort evidence supports MIF-CD74 as a follow-up candidate while retaining SPP1/TREM2/PLA2G7 as a macrophage-state readout.")
    fig.subplots_adjust(left=0.16, right=0.98, top=0.86, bottom=0.16)
    return save_publication_figure(fig, args.figure_dir / "nature_figure2_axis_evidence_perturbation")


def plot_specificity_heatmap(ax: plt.Axes, heatmap_table: pd.DataFrame) -> None:
    celltypes = [celltype for celltype in CELLTYPE_ORDER if celltype in heatmap_table.columns]
    matrix = heatmap_table[celltypes].to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("specificity", ["#F7F7F7", "#DCEEEF", "#42949E"])
    image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(np.arange(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=32, ha="right", fontsize=5.7)
    ax.set_yticks(np.arange(len(heatmap_table)))
    ax.set_yticklabels(heatmap_table["signature_label"], fontsize=6)
    ax.set_title("Large-reference cell-type specificity", loc="left", fontsize=7.4, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if value >= 0.34:
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.0, color="white" if value > 0.68 else NATURE_PALETTE["neutral_black"])
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("relative mean score", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.2, length=2)
    hide_all_spines(ax)
    ax.tick_params(length=0)
    add_panel_label(ax, "a", x=-0.06)


def plot_specificity_audit(ax: plt.Axes, specificity: pd.DataFrame) -> None:
    display = specificity.copy().iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    left = np.zeros(len(display))
    parts = [
        ("expected_fraction", "expected", NATURE_PALETTE["macrophage_axis"]),
        ("off_target_fraction", "off-target", "#C95A49"),
        ("missing_fraction", "missing", "#C9CCD1"),
    ]
    for col, label, color in parts:
        values = display[col].to_numpy(dtype=float)
        ax.barh(y, values, left=left, height=0.66, color=color, edgecolor="white", linewidth=0.35, label=label)
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(display["signature_label"], fontsize=5.8)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("fraction of marker genes", fontsize=6)
    ax.set_title("Marker specificity audit", loc="left", fontsize=7.2, fontweight="bold")
    for yi, row in display.iterrows():
        ax.text(1.02, yi, f"{int(row['expected'])}/{int(row['total_genes'])} expected", va="center", fontsize=5.2)
    ax.legend(fontsize=5.2, ncol=3, loc="lower left", bbox_to_anchor=(0.0, -0.32))
    add_panel_label(ax, "b", x=-0.13)


def plot_gene_compartment(ax: plt.Axes, genes: pd.DataFrame) -> None:
    display = genes.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [CELLTYPE_COLORS.get(celltype, NATURE_PALETTE["neutral_mid"]) for celltype in display["Cell_type.refined"]]
    ax.barh(y, display["mean_expression"], height=0.58, color=colors, edgecolor="black", linewidth=0.25)
    ax.set_yticks(y)
    ax.set_yticklabels(display["gene"].astype(str), fontsize=6)
    ax.set_xlabel("top cell-type mean", fontsize=6)
    ax.set_title("Candidate gene compartment", loc="left", fontsize=7.2, fontweight="bold")
    xmax = max(0.5, float(display["mean_expression"].max()))
    ax.set_xlim(0, xmax * 1.25)
    for yi, value, celltype in zip(y, display["mean_expression"], display["Cell_type.refined"]):
        ax.text(value + xmax * 0.02, yi, str(celltype).replace(" cells", ""), va="center", fontsize=5.2, color=NATURE_PALETTE["neutral_dark"])
    add_panel_label(ax, "c", x=-0.14)


def plot_figure3(signature_summary: pd.DataFrame, specificity_summary: pd.DataFrame, gene_top: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    selected = [
        "epithelial_progenitor_like",
        "spp1_macrophage",
        "refined_spp1_macrophage",
        "refined_c1q_macrophage",
        "refined_inflammatory_macrophage",
    ]
    heatmap = prepare_signature_celltype_heatmap(signature_summary, selected)
    specificity = prepare_specificity_status_summary(specificity_summary)
    genes = prepare_candidate_gene_specificity(gene_top)
    fig = plt.figure(figsize=(7.25, 5.75), constrained_layout=False)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.08, 1.0], width_ratios=[1.25, 0.75], hspace=0.55, wspace=0.58)
    ax_heat = fig.add_subplot(gs[0, :])
    ax_status = fig.add_subplot(gs[1, 0])
    ax_genes = fig.add_subplot(gs[1, 1])
    plot_specificity_heatmap(ax_heat, heatmap)
    plot_specificity_audit(ax_status, specificity)
    plot_gene_compartment(ax_genes, genes)
    add_figure_title(fig, "Specificity audit reframes the original SPP1 niche", "The large LUAD scRNA reference separates a MIF epithelial signal from macrophage-state readouts and exposes off-target refined SPP1 markers.")
    fig.subplots_adjust(left=0.16, right=0.975, top=0.86, bottom=0.16)
    return save_publication_figure(fig, args.figure_dir / "nature_figure3_specificity_audit")


def plot_spatial_heatmap(ax: plt.Axes, heatmap_table: pd.DataFrame) -> None:
    stages = [stage for stage in SPATIAL_STAGE_ORDER if stage in heatmap_table.columns]
    matrix = heatmap_table[stages].to_numpy(dtype=float)
    vmax = max(0.24, float(np.nanmax(np.abs(matrix))))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_xticks(np.arange(len(stages)))
    ax.set_xticklabels(stages, fontsize=6)
    ax.set_yticks(np.arange(len(heatmap_table)))
    ax.set_yticklabels(heatmap_table["row_label"], fontsize=5.8)
    ax.set_title("Spatial enrichment delta across ordered lesions", loc="left", fontsize=7.4, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.0, color="white" if abs(value) > vmax * 0.55 else NATURE_PALETTE["neutral_black"])
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("observed - null adjacency", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.2, length=2)
    hide_all_spines(ax)
    ax.tick_params(length=0)
    add_panel_label(ax, "a", x=-0.06)


def focus_line_table(spatial: pd.DataFrame) -> pd.DataFrame:
    table = spatial[spatial["axis_id"].isin(["mif_cd74_cxcr4", "spp1_trem2_macrophage_epithelial"])].copy()
    table = table[table["stage"].isin(SPATIAL_STAGE_ORDER)].copy()
    table["stage"] = pd.Categorical(table["stage"], categories=SPATIAL_STAGE_ORDER, ordered=True)
    table["axis_short"] = table["axis_id"].map(lambda value: "MIF-CD74" if value == "mif_cd74_cxcr4" else "SPP1 readout")
    table["side"] = table["evidence_type"].map(lambda value: "source" if str(value).startswith("source") else "receptor")
    table["line_label"] = table["axis_short"] + " " + table["side"]
    return table.sort_values(["axis_id", "side", "stage"]).reset_index(drop=True)


def late_axis_summary(spatial: pd.DataFrame) -> pd.DataFrame:
    late = spatial[spatial["stage"].isin(["MIA", "LUAD"])].copy()
    late["enrichment_delta_mean"] = pd.to_numeric(late["enrichment_delta_mean"], errors="coerce")
    summary = (
        late.groupby(["axis_id", "evidence_type"], dropna=False)
        .agg(mean_late_delta=("enrichment_delta_mean", "mean"), n_samples=("n_samples", "sum"))
        .reset_index()
    )
    label_map = {
        "mif_cd74_cxcr4": "MIF",
        "spp1_trem2_macrophage_epithelial": "SPP1",
        "c1q_apoe_trem2_lgals3": "C1Q",
        "cxcl9_cxcl10_cxcr3": "CXCL9",
        "inflammatory_il1_tnf_cxcl8": "IL1B",
    }
    summary["side"] = summary["evidence_type"].map(lambda value: "source" if str(value).startswith("source") else "receptor")
    summary["bar_label"] = summary["axis_id"].map(label_map) + " " + summary["side"].map({"source": "src", "receptor": "rec"})
    return summary.sort_values("mean_late_delta", ascending=False).reset_index(drop=True)


def plot_focus_trend(ax: plt.Axes, line_table: pd.DataFrame) -> None:
    styles = {
        "MIF-CD74 source": (NATURE_PALETTE["mif_axis"], "o", "-"),
        "MIF-CD74 receptor": (NATURE_PALETTE["mif_axis"], "s", "--"),
        "SPP1 readout source": (NATURE_PALETTE["macrophage_axis"], "o", "-"),
        "SPP1 readout receptor": (NATURE_PALETTE["macrophage_axis"], "s", "--"),
    }
    stage_to_x = {stage: i for i, stage in enumerate(SPATIAL_STAGE_ORDER)}
    for label, group in line_table.groupby("line_label", sort=False):
        color, marker, linestyle = styles.get(label, (NATURE_PALETTE["neutral_mid"], "o", "-"))
        xs = [stage_to_x[str(stage)] for stage in group["stage"]]
        ax.plot(xs, group["enrichment_delta_mean"], marker=marker, linestyle=linestyle, color=color, linewidth=1.25, markersize=3.8, label=label)
    ax.axhline(0, color=NATURE_PALETTE["neutral_mid"], linewidth=0.65, linestyle=":")
    ax.set_xticks(np.arange(len(SPATIAL_STAGE_ORDER)))
    ax.set_xticklabels(SPATIAL_STAGE_ORDER, fontsize=5.8)
    ax.set_ylabel("spatial enrichment delta", fontsize=6)
    ax.set_title("Focused progression trends", loc="left", fontsize=7.2, fontweight="bold")
    ax.legend(fontsize=5.1, ncol=2, loc="upper left", handlelength=1.6)
    add_panel_label(ax, "b", x=-0.14)


def plot_late_summary(ax: plt.Axes, summary: pd.DataFrame) -> None:
    display = summary.head(7).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [semantic_axis_color(axis_id) for axis_id in display["axis_id"]]
    ax.barh(y, display["mean_late_delta"], height=0.58, color=colors, edgecolor="black", linewidth=0.25)
    ax.axvline(0, color=NATURE_PALETTE["neutral_mid"], linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(display["bar_label"], fontsize=5.6)
    ax.set_xlabel("mean MIA/LUAD delta", fontsize=6)
    ax.set_title("Late-stage axis ranking", loc="left", fontsize=7.2, fontweight="bold")
    xmax = max(0.05, float(display["mean_late_delta"].max()))
    ax.set_xlim(min(-0.02, float(display["mean_late_delta"].min()) * 1.1), xmax * 1.24)
    for yi, value in zip(y, display["mean_late_delta"]):
        ax.text(value + xmax * 0.02, yi, f"{value:.3f}", va="center", fontsize=5.1)
    add_panel_label(ax, "c", x=-0.14)


def plot_paired_callout(ax: plt.Axes, paired_stats: pd.DataFrame) -> None:
    ax.set_axis_off()
    row = paired_stats[
        paired_stats["axis_id"].eq("mif_cd74_cxcr4")
        & paired_stats["evidence_type"].eq("source_near_epithelial_progenitor")
    ].iloc[0]
    ax.add_patch(FancyBboxPatch((0.02, 0.08), 0.96, 0.78, transform=ax.transAxes, boxstyle="round,pad=0.012,rounding_size=0.02", facecolor="#F8E8E5", edgecolor=NATURE_PALETTE["mif_axis"], linewidth=0.8))
    ax.text(0.08, 0.70, "Paired-patient support", transform=ax.transAxes, fontsize=6.6, fontweight="bold", color=NATURE_PALETTE["neutral_black"])
    ax.text(0.08, 0.48, f"MIF source delta = {row['paired_difference_mean']:.3f}", transform=ax.transAxes, fontsize=6.1, color=NATURE_PALETTE["mif_axis"], fontweight="bold")
    ax.text(0.08, 0.31, f"95% CI {row['ci_95_low']:.3f} to {row['ci_95_high']:.3f}", transform=ax.transAxes, fontsize=5.7, color=NATURE_PALETTE["neutral_dark"])
    ax.text(0.08, 0.16, f"n={int(row['n_paired_patients'])} patients; BH q={row['wilcoxon_q_bh']:.2g}", transform=ax.transAxes, fontsize=5.7, color=NATURE_PALETTE["neutral_dark"])
    add_panel_label(ax, "d", x=-0.12, y=1.02)


def plot_figure4(spatial: pd.DataFrame, paired_stats: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    heatmap = prepare_axis_stage_heatmap(spatial)
    lines = focus_line_table(spatial)
    late = late_axis_summary(spatial)
    fig = plt.figure(figsize=(7.25, 6.0), constrained_layout=False)
    gs = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1.20, 0.95], width_ratios=[1.0, 0.92, 0.62], hspace=0.55, wspace=0.72)
    ax_heat = fig.add_subplot(gs[0, :])
    ax_trend = fig.add_subplot(gs[1, 0])
    ax_late = fig.add_subplot(gs[1, 1])
    ax_callout = fig.add_subplot(gs[1, 2])
    plot_spatial_heatmap(ax_heat, heatmap)
    plot_focus_trend(ax_trend, lines)
    plot_late_summary(ax_late, late)
    plot_paired_callout(ax_callout, paired_stats)
    add_figure_title(fig, "Spatial progression supports source-side MIF enrichment", "GSE307534 spot-neighborhood summaries prioritize MIF-CD74/CXCR4 while preserving the Visium resolution boundary.")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.86, bottom=0.13)
    return save_publication_figure(fig, args.figure_dir / "nature_figure4_spatial_axis_progression")


def plot_perturbation_model(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Score-level target-prioritization model", loc="left", fontsize=7.2, fontweight="bold", pad=2)
    rounded_box(ax, (0.04, 0.58), 0.22, 0.26, "Epithelial", "MIF source score", "#F8E8E5", NATURE_PALETTE["mif_axis"], fontsize=6.0)
    rounded_box(ax, (0.39, 0.58), 0.22, 0.26, "Neighbor myeloid", "CD74 target score", "#E3F1F2", NATURE_PALETTE["macrophage_axis"], fontsize=6.0)
    rounded_box(ax, (0.74, 0.58), 0.20, 0.26, "Coupling", "spatial score", "#F7F7F7", "#A8A8A8", fontsize=6.0)
    rounded_box(ax, (0.20, 0.16), 0.24, 0.25, "Score reduction", "x0.5 or x0", "#FFF4DF", "#E7B35A", fontsize=6.0)
    rounded_box(ax, (0.57, 0.16), 0.24, 0.25, "Priority", "coupling loss", "#F8E8E5", NATURE_PALETTE["mif_axis"], fontsize=6.0)
    arrow_between(ax, (0.26, 0.71), (0.39, 0.71))
    arrow_between(ax, (0.61, 0.71), (0.74, 0.71))
    arrow_between(ax, (0.44, 0.28), (0.57, 0.28))
    arrow_between(ax, (0.34, 0.41), (0.48, 0.58))
    ax.text(0.04, 0.02, "Target prioritization, not causal validation", transform=ax.transAxes, fontsize=5.4, color=NATURE_PALETTE["neutral_dark"])
    add_panel_label(ax, "a", x=-0.06, y=1.02)


def plot_retention(ax: plt.Axes, dose: pd.DataFrame) -> None:
    ax.set_title("Coupling retention after score reduction", loc="left", fontsize=7.2, fontweight="bold")
    label_y = {"MIF": 0.0, "CD74": 0.27, "CD44": 0.78, "SPP1": 0.88, "TREM2": 0.95, "CXCR4": 1.00}
    for _, group in dose.groupby("perturbed_genes", sort=False):
        gene = str(group["perturbed_genes"].iloc[0])
        group = group.sort_values("perturbation_factor", ascending=False)
        ax.plot(group["perturbation_factor"], group["coupling_remaining"], marker="o", linewidth=1.35, markersize=3.8, color=GENE_COLORS.get(gene, NATURE_PALETTE["neutral_mid"]))
        if gene in label_y:
            ax.text(-0.07, label_y[gene], gene, color=GENE_COLORS.get(gene, NATURE_PALETTE["neutral_mid"]), fontsize=5.7, va="center")
    ax.axhline(1.0, color=NATURE_PALETTE["neutral_mid"], linewidth=0.65, linestyle=":")
    ax.set_xlim(1.05, -0.18)
    ax.set_ylim(-0.04, 1.08)
    ax.set_xticks([1.0, 0.5, 0.0])
    ax.set_xticklabels(["baseline", "x0.5", "x0"], fontsize=5.8)
    ax.set_ylabel("retained coupling", fontsize=6)
    ax.set_xlabel("retained score factor", fontsize=6)
    add_panel_label(ax, "b", x=-0.13)


def plot_loss_heatmap(ax: plt.Axes, stage_loss: pd.DataFrame) -> None:
    stages = [stage for stage in SPATIAL_STAGE_ORDER if stage in set(stage_loss["stage"].astype(str))]
    row_order = list(dict.fromkeys(stage_loss["row_label"].tolist()))
    matrix = (
        stage_loss.pivot_table(index="row_label", columns="stage", values="coupling_loss", aggfunc="mean", fill_value=0.0, observed=False)
        .reindex(index=row_order, columns=stages, fill_value=0.0)
    )
    values = matrix.to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("loss", ["#F7F7F7", "#F0C0CC", NATURE_PALETTE["mif_axis"]])
    image = ax.imshow(values, aspect="auto", vmin=0, vmax=max(1.0, float(np.nanmax(values))), cmap=cmap)
    ax.set_xticks(np.arange(len(stages)))
    ax.set_xticklabels(stages, fontsize=5.8)
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=5.8)
    ax.set_title("Full-dropout loss across stages", loc="left", fontsize=7.2, fontweight="bold")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.0, color="white" if value > 0.62 else NATURE_PALETTE["neutral_black"])
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("relative coupling loss", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.2, length=2)
    hide_all_spines(ax)
    ax.tick_params(length=0)
    add_panel_label(ax, "c", x=-0.08)


def plot_concordance(ax: plt.Axes, concordance: pd.DataFrame) -> None:
    ax.set_title("Cross-method concordance", loc="left", fontsize=7.2, fontweight="bold")
    offsets = {"MIF": (-0.07, 0.0, "right"), "CD74": (0.03, 0.0, "left"), "CD44": (0.03, 0.04, "left"), "SPP1": (0.03, 0.07, "left"), "TREM2": (0.04, 0.02, "left"), "CXCR4": (0.05, -0.02, "left"), "PLA2G7": (0.03, -0.02, "left")}
    for _, row in concordance.iterrows():
        gene = str(row["perturbed_genes"])
        ax.scatter(row["continuous_coupling_loss"], row["dropout_priority_score"], s=38 if gene in {"MIF", "CD74"} else 24, color=GENE_COLORS.get(gene, NATURE_PALETTE["neutral_mid"]), edgecolor="black", linewidth=0.3, zorder=3)
        dx, dy, ha = offsets.get(gene, (0.02, 0.0, "left"))
        ax.text(row["continuous_coupling_loss"] + dx, row["dropout_priority_score"] + dy, gene, fontsize=5.6, ha=ha, va="center", color=GENE_COLORS.get(gene, NATURE_PALETTE["neutral_mid"]))
    ax.set_xlabel("continuous coupling loss", fontsize=6)
    ax.set_ylabel("top-quantile priority", fontsize=6)
    ax.set_xlim(-0.03, max(1.04, float(concordance["continuous_coupling_loss"].max()) * 1.05))
    ax.set_ylim(-0.03, max(1.62, float(concordance["dropout_priority_score"].max()) * 1.05))
    ax.grid(axis="both", color="#E6E6E6", linewidth=0.5)
    ax.text(0.98, 0.04, "Receptor-side CD74 remains prioritized", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4, color=NATURE_PALETTE["neutral_dark"])
    add_panel_label(ax, "d", x=-0.13)


def plot_figure5(continuous_ranking: pd.DataFrame, continuous_effects: pd.DataFrame, virtual_ranking: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    genes = ("MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7")
    dose = prepare_perturbation_dose_response(continuous_ranking, genes=genes)
    stage_loss = prepare_perturbation_stage_loss(continuous_effects, continuous_ranking, genes=genes)
    concordance = prepare_perturbation_method_concordance(continuous_ranking, virtual_ranking, genes=genes)
    fig = plt.figure(figsize=(7.25, 5.85), constrained_layout=False)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[0.90, 1.10], width_ratios=[1.0, 1.0], hspace=0.58, wspace=0.58)
    ax_model = fig.add_subplot(gs[0, 0])
    ax_retention = fig.add_subplot(gs[0, 1])
    ax_loss = fig.add_subplot(gs[1, 0])
    ax_concordance = fig.add_subplot(gs[1, 1])
    plot_perturbation_model(ax_model)
    plot_retention(ax_retention, dose)
    plot_loss_heatmap(ax_loss, stage_loss)
    plot_concordance(ax_concordance, concordance)
    add_figure_title(fig, "Score-level prioritization ranks receptor-side CD74", "Score reduction is used as target prioritization, with explicit boundaries against causal interpretation.")
    fig.subplots_adjust(left=0.14, right=0.98, top=0.86, bottom=0.13)
    return save_publication_figure(fig, args.figure_dir / "nature_figure5_virtual_perturbation_priority")


def main() -> int:
    apply_nature_style(font_size=7.0)
    args = parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.metadata)
    evidence_matrix = pd.read_csv(args.evidence_matrix)
    perturbation_ranking = pd.read_csv(args.perturbation_ranking)
    continuous_effects = pd.read_csv(args.continuous_perturbation_effects)
    virtual_ranking = pd.read_csv(args.virtual_perturbation_ranking)
    signature_summary = pd.read_csv(args.signature_celltype_summary)
    specificity_summary = pd.read_csv(args.specificity_summary)
    gene_top = pd.read_csv(args.gene_top_celltype)
    spatial = pd.read_csv(args.axis_spatial_by_stage)
    paired_stats = pd.read_csv(args.paired_stats)

    saved: list[Path] = []
    saved.extend(plot_figure1(metadata, args))
    saved.extend(plot_figure2(evidence_matrix, perturbation_ranking, args))
    saved.extend(plot_figure3(signature_summary, specificity_summary, gene_top, args))
    saved.extend(plot_figure4(spatial, paired_stats, args))
    saved.extend(plot_figure5(perturbation_ranking, continuous_effects, virtual_ranking, args))

    print("Wrote Nature-style figure exports:")
    for path in saved:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
