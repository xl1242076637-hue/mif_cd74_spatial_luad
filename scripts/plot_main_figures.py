#!/usr/bin/env python
"""Generate manuscript draft figures for the early LUAD niche project."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
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


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


STAGE_COLORS = {
    "Normal": "#C9CCD1",
    "AAH": "#A6CEE3",
    "AIS": "#55B7B3",
    "MIA": "#E7B35A",
    "LUAD": "#C95A49",
    "IAC": "#8F3B2E",
    "Adjacent": "#8CB369",
    "Tumor": "#7C6CCF",
    "Primary tumor": "#9A4D8E",
    "Metastasis/effusion": "#767676",
    "Normal lymph node": "#D8D8D8",
    "LUSC": "#4D4D4D",
    "Unknown": "#EFEFEF",
}

GRADE_COLORS = {
    "lead": "#B64342",
    "strong": "#E28E2C",
    "supporting": "#5B8FD6",
    "benchmark_or_secondary": "#8F8F8F",
}

AXIS_COLORS = {
    "mif_cd74_cxcr4": "#B64342",
    "spp1_trem2_macrophage_epithelial": "#42949E",
    "cxcl9_cxcl10_cxcr3": "#5B8FD6",
    "c1q_apoe_trem2_lgals3": "#7BAA5B",
    "inflammatory_il1_tnf_cxcl8": "#8F8F8F",
}

CELLTYPE_COLORS = {
    "Epithelial cells": "#B64342",
    "Myeloid cells": "#42949E",
    "Fibroblasts": "#7BAA5B",
    "Endothelial cells": "#5B8FD6",
    "MAST cells": "#E28E2C",
    "B lymphocytes": "#7C6CCF",
    "T/NK cells": "#8F8F8F",
    "Unassigned": "#C9CCD1",
}

SPECIFICITY_COLORS = {
    "expected_fraction": "#42949E",
    "off_target_fraction": "#C95A49",
    "missing_fraction": "#C9CCD1",
}

GENE_COLORS = {
    "MIF": "#B64342",
    "CD74": "#C95A49",
    "CD44": "#E28E2C",
    "CXCR4": "#7C6CCF",
    "SPP1": "#42949E",
    "TREM2": "#7BAA5B",
    "PLA2G7": "#5B8FD6",
}

ROLE_TABLE_DISPLAY = {
    "GSE307534": ("Visium", "spatial axis + target prioritization"),
    "GSE308103": ("snRNA-seq", "cell-state progression"),
    "GSE131907": ("scRNA ref.", "specificity filtering"),
    "GSE282617": ("bulk RNA-seq", "marker trends"),
    "GSE189357": ("scRNA-seq", "state refinement"),
    "GSE189487": ("Visium", "pilot spatial niche"),
    "GSE164789": ("scRNA-seq", "tumor-adjacent contrast"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
    )
    parser.add_argument(
        "--evidence-matrix",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "main_axis_evidence_matrix.csv",
    )
    parser.add_argument(
        "--perturbation-ranking",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_continuous_perturbation_mia_luad_ranking.csv",
    )
    parser.add_argument(
        "--continuous-perturbation-effects",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_continuous_perturbation_effects.csv",
    )
    parser.add_argument(
        "--virtual-perturbation-ranking",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_virtual_perturbation_mia_luad_ranking.csv",
    )
    parser.add_argument(
        "--signature-celltype-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse131907_selected_signature_celltype_summary.csv",
    )
    parser.add_argument(
        "--specificity-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_signature_gse131907_specificity_summary.csv",
    )
    parser.add_argument(
        "--gene-top-celltype",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse131907_selected_gene_top_celltype.csv",
    )
    parser.add_argument(
        "--axis-spatial-by-stage",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_candidate_axis_spatial_by_stage.csv",
    )
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    return parser.parse_args()


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.05, y: float = 1.02) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
    )


def save_figure(fig: plt.Figure, base_path: Path) -> list[Path]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    for suffix, kwargs in {
        ".svg": {},
        ".pdf": {},
        ".tiff": {"dpi": 600},
        ".png": {"dpi": 300},
    }.items():
        out_path = base_path.with_suffix(suffix)
        fig.savefig(out_path, bbox_inches="tight", **kwargs)
        saved.append(out_path)
    plt.close(fig)
    return saved


def draw_workflow(ax: plt.Axes) -> None:
    ax.set_axis_off()
    steps = [
        ("Public cohorts", "spatial, snRNA,\nscRNA and bulk"),
        ("Specificity audit", "GSE131907\nreference cells"),
        ("Spatial scoring", "epithelial\nmyeloid niche"),
        ("Evidence ranking", "multi-cohort\naxis priority"),
        ("Perturbation", "score-level\ncoupling loss"),
        ("Candidate", "MIF-CD74\nplus SPP1 readout"),
    ]
    x_positions = np.linspace(0.03, 0.84, len(steps))
    width = 0.13
    y = 0.38
    height = 0.42
    for i, ((title, subtitle), x) in enumerate(zip(steps, x_positions)):
        color = "#F7F7F7" if i < len(steps) - 1 else "#F8E8E5"
        edge = "#B64342" if i == len(steps) - 1 else "#4D4D4D"
        ax.add_patch(Rectangle((x, y), width, height, facecolor=color, edgecolor=edge, linewidth=0.8))
        ax.text(x + width / 2, y + height * 0.65, title, ha="center", va="center", fontweight="bold", fontsize=7.2)
        ax.text(x + width / 2, y + height * 0.32, subtitle, ha="center", va="center", fontsize=6.2, color="#4D4D4D")
        if i < len(steps) - 1:
            ax.annotate(
                "",
                xy=(x_positions[i + 1] - 0.008, y + height / 2),
                xytext=(x + width + 0.008, y + height / 2),
                arrowprops={"arrowstyle": "-|>", "lw": 0.8, "color": "#4D4D4D"},
            )
    ax.text(
        0.03,
        0.12,
        "Progression frame: Normal -> AAH -> AIS -> MIA -> LUAD, with IAC/bulk and tumor-adjacent cohorts as support.",
        ha="left",
        va="center",
        fontsize=6.6,
        color="#4D4D4D",
    )
    add_panel_label(ax, "a", x=0.0, y=0.98)


def plot_stage_composition(ax: plt.Axes, composition: pd.DataFrame) -> None:
    stages = stage_columns(composition)
    y = np.arange(len(composition))
    left = np.zeros(len(composition))
    for stage in stages:
        values = composition[stage].to_numpy(dtype=float)
        if values.sum() == 0:
            continue
        ax.barh(
            y,
            values,
            left=left,
            color=STAGE_COLORS.get(stage, "#D8D8D8"),
            edgecolor="white",
            linewidth=0.3,
            label=stage,
            height=0.72,
        )
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(composition["series_accession"])
    ax.invert_yaxis()
    ax.set_xlabel("metadata samples")
    ax.set_title("Public cohort stage composition", loc="left", fontsize=7.5, fontweight="bold")
    max_count = float(composition["n_samples"].max())
    ax.set_xlim(0, max_count * 1.16)
    for yi, total, extent in zip(y, composition["n_samples"], composition["data_extent"]):
        ax.text(float(total) + max_count * 0.015, yi, str(extent), va="center", ha="left", fontsize=5.8, color="#4D4D4D")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=4, fontsize=5.8, handlelength=1.0)
    add_panel_label(ax, "b")


def plot_dataset_roles(ax: plt.Axes, composition: pd.DataFrame) -> None:
    ax.set_axis_off()
    ax.set_title("Evidence roles", loc="left", fontsize=7.5, fontweight="bold", pad=2)
    headers = ["Dataset", "Modality", "Use"]
    x_positions = [0.02, 0.31, 0.66]
    for x, header in zip(x_positions, headers):
        ax.text(x, 0.95, header, transform=ax.transAxes, ha="left", va="top", fontsize=6.5, fontweight="bold")
    y = 0.84
    step = 0.105
    for i, row in composition.iterrows():
        modality, analysis_use = ROLE_TABLE_DISPLAY.get(
            row["series_accession"],
            (row["modality"], row["analysis_use"]),
        )
        if i % 2 == 0:
            ax.add_patch(Rectangle((0.0, y - 0.075), 1.0, 0.092, transform=ax.transAxes, facecolor="#F7F7F7", edgecolor="none"))
        ax.text(0.02, y, row["series_accession"], transform=ax.transAxes, ha="left", va="top", fontsize=6.2, fontweight="bold")
        ax.text(0.31, y, modality, transform=ax.transAxes, ha="left", va="top", fontsize=5.8, color="#4D4D4D")
        ax.text(0.66, y, analysis_use, transform=ax.transAxes, ha="left", va="top", fontsize=5.6, color="#4D4D4D")
        y -= step
    add_panel_label(ax, "c", x=0.0, y=1.02)


def plot_figure1(metadata: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    composition = summarize_dataset_composition(metadata)
    composition.to_csv(args.table_dir / "figure1_dataset_composition_source.csv", index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 5.1), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[0.9, 1.55], width_ratios=[1.2, 1.0])
    ax_workflow = fig.add_subplot(gs[0, :])
    ax_stage = fig.add_subplot(gs[1, 0])
    ax_roles = fig.add_subplot(gs[1, 1])

    draw_workflow(ax_workflow)
    plot_stage_composition(ax_stage, composition)
    plot_dataset_roles(ax_roles, composition)
    fig.suptitle("Public-data framework for early LUAD epithelial-myeloid niche mining", x=0.02, ha="left", fontsize=9, fontweight="bold")
    return save_figure(fig, args.figure_dir / "figure1_workflow_dataset_composition")


def plot_priority(ax: plt.Axes, priority: pd.DataFrame) -> None:
    display = priority.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [GRADE_COLORS.get(grade, "#8F8F8F") for grade in display["evidence_grade"]]
    ax.barh(y, display["priority_score"], color=colors, edgecolor="black", linewidth=0.3, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(display["axis_short_label"])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("integrated priority score")
    ax.set_title("Axis ranking", loc="left", fontsize=7.5, fontweight="bold")
    for yi, value in zip(y, display["priority_score"]):
        ax.text(value + 0.015, yi, f"{value:.3f}", va="center", ha="left", fontsize=6.2)
    add_panel_label(ax, "a")


def plot_evidence_heatmap(ax: plt.Axes, heatmap_table: pd.DataFrame) -> None:
    component_columns = [
        "Spatial niche",
        "Specificity audit",
        "Bulk trend",
        "snRNA program",
        "Tumor-adjacent scRNA",
        "Target prioritization",
    ]
    matrix = heatmap_table[component_columns].to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("evidence", ["#F7F7F7", "#B4C0E4", "#0F4D92"])
    image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(np.arange(len(component_columns)))
    ax.set_xticklabels(component_columns, rotation=35, ha="right", fontsize=5.8)
    ax.set_yticks(np.arange(len(heatmap_table)))
    ax.set_yticklabels(heatmap_table["axis_short_label"], fontsize=6.0)
    ax.set_title("Normalized evidence components", loc="left", fontsize=7.5, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            color = "white" if value > 0.65 else "#272727"
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.4, color=color)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("0-1 display score", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    add_panel_label(ax, "b")


def plot_perturbation(ax: plt.Axes, perturbation: pd.DataFrame) -> None:
    display = perturbation.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [AXIS_COLORS.get(axis_id, "#8F8F8F") for axis_id in display["axis_id"]]
    ax.barh(y, display["coupling_loss"], color=colors, edgecolor="black", linewidth=0.3, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(display["perturbation_label"], fontsize=6.2)
    ax.set_xlim(0, max(1.02, float(display["coupling_loss"].max()) * 1.12))
    ax.set_xlabel("relative coupling loss after score dropout")
    ax.set_title("Continuous score-level target-prioritization effects", loc="left", fontsize=7.5, fontweight="bold")
    ax.axvline(0, color="#4D4D4D", linewidth=0.7)
    for yi, value, n_samples in zip(y, display["coupling_loss"], display["n_samples"]):
        ax.text(value + 0.015, yi, f"{value:.2f}; n={int(n_samples)}", va="center", ha="left", fontsize=5.8)
    ax.text(
        0.99,
        -0.2,
        "Interpretation: target-prioritization score, not causal proof",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.8,
        color="#4D4D4D",
    )
    add_panel_label(ax, "c")


def plot_figure2(evidence_matrix: pd.DataFrame, perturbation_ranking: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    priority = prepare_axis_priority(evidence_matrix)
    heatmap_table = prepare_evidence_heatmap(evidence_matrix)
    perturbation = prepare_top_perturbation_effects(perturbation_ranking, top_n=8)

    priority.to_csv(args.table_dir / "figure2_priority_source.csv", index=False, encoding="utf-8-sig")
    heatmap_table.to_csv(args.table_dir / "figure2_evidence_heatmap_source.csv", index=False, encoding="utf-8-sig")
    perturbation.to_csv(args.table_dir / "figure2_perturbation_source.csv", index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 5.6), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05], width_ratios=[0.82, 1.18])
    ax_priority = fig.add_subplot(gs[0, 0])
    ax_heatmap = fig.add_subplot(gs[0, 1])
    ax_perturb = fig.add_subplot(gs[1, :])

    plot_priority(ax_priority, priority)
    plot_evidence_heatmap(ax_heatmap, heatmap_table)
    plot_perturbation(ax_perturb, perturbation)
    fig.suptitle("Integrated evidence nominates MIF-CD74 for early LUAD follow-up", x=0.02, ha="left", fontsize=9, fontweight="bold")
    return save_figure(fig, args.figure_dir / "figure2_axis_evidence_perturbation")


def plot_signature_celltype_heatmap(ax: plt.Axes, heatmap_table: pd.DataFrame) -> None:
    celltypes = [celltype for celltype in CELLTYPE_ORDER if celltype in heatmap_table.columns]
    matrix = heatmap_table[celltypes].to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("signature_specificity", ["#F7F7F7", "#B4C0E4", "#0F4D92"])
    image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_xticks(np.arange(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=35, ha="right", fontsize=5.8)
    ax.set_yticks(np.arange(len(heatmap_table)))
    ax.set_yticklabels(heatmap_table["signature_label"], fontsize=6.2)
    ax.set_title("GSE131907 cell-type specificity of selected signatures", loc="left", fontsize=7.5, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if value >= 0.35:
                color = "white" if value > 0.68 else "#272727"
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.1, color=color)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("relative mean score", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    add_panel_label(ax, "a")


def plot_specificity_status(ax: plt.Axes, specificity: pd.DataFrame) -> None:
    display = specificity.copy().iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    left = np.zeros(len(display))
    labels = {
        "expected_fraction": "expected top cell type",
        "off_target_fraction": "off-target top cell type",
        "missing_fraction": "missing in reference",
    }
    for column in ["expected_fraction", "off_target_fraction", "missing_fraction"]:
        values = display[column].to_numpy(dtype=float)
        ax.barh(
            y,
            values,
            left=left,
            color=SPECIFICITY_COLORS[column],
            edgecolor="white",
            linewidth=0.4,
            height=0.72,
            label=labels[column],
        )
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(display["signature_label"], fontsize=6.2)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of marker genes")
    ax.set_title(
        "Specificity audit (teal=expected, red=off-target, gray=missing)",
        loc="left",
        fontsize=7.2,
        fontweight="bold",
    )
    for yi, expected, off_target, total in zip(y, display["expected"], display["off_target"], display["total_genes"]):
        ax.text(1.02, yi, f"{int(expected)}/{int(total)} expected; {int(off_target)} off-target", va="center", fontsize=5.6)
    add_panel_label(ax, "b", x=-0.12)


def plot_candidate_gene_specificity(ax: plt.Axes, genes: pd.DataFrame) -> None:
    display = genes.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [CELLTYPE_COLORS.get(celltype, "#8F8F8F") for celltype in display["Cell_type.refined"]]
    ax.barh(y, display["mean_expression"], color=colors, edgecolor="black", linewidth=0.3, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(display["gene"].astype(str), fontsize=6.5)
    ax.set_xlabel("top cell-type mean expression")
    ax.set_title("Top-expressing cell type for candidate genes", loc="left", fontsize=7.5, fontweight="bold")
    xmax = max(0.5, float(display["mean_expression"].max()))
    ax.set_xlim(0, xmax * 1.24)
    for yi, value, celltype in zip(y, display["mean_expression"], display["Cell_type.refined"]):
        ax.text(value + xmax * 0.025, yi, str(celltype).replace(" cells", ""), va="center", fontsize=5.8, color="#4D4D4D")
    add_panel_label(ax, "c", x=-0.12)


def plot_figure3(
    signature_celltype_summary: pd.DataFrame,
    specificity_summary: pd.DataFrame,
    gene_top_celltype: pd.DataFrame,
    args: argparse.Namespace,
) -> list[Path]:
    selected_signatures = [
        "epithelial_progenitor_like",
        "spp1_macrophage",
        "refined_spp1_macrophage",
        "refined_c1q_macrophage",
        "refined_inflammatory_macrophage",
    ]
    heatmap_table = prepare_signature_celltype_heatmap(signature_celltype_summary, selected_signatures)
    specificity = prepare_specificity_status_summary(specificity_summary)
    genes = prepare_candidate_gene_specificity(gene_top_celltype)

    heatmap_table.to_csv(args.table_dir / "figure3_signature_celltype_heatmap_source.csv", index=False, encoding="utf-8-sig")
    specificity.to_csv(args.table_dir / "figure3_specificity_status_source.csv", index=False, encoding="utf-8-sig")
    genes.to_csv(args.table_dir / "figure3_candidate_gene_specificity_source.csv", index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 5.9), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05], width_ratios=[1.25, 0.85])
    ax_heatmap = fig.add_subplot(gs[0, :])
    ax_status = fig.add_subplot(gs[1, 0])
    ax_genes = fig.add_subplot(gs[1, 1])

    plot_signature_celltype_heatmap(ax_heatmap, heatmap_table)
    plot_specificity_status(ax_status, specificity)
    plot_candidate_gene_specificity(ax_genes, genes)
    fig.suptitle("Specificity audit reframes the original SPP1-niche hypothesis", x=0.02, ha="left", fontsize=9, fontweight="bold")
    return save_figure(fig, args.figure_dir / "figure3_specificity_audit")


def plot_axis_stage_heatmap(ax: plt.Axes, heatmap_table: pd.DataFrame) -> None:
    stages = [stage for stage in SPATIAL_STAGE_ORDER if stage in heatmap_table.columns]
    matrix = heatmap_table[stages].to_numpy(dtype=float)
    vmax = max(0.2, float(np.nanmax(np.abs(matrix))))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_xticks(np.arange(len(stages)))
    ax.set_xticklabels(stages, fontsize=6.2)
    ax.set_yticks(np.arange(len(heatmap_table)))
    ax.set_yticklabels(heatmap_table["row_label"], fontsize=5.8)
    ax.set_title("GSE307534 spatial enrichment delta by stage", loc="left", fontsize=7.5, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            color = "white" if abs(value) > vmax * 0.55 else "#272727"
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.1, color=color)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("observed - null adjacency", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    add_panel_label(ax, "a")


def _line_table_for_axes(spatial_by_stage: pd.DataFrame) -> pd.DataFrame:
    focus_axes = ["mif_cd74_cxcr4", "spp1_trem2_macrophage_epithelial"]
    table = spatial_by_stage[spatial_by_stage["axis_id"].isin(focus_axes)].copy()
    table = table[table["stage"].isin(SPATIAL_STAGE_ORDER)].copy()
    table["stage"] = pd.Categorical(table["stage"], categories=SPATIAL_STAGE_ORDER, ordered=True)
    table["axis_short_label"] = table["axis_id"].map(lambda value: "MIF-CD74" if value == "mif_cd74_cxcr4" else "SPP1 readout")
    table["evidence_short"] = table["evidence_type"].map(lambda value: "source" if str(value).startswith("source") else "receptor")
    table["line_label"] = table["axis_short_label"] + " " + table["evidence_short"]
    return table.sort_values(["axis_id", "evidence_short", "stage"]).reset_index(drop=True)


def plot_focus_axis_trends(ax: plt.Axes, line_table: pd.DataFrame) -> None:
    line_styles = {
        "MIF-CD74 source": ("#B64342", "o", "-"),
        "MIF-CD74 receptor": ("#B64342", "s", "--"),
        "SPP1 readout source": ("#42949E", "o", "-"),
        "SPP1 readout receptor": ("#42949E", "s", "--"),
    }
    stage_to_x = {stage: i for i, stage in enumerate(SPATIAL_STAGE_ORDER)}
    for label, group in line_table.groupby("line_label", sort=False):
        color, marker, linestyle = line_styles.get(label, ("#767676", "o", "-"))
        xs = [stage_to_x[str(stage)] for stage in group["stage"]]
        ys = group["enrichment_delta_mean"].to_numpy(dtype=float)
        ax.plot(xs, ys, color=color, marker=marker, linestyle=linestyle, linewidth=1.4, markersize=4.2, label=label)
    ax.axhline(0, color="#767676", linewidth=0.7, linestyle=":")
    ax.set_xticks(np.arange(len(SPATIAL_STAGE_ORDER)))
    ax.set_xticklabels(SPATIAL_STAGE_ORDER, fontsize=6.2)
    ax.set_ylabel("spatial enrichment delta")
    ax.set_title("Focus-axis stage trends", loc="left", fontsize=7.5, fontweight="bold")
    ax.legend(fontsize=5.6, ncol=2, loc="upper left")
    add_panel_label(ax, "b", x=-0.13)


def _late_axis_summary(spatial_by_stage: pd.DataFrame) -> pd.DataFrame:
    late = spatial_by_stage[spatial_by_stage["stage"].isin(["MIA", "LUAD"])].copy()
    late["enrichment_delta_mean"] = pd.to_numeric(late["enrichment_delta_mean"], errors="coerce")
    summary = (
        late.groupby(["axis_id", "evidence_type"], dropna=False)
        .agg(
            mean_late_delta=("enrichment_delta_mean", "mean"),
            n_stage_rows=("stage", "size"),
            n_samples=("n_samples", "sum"),
        )
        .reset_index()
    )
    summary["axis_short_label"] = summary["axis_id"].map(lambda value: AXIS_COLORS.get(value, "#8F8F8F"))
    summary["axis_label"] = summary["axis_id"].map(
        {
            "mif_cd74_cxcr4": "MIF-CD74/CXCR4",
            "spp1_trem2_macrophage_epithelial": "SPP1/TREM2/PLA2G7",
            "c1q_apoe_trem2_lgals3": "C1Q/APOE/TREM2",
            "cxcl9_cxcl10_cxcr3": "CXCL9/10-CXCR3",
            "inflammatory_il1_tnf_cxcl8": "IL1B/TNF/CXCL8",
        }
    )
    summary["evidence_short"] = summary["evidence_type"].map(lambda value: "source" if str(value).startswith("source") else "receptor")
    summary["bar_label"] = summary["axis_label"] + " " + summary["evidence_short"]
    return summary.sort_values("mean_late_delta", ascending=False).reset_index(drop=True)


def plot_late_axis_summary(ax: plt.Axes, late_summary: pd.DataFrame) -> None:
    display = late_summary.head(8).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(display))
    colors = [AXIS_COLORS.get(axis_id, "#8F8F8F") for axis_id in display["axis_id"]]
    ax.barh(y, display["mean_late_delta"], color=colors, edgecolor="black", linewidth=0.3, height=0.68)
    ax.axvline(0, color="#767676", linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(display["bar_label"], fontsize=5.9)
    ax.set_xlabel("mean MIA/LUAD enrichment delta")
    ax.set_title("Late-stage spatial-axis summary", loc="left", fontsize=7.5, fontweight="bold")
    xmax = max(0.05, float(display["mean_late_delta"].max()))
    ax.set_xlim(min(-0.02, float(display["mean_late_delta"].min()) * 1.1), xmax * 1.22)
    for yi, value in zip(y, display["mean_late_delta"]):
        ax.text(value + xmax * 0.025, yi, f"{value:.3f}", va="center", fontsize=5.6)
    add_panel_label(ax, "c", x=-0.13)


def plot_figure4(spatial_by_stage: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    heatmap_table = prepare_axis_stage_heatmap(spatial_by_stage)
    line_table = _line_table_for_axes(spatial_by_stage)
    late_summary = _late_axis_summary(spatial_by_stage)

    heatmap_table.to_csv(args.table_dir / "figure4_axis_stage_heatmap_source.csv", index=False, encoding="utf-8-sig")
    line_table.to_csv(args.table_dir / "figure4_focus_axis_trend_source.csv", index=False, encoding="utf-8-sig")
    late_summary.to_csv(args.table_dir / "figure4_late_axis_summary_source.csv", index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.2, 1.0], width_ratios=[1.0, 1.0])
    ax_heatmap = fig.add_subplot(gs[0, :])
    ax_trend = fig.add_subplot(gs[1, 0])
    ax_late = fig.add_subplot(gs[1, 1])

    plot_axis_stage_heatmap(ax_heatmap, heatmap_table)
    plot_focus_axis_trends(ax_trend, line_table)
    plot_late_axis_summary(ax_late, late_summary)
    fig.suptitle("Spatial progression evidence supports MIF-CD74 epithelial-myeloid coupling", x=0.02, ha="left", fontsize=9, fontweight="bold")
    return save_figure(fig, args.figure_dir / "figure4_spatial_axis_progression")


def plot_perturbation_schematic(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Score-level target-prioritization model", loc="left", fontsize=7.5, fontweight="bold", pad=2)
    boxes = [
        (0.04, 0.58, 0.23, 0.24, "Epithelial\nprogenitor-like\nspot score", "#F8E8E5"),
        (0.38, 0.58, 0.23, 0.24, "Neighbor\nmyeloid target\nscore", "#E3F1F2"),
        (0.72, 0.58, 0.23, 0.24, "Continuous\ncoupling\nscore", "#F7F7F7"),
        (0.21, 0.18, 0.23, 0.22, "Score\nx0.5 or x0\nreduction", "#FFF4DF"),
        (0.56, 0.18, 0.23, 0.22, "Coupling\nloss and target\npriority", "#F7F7F7"),
    ]
    for x, y, w, h, label, color in boxes:
        ax.add_patch(Rectangle((x, y), w, h, transform=ax.transAxes, facecolor=color, edgecolor="#4D4D4D", linewidth=0.8))
        ax.text(x + w / 2, y + h / 2, label, transform=ax.transAxes, ha="center", va="center", fontsize=6.3)
    arrows = [
        ((0.27, 0.70), (0.38, 0.70)),
        ((0.61, 0.70), (0.72, 0.70)),
        ((0.44, 0.29), (0.56, 0.29)),
        ((0.35, 0.40), (0.46, 0.58)),
    ]
    for start, end in arrows:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            xycoords=ax.transAxes,
            textcoords=ax.transAxes,
            arrowprops={"arrowstyle": "-|>", "lw": 0.8, "color": "#4D4D4D"},
        )
    ax.text(
        0.04,
        0.03,
        "Interpretation boundary: score-level in-silico target prioritization, not causal wet-lab validation.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.8,
        color="#4D4D4D",
    )
    add_panel_label(ax, "a", x=-0.06, y=1.02)


def plot_dose_response(ax: plt.Axes, dose: pd.DataFrame) -> None:
    ax.set_title("Dose-response-like coupling retention", loc="left", fontsize=7.5, fontweight="bold")
    label_y = {
        "PLA2G7": 1.035,
        "CXCR4": 0.995,
        "TREM2": 0.955,
        "SPP1": 0.885,
        "CD44": 0.73,
        "CD74": 0.26,
        "MIF": 0.0,
    }
    for label, group in dose.groupby("line_label", sort=False):
        gene = str(group["perturbed_genes"].iloc[0])
        group = group.sort_values("perturbation_factor", ascending=False)
        ax.plot(
            group["perturbation_factor"],
            group["coupling_remaining"],
            marker="o",
            linewidth=1.5,
            markersize=4,
            color=GENE_COLORS.get(gene, "#767676"),
            label=label,
        )
        last = group.sort_values("perturbation_factor").iloc[0]
        if gene in {"MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2"}:
            ax.text(
                -0.04,
                label_y.get(gene, float(last["coupling_remaining"])),
                gene,
                ha="left",
                va="center",
                fontsize=5.7,
                color=GENE_COLORS.get(gene, "#767676"),
            )
    ax.axhline(1.0, color="#767676", linewidth=0.7, linestyle=":")
    ax.set_xlim(1.05, -0.18)
    ax.set_ylim(-0.04, 1.10)
    ax.set_xticks([1.0, 0.5, 0.0])
    ax.set_xticklabels(["baseline", "x0.5", "x0"])
    ax.set_ylabel("coupling retained vs baseline")
    ax.set_xlabel("retained score factor")
    add_panel_label(ax, "b", x=-0.12)


def plot_stage_loss_heatmap(ax: plt.Axes, stage_loss: pd.DataFrame) -> None:
    ax.set_title("Full-dropout coupling loss by stage", loc="left", fontsize=7.5, fontweight="bold")
    stages = [stage for stage in SPATIAL_STAGE_ORDER if stage in set(stage_loss["stage"].astype(str))]
    row_order = list(dict.fromkeys(stage_loss["row_label"].tolist()))
    matrix = (
        stage_loss.pivot_table(index="row_label", columns="stage", values="coupling_loss", aggfunc="mean", fill_value=0.0, observed=False)
        .reindex(index=row_order, columns=stages, fill_value=0.0)
    )
    values = matrix.to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("loss", ["#F7F7F7", "#F0C0CC", "#B64342"])
    image = ax.imshow(values, aspect="auto", vmin=0, vmax=max(1.0, float(np.nanmax(values))), cmap=cmap)
    ax.set_xticks(np.arange(len(stages)))
    ax.set_xticklabels(stages, fontsize=6.2)
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=5.9)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            color = "white" if value > 0.62 else "#272727"
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5.1, color=color)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("relative coupling loss", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    add_panel_label(ax, "c", x=-0.08)


def plot_method_concordance(ax: plt.Axes, concordance: pd.DataFrame) -> None:
    ax.set_title("Continuous vs top-quantile score dropout", loc="left", fontsize=7.5, fontweight="bold")
    label_offsets = {
        "MIF": (-0.13, 0.0, "right"),
        "CD74": (0.025, 0.0, "left"),
        "CD44": (0.035, 0.04, "left"),
        "SPP1": (0.035, 0.08, "left"),
        "TREM2": (0.08, 0.03, "left"),
        "CXCR4": (0.08, -0.02, "left"),
        "PLA2G7": (0.035, 0.0, "left"),
    }
    for _, row in concordance.iterrows():
        gene = str(row["perturbed_genes"])
        ax.scatter(
            row["continuous_coupling_loss"],
            row["dropout_priority_score"],
            s=38 if gene in {"MIF", "CD74"} else 26,
            color=GENE_COLORS.get(gene, "#767676"),
            edgecolor="black",
            linewidth=0.3,
            zorder=3,
        )
        dx, dy, ha = label_offsets.get(gene, (0.025, 0.0, "left"))
        ax.text(
            row["continuous_coupling_loss"] + dx,
            row["dropout_priority_score"] + dy,
            gene,
            fontsize=5.8,
            ha=ha,
            va="center",
            color=GENE_COLORS.get(gene, "#767676"),
        )
    ax.set_xlabel("continuous coupling loss")
    ax.set_ylabel("top-quantile dropout priority")
    ax.set_xlim(-0.02, max(1.05, float(concordance["continuous_coupling_loss"].max()) * 1.08))
    ax.set_ylim(-0.02, max(1.65, float(concordance["dropout_priority_score"].max()) * 1.08))
    ax.grid(axis="both", color="#E6E6E6", linewidth=0.5)
    ax.text(
        0.98,
        0.04,
        "MIF/CD74 remain prioritized across methods",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.8,
        color="#4D4D4D",
    )
    add_panel_label(ax, "d", x=-0.12)


def plot_figure5(
    continuous_ranking: pd.DataFrame,
    continuous_effects: pd.DataFrame,
    virtual_ranking: pd.DataFrame,
    args: argparse.Namespace,
) -> list[Path]:
    genes = ("MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7")
    dose = prepare_perturbation_dose_response(continuous_ranking, genes=genes)
    stage_loss = prepare_perturbation_stage_loss(continuous_effects, continuous_ranking, genes=genes)
    concordance = prepare_perturbation_method_concordance(continuous_ranking, virtual_ranking, genes=genes)

    dose.to_csv(args.table_dir / "figure5_dose_response_source.csv", index=False, encoding="utf-8-sig")
    stage_loss.to_csv(args.table_dir / "figure5_stage_loss_source.csv", index=False, encoding="utf-8-sig")
    concordance.to_csv(args.table_dir / "figure5_method_concordance_source.csv", index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[0.92, 1.08], width_ratios=[1.0, 1.0])
    ax_model = fig.add_subplot(gs[0, 0])
    ax_dose = fig.add_subplot(gs[0, 1])
    ax_stage = fig.add_subplot(gs[1, 0])
    ax_concordance = fig.add_subplot(gs[1, 1])

    plot_perturbation_schematic(ax_model)
    plot_dose_response(ax_dose, dose)
    plot_stage_loss_heatmap(ax_stage, stage_loss)
    plot_method_concordance(ax_concordance, concordance)
    fig.suptitle("Score-level in-silico target prioritization ranks receptor-side CD74", x=0.02, ha="left", fontsize=9, fontweight="bold")
    return save_figure(fig, args.figure_dir / "figure5_virtual_perturbation_priority")


def main() -> int:
    args = parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.table_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata)
    evidence_matrix = pd.read_csv(args.evidence_matrix)
    perturbation_ranking = pd.read_csv(args.perturbation_ranking)
    continuous_effects = pd.read_csv(args.continuous_perturbation_effects)
    virtual_ranking = pd.read_csv(args.virtual_perturbation_ranking)
    signature_celltype_summary = pd.read_csv(args.signature_celltype_summary)
    specificity_summary = pd.read_csv(args.specificity_summary)
    gene_top_celltype = pd.read_csv(args.gene_top_celltype)
    axis_spatial_by_stage = pd.read_csv(args.axis_spatial_by_stage)

    saved = []
    saved.extend(plot_figure1(metadata, args))
    saved.extend(plot_figure2(evidence_matrix, perturbation_ranking, args))
    saved.extend(plot_figure3(signature_celltype_summary, specificity_summary, gene_top_celltype, args))
    saved.extend(plot_figure4(axis_spatial_by_stage, args))
    saved.extend(plot_figure5(perturbation_ranking, continuous_effects, virtual_ranking, args))
    print("Wrote figure exports:")
    for path in saved:
        print(f"- {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
