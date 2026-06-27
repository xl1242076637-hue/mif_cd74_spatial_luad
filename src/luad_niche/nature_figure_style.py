"""Shared visual contract for Nature-style manuscript figure redraws."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt


EXPORT_SUFFIXES = (".svg", ".pdf", ".tiff", ".png")

NATURE_PALETTE = {
    "mif_axis": "#B94A45",
    "mif_soft": "#F3D3CF",
    "mif_pale": "#FAECEA",
    "macrophage_axis": "#2F8C8C",
    "macrophage_soft": "#DCEEEF",
    "macrophage_pale": "#EEF7F7",
    "immune_blue": "#5B7FCA",
    "immune_soft": "#E6ECF8",
    "green_support": "#6F9F6A",
    "green_soft": "#E8F0E3",
    "gold_support": "#D9A84E",
    "gold_soft": "#F7E9C8",
    "violet_support": "#8B75B8",
    "violet_soft": "#ECE7F5",
    "neutral_light": "#D8D8D8",
    "neutral_mid": "#8F8F8F",
    "neutral_dark": "#4A4A4A",
    "neutral_black": "#242424",
    "neutral_pale": "#F7F7F7",
    "panel_band": "#F4F4F4",
}

AXIS_COLOR_MAP = {
    "mif_cd74_cxcr4": NATURE_PALETTE["mif_axis"],
    "spp1_trem2_macrophage_epithelial": NATURE_PALETTE["macrophage_axis"],
    "cxcl9_cxcl10_cxcr3": NATURE_PALETTE["immune_blue"],
    "c1q_apoe_trem2_lgals3": NATURE_PALETTE["green_support"],
    "inflammatory_il1_tnf_cxcl8": NATURE_PALETTE["neutral_mid"],
}


def apply_nature_style(font_size: float = 7.0) -> None:
    """Apply compact editable-text matplotlib defaults for journal figures."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.linewidth"] = 0.75
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["legend.frameon"] = False


def semantic_axis_color(axis_id: str) -> str:
    """Return the cross-figure semantic color for a candidate axis."""
    return AXIS_COLOR_MAP.get(axis_id, NATURE_PALETTE["neutral_mid"])


def compact_panel_title(title: str, max_words: int = 3) -> str:
    """Trim analytic titles to the short labels used inside dense panels."""
    words = title.split()
    if len(words) <= max_words:
        return title
    return " ".join(words[:max_words])


def add_panel_label(ax, label: str, x: float = -0.08, y: float = 1.02) -> None:
    """Place a small Nature-style panel label."""
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color=NATURE_PALETTE["neutral_black"],
    )


def save_publication_figure(fig, base_path: Path, suffixes: Iterable[str] = EXPORT_SUFFIXES) -> list[Path]:
    """Save editable vector and high-resolution raster exports."""
    base_path.parent.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    for suffix in suffixes:
        output = base_path.with_suffix(suffix)
        kwargs = {"dpi": 600} if suffix == ".tiff" else {"dpi": 300} if suffix == ".png" else {}
        fig.savefig(output, bbox_inches="tight", **kwargs)
        saved.append(output)
    plt.close(fig)
    return saved
