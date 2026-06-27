#!/usr/bin/env python
"""Score GSE308103 snRNA-seq samples for LUAD precursor niche signatures."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.cell_states import add_within_group_high_flags, assign_dominant_labels  # noqa: E402
from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.sn_matrix import parse_gse308103_filename, read_tabular_selected_genes_with_totals  # noqa: E402


STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]
REFINED_PANEL_NAMES = {
    "epithelial_progenitor_like_vs_other_epithelial": "refined_epithelial_progenitor_like",
    "proliferating_epithelial_vs_other_epithelial": "refined_proliferating_epithelial",
    "spp1_macrophage_vs_other_macrophage": "refined_spp1_macrophage",
    "c1q_macrophage_vs_other_macrophage": "refined_c1q_macrophage",
    "inflammatory_macrophage_vs_other_macrophage": "refined_inflammatory_macrophage",
}
KEY_SCORE_COLUMNS = [
    "epithelial_score",
    "macrophage_score",
    "epithelial_progenitor_like_score",
    "spp1_macrophage_score",
    "c1q_macrophage_score",
    "inflammatory_macrophage_score",
    "refined_epithelial_progenitor_like_score",
    "refined_spp1_macrophage_score",
    "refined_c1q_macrophage_score",
    "refined_inflammatory_macrophage_score",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE308103" / "raw_counts",
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
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    parser.add_argument("--scale-factor", type=float, default=10_000.0)
    parser.add_argument("--class-min-score", type=float, default=0.05)
    parser.add_argument("--class-min-margin", type=float, default=0.05)
    parser.add_argument("--state-quantile", type=float, default=0.75)
    parser.add_argument("--limit-samples", type=int, default=None)
    return parser.parse_args()


def load_marker_config(path: Path) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    return (
        {name: list(genes) for name, genes in config["broad_classes"].items()},
        {name: list(genes) for name, genes in config["state_panels"].items()},
    )


def load_refined_panels(path: Path) -> dict[str, list[str]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        output_name: list(data[contrast])
        for contrast, output_name in REFINED_PANEL_NAMES.items()
        if contrast in data
    }


def flatten_genes(*panel_groups: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panels in panel_groups:
        for panel_genes in panels.values():
            genes.extend(gene for gene in panel_genes if gene not in genes)
    return genes


def sample_metadata(metadata: pd.DataFrame, sample_accession: str) -> dict[str, object]:
    rows = metadata[metadata["sample_accession"] == sample_accession]
    if rows.empty:
        return {"stage": "", "title": "", "source_name": ""}
    row = rows.iloc[0]
    return {
        "stage": row.get("interpreted_stage", ""),
        "title": row.get("title", ""),
        "source_name": row.get("source_name", ""),
    }


def assign_epithelial_states(table: pd.DataFrame, quantile: float) -> pd.DataFrame:
    result = add_within_group_high_flags(
        table,
        group_column="broad_class",
        score_column="epithelial_progenitor_like_score",
        group_value="epithelial",
        quantile=quantile,
        output_column="epithelial_progenitor_like_high",
    )
    result = add_within_group_high_flags(
        result,
        group_column="broad_class",
        score_column="at2_mature_score",
        group_value="epithelial",
        quantile=quantile,
        output_column="at2_mature_high",
    )
    result = add_within_group_high_flags(
        result,
        group_column="broad_class",
        score_column="proliferating_epithelial_score",
        group_value="epithelial",
        quantile=quantile,
        output_column="proliferating_epithelial_high",
    )
    result["epithelial_state"] = "not_epithelial"
    epithelial = result["broad_class"] == "epithelial"
    result.loc[epithelial, "epithelial_state"] = "other_epithelial"
    result.loc[epithelial & result["at2_mature_high"], "epithelial_state"] = "at2_mature_like"
    progenitor = epithelial & result["epithelial_progenitor_like_high"] & ~result["at2_mature_high"]
    result.loc[progenitor, "epithelial_state"] = "progenitor_like"
    proliferating = epithelial & result["proliferating_epithelial_high"]
    result.loc[proliferating, "epithelial_state"] = "proliferating_epithelial"
    return result


def assign_macrophage_states(table: pd.DataFrame, min_score: float, min_margin: float) -> pd.DataFrame:
    macrophage_columns = {
        "spp1_macrophage": "spp1_macrophage_score",
        "c1q_macrophage": "c1q_macrophage_score",
        "inflammatory_macrophage": "inflammatory_macrophage_score",
        "resident_macrophage": "resident_macrophage_score",
    }
    result = table.copy()
    subtype = assign_dominant_labels(
        result,
        macrophage_columns,
        min_score=min_score,
        min_margin=min_margin,
    )
    result["macrophage_state"] = "not_macrophage"
    macrophage = result["broad_class"] == "macrophage"
    result.loc[macrophage, "macrophage_state"] = subtype.loc[macrophage]
    result.loc[result["macrophage_state"].eq("unassigned"), "macrophage_state"] = "other_macrophage"
    result.loc[result["macrophage_state"].eq("ambiguous"), "macrophage_state"] = "ambiguous_macrophage"
    return result


def summarize_state_fractions(
    table: pd.DataFrame,
    state_column: str,
    denominator_filter: pd.Series | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if denominator_filter is None:
        denominator_filter = pd.Series(True, index=table.index)
    working = table.loc[denominator_filter].copy()
    sample_totals = working.groupby(["sample", "sample_name", "stage"]).size().reset_index(name="denominator")
    sample_counts = (
        working.groupby(["sample", "sample_name", "stage", state_column])
        .size()
        .reset_index(name="n_cells")
    )
    sample_summary = sample_counts.merge(sample_totals, on=["sample", "sample_name", "stage"], how="left")
    sample_summary["fraction"] = sample_summary["n_cells"] / sample_summary["denominator"]
    stage_summary = (
        sample_summary.groupby(["stage", state_column])
        .agg(
            n_samples=("sample", "nunique"),
            mean_fraction=("fraction", "mean"),
            sd_fraction=("fraction", "std"),
            total_cells=("n_cells", "sum"),
        )
        .reset_index()
    )
    return sample_summary, stage_summary


def summarize_scores(table: pd.DataFrame, score_columns: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    contexts = {
        "all_cells": pd.Series(True, index=table.index),
        "epithelial": table["broad_class"].eq("epithelial"),
        "macrophage": table["broad_class"].eq("macrophage"),
    }
    for context, mask in contexts.items():
        subset = table.loc[mask]
        for score_column in score_columns:
            records.append(
                {
                    "sample": table["sample"].iloc[0],
                    "sample_name": table["sample_name"].iloc[0],
                    "stage": table["stage"].iloc[0],
                    "context": context,
                    "score": score_column.removesuffix("_score"),
                    "n_cells": len(subset),
                    "mean_score": subset[score_column].mean() if len(subset) else None,
                }
            )
    sample_summary = pd.DataFrame(records)
    stage_summary = (
        sample_summary.groupby(["stage", "context", "score"])
        .agg(
            n_samples=("sample", "nunique"),
            mean_score=("mean_score", "mean"),
            sd_score=("mean_score", "std"),
            total_cells=("n_cells", "sum"),
        )
        .reset_index()
    )
    return sample_summary, stage_summary


def plot_stage_fractions(stage_summary: pd.DataFrame, state_column: str, output: Path, title: str) -> None:
    pivot = stage_summary.pivot_table(index="stage", columns=state_column, values="mean_fraction", fill_value=0)
    order = [stage for stage in STAGE_ORDER if stage in pivot.index]
    pivot = pivot.loc[order]
    ax = pivot.plot(kind="bar", stacked=True, figsize=(9, 4.8), width=0.8)
    ax.set_ylabel("Mean sample fraction")
    ax.set_xlabel("")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    ax.figure.savefig(output, dpi=220)
    plt.close(ax.figure)


def plot_refined_scores(stage_summary: pd.DataFrame, output: Path) -> None:
    subset = stage_summary[
        stage_summary["context"].isin(["epithelial", "macrophage"])
        & stage_summary["score"].isin(
            [
                "refined_epithelial_progenitor_like",
                "refined_spp1_macrophage",
                "refined_c1q_macrophage",
                "refined_inflammatory_macrophage",
            ]
        )
    ].copy()
    subset["label"] = subset["context"] + ":" + subset["score"]
    pivot = subset.pivot_table(index="stage", columns="label", values="mean_score")
    order = [stage for stage in STAGE_ORDER if stage in pivot.index]
    pivot = pivot.loc[order]
    ax = pivot.plot(kind="line", marker="o", figsize=(8, 4.8))
    ax.set_ylabel("Mean normalized score")
    ax.set_xlabel("")
    ax.set_title("GSE308103 refined discovery-signature scores")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    ax.figure.savefig(output, dpi=220)
    plt.close(ax.figure)


def main() -> int:
    args = parse_args()
    broad_classes, state_panels = load_marker_config(args.markers)
    refined_panels = load_refined_panels(args.refined_signatures)
    panels = {**broad_classes, **state_panels, **refined_panels}
    all_genes = flatten_genes(broad_classes, state_panels, refined_panels)
    metadata = pd.read_csv(args.metadata)
    files = sorted(args.input_dir.glob("GSM*.raw_counts.mtx.txt.gz"))
    if args.limit_samples:
        files = files[: args.limit_samples]
    if not files:
        raise SystemExit(f"No GSE308103 raw-count files found under {args.input_dir}")

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    compact_tables = []
    score_sample_tables = []
    genes_used_by_sample = {}

    for path in files:
        parsed = parse_gse308103_filename(path)
        sample_accession = parsed["sample_accession"]
        sample_name = parsed["sample_name"]
        meta = sample_metadata(metadata, sample_accession)
        counts, total_counts = read_tabular_selected_genes_with_totals(path, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        scores = compute_panel_scores(normalized, panels)
        genes_used_by_sample[sample_accession] = scores.attrs["panel_genes_used"]

        table = pd.DataFrame(
            {
                "sample": sample_accession,
                "sample_name": sample_name,
                "stage": meta["stage"],
                "cell_barcode": counts.index,
                "total_counts": total_counts.reindex(counts.index).to_numpy(),
            }
        )
        table = pd.concat([table, scores.reset_index(drop=True)], axis=1)
        class_columns = {name: f"{name}_score" for name in broad_classes}
        table["broad_class"] = assign_dominant_labels(
            table,
            class_columns,
            min_score=args.class_min_score,
            min_margin=args.class_min_margin,
        )
        table = assign_epithelial_states(table, quantile=args.state_quantile)
        table = assign_macrophage_states(
            table,
            min_score=args.class_min_score,
            min_margin=args.class_min_margin / 2,
        )

        key_score_columns = KEY_SCORE_COLUMNS
        present_key_scores = [column for column in key_score_columns if column in table.columns]
        compact_tables.append(
            table[
                [
                    "sample",
                    "sample_name",
                    "stage",
                    "cell_barcode",
                    "total_counts",
                    "broad_class",
                    "epithelial_state",
                    "macrophage_state",
                    *present_key_scores,
                ]
            ].copy()
        )
        score_sample, score_stage = summarize_scores(table, present_key_scores)
        score_sample_tables.append(score_sample)
        counts_by_class = table["broad_class"].value_counts().to_dict()
        print(f"{sample_accession}: stage={meta['stage']}; cells={len(table)}; classes={counts_by_class}")

    cells = pd.concat(compact_tables, ignore_index=True)
    cell_output = args.table_dir / "gse308103_snrna_cell_state_assignments.csv"
    cells.to_csv(cell_output, index=False, encoding="utf-8-sig")

    broad_sample, broad_stage = summarize_state_fractions(cells, "broad_class", denominator_filter=None)
    epithelial_sample, epithelial_stage = summarize_state_fractions(
        cells,
        "epithelial_state",
        denominator_filter=cells["broad_class"] == "epithelial",
    )
    macrophage_sample, macrophage_stage = summarize_state_fractions(
        cells,
        "macrophage_state",
        denominator_filter=cells["broad_class"] == "macrophage",
    )
    score_sample = pd.concat(score_sample_tables, ignore_index=True)
    score_stage = (
        score_sample.groupby(["stage", "context", "score"])
        .agg(
            n_samples=("sample", "nunique"),
            mean_score=("mean_score", "mean"),
            sd_score=("mean_score", "std"),
            total_cells=("n_cells", "sum"),
        )
        .reset_index()
    )

    outputs = {
        "gse308103_snrna_broad_class_sample_summary.csv": broad_sample,
        "gse308103_snrna_broad_class_stage_summary.csv": broad_stage,
        "gse308103_snrna_epithelial_state_sample_summary.csv": epithelial_sample,
        "gse308103_snrna_epithelial_state_stage_summary.csv": epithelial_stage,
        "gse308103_snrna_macrophage_state_sample_summary.csv": macrophage_sample,
        "gse308103_snrna_macrophage_state_stage_summary.csv": macrophage_stage,
        "gse308103_snrna_refined_signature_sample_summary.csv": score_sample,
        "gse308103_snrna_refined_signature_stage_summary.csv": score_stage,
    }
    for filename, frame in outputs.items():
        frame.to_csv(args.table_dir / filename, index=False, encoding="utf-8-sig")

    genes_used_output = args.table_dir / "gse308103_snrna_genes_used.json"
    genes_used_output.write_text(
        json.dumps(genes_used_by_sample, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    plot_stage_fractions(
        broad_stage,
        "broad_class",
        args.figure_dir / "gse308103_snrna_broad_class_stage_fractions.png",
        "GSE308103 broad marker-score classes",
    )
    plot_stage_fractions(
        epithelial_stage,
        "epithelial_state",
        args.figure_dir / "gse308103_snrna_epithelial_state_stage_fractions.png",
        "GSE308103 epithelial marker-score states",
    )
    plot_stage_fractions(
        macrophage_stage,
        "macrophage_state",
        args.figure_dir / "gse308103_snrna_macrophage_state_stage_fractions.png",
        "GSE308103 macrophage marker-score states",
    )
    plot_refined_scores(score_stage, args.figure_dir / "gse308103_snrna_refined_signature_stage_scores.png")
    print(f"Wrote compact cell assignments: {cell_output}")
    print(f"Wrote genes used: {genes_used_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
