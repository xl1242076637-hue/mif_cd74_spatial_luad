#!/usr/bin/env python
"""Score candidate gene panels in GSE189357 raw scRNA-seq 10x matrices."""

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

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE189357" / "raw_10x",
        help="Directory containing extracted GSE189357 10x files.",
    )
    parser.add_argument(
        "--genes",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_genes.yaml",
        help="Candidate marker gene YAML.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
        help="Annotated GEO sample metadata.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output table directory.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "figures",
        help="Output figure directory.",
    )
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=10_000.0,
        help="Per-cell total-count scaling factor before log1p.",
    )
    return parser.parse_args()


def load_panels(path: Path) -> dict[str, list[str]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    return {name: list(genes) for name, genes in config["panels"].items()}


def flatten_panel_genes(panels: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panel_genes in panels.values():
        genes.extend(gene for gene in panel_genes if gene not in genes)
    return genes


def sample_metadata(metadata: pd.DataFrame, sample_accession: str) -> dict[str, object]:
    rows = metadata[metadata["sample_accession"] == sample_accession]
    if rows.empty:
        return {
            "stage": "",
            "histological_type": "",
            "radiological_type": "",
            "gender": "",
            "title": "",
        }
    row = rows.iloc[0]
    return {
        "stage": row.get("interpreted_stage", ""),
        "histological_type": row.get("histological_type", ""),
        "radiological_type": row.get("radiological_type", ""),
        "gender": row.get("gender", ""),
        "title": row.get("title", ""),
    }


def plot_stage_panel_means(stage_summary: pd.DataFrame, output: Path) -> None:
    score_columns = [column for column in stage_summary.columns if column.endswith("_score_mean")]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    x = range(len(stage_summary))
    width = 0.18
    offsets = [(-1.5 + index) * width for index in range(len(score_columns))]
    for offset, column in zip(offsets, score_columns, strict=False):
        ax.bar(
            [value + offset for value in x],
            stage_summary[column],
            width=width,
            label=column.removesuffix("_score_mean"),
        )
    ax.set_xticks(list(x))
    ax.set_xticklabels(stage_summary["stage"])
    ax.set_ylabel("Mean log-normalized panel score")
    ax.set_title("GSE189357 scRNA panel scores by stage")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    panels = load_panels(args.genes)
    all_genes = flatten_panel_genes(panels)
    metadata = pd.read_csv(args.metadata)
    samples = discover_10x_samples(args.input_dir, require_tissue_positions=False)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    cell_tables = []
    sample_records = []
    genes_used_by_sample = {}
    for sample in samples.values():
        meta = sample_metadata(metadata, sample.sample_accession)
        obs = read_10x_obs(sample).rename(columns={"spot": "cell_barcode"})
        counts, total_counts = read_10x_selected_genes_with_totals(sample, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        scores = compute_panel_scores(normalized, panels)
        genes_used_by_sample[sample.sample_accession] = scores.attrs["panel_genes_used"]

        table = obs.merge(
            total_counts.rename_axis("cell_barcode").reset_index(),
            on="cell_barcode",
            how="left",
        )
        table = table.merge(scores.reset_index(names="cell_barcode"), on="cell_barcode", how="left")
        table.insert(0, "sample_name", sample.sample_name)
        table.insert(0, "sample", sample.sample_accession)
        table["cell_id"] = table["sample"] + "_" + table["cell_barcode"].astype(str)
        for key, value in meta.items():
            table[key] = value
        cell_tables.append(table)

        score_columns = [column for column in table.columns if column.endswith("_score")]
        record = {
            "sample": sample.sample_accession,
            "sample_name": sample.sample_name,
            "stage": meta["stage"],
            "histological_type": meta["histological_type"],
            "radiological_type": meta["radiological_type"],
            "gender": meta["gender"],
            "n_cells": len(table),
            "mean_total_counts": table["total_counts"].mean(),
            "median_total_counts": table["total_counts"].median(),
        }
        for column in score_columns:
            record[f"{column}_mean"] = table[column].mean()
            record[f"{column}_median"] = table[column].median()
        sample_records.append(record)
        print(
            f"{sample.sample_accession}: stage={meta['stage']}; cells={len(table)}; "
            f"median_counts={record['median_total_counts']:.0f}"
        )

    cells = pd.concat(cell_tables, ignore_index=True)
    cell_output = args.table_dir / "gse189357_scrna_panel_scores.csv"
    cells.to_csv(cell_output, index=False, encoding="utf-8-sig")

    sample_summary = pd.DataFrame(sample_records)
    sample_output = args.table_dir / "gse189357_scrna_panel_score_sample_summary.csv"
    sample_summary.to_csv(sample_output, index=False, encoding="utf-8-sig")

    stage_summary = (
        sample_summary.groupby("stage")
        .agg(
            n_samples=("sample", "nunique"),
            total_cells=("n_cells", "sum"),
            mean_total_counts=("mean_total_counts", "mean"),
            median_total_counts=("median_total_counts", "median"),
            epithelial_progenitor_score_mean=("epithelial_progenitor_score_mean", "mean"),
            proliferation_emt_score_mean=("proliferation_emt_score_mean", "mean"),
            macrophage_niche_score_mean=("macrophage_niche_score_mean", "mean"),
            ligand_receptor_axes_score_mean=("ligand_receptor_axes_score_mean", "mean"),
        )
        .reset_index()
    )
    order = {"AIS": 0, "MIA": 1, "IAC": 2}
    stage_summary["_stage_order"] = stage_summary["stage"].map(order).fillna(len(order))
    stage_summary = stage_summary.sort_values(["_stage_order", "stage"]).drop(columns="_stage_order")
    stage_output = args.table_dir / "gse189357_scrna_panel_score_stage_summary.csv"
    stage_summary.to_csv(stage_output, index=False, encoding="utf-8-sig")

    genes_used_output = args.table_dir / "gse189357_scrna_panel_genes_used.json"
    genes_used_output.write_text(
        json.dumps(genes_used_by_sample, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    plot_stage_panel_means(
        stage_summary,
        args.figure_dir / "gse189357_scrna_panel_score_stage_summary.png",
    )

    print(f"Wrote cell scores: {cell_output}")
    print(f"Wrote sample summary: {sample_output}")
    print(f"Wrote stage summary: {stage_output}")
    print(f"Wrote genes used: {genes_used_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
