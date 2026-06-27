#!/usr/bin/env python
"""Summarize candidate niche-marker expression across GSE282617 stages."""

from __future__ import annotations

import argparse
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

from luad_niche.expression import (  # noqa: E402
    STAGE_ORDER,
    expression_by_gene,
    marker_stage_table,
    marker_trend_summary,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--expression",
        type=Path,
        default=PROJECT_ROOT / "data" / "raw" / "GSE282617" / "GSE282617_processed_data.csv.gz",
        help="GSE282617 processed expression CSV.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
        help="Annotated GEO sample metadata CSV.",
    )
    parser.add_argument(
        "--genes",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_genes.yaml",
        help="Candidate marker gene YAML.",
    )
    parser.add_argument(
        "--table-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse282617_candidate_marker_stage_means.csv",
        help="Output CSV table.",
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "figures" / "gse282617_candidate_marker_stage_means.png",
        help="Output PNG figure.",
    )
    parser.add_argument(
        "--trend-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse282617_candidate_marker_trends.csv",
        help="Output trend summary CSV.",
    )
    return parser.parse_args()


def flatten_gene_panels(path: Path) -> list[str]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    genes: list[str] = []
    for panel_genes in config["panels"].values():
        for gene in panel_genes:
            if gene not in genes:
                genes.append(gene)
    return genes


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.expression)
    metadata = pd.read_csv(args.metadata)
    metadata = metadata[metadata["series_accession"] == "GSE282617"].copy()
    expr = expression_by_gene(df, value_prefix="FPKM.")
    markers = flatten_gene_panels(args.genes)
    table = marker_stage_table(expr, metadata, markers)

    args.table_output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.table_output, index=False, encoding="utf-8-sig")
    trends = marker_trend_summary(table)
    trends.to_csv(args.trend_output, index=False, encoding="utf-8-sig")

    plot_df = table[table["status"] == "present"].copy()
    pivot = plot_df.pivot(index="gene", columns="stage", values="mean_expression")
    present_stages = [stage for stage in STAGE_ORDER if stage in pivot.columns]
    pivot = pivot[present_stages]
    z = pivot.sub(pivot.mean(axis=1), axis=0).div(pivot.std(axis=1).replace(0, 1), axis=0)
    z = z.fillna(0)

    fig_height = max(6, 0.25 * len(z))
    fig, ax = plt.subplots(figsize=(7, fig_height))
    im = ax.imshow(z.values, aspect="auto", cmap="vlag" if "vlag" in plt.colormaps() else "coolwarm")
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels(z.columns)
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=8)
    ax.set_title("GSE282617 candidate marker stage means (row z-score)")
    fig.colorbar(im, ax=ax, label="row z-score")
    fig.tight_layout()
    args.figure_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure_output, dpi=200)
    plt.close(fig)

    missing = table[table["status"] == "missing"]["gene"].tolist()
    print(f"Wrote table: {args.table_output}")
    print(f"Wrote trends: {args.trend_output}")
    print(f"Wrote figure: {args.figure_output}")
    print(f"Present markers: {len(set(markers) - set(missing))}; missing markers: {len(missing)}")
    if missing:
        print("Missing:", ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
