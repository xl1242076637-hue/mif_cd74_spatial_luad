#!/usr/bin/env python
"""Score candidate gene panels in a STOMICS/AnnData h5ad spatial sample."""

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

from luad_niche.h5ad import compute_panel_scores, read_h5ad_obs, read_h5ad_selected_genes  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=PROJECT_ROOT
        / "data"
        / "raw"
        / "STDS0000125"
        / "GSM5702474"
        / "GSM5702474_10x_Visium_processed.h5ad",
        help="Input h5ad file.",
    )
    parser.add_argument(
        "--genes",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_genes.yaml",
        help="Candidate marker gene YAML.",
    )
    parser.add_argument(
        "--sample",
        default="GSM5702474",
        help="Sample identifier used in output filenames.",
    )
    parser.add_argument(
        "--table-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gsm5702474_spatial_panel_scores.csv",
        help="Output per-spot score table.",
    )
    parser.add_argument(
        "--genes-used-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gsm5702474_panel_genes_used.json",
        help="Output JSON file listing genes used for each panel.",
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "figures" / "gsm5702474_spatial_panel_scores.png",
        help="Output spatial score figure.",
    )
    return parser.parse_args()


def load_panels(path: Path) -> dict[str, list[str]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    return {name: list(genes) for name, genes in config["panels"].items()}


def main() -> int:
    args = parse_args()
    panels = load_panels(args.genes)
    all_genes = []
    for genes in panels.values():
        all_genes.extend(gene for gene in genes if gene not in all_genes)

    obs = read_h5ad_obs(args.h5ad)
    expr = read_h5ad_selected_genes(args.h5ad, all_genes)
    scores = compute_panel_scores(expr, panels)
    table = obs.merge(scores.reset_index(names="spot"), on="spot", how="left")
    table.insert(0, "sample", args.sample)

    args.table_output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.table_output, index=False, encoding="utf-8-sig")
    args.genes_used_output.parent.mkdir(parents=True, exist_ok=True)
    args.genes_used_output.write_text(
        json.dumps(scores.attrs["panel_genes_used"], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    score_columns = [column for column in table.columns if column.endswith("_score")]
    ncols = 2
    nrows = (len(score_columns) + 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 4.5 * nrows), squeeze=False)
    for ax, score_column in zip(axes.ravel(), score_columns, strict=False):
        sc = ax.scatter(
            table["x"],
            -table["y"],
            c=table[score_column],
            s=8,
            cmap="viridis",
            linewidths=0,
        )
        ax.set_title(score_column)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    for ax in axes.ravel()[len(score_columns) :]:
        ax.axis("off")
    fig.suptitle(f"{args.sample} spatial candidate panel scores", y=0.995)
    fig.tight_layout()
    args.figure_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure_output, dpi=220)
    plt.close(fig)

    print(f"Wrote table: {args.table_output}")
    print(f"Wrote genes-used JSON: {args.genes_used_output}")
    print(f"Wrote figure: {args.figure_output}")
    print(f"Spots: {len(table)}; panels: {len(score_columns)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

