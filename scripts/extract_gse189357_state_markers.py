#!/usr/bin/env python
"""Extract aggregate marker signatures from GSE189357 marker-score cell states."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.io import mmread


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.differential import aggregate_contrast, rank_markers  # noqa: E402
from luad_niche.tenx import discover_10x_samples, read_10x_barcodes, read_10x_features  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE189357" / "raw_10x",
        help="Directory containing extracted GSE189357 10x files.",
    )
    parser.add_argument(
        "--assignments",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_scrna_cell_state_assignments.csv",
        help="Cell-state assignment table from refine_gse189357_cell_states.py.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output table directory.",
    )
    parser.add_argument("--top-n", type=int, default=30, help="Top genes to export per contrast.")
    return parser.parse_args()


def contrast_masks(assignments: pd.DataFrame) -> dict[str, tuple[pd.Series, pd.Series]]:
    epithelial_background = assignments["broad_class"].eq("epithelial")
    macrophage_background = assignments["broad_class"].eq("macrophage")
    return {
        "epithelial_progenitor_like_vs_other_epithelial": (
            assignments["epithelial_state"].eq("progenitor_like"),
            epithelial_background & ~assignments["epithelial_state"].eq("progenitor_like"),
        ),
        "proliferating_epithelial_vs_other_epithelial": (
            assignments["epithelial_state"].eq("proliferating_epithelial"),
            epithelial_background & ~assignments["epithelial_state"].eq("proliferating_epithelial"),
        ),
        "spp1_macrophage_vs_other_macrophage": (
            assignments["macrophage_state"].eq("spp1_macrophage"),
            macrophage_background & ~assignments["macrophage_state"].eq("spp1_macrophage"),
        ),
        "c1q_macrophage_vs_other_macrophage": (
            assignments["macrophage_state"].eq("c1q_macrophage"),
            macrophage_background & ~assignments["macrophage_state"].eq("c1q_macrophage"),
        ),
        "inflammatory_macrophage_vs_other_macrophage": (
            assignments["macrophage_state"].eq("inflammatory_macrophage"),
            macrophage_background & ~assignments["macrophage_state"].eq("inflammatory_macrophage"),
        ),
        "resident_macrophage_vs_other_macrophage": (
            assignments["macrophage_state"].eq("resident_macrophage"),
            macrophage_background & ~assignments["macrophage_state"].eq("resident_macrophage"),
        ),
    }


def main() -> int:
    args = parse_args()
    samples = discover_10x_samples(args.input_dir, require_tissue_positions=False)
    assignments = pd.read_csv(args.assignments)
    accumulators: dict[str, list[pd.DataFrame]] = {}
    totals: dict[str, dict[str, int]] = {}

    for sample in samples.values():
        sample_assignments = assignments[assignments["sample"].eq(sample.sample_accession)].copy()
        barcodes = read_10x_barcodes(sample)
        sample_assignments = sample_assignments.set_index("cell_barcode").loc[barcodes].reset_index()
        features = read_10x_features(sample)
        matrix = mmread(sample.matrix).tocsr()
        masks = contrast_masks(sample_assignments)

        for contrast, (group_mask, background_mask) in masks.items():
            n_group = int(group_mask.sum())
            n_background = int(background_mask.sum())
            if n_group == 0 or n_background == 0:
                continue
            aggregate = aggregate_contrast(matrix, group_mask.to_numpy(), background_mask.to_numpy())
            aggregate["gene"] = features["gene"].astype(str).to_numpy()
            aggregate = aggregate[
                ["gene", "group_sum", "background_sum", "group_detected", "background_detected"]
            ]
            accumulators.setdefault(contrast, []).append(aggregate)
            totals.setdefault(contrast, {"n_group": 0, "n_background": 0})
            totals[contrast]["n_group"] += n_group
            totals[contrast]["n_background"] += n_background
        print(f"Aggregated full matrix markers for {sample.sample_accession}")

    ranked_frames = []
    signatures: dict[str, list[str]] = {}
    for contrast, frames in accumulators.items():
        combined = (
            pd.concat(frames, ignore_index=True)
            .groupby("gene", as_index=False)
            .agg(
                group_sum=("group_sum", "sum"),
                background_sum=("background_sum", "sum"),
                group_detected=("group_detected", "sum"),
                background_detected=("background_detected", "sum"),
            )
        )
        combined["n_group"] = totals[contrast]["n_group"]
        combined["n_background"] = totals[contrast]["n_background"]
        ranked = rank_markers(combined)
        ranked.insert(0, "contrast", contrast)
        ranked_frames.append(ranked)
        filtered = ranked[(ranked["log2fc"] > 0) & (ranked["pct_group"] >= 0.05)]
        signatures[contrast] = filtered.head(args.top_n)["gene"].tolist()

    all_ranked = pd.concat(ranked_frames, ignore_index=True)
    top_ranked = (
        all_ranked[(all_ranked["log2fc"] > 0) & (all_ranked["pct_group"] >= 0.05)]
        .groupby("contrast", group_keys=False)
        .head(args.top_n)
    )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    ranked_output = args.table_dir / "gse189357_refined_state_markers.csv"
    top_output = args.table_dir / "gse189357_refined_state_top_markers.csv"
    json_output = args.table_dir / "gse189357_refined_state_signature_genes.json"
    all_ranked.to_csv(ranked_output, index=False, encoding="utf-8-sig")
    top_ranked.to_csv(top_output, index=False, encoding="utf-8-sig")
    json_output.write_text(json.dumps(signatures, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Wrote ranked markers: {ranked_output}")
    print(f"Wrote top markers: {top_output}")
    print(f"Wrote signature JSON: {json_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
