#!/usr/bin/env python
"""Summarize original scTenifoldKnk sensitivity output against the main GRN run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]

EPITHELIAL_MARKERS = {
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT19",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "CD24",
    "PERP",
    "MDK",
}
MACROPHAGE_MARKERS = {
    "APOE",
    "APOC1",
    "C1QA",
    "C1QB",
    "C1QC",
    "TREM2",
    "PLA2G7",
    "LPL",
    "LRP1",
    "LGALS3",
    "CHIT1",
    "CYP27A1",
    "FABP3",
    "RAB42",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--original-summary",
        type=Path,
        default=PROJECT_ROOT
        / "results"
        / "tables"
        / "sctenifoldknk_original_expanded"
        / "sctenifoldknk_original_expanded_summary.csv",
    )
    parser.add_argument(
        "--main-grn-ranking",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_grn_virtual_perturbation_target_ranking.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    return parser.parse_args()


def split_genes(value: object) -> list[str]:
    return [gene.strip() for gene in str(value).split(",") if gene.strip()]


def classify_boundary(row: pd.Series) -> str:
    genes = split_genes(row["original_top_10_non_target_genes"])
    p_adj = pd.to_numeric(row["original_top_non_target_p_adj"], errors="coerce")
    epithelial_hits = set(genes).intersection(EPITHELIAL_MARKERS)
    macrophage_hits = set(genes).intersection(MACROPHAGE_MARKERS)

    if row["broad_class"] == "epithelial":
        if epithelial_hits and pd.notna(p_adj) and p_adj < 0.05:
            return "supports epithelial network-context sensitivity"
        if epithelial_hits:
            return "directionally epithelial but not FDR-supported"
        return "technical sensitivity only"

    if pd.notna(p_adj) and p_adj < 0.05 and macrophage_hits and len(epithelial_hits) <= 1:
        return "supports macrophage network-context sensitivity"
    if epithelial_hits:
        return "reduced-panel boundary: off-class epithelial genes among top hits"
    if macrophage_hits:
        return "directionally macrophage but not FDR-supported"
    return "technical sensitivity only"


def main() -> int:
    args = parse_args()
    args.table_dir.mkdir(parents=True, exist_ok=True)

    original = pd.read_csv(args.original_summary)
    main_grn = pd.read_csv(args.main_grn_ranking)

    original = original.rename(
        columns={
            "top_non_target_gene_by_p_adj": "original_top_non_target_gene_by_p_adj",
            "top_non_target_p_adj": "original_top_non_target_p_adj",
            "top_10_non_target_genes": "original_top_10_non_target_genes",
        }
    )
    keep_original = [
        "broad_class",
        "target_gene",
        "n_input_genes",
        "n_input_cells",
        "original_top_non_target_gene_by_p_adj",
        "original_top_non_target_p_adj",
        "original_top_10_non_target_genes",
        "sctenifoldknk_version",
        "nc_nNet",
        "nc_nCells",
        "td_maxIter",
    ]
    original = original[keep_original].copy()

    keep_main = [
        "broad_class",
        "target_gene",
        "top_impacted_signature",
        "top_signature_mean_impact",
        "top_impacted_genes",
    ]
    merged = original.merge(main_grn[keep_main], on=["broad_class", "target_gene"], how="left")
    merged["original_gene_overlap_with_main_grn_top_genes"] = merged.apply(
        lambda row: ",".join(
            sorted(
                set(split_genes(row["original_top_10_non_target_genes"])).intersection(
                    split_genes(row["top_impacted_genes"])
                )
            )
        ),
        axis=1,
    )
    merged["epithelial_marker_hits_in_original_top10"] = merged["original_top_10_non_target_genes"].map(
        lambda value: ",".join(sorted(set(split_genes(value)).intersection(EPITHELIAL_MARKERS)))
    )
    merged["macrophage_marker_hits_in_original_top10"] = merged["original_top_10_non_target_genes"].map(
        lambda value: ",".join(sorted(set(split_genes(value)).intersection(MACROPHAGE_MARKERS)))
    )
    merged["interpretation_boundary"] = merged.apply(classify_boundary, axis=1)
    merged["manuscript_use"] = (
        "Report as original scTenifoldKnk sensitivity/boundary analysis; do not use as causal knockout evidence."
    )

    output = args.table_dir / "sctenifoldknk_original_expanded_interpretation.csv"
    merged.to_csv(output, index=False, encoding="utf-8-sig")
    print(f"Wrote {len(merged)} original scTenifoldKnk interpretation rows: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
