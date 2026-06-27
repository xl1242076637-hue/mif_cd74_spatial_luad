#!/usr/bin/env python
"""Run scTenifoldKnk-inspired GRN-level virtual perturbation on GSE308103."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.grn_perturbation import (  # noqa: E402
    build_correlation_network,
    summarize_signature_impacts,
    virtual_outgoing_knockout,
)
from luad_niche.sn_matrix import parse_gse308103_filename, read_tabular_selected_genes  # noqa: E402


TARGETS_BY_CLASS = {
    "epithelial": ["MIF"],
    "macrophage": ["CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7"],
}

SIGNATURES_BY_CLASS = {
    "epithelial": [
        "epithelial_progenitor_like_vs_other_epithelial",
        "proliferating_epithelial_vs_other_epithelial",
    ],
    "macrophage": [
        "spp1_macrophage_vs_other_macrophage",
        "c1q_macrophage_vs_other_macrophage",
        "inflammatory_macrophage_vs_other_macrophage",
        "resident_macrophage_vs_other_macrophage",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE308103" / "raw_counts",
    )
    parser.add_argument(
        "--assignments",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_snrna_cell_state_assignments.csv",
    )
    parser.add_argument(
        "--mechanisms",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_mechanisms.yaml",
    )
    parser.add_argument(
        "--signatures",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_signature_genes_gse131907_specificity_filtered.json",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--stages", default="MIA,LUAD")
    parser.add_argument("--classes", default="epithelial,macrophage")
    parser.add_argument("--max-cells-per-class", type=int, default=6000)
    parser.add_argument("--min-abs-correlation", type=float, default=0.05)
    parser.add_argument("--n-steps", type=int, default=3)
    parser.add_argument("--restart", type=float, default=0.15)
    parser.add_argument("--decay", type=float, default=0.5)
    parser.add_argument("--random-seed", type=int, default=7)
    return parser.parse_args()


def parse_csv_values(text: str) -> list[str]:
    return [value.strip() for value in text.split(",") if value.strip()]


def load_axes_genes(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    genes: list[str] = []
    for axis in config.get("axes", []):
        for field in ("source_genes", "target_genes", "bulk_genes", "perturbation_genes"):
            genes.extend(axis.get(field, []) or [])
    return list(dict.fromkeys(genes))


def load_signature_genes(path: Path) -> dict[str, list[str]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {key: list(dict.fromkeys(value)) for key, value in data.items()}


def raw_count_files_by_accession(raw_dir: Path) -> dict[str, Path]:
    files = {}
    for path in sorted(raw_dir.glob("*.raw_counts.mtx.txt.gz")):
        parsed = parse_gse308103_filename(path)
        files[parsed["sample_accession"]] = path
    return files


def load_class_expression(
    assignments: pd.DataFrame,
    raw_files: dict[str, Path],
    genes: list[str],
    *,
    scale_factor: float = 10_000.0,
) -> pd.DataFrame:
    matrices: list[pd.DataFrame] = []
    for sample, sample_cells in assignments.groupby("sample", sort=True):
        path = raw_files.get(str(sample))
        if path is None:
            continue
        barcodes = sample_cells["cell_barcode"].astype(str).tolist()
        counts = read_tabular_selected_genes(path, genes)
        present = [barcode for barcode in barcodes if barcode in counts.index]
        if not present:
            continue
        counts = counts.loc[present]
        totals = sample_cells.set_index("cell_barcode").loc[present, "total_counts"].astype(float)
        normalized = normalize_log1p_counts_by_totals(counts, totals, scale_factor=scale_factor)
        normalized.index = [f"{sample}|{barcode}" for barcode in normalized.index]
        matrices.append(normalized)
    if not matrices:
        return pd.DataFrame()
    return pd.concat(matrices, axis=0).fillna(0.0)


def class_signatures(signatures: dict[str, list[str]], broad_class: str) -> dict[str, list[str]]:
    return {
        name: signatures[name]
        for name in SIGNATURES_BY_CLASS.get(broad_class, [])
        if name in signatures
    }


def main() -> int:
    args = parse_args()
    stages = set(parse_csv_values(args.stages))
    classes = parse_csv_values(args.classes)
    signatures = load_signature_genes(args.signatures)
    axis_genes = load_axes_genes(args.mechanisms)
    signature_genes = [gene for genes in signatures.values() for gene in genes]
    target_genes = [gene for genes in TARGETS_BY_CLASS.values() for gene in genes]
    gene_universe = list(dict.fromkeys([*axis_genes, *signature_genes, *target_genes]))

    usecols = ["sample", "stage", "cell_barcode", "total_counts", "broad_class"]
    assignments = pd.read_csv(args.assignments, usecols=usecols)
    assignments = assignments[assignments["stage"].isin(stages)].copy()
    raw_files = raw_count_files_by_accession(args.raw_dir)

    gene_effect_rows: list[pd.DataFrame] = []
    signature_rows: list[pd.DataFrame] = []
    network_rows: list[dict[str, object]] = []
    for broad_class in classes:
        subset = assignments[assignments["broad_class"].eq(broad_class)].copy()
        if subset.empty:
            continue
        if args.max_cells_per_class and len(subset) > args.max_cells_per_class:
            subset = subset.sample(args.max_cells_per_class, random_state=args.random_seed)
        expr = load_class_expression(subset, raw_files, gene_universe)
        network = build_correlation_network(expr, min_abs_correlation=args.min_abs_correlation)
        targets = [gene for gene in TARGETS_BY_CLASS.get(broad_class, []) if gene in network.index]
        network_rows.append(
            {
                "broad_class": broad_class,
                "stages": ",".join(sorted(stages)),
                "n_cells": int(expr.shape[0]),
                "n_input_genes": int(expr.shape[1]),
                "n_network_genes": int(network.shape[0]),
                "n_edges": int((network > 0).sum().sum()),
                "tested_targets": ",".join(targets),
                "missing_targets": ",".join(
                    gene for gene in TARGETS_BY_CLASS.get(broad_class, []) if gene not in network.index
                ),
            }
        )
        for target in targets:
            effects = virtual_outgoing_knockout(
                network,
                [target],
                n_steps=args.n_steps,
                restart=args.restart,
                decay=args.decay,
            )
            effects.insert(0, "target_gene", target)
            effects.insert(0, "broad_class", broad_class)
            effects.insert(0, "stages", ",".join(sorted(stages)))
            gene_effect_rows.append(effects)

            sig_summary = summarize_signature_impacts(effects, class_signatures(signatures, broad_class))
            sig_summary.insert(0, "target_gene", target)
            sig_summary.insert(0, "broad_class", broad_class)
            sig_summary.insert(0, "stages", ",".join(sorted(stages)))
            signature_rows.append(sig_summary)

    args.table_dir.mkdir(parents=True, exist_ok=True)
    gene_effects = pd.concat(gene_effect_rows, axis=0, ignore_index=True) if gene_effect_rows else pd.DataFrame()
    signature_impacts = (
        pd.concat(signature_rows, axis=0, ignore_index=True) if signature_rows else pd.DataFrame()
    )
    if signature_impacts.empty:
        target_ranking = pd.DataFrame()
    else:
        target_ranking = (
            signature_impacts.sort_values(
                ["broad_class", "target_gene", "mean_impact_score", "max_impact_score"],
                ascending=[True, True, False, False],
            )
            .groupby(["broad_class", "target_gene"], as_index=False)
            .head(1)
            .rename(
                columns={
                    "signature": "top_impacted_signature",
                    "mean_impact_score": "top_signature_mean_impact",
                    "max_impact_score": "top_signature_max_impact",
                    "sum_impact_score": "top_signature_sum_impact",
                }
            )
        )
        top_genes = (
            gene_effects.sort_values(["broad_class", "target_gene", "impact_score"], ascending=[True, True, False])
            .groupby(["broad_class", "target_gene"], as_index=False)
            .head(5)
            .groupby(["broad_class", "target_gene"])["gene"]
            .apply(lambda values: ",".join(values))
            .reset_index(name="top_impacted_genes")
        )
        target_ranking = target_ranking.merge(top_genes, on=["broad_class", "target_gene"], how="left")
        target_ranking = target_ranking.sort_values(
            ["top_signature_mean_impact", "top_signature_max_impact"],
            ascending=[False, False],
        ).reset_index(drop=True)
    network_summary = pd.DataFrame(network_rows)
    gene_effects.to_csv(
        args.table_dir / "gse308103_grn_virtual_perturbation_gene_effects.csv",
        index=False,
        encoding="utf-8-sig",
    )
    signature_impacts.to_csv(
        args.table_dir / "gse308103_grn_virtual_perturbation_signature_impacts.csv",
        index=False,
        encoding="utf-8-sig",
    )
    target_ranking.to_csv(
        args.table_dir / "gse308103_grn_virtual_perturbation_target_ranking.csv",
        index=False,
        encoding="utf-8-sig",
    )
    network_summary.to_csv(
        args.table_dir / "gse308103_grn_virtual_perturbation_network_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    config = {
        "method": "scTenifoldKnk-inspired outgoing-edge virtual perturbation on positive correlation GRN",
        "stages": sorted(stages),
        "classes": classes,
        "max_cells_per_class": args.max_cells_per_class,
        "min_abs_correlation": args.min_abs_correlation,
        "n_steps": args.n_steps,
        "restart": args.restart,
        "decay": args.decay,
        "targets_by_class": TARGETS_BY_CLASS,
    }
    (args.table_dir / "gse308103_grn_virtual_perturbation_config.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(
        f"Wrote {len(gene_effects)} gene-effect rows, "
        f"{len(signature_impacts)} signature-impact rows and {len(target_ranking)} target-ranking rows."
    )
    print(network_summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
