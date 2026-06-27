#!/usr/bin/env python
"""Run parameter/subsampling robustness checks for GSE308103 GRN perturbation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_DIR, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from luad_niche.grn_perturbation import (  # noqa: E402
    build_correlation_network,
    summarize_signature_impacts,
    summarize_target_ranking_stability,
    virtual_outgoing_knockout,
)
from run_gse308103_grn_virtual_perturbation import (  # noqa: E402
    SIGNATURES_BY_CLASS,
    TARGETS_BY_CLASS,
    class_signatures,
    load_axes_genes,
    load_class_expression,
    load_signature_genes,
    parse_csv_values,
    raw_count_files_by_accession,
)


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
        default=PROJECT_ROOT
        / "results"
        / "tables"
        / "gse189357_refined_signature_genes_gse131907_specificity_filtered.json",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--stages", default="MIA,LUAD")
    parser.add_argument("--classes", default="epithelial,macrophage")
    parser.add_argument("--max-cells-per-class", type=int, default=6000)
    parser.add_argument("--thresholds", default="0.03,0.05,0.08")
    parser.add_argument("--random-seeds", default="7,11,23")
    parser.add_argument("--n-steps", type=int, default=3)
    parser.add_argument("--restart", type=float, default=0.15)
    parser.add_argument("--decay", type=float, default=0.5)
    return parser.parse_args()


def parse_float_values(text: str) -> list[float]:
    return [float(value.strip()) for value in text.split(",") if value.strip()]


def parse_int_values(text: str) -> list[int]:
    return [int(value.strip()) for value in text.split(",") if value.strip()]


def rank_targets(gene_effects: pd.DataFrame, signature_impacts: pd.DataFrame) -> pd.DataFrame:
    """Collapse gene and signature impact rows to a target-level ranking."""
    if signature_impacts.empty:
        return pd.DataFrame()
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
        ["top_signature_mean_impact", "top_signature_max_impact", "target_gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    target_ranking.insert(0, "display_rank", range(1, len(target_ranking) + 1))
    return target_ranking


def load_expression_by_class(
    assignments: pd.DataFrame,
    raw_files: dict[str, Path],
    gene_universe: list[str],
    *,
    classes: list[str],
    random_seed: int,
    max_cells_per_class: int,
) -> dict[str, pd.DataFrame]:
    """Load selected-gene expression matrices once per class and seed."""
    expression_by_class: dict[str, pd.DataFrame] = {}
    for broad_class in classes:
        subset = assignments[assignments["broad_class"].eq(broad_class)].copy()
        if subset.empty:
            continue
        if max_cells_per_class and len(subset) > max_cells_per_class:
            subset = subset.sample(max_cells_per_class, random_state=random_seed)
        expression_by_class[broad_class] = load_class_expression(subset, raw_files, gene_universe)
    return expression_by_class


def run_threshold_from_expressions(
    expression_by_class: dict[str, pd.DataFrame],
    signatures: dict[str, list[str]],
    *,
    stages: set[str],
    classes: list[str],
    threshold: float,
    random_seed: int,
    n_steps: int,
    restart: float,
    decay: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_effect_rows: list[pd.DataFrame] = []
    signature_rows: list[pd.DataFrame] = []
    network_rows: list[dict[str, object]] = []
    for broad_class in classes:
        expr = expression_by_class.get(broad_class, pd.DataFrame())
        if expr.empty:
            continue
        network = build_correlation_network(expr, min_abs_correlation=threshold)
        targets = [gene for gene in TARGETS_BY_CLASS.get(broad_class, []) if gene in network.index]
        network_rows.append(
            {
                "broad_class": broad_class,
                "n_cells": int(expr.shape[0]),
                "n_input_genes": int(expr.shape[1]),
                "n_network_genes": int(network.shape[0]),
                "n_edges": int((network > 0).sum().sum()),
                "tested_targets": ",".join(targets),
            }
        )
        for target in targets:
            effects = virtual_outgoing_knockout(
                network,
                [target],
                n_steps=n_steps,
                restart=restart,
                decay=decay,
            )
            effects.insert(0, "target_gene", target)
            effects.insert(0, "broad_class", broad_class)
            gene_effect_rows.append(effects)

            sig_summary = summarize_signature_impacts(effects, class_signatures(signatures, broad_class))
            sig_summary.insert(0, "target_gene", target)
            sig_summary.insert(0, "broad_class", broad_class)
            signature_rows.append(sig_summary)
    gene_effects = pd.concat(gene_effect_rows, axis=0, ignore_index=True) if gene_effect_rows else pd.DataFrame()
    signature_impacts = pd.concat(signature_rows, axis=0, ignore_index=True) if signature_rows else pd.DataFrame()
    target_ranking = rank_targets(gene_effects, signature_impacts)
    if not target_ranking.empty:
        target_ranking.insert(0, "min_abs_correlation", threshold)
        target_ranking.insert(0, "random_seed", random_seed)
        target_ranking.insert(0, "run_id", f"seed{random_seed}_corr{threshold:g}")
        target_ranking.insert(0, "stages", ",".join(sorted(stages)))
    network_summary = pd.DataFrame(network_rows)
    if not network_summary.empty:
        network_summary.insert(0, "min_abs_correlation", threshold)
        network_summary.insert(0, "random_seed", random_seed)
        network_summary.insert(0, "run_id", f"seed{random_seed}_corr{threshold:g}")
        network_summary.insert(0, "stages", ",".join(sorted(stages)))
    return target_ranking, network_summary


def main() -> int:
    args = parse_args()
    stages = set(parse_csv_values(args.stages))
    classes = parse_csv_values(args.classes)
    thresholds = parse_float_values(args.thresholds)
    random_seeds = parse_int_values(args.random_seeds)
    signatures = load_signature_genes(args.signatures)
    axis_genes = load_axes_genes(args.mechanisms)
    signature_genes = [gene for genes in signatures.values() for gene in genes]
    target_genes = [gene for genes in TARGETS_BY_CLASS.values() for gene in genes]
    gene_universe = list(dict.fromkeys([*axis_genes, *signature_genes, *target_genes]))

    usecols = ["sample", "stage", "cell_barcode", "total_counts", "broad_class"]
    assignments = pd.read_csv(args.assignments, usecols=usecols)
    assignments = assignments[assignments["stage"].isin(stages)].copy()
    raw_files = raw_count_files_by_accession(args.raw_dir)

    detail_rows: list[pd.DataFrame] = []
    network_rows: list[pd.DataFrame] = []
    for random_seed in random_seeds:
        expression_by_class = load_expression_by_class(
            assignments,
            raw_files,
            gene_universe,
            classes=classes,
            random_seed=random_seed,
            max_cells_per_class=args.max_cells_per_class,
        )
        for threshold in thresholds:
            detail, network = run_threshold_from_expressions(
                expression_by_class,
                signatures,
                stages=stages,
                classes=classes,
                threshold=threshold,
                random_seed=random_seed,
                n_steps=args.n_steps,
                restart=args.restart,
                decay=args.decay,
            )
            detail_rows.append(detail)
            network_rows.append(network)
            print(f"Finished seed={random_seed}, corr={threshold:g}: {len(detail)} target rows")

    detail = pd.concat(detail_rows, axis=0, ignore_index=True) if detail_rows else pd.DataFrame()
    network_summary = pd.concat(network_rows, axis=0, ignore_index=True) if network_rows else pd.DataFrame()
    summary = summarize_target_ranking_stability(detail) if not detail.empty else pd.DataFrame()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.table_dir / "gse308103_grn_virtual_perturbation_robustness_detail.csv"
    summary_path = args.table_dir / "gse308103_grn_virtual_perturbation_robustness_summary.csv"
    network_path = args.table_dir / "gse308103_grn_virtual_perturbation_robustness_network_summary.csv"
    config_path = args.table_dir / "gse308103_grn_virtual_perturbation_robustness_config.json"
    detail.to_csv(detail_path, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    network_summary.to_csv(network_path, index=False, encoding="utf-8-sig")
    config = {
        "method": "parameter/subsampling robustness for scTenifoldKnk-inspired outgoing-edge virtual perturbation",
        "stages": sorted(stages),
        "classes": classes,
        "thresholds": thresholds,
        "random_seeds": random_seeds,
        "max_cells_per_class": args.max_cells_per_class,
        "n_steps": args.n_steps,
        "restart": args.restart,
        "decay": args.decay,
        "signature_classes": SIGNATURES_BY_CLASS,
    }
    config_path.write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote robustness detail: {detail_path}")
    print(f"Wrote robustness summary: {summary_path}")
    print(f"Wrote robustness network summary: {network_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
