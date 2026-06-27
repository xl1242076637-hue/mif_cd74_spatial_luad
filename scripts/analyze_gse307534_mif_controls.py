#!/usr/bin/env python
"""Run expression-matched random-gene and density controls for source-side MIF."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.io import mmread


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_DIR, SCRIPTS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.spatial_controls import (  # noqa: E402
    select_expression_matched_controls,
    spearman_association,
    summarize_random_control_distribution,
)
from luad_niche.spatial_statistics import (  # noqa: E402
    add_patient_phase,
    build_paired_patient_differences,
    summarize_paired_patient_differences,
)
from luad_niche.tenx import (  # noqa: E402
    TenXSampleFiles,
    discover_10x_samples,
    read_10x_features,
    read_10x_obs,
    read_10x_selected_genes_with_totals,
)
from score_gse307534_refined_signatures import (  # noqa: E402
    adjacency_for_target,
    flatten_genes,
    load_signatures,
    sample_accessions_from_filelist,
    sample_metadata,
)


STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE307534" / "raw_visium",
    )
    parser.add_argument(
        "--filelist",
        type=Path,
        default=PROJECT_ROOT / "data" / "raw" / "GSE307534" / "filelist.txt",
    )
    parser.add_argument(
        "--signatures",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_signature_genes_gse131907_specificity_filtered.json",
    )
    parser.add_argument(
        "--mechanisms",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_mechanisms.yaml",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
    )
    parser.add_argument(
        "--axis-adjacency",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_candidate_axis_spatial_adjacency.csv",
    )
    parser.add_argument(
        "--broad-adjacency",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_refined_signature_adjacency.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    parser.add_argument("--docs-dir", type=Path, default=PROJECT_ROOT / "docs")
    parser.add_argument("--target-gene", default="MIF")
    parser.add_argument("--n-controls", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=25)
    parser.add_argument("--bootstrap-iterations", type=int, default=2_000)
    parser.add_argument("--scale-factor", type=float, default=10_000.0)
    parser.add_argument("--quantile", type=float, default=0.75)
    parser.add_argument("--radius-multiplier", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=307534)
    return parser.parse_args()


def load_excluded_genes(mechanisms_path: Path, signatures: dict[str, list[str]]) -> set[str]:
    """Exclude candidate and signature genes from random-control selection."""
    config = yaml.safe_load(Path(mechanisms_path).read_text(encoding="utf-8")) or {}
    excluded = set(flatten_genes(signatures))
    for axis in config.get("axes", []):
        for field in ("source_genes", "target_genes", "bulk_genes", "perturbation_genes"):
            excluded.update(axis.get(field, []) or [])
    return excluded


def choose_reference_samples(
    samples: list[TenXSampleFiles],
    metadata: pd.DataFrame,
    accession_by_sample: dict[str, str],
) -> list[tuple[TenXSampleFiles, str, dict[str, object]]]:
    """Choose one deterministic reference section per ordered stage."""
    selected = []
    seen = set()
    for sample in sorted(samples, key=lambda item: item.sample_name):
        accession = accession_by_sample.get(sample.sample_name, sample.sample_accession)
        meta = sample_metadata(metadata, accession)
        stage = str(meta["stage"])
        if stage in STAGE_ORDER and stage not in seen:
            selected.append((sample, accession, meta))
            seen.add(stage)
    return selected


def reference_gene_summary(
    reference_samples: list[tuple[TenXSampleFiles, str, dict[str, object]]],
    *,
    scale_factor: float,
) -> pd.DataFrame:
    """Estimate gene-level mean normalized expression in stage-balanced references."""
    tables = []
    for sample, accession, meta in reference_samples:
        obs = read_10x_obs(sample)
        in_tissue = obs["in_tissue"].eq(1).fillna(False).to_numpy()
        features = read_10x_features(sample)
        matrix = mmread(sample.matrix).tocsr()[:, in_tissue]
        totals = np.asarray(matrix.sum(axis=0)).ravel()
        factors = scale_factor / np.where(totals > 0, totals, 1.0)
        normalized = matrix.multiply(factors)
        normalized.data = np.log1p(normalized.data)
        means = np.asarray(normalized.mean(axis=1)).ravel()
        table = pd.DataFrame({"gene": features["gene"].astype(str), "mean_expression": means})
        table = table.groupby("gene", as_index=False).agg(mean_expression=("mean_expression", "sum"))
        table["sample"] = accession
        table["stage"] = meta["stage"]
        tables.append(table)
        print(f"Reference expression: {accession} {meta['stage']}")
    combined = pd.concat(tables, ignore_index=True)
    return (
        combined.groupby("gene", as_index=False)
        .agg(
            mean_expression=("mean_expression", "mean"),
            n_reference_samples=("sample", "nunique"),
        )
        .sort_values("gene")
        .reset_index(drop=True)
    )


def run_random_control_adjacency(
    samples: list[TenXSampleFiles],
    metadata: pd.DataFrame,
    accession_by_sample: dict[str, str],
    signatures: dict[str, list[str]],
    control_genes: list[str],
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Calculate source-epithelial adjacency for MIF and matched single-gene controls."""
    panels = {"epithelial_progenitor_like": signatures["epithelial_progenitor_like"]}
    panels.update({f"control_{gene}": [gene] for gene in control_genes})
    genes = flatten_genes(panels)
    rows = []
    for sample_index, sample in enumerate(sorted(samples, key=lambda item: item.sample_name)):
        accession = accession_by_sample.get(sample.sample_name, sample.sample_accession)
        meta = sample_metadata(metadata, accession)
        obs = read_10x_obs(sample)
        counts, total_counts = read_10x_selected_genes_with_totals(sample, genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        scores = compute_panel_scores(normalized, panels)
        table = obs.merge(scores.reset_index(names="spot"), on="spot", how="left")
        table = table[table["in_tissue"].eq(1)].copy()
        table.insert(0, "sample_name", sample.sample_name)
        table.insert(0, "sample", accession)
        table["stage"] = meta["stage"]
        for gene_index, gene in enumerate(control_genes):
            stats = adjacency_for_target(
                table,
                source_column="epithelial_progenitor_like_score",
                target_column=f"control_{gene}_score",
                quantile=args.quantile,
                radius_multiplier=args.radius_multiplier,
                permutations=args.permutations,
                seed=args.seed + sample_index * 100 + gene_index,
            )
            stats.update(
                {
                    "axis_id": gene,
                    "axis_label": f"{gene} single-gene control",
                    "evidence_type": "source_near_epithelial_progenitor",
                    "control_gene": gene,
                    "is_target_gene": gene == args.target_gene,
                }
            )
            rows.append(stats)
        print(f"Control adjacency {sample_index + 1}/{len(samples)}: {accession} {meta['stage']}")
    return pd.DataFrame(rows)


def summarize_broad_controls(path: Path, *, bootstrap_iterations: int, seed: int) -> pd.DataFrame:
    """Summarize existing broad macrophage signature adjacency as paired controls."""
    table = pd.read_csv(path)
    table["axis_id"] = table["target"].astype(str)
    table["axis_label"] = table["target"].astype(str)
    table["evidence_type"] = "broad_signature_near_epithelial_progenitor"
    paired = build_paired_patient_differences(table)
    return summarize_paired_patient_differences(
        paired,
        bootstrap_iterations=bootstrap_iterations,
        seed=seed,
    )


def summarize_density_controls(path: Path) -> pd.DataFrame:
    """Measure whether source-side MIF enrichment tracks basic section geometry."""
    table = add_patient_phase(pd.read_csv(path))
    table = table[
        table["axis_id"].eq("mif_cd74_cxcr4")
        & table["evidence_type"].eq("source_near_epithelial_progenitor")
        & table["status"].eq("ok")
    ].copy()
    records = []
    for column in ("n_spots", "null_mean", "radius"):
        record = spearman_association(table, "enrichment_delta", column)
        record["analysis_level"] = "sample"
        records.append(record)
    patient_phase = (
        table[table["phase"].isin(["precursor", "late"])]
        .groupby(["patient_id", "phase"], as_index=False)
        .agg(enrichment_delta=("enrichment_delta", "mean"), n_spots=("n_spots", "mean"))
    )
    paired = patient_phase.pivot(index="patient_id", columns="phase", values=["enrichment_delta", "n_spots"])
    paired.columns = [f"{metric}_{phase}" for metric, phase in paired.columns]
    paired = paired.dropna().reset_index()
    paired["enrichment_delta_change"] = paired["enrichment_delta_late"] - paired["enrichment_delta_precursor"]
    paired["n_spots_change"] = paired["n_spots_late"] - paired["n_spots_precursor"]
    record = spearman_association(paired, "enrichment_delta_change", "n_spots_change")
    record["analysis_level"] = "paired_patient_change"
    records.append(record)
    return pd.DataFrame(records)


def save_control_figure(
    output_stem: Path,
    random_paired: pd.DataFrame,
    random_summary: pd.DataFrame,
    broad_paired: pd.DataFrame,
    density_table: pd.DataFrame,
) -> None:
    """Render a compact supplementary control figure."""
    target_value = float(random_summary.iloc[0]["target_paired_difference_mean"])
    random_only = random_paired[~random_paired["axis_id"].eq("MIF")].copy()
    density = add_patient_phase(density_table)
    density = density[
        density["axis_id"].eq("mif_cd74_cxcr4")
        & density["evidence_type"].eq("source_near_epithelial_progenitor")
        & density["status"].eq("ok")
    ].copy()

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    ax = axes[0]
    ax.scatter(np.arange(len(random_only)), random_only["paired_difference_mean"], color="#8F8F8F", s=22)
    ax.axhline(target_value, color="#B64342", linewidth=1.6, label=f"MIF = {target_value:.3f}")
    ax.set_title("Expression-matched single-gene controls", loc="left", fontsize=8, fontweight="bold")
    ax.set_xlabel("matched control gene")
    ax.set_ylabel("paired late-minus-precursor delta")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    broad = broad_paired.sort_values("paired_difference_mean", ascending=True)
    ax.barh(broad["axis_id"], broad["paired_difference_mean"], color="#42949E")
    ax.axvline(target_value, color="#B64342", linewidth=1.4, linestyle="--", label="MIF source-side")
    ax.set_title("Broad macrophage-signature controls", loc="left", fontsize=8, fontweight="bold")
    ax.set_xlabel("paired late-minus-precursor delta")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[2]
    ax.scatter(density["n_spots"], density["enrichment_delta"], color="#5B8FD6", s=22, alpha=0.8)
    association = spearman_association(density, "enrichment_delta", "n_spots")
    ax.set_title("MIF enrichment vs tissue spots", loc="left", fontsize=8, fontweight="bold")
    ax.set_xlabel("in-tissue spots")
    ax.set_ylabel("MIF source-side enrichment")
    ax.text(
        0.04,
        0.96,
        f"Spearman r={association['spearman_r']:.2f}\np={association['spearman_p']:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=7)
    fig.suptitle("GSE307534 source-side MIF spatial controls", x=0.01, ha="left", fontsize=10, fontweight="bold")
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".svg"))
    fig.savefig(output_stem.with_suffix(".pdf"))
    fig.savefig(output_stem.with_suffix(".png"), dpi=220)
    fig.savefig(output_stem.with_suffix(".tiff"), dpi=300)
    plt.close(fig)


def write_summary(
    path: Path,
    references: list[tuple[TenXSampleFiles, str, dict[str, object]]],
    matched: pd.DataFrame,
    random_summary: pd.DataFrame,
    broad_summary: pd.DataFrame,
    density_summary: pd.DataFrame,
) -> None:
    result = random_summary.iloc[0]
    lines = [
        "# GSE307534 MIF Spatial Controls",
        "",
        "Date: 2026-05-30",
        "",
        "## Design",
        "",
        "- Selected 20 expression-matched single-gene controls using one reference Visium section per ordered stage.",
        "- Excluded candidate-axis genes, discovery-signature genes, mitochondrial genes, and ribosomal genes.",
        "- Recomputed epithelial-progenitor-neighborhood adjacency for `MIF` and matched genes across all 56 GSE307534 sections.",
        "- Reused broad macrophage-signature adjacency and basic tissue-geometry metrics as additional controls.",
        "",
        "Reference sections:",
        "",
    ]
    lines.extend(f"- `{accession}`: `{meta['stage']}` (`{sample.sample_name}`)" for sample, accession, meta in references)
    lines.extend(
        [
            "",
            "## Random-Gene Control Result",
            "",
            f"- `MIF` paired late-minus-precursor enrichment delta: {result['target_paired_difference_mean']:.3f}.",
            f"- Matched control genes: {int(result['n_control_genes'])}.",
            f"- `MIF` percentile versus matched controls: {result['target_percentile']:.3f}.",
            f"- Empirical upper-tail p-value: {result['empirical_upper_p']:.3f}.",
            "",
            "Matched genes:",
            "",
            ", ".join(f"`{gene}`" for gene in matched["gene"]),
            "",
            "## Tissue-Density Controls",
            "",
            density_summary.to_markdown(index=False),
            "",
            "## Broad Macrophage-Signature Controls",
            "",
            broad_summary[["axis_id", "n_paired_patients", "paired_difference_mean", "ci_95_low", "ci_95_high", "wilcoxon_q_bh"]].to_markdown(index=False),
            "",
            "## Interpretation Boundary",
            "",
            "These controls evaluate whether the source-side MIF spatial result behaves like generic tissue-density, broad macrophage-signature, or expression-matched single-gene effects. They do not prove causal MIF signaling.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.docs_dir.mkdir(parents=True, exist_ok=True)
    signatures = load_signatures(args.signatures)
    metadata = pd.read_csv(args.metadata)
    accession_by_sample = sample_accessions_from_filelist(args.filelist)
    samples = list(discover_10x_samples(args.input_dir, require_tissue_positions=True).values())
    references = choose_reference_samples(samples, metadata, accession_by_sample)
    gene_summary = reference_gene_summary(references, scale_factor=args.scale_factor)
    excluded = load_excluded_genes(args.mechanisms, signatures)
    matched = select_expression_matched_controls(
        gene_summary,
        target_gene=args.target_gene,
        n_controls=args.n_controls,
        min_reference_samples=len(references),
        excluded_genes=excluded,
    )
    selected_genes = [args.target_gene] + matched["gene"].tolist()
    adjacency = run_random_control_adjacency(
        samples,
        metadata,
        accession_by_sample,
        signatures,
        selected_genes,
        args,
    )
    paired_differences = build_paired_patient_differences(adjacency)
    paired_stats = summarize_paired_patient_differences(
        paired_differences,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    random_summary = summarize_random_control_distribution(paired_stats, target_gene=args.target_gene)
    broad_summary = summarize_broad_controls(
        args.broad_adjacency,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    density_summary = summarize_density_controls(args.axis_adjacency)

    gene_summary.to_csv(args.table_dir / "gse307534_mif_control_reference_gene_summary.csv", index=False)
    matched.to_csv(args.table_dir / "gse307534_mif_expression_matched_controls.csv", index=False)
    adjacency.to_csv(args.table_dir / "gse307534_mif_random_control_adjacency.csv", index=False)
    paired_differences.to_csv(args.table_dir / "gse307534_mif_random_control_paired_differences.csv", index=False)
    paired_stats.to_csv(args.table_dir / "gse307534_mif_random_control_paired_stats.csv", index=False)
    random_summary.to_csv(args.table_dir / "gse307534_mif_random_control_summary.csv", index=False)
    broad_summary.to_csv(args.table_dir / "gse307534_mif_broad_signature_control_paired_stats.csv", index=False)
    density_summary.to_csv(args.table_dir / "gse307534_mif_density_control_summary.csv", index=False)
    (args.table_dir / "gse307534_mif_control_selected_genes.json").write_text(
        json.dumps({"target_gene": args.target_gene, "matched_control_genes": matched["gene"].tolist()}, indent=2),
        encoding="utf-8",
    )
    save_control_figure(
        args.figure_dir / "supplementary_figure_mif_spatial_controls",
        paired_stats,
        random_summary,
        broad_summary,
        pd.read_csv(args.axis_adjacency),
    )
    write_summary(
        args.docs_dir / "gse307534_mif_spatial_controls.md",
        references,
        matched,
        random_summary,
        broad_summary,
        density_summary,
    )
    print(f"Wrote MIF spatial controls using {len(matched)} matched genes across {len(samples)} samples.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

