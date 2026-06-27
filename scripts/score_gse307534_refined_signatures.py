#!/usr/bin/env python
"""Map refined GSE189357 signatures onto GSE307534 LungPCA Visium samples."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.spatial_niche import (  # noqa: E402
    finite_coordinate_mask,
    median_nearest_neighbor_distance,
    nearest_target_fraction,
    permutation_adjacency_test,
    top_quantile_mask,
)
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402


STAGE_ORDER = ["Normal", "AAH", "AIS", "MIA", "LUAD"]
FILELIST_SAMPLE_RE = re.compile(r"^(GSM\d+)_(.+)\.tar\.gz$")
DUPLICATE_SUFFIX_RE = re.compile(r"^(.+)-(\d+)$")


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
    parser.add_argument("--quantile", type=float, default=0.75)
    parser.add_argument("--radius-multiplier", type=float, default=1.0)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=307534)
    parser.add_argument("--limit-samples", type=int, default=None)
    return parser.parse_args()


def clean_signature_name(name: str) -> str:
    return name.replace("_vs_other_epithelial", "").replace("_vs_other_macrophage", "")


def load_signatures(path: Path) -> dict[str, list[str]]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return {clean_signature_name(name): list(genes) for name, genes in raw.items()}


def flatten_genes(panels: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panel_genes in panels.values():
        genes.extend(gene for gene in panel_genes if gene not in genes)
    return genes


def sample_accessions_from_filelist(path: Path) -> dict[str, str]:
    entries: list[tuple[str, str]] = []
    mapping: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 2 or parts[0] != "File":
                continue
            match = FILELIST_SAMPLE_RE.match(parts[1])
            if match:
                accession = match.group(1)
                sample_name = match.group(2)
                entries.append((sample_name, accession))
                mapping[sample_name] = accession
    duplicate_groups: dict[str, list[tuple[str, str, int | None]]] = {}
    for sample_name, accession in entries:
        duplicate_match = DUPLICATE_SUFFIX_RE.match(sample_name)
        if duplicate_match:
            base_name = duplicate_match.group(1)
            duplicate_index = int(duplicate_match.group(2))
        else:
            base_name = sample_name
            duplicate_index = None
        duplicate_groups.setdefault(base_name, []).append((sample_name, accession, duplicate_index))
    for base_name, records in duplicate_groups.items():
        if not any(duplicate_index is not None for _, _, duplicate_index in records):
            continue
        for sample_name, accession, duplicate_index in records:
            alias_index = 1 if duplicate_index is None else duplicate_index + 1
            mapping[f"{base_name}{alias_index}"] = accession
    return mapping


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


def adjacency_for_target(
    table: pd.DataFrame,
    source_column: str,
    target_column: str,
    quantile: float,
    radius_multiplier: float,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    usable = table[table[["x", "y", source_column, target_column]].notna().all(axis=1)].copy()
    finite_mask = finite_coordinate_mask(usable[["x", "y"]].to_numpy(dtype=float))
    usable = usable.loc[finite_mask].reset_index(drop=True)
    context = usable.iloc[0] if len(usable) else table.iloc[0]
    coords = usable[["x", "y"]].to_numpy(dtype=float)
    if len(usable) >= 2:
        radius = radius_multiplier * median_nearest_neighbor_distance(coords)
        source_high = top_quantile_mask(usable[source_column], quantile).to_numpy()
        target_high = top_quantile_mask(usable[target_column], quantile).to_numpy()
    else:
        radius = float("nan")
        source_high = pd.Series([], dtype=bool).to_numpy()
        target_high = pd.Series([], dtype=bool).to_numpy()
    overlap = source_high & target_high
    n_source_high = int(source_high.sum())
    n_target_high = int(target_high.sum())
    if len(usable) < 2:
        stats = invalid_adjacency_stats("insufficient_spots", radius, permutations)
    elif n_source_high == 0 or n_target_high == 0:
        stats = invalid_adjacency_stats("insufficient_high_spots", radius, permutations)
    else:
        stats = permutation_adjacency_test(
            coords,
            source_high,
            target_high,
            radius=radius,
            n_permutations=permutations,
            seed=seed,
        )
        stats["status"] = "ok"
    stats.update(
        {
            "sample": context["sample"],
            "sample_name": context["sample_name"],
            "stage": context["stage"],
            "source": source_column.removesuffix("_score"),
            "target": target_column.removesuffix("_score"),
            "quantile": quantile,
            "radius_multiplier": radius_multiplier,
            "n_spots": len(usable),
            "n_source_high": n_source_high,
            "n_target_high": n_target_high,
            "n_overlap_high": int(overlap.sum()),
            "overlap_fraction_of_source": (
                float(overlap.sum() / n_source_high) if n_source_high else float("nan")
            ),
            "target_to_source_fraction": (
                nearest_target_fraction(coords, target_high, source_high, radius=radius)
                if n_source_high and n_target_high
                else float("nan")
            ),
        }
    )
    return stats


def invalid_adjacency_stats(status: str, radius: float, n_permutations: int) -> dict[str, object]:
    return {
        "observed_fraction": float("nan"),
        "null_mean": float("nan"),
        "null_sd": float("nan"),
        "enrichment_delta": float("nan"),
        "empirical_p_greater": float("nan"),
        "empirical_p_less": float("nan"),
        "n_permutations": n_permutations,
        "radius": radius,
        "status": status,
    }


def summarize_adjacency_by_target_stage(adjacency: pd.DataFrame) -> pd.DataFrame:
    summary_input = adjacency.copy()
    if "status" not in summary_input.columns:
        summary_input["status"] = "ok"
    group_columns = ["target", "stage"]
    all_counts = (
        summary_input.groupby(group_columns, dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            n_tests=("sample", "size"),
            n_invalid_tests=("status", lambda values: int((values != "ok").sum())),
        )
        .reset_index()
    )
    valid = summary_input[summary_input["status"] == "ok"].copy()
    if valid.empty:
        valid_summary = pd.DataFrame(columns=group_columns + [
            "n_valid_samples",
            "n_valid_tests",
            "observed_fraction_mean",
            "null_mean_mean",
            "enrichment_delta_mean",
            "empirical_p_greater_median",
            "empirical_p_less_median",
        ])
    else:
        valid_summary = (
            valid.groupby(group_columns, dropna=False)
            .agg(
                n_valid_samples=("sample", "nunique"),
                n_valid_tests=("sample", "size"),
                observed_fraction_mean=("observed_fraction", "mean"),
                null_mean_mean=("null_mean", "mean"),
                enrichment_delta_mean=("enrichment_delta", "mean"),
                empirical_p_greater_median=("empirical_p_greater", "median"),
                empirical_p_less_median=("empirical_p_less", "median"),
            )
            .reset_index()
        )
    summary = all_counts.merge(valid_summary, on=group_columns, how="left")
    for column in ["n_valid_samples", "n_valid_tests"]:
        summary[column] = summary[column].fillna(0).astype(int)
    stage_order = {stage: index for index, stage in enumerate(STAGE_ORDER)}
    summary["_stage_order"] = summary["stage"].map(stage_order).fillna(len(stage_order))
    summary = summary.sort_values(["target", "_stage_order", "stage"]).drop(columns="_stage_order")
    return summary


def summarize_score_by_stage(scores: pd.DataFrame) -> pd.DataFrame:
    score_columns = [column for column in scores.columns if column.endswith("_score")]
    records = []
    for column in score_columns:
        sample_means = (
            scores.groupby(["sample", "sample_name", "stage"])[column]
            .mean()
            .reset_index(name="mean_score")
        )
        sample_means["signature"] = column.removesuffix("_score")
        records.append(sample_means)
    sample_long = pd.concat(records, ignore_index=True)
    return (
        sample_long.groupby(["stage", "signature"])
        .agg(n_samples=("sample", "nunique"), mean_score=("mean_score", "mean"), sd_score=("mean_score", "std"))
        .reset_index()
    )


def plot_stage_target_summary(stage_summary: pd.DataFrame, output: Path) -> None:
    pivot = stage_summary.pivot_table(index="stage", columns="target", values="enrichment_delta_mean", fill_value=0)
    order = [stage for stage in STAGE_ORDER if stage in pivot.index]
    pivot = pivot.loc[order]
    ax = pivot.plot(kind="bar", figsize=(9, 4.8), width=0.8)
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_ylabel("Observed minus permuted-null adjacency")
    ax.set_xlabel("")
    ax.set_title("GSE307534 refined signature adjacency by stage")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    ax.figure.savefig(output, dpi=220)
    plt.close(ax.figure)


def main() -> int:
    args = parse_args()
    signatures = load_signatures(args.signatures)
    all_genes = flatten_genes(signatures)
    metadata = pd.read_csv(args.metadata)
    accession_by_sample = sample_accessions_from_filelist(args.filelist)
    samples = list(discover_10x_samples(args.input_dir, require_tissue_positions=True).values())
    samples = sorted(samples, key=lambda sample: sample.sample_name)
    if args.limit_samples:
        samples = samples[: args.limit_samples]

    score_tables = []
    adjacency_rows = []
    genes_used = {}
    for sample_index, sample in enumerate(samples):
        sample_accession = accession_by_sample.get(sample.sample_name, sample.sample_accession)
        meta = sample_metadata(metadata, sample_accession)
        obs = read_10x_obs(sample)
        counts, total_counts = read_10x_selected_genes_with_totals(sample, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        scores = compute_panel_scores(normalized, signatures)
        genes_used[sample_accession] = scores.attrs["panel_genes_used"]
        table = obs.merge(total_counts.rename_axis("spot").reset_index(), on="spot", how="left")
        table = table.merge(scores.reset_index(names="spot"), on="spot", how="left")
        table = table[table["in_tissue"] == 1].copy()
        table.insert(0, "sample_name", sample.sample_name)
        table.insert(0, "sample", sample_accession)
        for key, value in meta.items():
            table[key] = value
        score_tables.append(table)

        source_column = "epithelial_progenitor_like_score"
        target_columns = [
            "spp1_macrophage_score",
            "c1q_macrophage_score",
            "inflammatory_macrophage_score",
            "resident_macrophage_score",
        ]
        target_columns = [column for column in target_columns if column in table.columns]
        for target_index, target_column in enumerate(target_columns):
            adjacency = adjacency_for_target(
                table,
                source_column=source_column,
                target_column=target_column,
                quantile=args.quantile,
                radius_multiplier=args.radius_multiplier,
                permutations=args.permutations,
                seed=args.seed + sample_index * 10 + target_index,
            )
            adjacency_rows.append(adjacency)
            if adjacency["status"] == "ok":
                print(
                    f"{sample_accession} {sample.sample_name} {meta['stage']} {target_column}: "
                    f"delta={adjacency['enrichment_delta']:.3f}; "
                    f"p_greater={adjacency['empirical_p_greater']:.4f}"
                )
            else:
                print(
                    f"{sample_accession} {sample.sample_name} {meta['stage']} {target_column}: "
                    f"status={adjacency['status']}; "
                    f"n_source_high={adjacency['n_source_high']}; "
                    f"n_target_high={adjacency['n_target_high']}"
                )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.concat(score_tables, ignore_index=True)
    adjacency = pd.DataFrame(adjacency_rows)
    stage_summary = summarize_adjacency_by_target_stage(adjacency)
    score_stage = summarize_score_by_stage(scores)

    scores.to_csv(args.table_dir / "gse307534_refined_signature_spatial_scores.csv", index=False, encoding="utf-8-sig")
    adjacency.to_csv(args.table_dir / "gse307534_refined_signature_adjacency.csv", index=False, encoding="utf-8-sig")
    stage_summary.to_csv(
        args.table_dir / "gse307534_refined_signature_adjacency_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    score_stage.to_csv(
        args.table_dir / "gse307534_refined_signature_score_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.table_dir / "gse307534_refined_signature_genes_used.json").write_text(
        json.dumps(genes_used, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    plot_stage_target_summary(
        stage_summary,
        args.figure_dir / "gse307534_refined_signature_adjacency_by_stage.png",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
