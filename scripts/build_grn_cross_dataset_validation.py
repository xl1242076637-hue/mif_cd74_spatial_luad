#!/usr/bin/env python
"""Build cross-dataset validation summaries for GRN-prioritized signatures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.grn_validation import (  # noqa: E402
    canonical_signature_state,
    expected_celltype_for_state,
    summarize_group_delta,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--grn-target-ranking",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_grn_virtual_perturbation_target_ranking.csv",
    )
    parser.add_argument(
        "--gse307534-spatial",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_refined_signature_score_by_stage.csv",
    )
    parser.add_argument(
        "--gse189357-macrophage",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_scrna_macrophage_state_stage_summary.csv",
    )
    parser.add_argument(
        "--gse189357-epithelial",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_scrna_epithelial_state_stage_summary.csv",
    )
    parser.add_argument(
        "--gse164789-macrophage",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse164789_scrna_macrophage_state_stage_summary.csv",
    )
    parser.add_argument(
        "--gse164789-epithelial",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse164789_scrna_epithelial_state_stage_summary.csv",
    )
    parser.add_argument(
        "--gse131907-specificity",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse131907_selected_signature_celltype_summary.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    return parser.parse_args()


def trend_direction(delta: float) -> str:
    if pd.isna(delta):
        return "not_available"
    if delta > 0:
        return "higher_in_comparison"
    if delta < 0:
        return "lower_in_comparison"
    return "unchanged"


def state_table_for_expected_celltype(
    expected_celltype: str,
    macrophage_table: pd.DataFrame,
    epithelial_table: pd.DataFrame,
) -> tuple[pd.DataFrame, str]:
    if expected_celltype == "Myeloid cells":
        return macrophage_table, "macrophage_state"
    if expected_celltype == "Epithelial cells":
        return epithelial_table, "epithelial_state"
    return pd.DataFrame(), "state"


def add_delta_row(
    rows: list[dict[str, object]],
    *,
    target: pd.Series,
    dataset: str,
    validation_layer: str,
    table: pd.DataFrame,
    state_column: str,
    value_column: str,
    baseline_groups: list[str],
    comparison_groups: list[str],
    metric_name: str,
) -> None:
    state = target["canonical_state"]
    summary = summarize_group_delta(
        table,
        state_column=state_column,
        state_value=state,
        group_column="stage",
        value_column=value_column,
        baseline_groups=baseline_groups,
        comparison_groups=comparison_groups,
    )
    rows.append(
        {
            "target_gene": target["target_gene"],
            "broad_class": target["broad_class"],
            "grn_top_signature": target["top_impacted_signature"],
            "canonical_state": state,
            "dataset": dataset,
            "validation_layer": validation_layer,
            "metric_name": metric_name,
            "baseline_groups": summary["baseline_groups"],
            "comparison_groups": summary["comparison_groups"],
            "baseline_mean": summary["baseline_mean"],
            "comparison_mean": summary["comparison_mean"],
            "delta": summary["delta"],
            "direction": trend_direction(summary["delta"]),
            "baseline_n_samples": summary["baseline_n_samples"],
            "comparison_n_samples": summary["comparison_n_samples"],
            "expected_celltype": target["expected_celltype"],
            "top_celltype": "",
            "top_celltype_mean_score": "",
            "expected_celltype_mean_score": "",
        }
    )


def add_specificity_row(
    rows: list[dict[str, object]],
    *,
    target: pd.Series,
    specificity: pd.DataFrame,
) -> None:
    state = target["canonical_state"]
    expected = target["expected_celltype"]
    subset = specificity[specificity["score"].astype(str).eq(str(state))].copy()
    subset["mean_score"] = pd.to_numeric(subset["mean_score"], errors="coerce").fillna(0.0)
    if subset.empty:
        top_celltype = ""
        top_score = float("nan")
        expected_score = float("nan")
        direction = "not_available"
    else:
        top = subset.sort_values("mean_score", ascending=False).iloc[0]
        top_celltype = str(top["Cell_type.refined"])
        top_score = float(top["mean_score"])
        expected_rows = subset[subset["Cell_type.refined"].astype(str).eq(expected)]
        expected_score = float(expected_rows["mean_score"].iloc[0]) if not expected_rows.empty else float("nan")
        direction = "expected_top_celltype" if top_celltype == expected else "off_target_top_celltype"
    rows.append(
        {
            "target_gene": target["target_gene"],
            "broad_class": target["broad_class"],
            "grn_top_signature": target["top_impacted_signature"],
            "canonical_state": state,
            "dataset": "GSE131907",
            "validation_layer": "celltype_specificity",
            "metric_name": "top_celltype_mean_score",
            "baseline_groups": "",
            "comparison_groups": "",
            "baseline_mean": "",
            "comparison_mean": "",
            "delta": "",
            "direction": direction,
            "baseline_n_samples": "",
            "comparison_n_samples": "",
            "expected_celltype": expected,
            "top_celltype": top_celltype,
            "top_celltype_mean_score": top_score,
            "expected_celltype_mean_score": expected_score,
        }
    )


def summarize_validation(validation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (target_gene, state), group in validation.groupby(["target_gene", "canonical_state"], sort=False):
        rows.append(
            {
                "target_gene": target_gene,
                "canonical_state": state,
                "n_validation_rows": len(group),
                "n_expected_specificity_rows": int(group["direction"].eq("expected_top_celltype").sum()),
                "n_higher_in_comparison_rows": int(group["direction"].eq("higher_in_comparison").sum()),
                "n_lower_in_comparison_rows": int(group["direction"].eq("lower_in_comparison").sum()),
                "datasets": ",".join(group["dataset"].astype(str).drop_duplicates()),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    target_ranking = pd.read_csv(args.grn_target_ranking)
    target_ranking["canonical_state"] = target_ranking["top_impacted_signature"].map(canonical_signature_state)
    target_ranking["expected_celltype"] = target_ranking["canonical_state"].map(expected_celltype_for_state)

    gse307534 = pd.read_csv(args.gse307534_spatial)
    gse189357_mac = pd.read_csv(args.gse189357_macrophage)
    gse189357_epi = pd.read_csv(args.gse189357_epithelial)
    gse164789_mac = pd.read_csv(args.gse164789_macrophage)
    gse164789_epi = pd.read_csv(args.gse164789_epithelial)
    gse131907 = pd.read_csv(args.gse131907_specificity)

    rows: list[dict[str, object]] = []
    for _, target in target_ranking.iterrows():
        add_delta_row(
            rows,
            target=target,
            dataset="GSE307534",
            validation_layer="visium_spatial_signature_score",
            table=gse307534,
            state_column="signature",
            value_column="mean_score",
            baseline_groups=["AAH", "AIS"],
            comparison_groups=["MIA", "LUAD"],
            metric_name="late_minus_precursor_mean_score",
        )
        expected = target["expected_celltype"]
        gse189357_table, gse189357_state_column = state_table_for_expected_celltype(
            expected,
            gse189357_mac,
            gse189357_epi,
        )
        add_delta_row(
            rows,
            target=target,
            dataset="GSE189357",
            validation_layer="scrna_state_fraction",
            table=gse189357_table,
            state_column=gse189357_state_column,
            value_column="mean_fraction",
            baseline_groups=["AIS"],
            comparison_groups=["IAC"],
            metric_name="iac_minus_ais_fraction",
        )
        gse164789_table, gse164789_state_column = state_table_for_expected_celltype(
            expected,
            gse164789_mac,
            gse164789_epi,
        )
        add_delta_row(
            rows,
            target=target,
            dataset="GSE164789",
            validation_layer="scrna_tumor_adjacent_state_fraction",
            table=gse164789_table,
            state_column=gse164789_state_column,
            value_column="mean_fraction",
            baseline_groups=["Adjacent"],
            comparison_groups=["Tumor"],
            metric_name="tumor_minus_adjacent_fraction",
        )
        add_specificity_row(rows, target=target, specificity=gse131907)

    validation = pd.DataFrame(rows)
    summary = summarize_validation(validation)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    validation_path = args.table_dir / "grn_cross_dataset_signature_validation.csv"
    summary_path = args.table_dir / "grn_cross_dataset_signature_validation_summary.csv"
    validation.to_csv(validation_path, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    print(f"Wrote cross-dataset validation: {validation_path}")
    print(f"Wrote cross-dataset validation summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
