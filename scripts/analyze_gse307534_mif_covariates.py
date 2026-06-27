#!/usr/bin/env python
"""Run covariate-adjusted sensitivity models for source-side MIF spatial enrichment."""

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

from luad_niche.spatial_controls import build_paired_change_table, fit_ols_sensitivity  # noqa: E402
from luad_niche.spatial_statistics import add_patient_phase  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mif-adjacency",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_mif_random_control_adjacency.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--docs-dir", type=Path, default=PROJECT_ROOT / "docs")
    return parser.parse_args()


def source_mif_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    table = add_patient_phase(table)
    table = table[
        table["axis_id"].eq("MIF")
        & table["evidence_type"].eq("source_near_epithelial_progenitor")
        & table["status"].eq("ok")
    ].copy()
    table = table[table["phase"].isin(["precursor", "late"])].copy()
    table["late_indicator"] = table["phase"].eq("late").astype(float)
    return table


def run_models(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    covariate_sets = [
        ("unadjusted", []),
        ("adjusted_n_spots", ["n_spots"]),
        ("adjusted_n_spots_null_mean", ["n_spots", "null_mean"]),
        ("adjusted_full_geometry", ["n_spots", "null_mean", "radius"]),
    ]
    sample_records = [
        fit_ols_sensitivity(
            table,
            outcome_column="enrichment_delta",
            effect_column="late_indicator",
            covariate_columns=covariates,
            model_id=f"sample_{name}",
        )
        for name, covariates in covariate_sets
    ]
    paired = build_paired_change_table(
        table,
        value_column="enrichment_delta",
        covariate_columns=["n_spots", "null_mean", "radius"],
    )
    paired_records = [
        fit_ols_sensitivity(
            paired,
            outcome_column="enrichment_delta_change",
            effect_column="intercept",
            covariate_columns=[f"{column}_change" for column in covariates],
            model_id=f"paired_change_{name}",
        )
        for name, covariates in covariate_sets
    ]
    return pd.DataFrame(sample_records + paired_records), paired


def write_summary(path: Path, models: pd.DataFrame, paired: pd.DataFrame) -> None:
    lines = [
        "# GSE307534 MIF Covariate Sensitivity",
        "",
        "Date: 2026-05-31",
        "",
        "## Design",
        "",
        "- Outcome: source-side `MIF` epithelial-progenitor-neighborhood enrichment delta.",
        "- Sample-level models use AAH/AIS as precursor and MIA/LUAD as late lesions.",
        "- Paired-change models use patient-level late-minus-precursor changes.",
        "- Covariates tested: in-tissue spot count, permutation-null mean, and neighborhood radius.",
        "",
        "## Model Summary",
        "",
        models.to_markdown(index=False),
        "",
        "## Interpretation Boundary",
        "",
        "These are linear sensitivity models for robustness assessment. They do not prove causal MIF signaling and do not replace the paired-patient spatial statistics.",
        "",
        f"Paired patients retained in change models: {len(paired)}.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    table = source_mif_table(args.mif_adjacency)
    models, paired = run_models(table)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.docs_dir.mkdir(parents=True, exist_ok=True)
    models_output = args.table_dir / "gse307534_mif_covariate_sensitivity_models.csv"
    paired_output = args.table_dir / "gse307534_mif_covariate_paired_changes.csv"
    models.to_csv(models_output, index=False, encoding="utf-8-sig")
    paired.to_csv(paired_output, index=False, encoding="utf-8-sig")
    write_summary(args.docs_dir / "gse307534_mif_covariate_sensitivity.md", models, paired)
    print(f"Wrote models: {models_output}")
    print(f"Wrote paired changes: {paired_output}")
    print(f"Paired patients: {len(paired)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
