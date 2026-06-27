import math

import pandas as pd

from luad_niche.spatial_controls import (
    build_paired_change_table,
    fit_ols_sensitivity,
    select_expression_matched_controls,
    spearman_association,
    summarize_random_control_distribution,
)


def test_select_expression_matched_controls_excludes_candidates_and_ribosomal_genes():
    summary = pd.DataFrame(
        {
            "gene": ["MIF", "A", "B", "RPL1", "C", "EXCLUDED"],
            "mean_expression": [2.0, 2.1, 1.9, 2.0, 5.0, 2.01],
            "n_reference_samples": [5, 5, 5, 5, 5, 5],
        }
    )

    controls = select_expression_matched_controls(
        summary,
        target_gene="MIF",
        n_controls=2,
        min_reference_samples=5,
        excluded_genes=["EXCLUDED"],
    )

    assert controls["gene"].tolist() == ["A", "B"]


def test_summarize_random_control_distribution_reports_empirical_upper_tail():
    paired = pd.DataFrame(
        {
            "axis_id": ["MIF", "A", "B", "C"],
            "paired_difference_mean": [0.5, 0.1, 0.2, 0.6],
        }
    )

    summary = summarize_random_control_distribution(paired, target_gene="MIF").iloc[0]

    assert summary["n_control_genes"] == 3
    assert math.isclose(summary["target_percentile"], 2 / 3)
    assert math.isclose(summary["empirical_upper_p"], 0.5)


def test_spearman_association_reports_monotonic_relationship():
    table = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})

    result = spearman_association(table, "x", "y")

    assert result["n"] == 4
    assert math.isclose(result["spearman_r"], 1.0)


def test_fit_ols_sensitivity_reports_adjusted_binary_effect():
    table = pd.DataFrame(
        {
            "late": [0, 0, 1, 1],
            "covariate": [0.0, 1.0, 0.0, 1.0],
            "outcome": [1.0, 2.0, 4.0, 5.0],
        }
    )

    result = fit_ols_sensitivity(
        table,
        outcome_column="outcome",
        effect_column="late",
        covariate_columns=["covariate"],
        model_id="toy",
    )

    assert result["model_id"] == "toy"
    assert result["n_observations"] == 4
    assert math.isclose(result["effect_estimate"], 3.0)
    assert result["n_covariates"] == 1


def test_build_paired_change_table_returns_late_minus_precursor_covariate_changes():
    table = pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P2", "P2", "P3"],
            "phase": ["precursor", "late", "precursor", "late", "late"],
            "enrichment_delta": [0.2, 0.5, 0.1, 0.0, 0.8],
            "n_spots": [10, 15, 5, 7, 20],
            "null_mean": [0.3, 0.4, 0.2, 0.2, 0.5],
        }
    )

    result = build_paired_change_table(
        table,
        value_column="enrichment_delta",
        covariate_columns=["n_spots", "null_mean"],
    )

    assert result["patient_id"].tolist() == ["P1", "P2"]
    assert result["enrichment_delta_change"].tolist() == [0.3, -0.1]
    assert result["n_spots_change"].tolist() == [5.0, 2.0]
