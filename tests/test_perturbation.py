import pandas as pd

from luad_niche.perturbation import apply_expression_perturbation, summarize_perturbation_effects


def test_apply_expression_perturbation_scales_present_genes_without_mutating_input():
    expr = pd.DataFrame(
        {
            "MIF": [2.0, 4.0],
            "CD74": [1.0, 3.0],
            "EPCAM": [5.0, 6.0],
        },
        index=["spot1", "spot2"],
    )

    perturbed, metadata = apply_expression_perturbation(expr, ["MIF", "CXCR4"], factor=0.5)

    assert perturbed["MIF"].tolist() == [1.0, 2.0]
    assert perturbed["CD74"].tolist() == [1.0, 3.0]
    assert perturbed["EPCAM"].tolist() == [5.0, 6.0]
    assert expr["MIF"].tolist() == [2.0, 4.0]
    assert metadata == {"present_genes": ["MIF"], "missing_genes": ["CXCR4"], "factor": 0.5}


def test_summarize_perturbation_effects_reports_mean_delta_by_stage():
    effects = pd.DataFrame(
        {
            "perturbation_id": ["MIF_ko", "MIF_ko", "MIF_ko"],
            "axis_id": ["mif_axis", "mif_axis", "mif_axis"],
            "evidence_type": ["source", "source", "target"],
            "stage": ["MIA", "MIA", "MIA"],
            "sample": ["A", "B", "A"],
            "baseline_panel_mean": [2.0, 4.0, 1.0],
            "perturbed_panel_mean": [1.0, 2.0, 0.8],
            "panel_mean_delta": [-1.0, -2.0, -0.2],
            "baseline_observed_fraction": [0.5, 0.7, 0.3],
            "perturbed_observed_fraction": [0.25, 0.5, 0.2],
            "observed_fraction_delta": [-0.25, -0.2, -0.1],
        }
    )

    summary = summarize_perturbation_effects(effects)

    source = summary[summary["evidence_type"].eq("source")].iloc[0]
    assert source["n_samples"] == 2
    assert source["panel_mean_delta_mean"] == -1.5
    assert source["observed_fraction_delta_mean"] == -0.225

    target = summary[summary["evidence_type"].eq("target")].iloc[0]
    assert target["n_samples"] == 1
    assert target["panel_mean_delta_mean"] == -0.2
