import pandas as pd

from luad_niche.spatial_coupling import continuous_spatial_coupling, summarize_continuous_effects


def test_continuous_spatial_coupling_uses_neighbor_target_scores_without_self():
    table = pd.DataFrame(
        {
            "x": [0.0, 1.0, 4.0],
            "y": [0.0, 0.0, 0.0],
            "source_score": [2.0, 0.0, 1.0],
            "target_score": [10.0, 4.0, 20.0],
        }
    )

    result = continuous_spatial_coupling(table, "source_score", "target_score", radius=1.5)

    assert result["n_spots"] == 3
    assert result["n_spots_with_neighbors"] == 2
    assert result["source_mean"] == 1.0
    assert result["target_mean"] == 34.0 / 3.0
    assert result["neighbor_target_mean"] == 7.0
    assert result["coupling_score"] == 4.0
    assert result["source_weighted_neighbor_target_mean"] == 4.0


def test_summarize_continuous_effects_reports_mia_luad_mean_deltas():
    effects = pd.DataFrame(
        {
            "perturbation_id": ["MIF_x0", "MIF_x0", "MIF_x0"],
            "axis_id": ["mif_axis", "mif_axis", "mif_axis"],
            "evidence_type": ["source", "source", "source"],
            "stage": ["MIA", "LUAD", "AIS"],
            "sample": ["A", "B", "C"],
            "baseline_coupling_score": [2.0, 4.0, 8.0],
            "perturbed_coupling_score": [1.0, 2.0, 7.0],
            "coupling_delta": [-1.0, -2.0, -1.0],
            "coupling_relative_delta": [-0.5, -0.5, -0.125],
        }
    )

    summary = summarize_continuous_effects(effects)

    row = summary.iloc[0]
    assert row["n_samples"] == 2
    assert row["coupling_delta_mean"] == -1.5
    assert row["coupling_relative_delta_mean"] == -0.5
