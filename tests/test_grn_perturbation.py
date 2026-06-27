import pandas as pd

from luad_niche.grn_perturbation import (
    build_correlation_network,
    summarize_target_ranking_stability,
    summarize_signature_impacts,
    virtual_outgoing_knockout,
)


def test_build_correlation_network_returns_row_normalized_positive_network():
    expr = pd.DataFrame(
        {
            "A": [0, 1, 2, 3, 4],
            "B": [0, 2, 4, 6, 8],
            "C": [4, 3, 2, 1, 0],
            "D": [1, 1, 1, 1, 1],
        }
    )

    network = build_correlation_network(expr, min_abs_correlation=0.1)

    assert list(network.index) == ["A", "B"]
    assert list(network.columns) == ["A", "B"]
    assert network.loc["A", "A"] == 0
    assert network.loc["A", "B"] > 0
    assert network.sum(axis=1).round(6).tolist() == [1.0, 1.0]


def test_virtual_outgoing_knockout_prioritizes_downstream_neighbors():
    network = pd.DataFrame(
        {
            "A": [0.0, 0.0, 0.0],
            "B": [0.8, 0.0, 0.0],
            "C": [0.2, 1.0, 0.0],
        },
        index=["A", "B", "C"],
    )

    effects = virtual_outgoing_knockout(network, ["A"], restart=0.0, n_steps=2)

    assert effects.loc[effects["gene"].eq("B"), "impact_score"].iloc[0] > 0
    assert effects.loc[effects["gene"].eq("C"), "impact_score"].iloc[0] > 0
    assert effects.sort_values("impact_score", ascending=False)["gene"].iloc[0] == "B"
    assert "A" not in effects["gene"].tolist()


def test_summarize_signature_impacts_uses_present_genes_only():
    effects = pd.DataFrame(
        {
            "gene": ["B", "C", "D"],
            "impact_score": [0.6, 0.2, 0.0],
        }
    )

    summary = summarize_signature_impacts(
        effects,
        {"state_one": ["B", "C", "MISSING"], "state_two": ["D"]},
    )

    row = summary[summary["signature"].eq("state_one")].iloc[0]
    assert row["n_present_genes"] == 2
    assert row["mean_impact_score"] == 0.4
    assert row["max_impact_score"] == 0.6


def test_summarize_target_ranking_stability_tracks_rank_and_signature_mode():
    detail = pd.DataFrame(
        {
            "run_id": ["r1", "r1", "r2", "r2", "r3", "r3"],
            "broad_class": ["macrophage"] * 6,
            "target_gene": ["CD74", "CXCR4", "CD74", "CXCR4", "CD74", "CXCR4"],
            "display_rank": [1, 2, 2, 1, 1, 2],
            "top_impacted_signature": [
                "c1q_macrophage_vs_other_macrophage",
                "spp1_macrophage_vs_other_macrophage",
                "c1q_macrophage_vs_other_macrophage",
                "spp1_macrophage_vs_other_macrophage",
                "resident_macrophage_vs_other_macrophage",
                "spp1_macrophage_vs_other_macrophage",
            ],
            "top_signature_mean_impact": [0.2, 0.1, 0.15, 0.16, 0.18, 0.11],
        }
    )

    summary = summarize_target_ranking_stability(detail)

    cd74 = summary[summary["target_gene"].eq("CD74")].iloc[0]
    assert cd74["n_runs"] == 3
    assert cd74["median_rank"] == 1.0
    assert cd74["rank_min"] == 1
    assert cd74["rank_max"] == 2
    assert cd74["top_signature_mode"] == "c1q_macrophage_vs_other_macrophage"
    assert cd74["top_signature_mode_fraction"] == 2 / 3
