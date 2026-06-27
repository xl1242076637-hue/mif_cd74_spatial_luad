import math

import pandas as pd

from luad_niche.spatial_axis import build_axis_gene_panels, summarize_axis_adjacency_by_stage


def test_build_axis_gene_panels_deduplicates_source_and_target_genes():
    axes = [
        {
            "id": "spp1_axis",
            "source_genes": ["SPP1", "TREM2", "SPP1"],
            "target_genes": ["ITGAV", "ITGB1"],
        },
        {
            "id": "empty_axis",
            "source_genes": [],
            "target_genes": [],
        },
    ]

    panels = build_axis_gene_panels(axes)

    assert panels == {
        "spp1_axis_source": ["SPP1", "TREM2"],
        "spp1_axis_target": ["ITGAV", "ITGB1"],
    }


def test_summarize_axis_adjacency_by_stage_excludes_invalid_effects_from_means():
    adjacency = pd.DataFrame(
        {
            "axis_id": ["spp1_axis", "spp1_axis", "spp1_axis"],
            "evidence_type": ["source_near_epithelial", "source_near_epithelial", "target_near_epithelial"],
            "stage": ["MIA", "MIA", "MIA"],
            "sample": ["A", "B", "A"],
            "status": ["ok", "insufficient_high_spots", "ok"],
            "enrichment_delta": [0.2, math.nan, -0.1],
            "observed_fraction": [0.5, math.nan, 0.1],
            "null_mean": [0.3, math.nan, 0.2],
            "empirical_p_greater": [0.03, math.nan, 0.9],
            "empirical_p_less": [0.98, math.nan, 0.05],
        }
    )

    summary = summarize_axis_adjacency_by_stage(adjacency)

    source = summary[
        summary["axis_id"].eq("spp1_axis") & summary["evidence_type"].eq("source_near_epithelial")
    ].iloc[0]
    assert source["n_tests"] == 2
    assert source["n_valid_tests"] == 1
    assert source["n_invalid_tests"] == 1
    assert source["enrichment_delta_mean"] == 0.2

    target = summary[
        summary["axis_id"].eq("spp1_axis") & summary["evidence_type"].eq("target_near_epithelial")
    ].iloc[0]
    assert target["enrichment_delta_mean"] == -0.1
