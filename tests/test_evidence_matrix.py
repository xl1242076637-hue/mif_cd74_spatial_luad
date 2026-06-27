import pandas as pd

from luad_niche.evidence_matrix import (
    best_continuous_perturbation_by_axis,
    mean_axis_spatial_by_evidence,
)


def test_mean_axis_spatial_by_evidence_averages_mia_and_luad_rows():
    spatial = pd.DataFrame(
        {
            "axis_id": ["mif_axis", "mif_axis", "mif_axis", "other"],
            "evidence_type": [
                "source_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
            ],
            "stage": ["MIA", "LUAD", "LUAD", "MIA"],
            "enrichment_delta_mean": [0.2, 0.4, 0.1, 1.0],
        }
    )

    result = mean_axis_spatial_by_evidence(spatial)

    source = result[result["axis_id"].eq("mif_axis") & result["evidence_type"].str.startswith("source")].iloc[0]
    assert source["mean_axis_spatial_delta"] == 0.3

    target = result[result["axis_id"].eq("mif_axis") & result["evidence_type"].str.startswith("target")].iloc[0]
    assert target["mean_axis_spatial_delta"] == 0.1


def test_best_continuous_perturbation_by_axis_selects_highest_priority_per_axis():
    ranking = pd.DataFrame(
        {
            "axis_id": ["mif_axis", "mif_axis", "spp1_axis"],
            "perturbed_genes": ["CXCR4", "CD74", "SPP1"],
            "perturbation_factor": [0.0, 0.0, 0.0],
            "evidence_type": ["target", "target", "source"],
            "coupling_relative_delta_mean": [-0.1, -0.7, -0.2],
            "continuous_priority_score": [0.1, 0.7, 0.2],
        }
    )

    result = best_continuous_perturbation_by_axis(ranking)

    mif = result[result["axis_id"].eq("mif_axis")].iloc[0]
    assert mif["top_perturbed_genes"] == "CD74"
    assert mif["top_continuous_priority_score"] == 0.7
    assert mif["top_coupling_relative_delta"] == -0.7

    spp1 = result[result["axis_id"].eq("spp1_axis")].iloc[0]
    assert spp1["top_perturbed_genes"] == "SPP1"
