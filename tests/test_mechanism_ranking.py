import pandas as pd

from luad_niche.mechanism_ranking import (
    gene_specificity_fraction,
    mean_axis_spatial_delta,
    rank_candidate_axes,
    summarize_perturbation_candidates,
)


def test_gene_specificity_fraction_counts_expected_present_genes_only():
    top = pd.DataFrame(
        {
            "gene": ["SPP1", "TREM2", "WFDC2"],
            "Cell_type.refined": ["Myeloid cells", "Myeloid cells", "Epithelial cells"],
        }
    )

    result = gene_specificity_fraction(["SPP1", "TREM2", "MISSING"], top, "Myeloid cells")

    assert result["n_present"] == 2
    assert result["n_expected"] == 2
    assert result["fraction"] == 1.0
    assert result["missing_genes"] == "MISSING"


def test_rank_candidate_axes_orders_axes_by_multicohort_evidence():
    axes = [
        {
            "id": "spp1_trem2",
            "label": "SPP1-TREM2 macrophage program",
            "source_genes": ["SPP1", "TREM2"],
            "source_expected_celltype": "Myeloid cells",
            "target_genes": ["EPCAM"],
            "target_expected_celltype": "Epithelial cells",
            "bulk_genes": ["SPP1", "TREM2"],
            "spatial_target": "spp1_macrophage",
            "snrna_context": "macrophage",
            "snrna_score": "spp1_macrophage",
            "scrna_context": "macrophage",
            "scrna_score": "spp1_macrophage",
            "perturbation_genes": ["SPP1", "TREM2"],
        },
        {
            "id": "weak_axis",
            "label": "Weak axis",
            "source_genes": ["WFDC2"],
            "source_expected_celltype": "Myeloid cells",
            "target_genes": ["LYZ"],
            "target_expected_celltype": "Epithelial cells",
            "bulk_genes": ["WFDC2"],
            "spatial_target": "resident_macrophage",
            "snrna_context": "macrophage",
            "snrna_score": "resident_macrophage",
            "scrna_context": "macrophage",
            "scrna_score": "resident_macrophage",
            "perturbation_genes": ["WFDC2"],
        },
    ]
    gene_top = pd.DataFrame(
        {
            "gene": ["SPP1", "TREM2", "EPCAM", "WFDC2", "LYZ"],
            "Cell_type.refined": ["Myeloid cells", "Myeloid cells", "Epithelial cells", "Epithelial cells", "Myeloid cells"],
        }
    )
    bulk = pd.DataFrame(
        {
            "gene": ["SPP1", "TREM2", "WFDC2"],
            "delta_iac_vs_normal": [50.0, 5.0, -10.0],
        }
    )
    spatial = pd.DataFrame(
        {
            "target": ["spp1_macrophage", "spp1_macrophage", "resident_macrophage"],
            "stage": ["MIA", "LUAD", "MIA"],
            "enrichment_delta_mean": [0.08, 0.09, -0.02],
        }
    )
    snrna = pd.DataFrame(
        {
            "stage": ["Normal", "LUAD", "Normal", "LUAD"],
            "context": ["macrophage", "macrophage", "macrophage", "macrophage"],
            "score": ["spp1_macrophage", "spp1_macrophage", "resident_macrophage", "resident_macrophage"],
            "mean_score": [0.1, 0.3, 0.2, 0.1],
        }
    )
    scrna = pd.DataFrame(
        {
            "stage": ["Adjacent", "Tumor", "Adjacent", "Tumor"],
            "context": ["macrophage", "macrophage", "macrophage", "macrophage"],
            "score": ["spp1_macrophage", "spp1_macrophage", "resident_macrophage", "resident_macrophage"],
            "mean_score": [0.1, 0.2, 0.2, 0.1],
        }
    )

    ranked = rank_candidate_axes(axes, gene_top, bulk, spatial, snrna, scrna)

    assert ranked.iloc[0]["axis_id"] == "spp1_trem2"
    assert ranked.iloc[0]["priority_score"] > ranked.iloc[1]["priority_score"]
    assert ranked.iloc[0]["source_specificity_fraction"] == 1.0
    assert ranked.iloc[0]["target_specificity_fraction"] == 1.0


def test_rank_candidate_axes_can_use_direct_axis_spatial_evidence():
    axes = [
        {
            "id": "spp1_trem2",
            "source_genes": ["SPP1"],
            "source_expected_celltype": "Myeloid cells",
            "target_genes": ["EPCAM"],
            "target_expected_celltype": "Epithelial cells",
            "bulk_genes": ["SPP1"],
        }
    ]
    gene_top = pd.DataFrame(
        {
            "gene": ["SPP1", "EPCAM"],
            "Cell_type.refined": ["Myeloid cells", "Epithelial cells"],
        }
    )
    bulk = pd.DataFrame({"gene": ["SPP1"], "delta_iac_vs_normal": [0.0]})
    spatial = pd.DataFrame({"target": [], "stage": [], "enrichment_delta_mean": []})
    axis_spatial = pd.DataFrame(
        {
            "axis_id": ["spp1_trem2", "spp1_trem2", "spp1_trem2"],
            "evidence_type": [
                "source_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
            ],
            "stage": ["MIA", "LUAD", "AIS"],
            "enrichment_delta_mean": [0.1, 0.2, 0.5],
        }
    )
    empty_scores = pd.DataFrame({"stage": [], "context": [], "score": [], "mean_score": []})

    ranked = rank_candidate_axes(axes, gene_top, bulk, spatial, empty_scores, empty_scores, axis_spatial)

    assert mean_axis_spatial_delta(axis_spatial, "spp1_trem2") == 0.15000000000000002
    assert ranked.iloc[0]["axis_spatial_mia_luad_delta"] == 0.15000000000000002
    assert ranked.iloc[0]["axis_spatial_score"] == 0.7500000000000001


def test_summarize_perturbation_candidates_inherits_axis_priority():
    axes = [
        {"id": "axis_a", "perturbation_genes": ["SPP1", "TREM2"]},
        {"id": "axis_b", "perturbation_genes": ["SPP1"]},
    ]
    ranked = pd.DataFrame(
        {
            "axis_id": ["axis_a", "axis_b"],
            "priority_score": [0.8, 0.4],
        }
    )

    perturbations = summarize_perturbation_candidates(axes, ranked)

    spp1 = perturbations[perturbations["gene"] == "SPP1"].iloc[0]
    assert spp1["n_axes"] == 2
    assert spp1["max_axis_priority"] == 0.8
    assert spp1["mean_axis_priority"] == 0.6
