import pandas as pd

from luad_niche.nature_figure_style import (
    EXPORT_SUFFIXES,
    NATURE_PALETTE,
    compact_panel_title,
    semantic_axis_color,
)
from luad_niche.figure_data import (
    prepare_axis_stage_heatmap,
    prepare_axis_priority,
    prepare_candidate_gene_specificity,
    prepare_evidence_heatmap,
    prepare_grn_target_ranking_source,
    prepare_grn_top_gene_source,
    prepare_perturbation_dose_response,
    prepare_perturbation_method_concordance,
    prepare_perturbation_stage_loss,
    prepare_signature_refinement_source,
    prepare_signature_refinement_status_summary,
    prepare_signature_celltype_heatmap,
    prepare_specificity_status_summary,
    prepare_top_perturbation_effects,
    stage_columns,
    summarize_dataset_composition,
)


def _is_hex_color(value: str) -> bool:
    return isinstance(value, str) and len(value) == 7 and value.startswith("#")


def test_nature_redesign_visual_contract_defines_editable_exports_and_semantic_colors():
    assert EXPORT_SUFFIXES == (".svg", ".pdf", ".tiff", ".png")
    assert _is_hex_color(NATURE_PALETTE["mif_axis"])
    assert _is_hex_color(NATURE_PALETTE["macrophage_axis"])
    assert NATURE_PALETTE["mif_axis"] != NATURE_PALETTE["macrophage_axis"]
    assert semantic_axis_color("mif_cd74_cxcr4") == NATURE_PALETTE["mif_axis"]
    assert semantic_axis_color("spp1_trem2_macrophage_epithelial") == NATURE_PALETTE["macrophage_axis"]
    assert semantic_axis_color("unknown_axis") == NATURE_PALETTE["neutral_mid"]


def test_compact_panel_title_keeps_panel_titles_short_without_losing_gene_symbols():
    title = compact_panel_title("Continuous score-level target-prioritization effects for MIF-CD74/CXCR4")

    assert title == "Continuous score-level target-prioritization"


def test_summarize_dataset_composition_counts_stages_and_progression_rows():
    metadata = pd.DataFrame(
        {
            "series_accession": ["GSE307534", "GSE307534", "GSE164789"],
            "sample_accession": ["A", "B", "C"],
            "interpreted_stage": ["Normal", "LUAD", "Tumor"],
            "include_in_luad_progression": [True, True, False],
        }
    )

    result = summarize_dataset_composition(metadata)

    spatial = result[result["series_accession"].eq("GSE307534")].iloc[0]
    assert spatial["Normal"] == 1
    assert spatial["LUAD"] == 1
    assert spatial["n_samples"] == 2
    assert spatial["n_progression_samples"] == 2
    assert spatial["modality"] == "Visium spatial"
    assert stage_columns(result)[:2] == ["Normal", "LUAD"]


def test_prepare_axis_priority_adds_short_labels_and_orders_scores():
    evidence = pd.DataFrame(
        {
            "axis_id": ["spp1_trem2_macrophage_epithelial", "mif_cd74_cxcr4"],
            "evidence_grade": ["strong", "lead"],
            "priority_score": [0.5, 0.9],
        }
    )

    result = prepare_axis_priority(evidence)

    assert result.iloc[0]["axis_id"] == "mif_cd74_cxcr4"
    assert result.iloc[0]["axis_short_label"] == "MIF-CD74/CXCR4"
    assert result.iloc[0]["display_rank"] == 1


def test_prepare_evidence_heatmap_scales_components_to_unit_interval():
    evidence = pd.DataFrame(
        {
            "axis_id": ["mif_cd74_cxcr4"],
            "evidence_grade": ["lead"],
            "priority_score": [0.8],
            "source_spatial_delta": [0.2],
            "target_spatial_delta": [0.1],
            "source_specificity_fraction": [1.0],
            "target_specificity_fraction": [0.5],
            "bulk_delta_mean": [80.0],
            "snrna_program_delta": [-0.1],
            "scrna_tumor_adjacent_delta": [0.05],
            "top_continuous_priority_score": [0.7],
        }
    )

    result = prepare_evidence_heatmap(evidence)

    row = result.iloc[0]
    assert row["Spatial niche"] == 0.75
    assert row["Specificity audit"] == 0.75
    assert row["Bulk trend"] == 1.0
    assert row["snRNA program"] == 0.0
    assert row["Tumor-adjacent scRNA"] == 0.5
    assert row["Target prioritization"] == 0.7


def test_prepare_top_perturbation_effects_returns_positive_coupling_losses():
    perturbation = pd.DataFrame(
        {
            "axis_id": ["mif_cd74_cxcr4", "mif_cd74_cxcr4", "other_axis"],
            "perturbation_type": ["gene", "gene", "axis"],
            "perturbed_genes": ["CD74", "CXCR4", "MIF,CD74"],
            "perturbation_factor": [0.0, 0.5, 0.0],
            "evidence_type": [
                "target_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
            ],
            "n_samples": [30, 30, 30],
            "coupling_relative_delta_mean": [-0.74, -0.03, -1.0],
            "continuous_priority_score": [0.74, 0.03, 1.0],
        }
    )

    result = prepare_top_perturbation_effects(perturbation, top_n=1)

    assert len(result) == 1
    assert result.iloc[0]["perturbed_genes"] == "CD74"
    assert result.iloc[0]["coupling_loss"] == 0.74
    assert result.iloc[0]["evidence_short"] == "receptor-side"


def test_prepare_signature_celltype_heatmap_normalizes_each_signature():
    summary = pd.DataFrame(
        {
            "Cell_type.refined": ["Epithelial cells", "Myeloid cells", "Epithelial cells", "Myeloid cells"],
            "score": ["epi", "epi", "macro", "macro"],
            "mean_score": [2.0, 1.0, 0.5, 1.0],
        }
    )

    result = prepare_signature_celltype_heatmap(
        summary,
        signatures=["epi", "macro"],
        celltypes=["Epithelial cells", "Myeloid cells"],
    )

    epi = result[result["signature"].eq("epi")].iloc[0]
    assert epi["Epithelial cells"] == 1.0
    assert epi["Myeloid cells"] == 0.5

    macro = result[result["signature"].eq("macro")].iloc[0]
    assert macro["Epithelial cells"] == 0.5
    assert macro["Myeloid cells"] == 1.0


def test_prepare_specificity_status_summary_pivots_counts_and_fractions():
    summary = pd.DataFrame(
        {
            "signature": ["spp1_macrophage_vs_other_macrophage"] * 3,
            "specificity_status": ["expected", "off_target", "missing"],
            "n_genes": [14, 15, 1],
        }
    )

    result = prepare_specificity_status_summary(summary)

    row = result.iloc[0]
    assert row["signature_label"] == "SPP1 macrophage"
    assert row["total_genes"] == 30
    assert row["expected_fraction"] == 14 / 30
    assert row["off_target_fraction"] == 15 / 30


def test_prepare_signature_refinement_source_marks_retained_expected_genes():
    audit = pd.DataFrame(
        {
            "signature": ["spp1_macrophage_vs_other_macrophage"] * 4,
            "rank": [3, 1, 4, 2],
            "gene": ["C", "A", "D", "B"],
            "expected_celltype": ["Myeloid cells"] * 4,
            "top_celltype": ["Fibroblasts", "Myeloid cells", "", "Epithelial cells"],
            "specificity_status": ["off_target", "expected", "missing", "off_target"],
        }
    )

    result = prepare_signature_refinement_source(audit, filtered_genes=["A"])

    assert result["gene"].tolist() == ["A", "B", "C", "D"]
    assert result["retained_after_audit"].tolist() == [True, False, False, False]
    assert result["retained_label"].tolist() == ["Retained", "Removed", "Removed", "Removed"]
    assert result.loc[result["gene"].eq("D"), "top_celltype"].iloc[0] == "Missing in reference"


def test_prepare_signature_refinement_status_summary_counts_original_and_retained_rows():
    source = pd.DataFrame(
        {
            "specificity_status": ["expected", "expected", "off_target", "missing"],
            "retained_after_audit": [True, False, False, False],
        }
    )

    result = prepare_signature_refinement_status_summary(source)

    assert result.loc[result["specificity_status"].eq("expected"), "original_n"].iloc[0] == 2
    assert result.loc[result["specificity_status"].eq("expected"), "retained_n"].iloc[0] == 1
    assert result.loc[result["specificity_status"].eq("off_target"), "retained_n"].iloc[0] == 0


def test_prepare_candidate_gene_specificity_keeps_requested_gene_order():
    gene_top = pd.DataFrame(
        {
            "gene": ["CD74", "MIF", "SPP1"],
            "Cell_type.refined": ["Myeloid cells", "Epithelial cells", "Myeloid cells"],
            "mean_expression": [3.0, 1.0, 2.0],
        }
    )

    result = prepare_candidate_gene_specificity(gene_top, genes=["MIF", "CD74", "SPP1"])

    assert result["gene"].astype(str).tolist() == ["MIF", "CD74", "SPP1"]
    assert result.iloc[0]["Cell_type.refined"] == "Epithelial cells"


def test_prepare_axis_stage_heatmap_creates_ordered_axis_evidence_rows():
    spatial = pd.DataFrame(
        {
            "axis_id": ["mif_cd74_cxcr4", "mif_cd74_cxcr4", "spp1_trem2_macrophage_epithelial"],
            "evidence_type": [
                "source_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
            ],
            "stage": ["MIA", "LUAD", "MIA"],
            "n_samples": [4, 26, 4],
            "enrichment_delta_mean": [0.2, 0.1, 0.05],
        }
    )

    result = prepare_axis_stage_heatmap(spatial, stages=["MIA", "LUAD"])

    assert result.iloc[0]["row_label"] == "MIF-CD74/CXCR4 (source-side)"
    assert result.iloc[0]["MIA"] == 0.2
    assert result.iloc[1]["row_label"] == "MIF-CD74/CXCR4 (receptor-side)"
    assert result.iloc[1]["LUAD"] == 0.1


def _continuous_ranking_fixture():
    return pd.DataFrame(
        {
            "perturbation_id": ["mif_x0", "mif_x05", "cd74_x0", "cd74_x05", "trem2_source_x0", "trem2_target_x0"],
            "perturbation_type": ["gene"] * 6,
            "perturbed_genes": ["MIF", "MIF", "CD74", "CD74", "TREM2", "TREM2"],
            "perturbation_factor": [0.0, 0.5, 0.0, 0.5, 0.0, 0.0],
            "axis_id": ["mif_axis", "mif_axis", "mif_axis", "mif_axis", "spp1_axis", "c1q_axis"],
            "evidence_type": [
                "source_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "source_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
            ],
            "n_samples": [30] * 6,
            "coupling_relative_delta_mean": [-1.0, -0.5, -0.7, -0.35, -0.02, -0.06],
            "continuous_priority_score": [1.0, 0.5, 0.7, 0.35, 0.02, 0.06],
        }
    )


def test_prepare_perturbation_dose_response_adds_baseline_and_best_duplicate_rows():
    result = prepare_perturbation_dose_response(_continuous_ranking_fixture(), genes=["MIF", "CD74", "TREM2"])

    mif = result[result["perturbed_genes"].eq("MIF")]
    assert mif["perturbation_factor"].tolist() == [1.0, 0.5, 0.0]
    assert mif["coupling_remaining"].tolist() == [1.0, 0.5, 0.0]

    trem2_full = result[result["perturbed_genes"].eq("TREM2") & result["perturbation_factor"].eq(0.0)].iloc[0]
    assert trem2_full["evidence_short"] == "receptor-side"
    assert trem2_full["coupling_loss"] == 0.06


def test_prepare_perturbation_stage_loss_uses_best_full_dropout_id_by_gene():
    effects = pd.DataFrame(
        {
            "perturbation_id": ["mif_x0", "mif_x0", "cd74_x0", "trem2_target_x0"],
            "stage": ["MIA", "LUAD", "MIA", "MIA"],
            "sample": ["A", "B", "A", "A"],
            "coupling_relative_delta": [-1.0, -0.8, -0.4, -0.1],
        }
    )

    result = prepare_perturbation_stage_loss(
        effects,
        _continuous_ranking_fixture(),
        genes=["MIF", "CD74", "TREM2"],
        stages=["MIA", "LUAD"],
    )

    mif_mia = result[result["row_label"].eq("MIF (source-side)") & result["stage"].astype(str).eq("MIA")].iloc[0]
    assert mif_mia["coupling_loss"] == 1.0

    trem2 = result[result["perturbed_genes"].eq("TREM2")].iloc[0]
    assert trem2["row_label"] == "TREM2 (receptor-side)"


def test_prepare_perturbation_method_concordance_merges_continuous_and_top_quantile_outputs():
    virtual = pd.DataFrame(
        {
            "perturbation_id": ["mif_x0", "cd74_x0", "trem2_target_x0"],
            "axis_id": ["mif_axis", "mif_axis", "c1q_axis"],
            "evidence_type": [
                "source_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
                "target_near_epithelial_progenitor",
            ],
            "panel_relative_delta": [-1.0, -0.6, -0.05],
            "observed_fraction_delta": [-0.5, 0.0, 0.01],
            "dropout_priority_score": [1.5, 0.6, 0.05],
        }
    )

    result = prepare_perturbation_method_concordance(
        _continuous_ranking_fixture(),
        virtual,
        genes=["MIF", "CD74", "TREM2"],
    )

    mif = result[result["perturbed_genes"].eq("MIF")].iloc[0]
    assert mif["continuous_coupling_loss"] == 1.0
    assert mif["top_quantile_panel_loss"] == 1.0
    assert mif["dropout_priority_score"] == 1.5


def _grn_ranking_fixture():
    return pd.DataFrame(
        {
            "broad_class": ["macrophage", "macrophage", "epithelial"],
            "target_gene": ["CD74", "CXCR4", "MIF"],
            "top_impacted_signature": [
                "c1q_macrophage_vs_other_macrophage",
                "spp1_macrophage_vs_other_macrophage",
                "proliferating_epithelial_vs_other_epithelial",
            ],
            "top_signature_mean_impact": [0.015, 0.012, 0.013],
            "top_signature_max_impact": [0.04, 0.03, 0.02],
            "top_signature_sum_impact": [0.30, 0.25, 0.20],
            "top_impacted_genes": [
                "C1QA,C1QB,C1QC,APOE,ACP5",
                "PRDM1,GPR183,ADAMDEC1,PLA2G7,APOE",
                "MDK,EPCAM,KRT8,CD24,PERP",
            ],
        }
    )


def test_prepare_grn_target_ranking_source_orders_targets_and_shortens_signatures():
    result = prepare_grn_target_ranking_source(_grn_ranking_fixture())

    assert result["target_gene"].tolist() == ["CD74", "MIF", "CXCR4"]
    assert result["display_rank"].tolist() == [1, 2, 3]
    assert result.loc[result["target_gene"].eq("CD74"), "signature_label"].iloc[0] == "C1Q macrophage"
    assert result.loc[result["target_gene"].eq("MIF"), "class_label"].iloc[0] == "Epithelial"


def test_prepare_grn_top_gene_source_splits_top_impacted_genes():
    ranking = prepare_grn_target_ranking_source(_grn_ranking_fixture())

    result = prepare_grn_top_gene_source(ranking, top_n=3)

    cd74 = result[result["target_gene"].eq("CD74")]
    assert cd74["impacted_gene"].tolist() == ["C1QA", "C1QB", "C1QC"]
    assert cd74["gene_position"].tolist() == [1, 2, 3]
    assert result.loc[result["target_gene"].eq("MIF"), "signature_label"].iloc[0] == "Proliferating epi"
