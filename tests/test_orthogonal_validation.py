import pandas as pd

from luad_niche.orthogonal_validation import (
    FOCUSED_GENE_CONTEXTS,
    prepare_bulk_focused_stage_source,
    prepare_snrna_focused_stage_source,
    summarize_compartment_gene_expression,
    summarize_stage_expression,
)


def test_summarize_compartment_gene_expression_reports_mean_and_detection_fraction():
    table = pd.DataFrame(
        {
            "broad_class": ["epithelial", "epithelial", "macrophage", "macrophage"],
            "MIF": [2.0, 0.0, 1.0, 0.0],
            "CD74": [0.0, 0.0, 3.0, 1.0],
        },
        index=["cell1", "cell2", "cell3", "cell4"],
    )

    result = summarize_compartment_gene_expression(
        table,
        sample="GSM1",
        sample_name="P1_Normal",
        stage="Normal",
        genes=["MIF", "CD74"],
    )

    mif_epithelial = result[
        result["gene"].eq("MIF") & result["context"].eq("epithelial")
    ].iloc[0]
    cd74_macrophage = result[
        result["gene"].eq("CD74") & result["context"].eq("macrophage")
    ].iloc[0]
    assert mif_epithelial["n_cells"] == 2
    assert mif_epithelial["mean_expression"] == 1.0
    assert mif_epithelial["detection_fraction"] == 0.5
    assert cd74_macrophage["mean_expression"] == 2.0
    assert cd74_macrophage["detection_fraction"] == 1.0


def test_summarize_stage_expression_averages_sample_level_values():
    sample_summary = pd.DataFrame(
        {
            "sample": ["A", "B"],
            "stage": ["Normal", "Normal"],
            "context": ["epithelial", "epithelial"],
            "gene": ["MIF", "MIF"],
            "n_cells": [10, 20],
            "mean_expression": [1.0, 3.0],
            "detection_fraction": [0.1, 0.3],
        }
    )

    result = summarize_stage_expression(sample_summary)

    row = result.iloc[0]
    assert row["n_samples"] == 2
    assert row["total_cells"] == 30
    assert row["mean_expression"] == 2.0
    assert row["mean_detection_fraction"] == 0.2


def test_prepare_snrna_focused_stage_source_keeps_gene_specific_contexts_and_zscores():
    stage_summary = pd.DataFrame(
        {
            "stage": ["Normal", "LUAD", "Normal", "LUAD", "Normal", "LUAD"],
            "context": [
                "epithelial",
                "epithelial",
                "macrophage",
                "macrophage",
                "all_cells",
                "all_cells",
            ],
            "gene": ["MIF", "MIF", "CD74", "CD74", "MIF", "MIF"],
            "n_samples": [3, 4, 3, 4, 3, 4],
            "mean_expression": [1.0, 3.0, 2.0, 6.0, 9.0, 10.0],
            "sd_expression": [0.1, 0.2, 0.1, 0.3, 0.2, 0.2],
            "mean_detection_fraction": [0.1, 0.3, 0.2, 0.6, 0.8, 0.9],
            "total_cells": [20, 30, 20, 30, 50, 60],
        }
    )

    result = prepare_snrna_focused_stage_source(
        stage_summary,
        gene_contexts={"MIF": "epithelial", "CD74": "macrophage"},
        stages=["Normal", "LUAD"],
    )

    assert set(result["context"]) == {"epithelial", "macrophage"}
    assert result[result["gene"].eq("MIF")]["mean_expression"].tolist() == [1.0, 3.0]
    assert result[result["gene"].eq("CD74")]["mean_expression"].tolist() == [2.0, 6.0]
    assert result.groupby("gene")["row_zscore"].mean().round(10).tolist() == [0.0, 0.0]


def test_prepare_bulk_focused_stage_source_filters_genes_and_orders_stages():
    marker_table = pd.DataFrame(
        {
            "gene": ["MIF", "MIF", "CD74", "CD74", "OTHER", "OTHER"],
            "stage": ["IAC", "Normal", "IAC", "Normal", "IAC", "Normal"],
            "mean_expression": [4.0, 2.0, 9.0, 3.0, 10.0, 2.0],
            "status": ["present"] * 6,
        }
    )

    result = prepare_bulk_focused_stage_source(
        marker_table,
        genes=["MIF", "CD74"],
        stages=["Normal", "IAC"],
    )

    assert result["gene"].drop_duplicates().tolist() == ["MIF", "CD74"]
    assert result[result["gene"].eq("MIF")]["stage"].tolist() == ["Normal", "IAC"]
    assert result[result["gene"].eq("CD74")]["delta_late_vs_normal"].unique().tolist() == [6.0]


def test_default_focused_gene_contexts_keep_mif_source_and_cd74_receptor_roles():
    assert FOCUSED_GENE_CONTEXTS["MIF"] == "epithelial"
    assert FOCUSED_GENE_CONTEXTS["CD74"] == "macrophage"
