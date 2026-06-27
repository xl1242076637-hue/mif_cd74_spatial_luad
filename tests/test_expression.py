import pandas as pd

from luad_niche.expression import (
    compute_stage_means,
    expression_by_gene,
    marker_stage_table,
    marker_trend_summary,
    normalize_log1p_counts,
    normalize_log1p_counts_by_totals,
)


def test_expression_by_gene_extracts_prefixed_sample_columns():
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2"],
            "gene_name": ["IL1B", "SRGN"],
            "FPKM.AIS1": [1.0, 2.0],
            "FPKM.IAC1": [3.0, 4.0],
            "count.AIS1": [10, 20],
        }
    )

    expr = expression_by_gene(df, value_prefix="FPKM.")

    assert list(expr.index) == ["IL1B", "SRGN"]
    assert list(expr.columns) == ["AIS1", "IAC1"]
    assert expr.loc["SRGN", "IAC1"] == 4.0


def test_compute_stage_means_uses_sample_metadata():
    expr = pd.DataFrame(
        {
            "S1": [1.0, 3.0],
            "S2": [2.0, 5.0],
            "S3": [7.0, 9.0],
        },
        index=["IL1B", "SRGN"],
    )
    metadata = pd.DataFrame(
        {
            "title": ["S1", "S2", "S3"],
            "interpreted_stage": ["AIS", "AIS", "IAC"],
            "include_in_luad_progression": [True, True, True],
        }
    )

    means = compute_stage_means(expr, metadata)

    assert means.loc["IL1B", "AIS"] == 1.5
    assert means.loc["SRGN", "IAC"] == 9.0


def test_marker_stage_table_reports_missing_markers():
    expr = pd.DataFrame({"S1": [1.0]}, index=["IL1B"])
    metadata = pd.DataFrame(
        {
            "title": ["S1"],
            "interpreted_stage": ["AIS"],
            "include_in_luad_progression": [True],
        }
    )

    table = marker_stage_table(expr, metadata, markers=["IL1B", "SPP1"])

    assert table.loc[table["gene"] == "IL1B", "status"].iloc[0] == "present"
    assert table.loc[table["gene"] == "SPP1", "status"].iloc[0] == "missing"


def test_marker_trend_summary_reports_delta_and_max_stage():
    table = pd.DataFrame(
        {
            "gene": ["IL1B", "IL1B", "SRGN", "SRGN"],
            "stage": ["Normal", "IAC", "Normal", "IAC"],
            "mean_expression": [1.0, 5.0, 8.0, 2.0],
            "status": ["present", "present", "present", "present"],
        }
    )

    trends = marker_trend_summary(table)

    il1b = trends[trends["gene"] == "IL1B"].iloc[0]
    srgn = trends[trends["gene"] == "SRGN"].iloc[0]
    assert il1b["delta_iac_vs_normal"] == 4.0
    assert il1b["max_stage"] == "IAC"
    assert srgn["delta_iac_vs_normal"] == -6.0
    assert srgn["max_stage"] == "Normal"


def test_normalize_log1p_counts_scales_each_spot_and_keeps_zero_rows_zero():
    expr = pd.DataFrame(
        {"EPCAM": [5.0, 50.0, 0.0], "SPP1": [5.0, 50.0, 0.0]},
        index=["spot_low_depth", "spot_high_depth", "spot_zero"],
    )

    normalized = normalize_log1p_counts(expr, scale_factor=10.0)

    assert normalized.loc["spot_low_depth", "EPCAM"] == normalized.loc["spot_high_depth", "EPCAM"]
    assert normalized.loc["spot_zero", "EPCAM"] == 0.0
    assert normalized.loc["spot_zero", "SPP1"] == 0.0


def test_normalize_log1p_counts_by_totals_uses_external_library_sizes():
    expr = pd.DataFrame({"EPCAM": [5.0, 50.0, 0.0]}, index=["low", "high", "zero"])
    totals = pd.Series({"low": 10.0, "high": 100.0, "zero": 0.0})

    normalized = normalize_log1p_counts_by_totals(expr, totals, scale_factor=10.0)

    assert normalized.loc["low", "EPCAM"] == normalized.loc["high", "EPCAM"]
    assert normalized.loc["zero", "EPCAM"] == 0.0
