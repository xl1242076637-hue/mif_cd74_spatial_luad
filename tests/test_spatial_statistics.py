import math

import pandas as pd

from luad_niche.spatial_statistics import (
    add_patient_phase,
    benjamini_hochberg,
    build_paired_patient_differences,
    extract_patient_id,
    summarize_late_vs_precursor,
    summarize_paired_patient_differences,
)


def _adjacency() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "axis_id": ["axis"] * 7,
            "axis_label": ["Axis"] * 7,
            "evidence_type": ["source"] * 7,
            "sample_name": ["P1_AAH", "P1_LUAD", "P2_AIS", "P2_LUAD", "P3_AIS1", "P3_AIS2", "P3_LUAD"],
            "stage": ["AAH", "LUAD", "AIS", "LUAD", "AIS", "AIS", "LUAD"],
            "status": ["ok"] * 7,
            "enrichment_delta": [0.1, 0.4, 0.2, 0.5, 0.1, 0.3, 0.7],
        }
    )


def test_extract_patient_id_and_add_phase():
    assert extract_patient_id("P21_AIS2") == "P21"
    assert extract_patient_id("unknown") is None

    annotated = add_patient_phase(_adjacency())

    assert annotated.loc[0, "patient_id"] == "P1"
    assert annotated.loc[0, "phase"] == "precursor"
    assert annotated.loc[1, "phase"] == "late"


def test_benjamini_hochberg_preserves_nan_and_adjusts_order():
    adjusted = benjamini_hochberg([0.01, 0.04, math.nan, 0.03])

    assert math.isclose(adjusted[0], 0.03)
    assert math.isclose(adjusted[1], 0.04)
    assert math.isnan(adjusted[2])
    assert math.isclose(adjusted[3], 0.04)


def test_patient_aggregation_averages_duplicate_lesions():
    sample_summary = summarize_late_vs_precursor(_adjacency(), aggregate_by_patient=False, bootstrap_iterations=50)
    patient_summary = summarize_late_vs_precursor(_adjacency(), aggregate_by_patient=True, bootstrap_iterations=50)

    assert sample_summary.iloc[0]["n_precursor"] == 4
    assert patient_summary.iloc[0]["n_precursor"] == 3
    assert math.isclose(patient_summary.iloc[0]["precursor_mean"], (0.1 + 0.2 + 0.2) / 3)


def test_paired_patient_summary_uses_late_minus_precursor_difference():
    paired = build_paired_patient_differences(_adjacency())
    summary = summarize_paired_patient_differences(paired, bootstrap_iterations=50)

    assert len(paired) == 3
    assert set(paired["late_minus_precursor"].round(3)) == {0.3, 0.5}
    assert summary.iloc[0]["n_paired_patients"] == 3
    assert math.isclose(summary.iloc[0]["paired_difference_mean"], (0.3 + 0.3 + 0.5) / 3)

