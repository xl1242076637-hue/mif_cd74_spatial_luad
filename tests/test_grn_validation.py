import pandas as pd

from luad_niche.grn_validation import (
    canonical_signature_state,
    expected_celltype_for_state,
    summarize_group_delta,
)


def test_canonical_signature_state_maps_grn_signature_names():
    assert canonical_signature_state("c1q_macrophage_vs_other_macrophage") == "c1q_macrophage"
    assert canonical_signature_state("spp1_macrophage_vs_other_macrophage") == "spp1_macrophage"
    assert canonical_signature_state("proliferating_epithelial_vs_other_epithelial") == "proliferating_epithelial"


def test_expected_celltype_for_state_uses_broad_compartment():
    assert expected_celltype_for_state("c1q_macrophage") == "Myeloid cells"
    assert expected_celltype_for_state("proliferating_epithelial") == "Epithelial cells"
    assert expected_celltype_for_state("unknown_state") == ""


def test_summarize_group_delta_reports_comparison_minus_baseline_mean():
    table = pd.DataFrame(
        {
            "stage": ["AIS", "MIA", "IAC", "IAC"],
            "state": ["spp1_macrophage", "spp1_macrophage", "spp1_macrophage", "c1q_macrophage"],
            "mean_fraction": [0.2, 0.4, 0.5, 0.9],
            "n_samples": [2, 2, 3, 3],
        }
    )

    result = summarize_group_delta(
        table,
        state_column="state",
        state_value="spp1_macrophage",
        group_column="stage",
        value_column="mean_fraction",
        baseline_groups=["AIS"],
        comparison_groups=["MIA", "IAC"],
    )

    assert result["baseline_mean"] == 0.2
    assert result["comparison_mean"] == 0.45
    assert result["delta"] == 0.25
    assert result["baseline_n_samples"] == 2
    assert result["comparison_n_samples"] == 5
