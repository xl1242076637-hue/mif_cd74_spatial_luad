import pandas as pd

from luad_niche.signatures import build_top_n_signature_panels


def test_build_top_n_signature_panels_filters_and_names_panels():
    ranked = pd.DataFrame(
        {
            "contrast": ["a_vs_b", "a_vs_b", "a_vs_b", "c_vs_d"],
            "gene": ["G1", "G2", "LOW", "G3"],
            "log2fc": [2.0, 1.0, -1.0, 3.0],
            "pct_group": [0.5, 0.1, 0.9, 0.7],
        }
    )

    panels = build_top_n_signature_panels(
        ranked,
        contrasts=["a_vs_b", "c_vs_d"],
        top_ns=[1, 2],
        min_pct_group=0.05,
    )

    assert panels["a_vs_b_top1"] == ["G1"]
    assert panels["a_vs_b_top2"] == ["G1", "G2"]
    assert panels["c_vs_d_top2"] == ["G3"]
