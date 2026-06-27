import pandas as pd

from luad_niche.cell_states import assign_dominant_labels, add_within_group_high_flags


def test_assign_dominant_labels_uses_highest_score_when_margin_is_clear():
    scores = pd.DataFrame(
        {
            "epithelial_score": [1.2, 0.1],
            "macrophage_score": [0.2, 1.1],
            "t_cell_score": [0.1, 0.3],
        },
        index=["cell1", "cell2"],
    )

    labels = assign_dominant_labels(
        scores,
        {
            "epithelial": "epithelial_score",
            "macrophage": "macrophage_score",
            "t_cell": "t_cell_score",
        },
        min_score=0.05,
        min_margin=0.05,
    )

    assert labels.to_dict() == {"cell1": "epithelial", "cell2": "macrophage"}


def test_assign_dominant_labels_marks_low_or_tied_cells():
    scores = pd.DataFrame(
        {
            "epithelial_score": [0.01, 0.50],
            "macrophage_score": [0.02, 0.48],
        },
        index=["low", "tied"],
    )

    labels = assign_dominant_labels(
        scores,
        {"epithelial": "epithelial_score", "macrophage": "macrophage_score"},
        min_score=0.05,
        min_margin=0.05,
    )

    assert labels.to_dict() == {"low": "unassigned", "tied": "ambiguous"}


def test_add_within_group_high_flags_marks_top_quantile_inside_group_only():
    df = pd.DataFrame(
        {
            "broad_class": ["epithelial", "epithelial", "macrophage", "epithelial"],
            "progenitor_score": [1.0, 2.0, 10.0, 3.0],
        },
        index=["a", "b", "c", "d"],
    )

    flagged = add_within_group_high_flags(
        df,
        group_column="broad_class",
        score_column="progenitor_score",
        group_value="epithelial",
        quantile=0.5,
        output_column="progenitor_high",
    )

    assert flagged["progenitor_high"].to_dict() == {"a": False, "b": False, "c": False, "d": True}
