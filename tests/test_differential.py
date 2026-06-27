import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from luad_niche.differential import aggregate_contrast, rank_markers


def test_aggregate_contrast_sums_counts_and_detection_by_group():
    matrix = csr_matrix(
        np.array(
            [
                [5, 4, 0, 0],
                [0, 1, 3, 4],
            ]
        )
    )
    group_mask = np.array([True, True, False, False])
    background_mask = np.array([False, False, True, True])

    aggregate = aggregate_contrast(matrix, group_mask, background_mask)

    assert aggregate["group_sum"].tolist() == [9, 1]
    assert aggregate["background_sum"].tolist() == [0, 7]
    assert aggregate["group_detected"].tolist() == [2, 1]
    assert aggregate["background_detected"].tolist() == [0, 2]
    assert aggregate.attrs["n_group"] == 2
    assert aggregate.attrs["n_background"] == 2


def test_rank_markers_prioritizes_high_group_expression_and_detection():
    aggregate = pd.DataFrame(
        {
            "gene": ["GENE_A", "GENE_B"],
            "group_sum": [9.0, 1.0],
            "background_sum": [0.0, 7.0],
            "group_detected": [2, 1],
            "background_detected": [0, 2],
            "n_group": [2, 2],
            "n_background": [2, 2],
        }
    )

    ranked = rank_markers(aggregate)

    assert ranked.iloc[0]["gene"] == "GENE_A"
    assert ranked.iloc[0]["log2fc"] > 0
    assert ranked.iloc[0]["pct_group"] == 1.0
