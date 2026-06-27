import numpy as np
import pandas as pd

from luad_niche.spatial_niche import (
    finite_coordinate_mask,
    median_nearest_neighbor_distance,
    nearest_target_fraction,
    permutation_adjacency_test,
    summarize_adjacency_by_stage,
    top_quantile_mask,
)


def test_top_quantile_mask_selects_high_values():
    values = pd.Series([1.0, 2.0, 3.0, 4.0], index=["a", "b", "c", "d"])

    mask = top_quantile_mask(values, quantile=0.75)

    assert mask.to_dict() == {"a": False, "b": False, "c": False, "d": True}


def test_nearest_target_fraction_counts_sources_near_targets():
    coords = np.array([[0.0, 0.0], [0.0, 1.0], [10.0, 10.0]])
    source_mask = np.array([True, False, True])
    target_mask = np.array([False, True, False])

    fraction = nearest_target_fraction(coords, source_mask, target_mask, radius=1.5)

    assert fraction == 0.5


def test_median_nearest_neighbor_distance_uses_first_nonself_neighbor():
    coords = np.array([[0.0, 0.0], [0.0, 2.0], [10.0, 10.0]])

    assert median_nearest_neighbor_distance(coords) == 2.0


def test_finite_coordinate_mask_flags_nan_rows():
    coords = np.array([[0.0, 0.0], [np.nan, np.nan], [2.0, 2.0]])

    assert finite_coordinate_mask(coords).tolist() == [True, False, True]


def test_permutation_adjacency_test_reports_both_tail_probabilities_and_delta():
    coords = np.array([[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]])
    source_mask = np.array([True, False, True, False])
    target_mask = np.array([False, True, False, True])

    stats = permutation_adjacency_test(
        coords,
        source_mask,
        target_mask,
        radius=1.5,
        n_permutations=10,
        seed=1,
    )

    assert "empirical_p_less" in stats
    assert "enrichment_delta" in stats
    assert stats["enrichment_delta"] == stats["observed_fraction"] - stats["null_mean"]


def test_summarize_adjacency_by_stage_reports_mean_effect_size():
    adjacency = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3"],
            "stage": ["AIS", "AIS", "IAC"],
            "observed_fraction": [0.5, 0.7, 0.9],
            "null_mean": [0.4, 0.5, 0.3],
            "empirical_p_greater": [0.2, 0.1, 0.01],
            "empirical_p_less": [0.8, 0.9, 0.99],
        }
    )

    summary = summarize_adjacency_by_stage(adjacency)

    ais = summary[summary["stage"] == "AIS"].iloc[0]
    assert ais["n_samples"] == 2
    assert ais["observed_fraction_mean"] == 0.6
    assert ais["null_mean_mean"] == 0.45
    assert ais["enrichment_delta_mean"] == 0.15
    assert ais["empirical_p_less_median"] == 0.85
