import math
import importlib.util
from pathlib import Path

import pandas as pd

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "score_gse307534_refined_signatures.py"
SPEC = importlib.util.spec_from_file_location("score_gse307534_refined_signatures", SCRIPT_PATH)
assert SPEC is not None
scoring = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(scoring)

adjacency_for_target = scoring.adjacency_for_target
sample_accessions_from_filelist = scoring.sample_accessions_from_filelist
summarize_adjacency_by_target_stage = (
    scoring.summarize_adjacency_by_target_stage
)


def test_sample_accessions_from_filelist_adds_aliases_for_duplicate_lesions(tmp_path):
    filelist = tmp_path / "filelist.txt"
    filelist.write_text(
        "\n".join(
            [
                "#Archive/File\tName\tTime\tSize\tType",
                "File\tGSM9226175_P4_AAH.tar.gz\t08/29/2025\t1\tTAR",
                "File\tGSM9226176_P4_AAH-1.tar.gz\t08/29/2025\t1\tTAR",
                "File\tGSM9226211_P21_AIS.tar.gz\t08/29/2025\t1\tTAR",
                "File\tGSM9226212_P21_AIS-1.tar.gz\t08/29/2025\t1\tTAR",
                "File\tGSM9226189_P10_MIA.tar.gz\t08/29/2025\t1\tTAR",
            ]
        ),
        encoding="utf-8",
    )

    mapping = sample_accessions_from_filelist(filelist)

    assert mapping["P4_AAH"] == "GSM9226175"
    assert mapping["P4_AAH1"] == "GSM9226175"
    assert mapping["P4_AAH-1"] == "GSM9226176"
    assert mapping["P4_AAH2"] == "GSM9226176"
    assert mapping["P21_AIS1"] == "GSM9226211"
    assert mapping["P21_AIS2"] == "GSM9226212"
    assert mapping["P10_MIA"] == "GSM9226189"
    assert "P10_MIA1" not in mapping


def test_adjacency_for_target_marks_empty_high_spot_sets_invalid():
    table = pd.DataFrame(
        {
            "sample": ["GSMTEST"] * 4,
            "sample_name": ["P1_AAH"] * 4,
            "stage": ["AAH"] * 4,
            "x": [0.0, 0.0, 1.0, 1.0],
            "y": [0.0, 1.0, 0.0, 1.0],
            "epithelial_progenitor_like_score": [0.0, 0.0, 0.0, 0.0],
            "spp1_macrophage_score": [0.0, 1.0, 2.0, 3.0],
        }
    )

    stats = adjacency_for_target(
        table,
        source_column="epithelial_progenitor_like_score",
        target_column="spp1_macrophage_score",
        quantile=0.75,
        radius_multiplier=1.0,
        permutations=5,
        seed=1,
    )

    assert stats["status"] == "insufficient_high_spots"
    assert stats["n_source_high"] == 0
    assert stats["n_target_high"] == 1
    assert math.isnan(stats["observed_fraction"])
    assert math.isnan(stats["null_mean"])
    assert math.isnan(stats["enrichment_delta"])
    assert math.isnan(stats["empirical_p_greater"])


def test_summarize_adjacency_by_target_stage_counts_valid_samples_only_for_means():
    adjacency = pd.DataFrame(
        {
            "target": ["spp1_macrophage", "spp1_macrophage"],
            "stage": ["AIS", "AIS"],
            "sample": ["s1", "s2"],
            "status": ["ok", "insufficient_high_spots"],
            "observed_fraction": [0.6, math.nan],
            "null_mean": [0.4, math.nan],
            "enrichment_delta": [0.2, math.nan],
            "empirical_p_greater": [0.05, math.nan],
            "empirical_p_less": [0.95, math.nan],
        }
    )

    summary = summarize_adjacency_by_target_stage(adjacency)

    row = summary.iloc[0]
    assert row["n_samples"] == 2
    assert row["n_valid_samples"] == 1
    assert row["n_invalid_tests"] == 1
    assert row["observed_fraction_mean"] == 0.6
    assert row["enrichment_delta_mean"] == 0.2
