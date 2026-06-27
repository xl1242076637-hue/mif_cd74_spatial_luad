"""Signature panel helpers."""

from __future__ import annotations

import pandas as pd


def build_top_n_signature_panels(
    ranked_markers: pd.DataFrame,
    contrasts: list[str],
    top_ns: list[int],
    min_pct_group: float = 0.05,
    min_log2fc: float = 0.0,
) -> dict[str, list[str]]:
    """Build top-N gene panels from ranked marker tables."""
    panels: dict[str, list[str]] = {}
    for contrast in contrasts:
        filtered = ranked_markers[
            ranked_markers["contrast"].eq(contrast)
            & (ranked_markers["log2fc"] > min_log2fc)
            & (ranked_markers["pct_group"] >= min_pct_group)
        ]
        genes = filtered["gene"].astype(str).drop_duplicates().tolist()
        for top_n in top_ns:
            panels[f"{contrast}_top{top_n}"] = genes[:top_n]
    return panels
