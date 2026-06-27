"""Data preparation helpers for manuscript-facing LUAD niche figures."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd


STAGE_ORDER = (
    "Normal",
    "AAH",
    "AIS",
    "MIA",
    "LUAD",
    "IAC",
    "Adjacent",
    "Tumor",
    "Primary tumor",
    "Metastasis/effusion",
    "Normal lymph node",
    "LUSC",
    "Unknown",
)

SPATIAL_STAGE_ORDER = ("Normal", "AAH", "AIS", "MIA", "LUAD")

CELLTYPE_ORDER = (
    "Epithelial cells",
    "Myeloid cells",
    "Fibroblasts",
    "Endothelial cells",
    "MAST cells",
    "B lymphocytes",
    "T/NK cells",
    "Unassigned",
)

PROGRESSION_STAGES = ("Normal", "AAH", "AIS", "MIA", "LUAD", "IAC")

DATASET_ORDER = (
    "GSE307534",
    "GSE308103",
    "GSE131907",
    "GSE282617",
    "GSE189357",
    "GSE189487",
    "GSE164789",
)

DATASET_ANNOTATIONS = {
    "GSE307534": {
        "modality": "Visium spatial",
        "role": "main spatial progression",
        "analysis_use": "axis localization and perturbation",
        "data_extent": "56 tissue sections",
    },
    "GSE308103": {
        "modality": "snRNA-seq",
        "role": "main single-nucleus progression",
        "analysis_use": "cell-state trajectory support",
        "data_extent": "75 samples",
    },
    "GSE131907": {
        "modality": "scRNA-seq reference",
        "role": "specificity audit",
        "analysis_use": "cell-type specificity filtering",
        "data_extent": "208,506 cells",
    },
    "GSE282617": {
        "modality": "bulk RNA-seq",
        "role": "bulk progression trend",
        "analysis_use": "marker trend support",
        "data_extent": "70 samples",
    },
    "GSE189357": {
        "modality": "scRNA-seq",
        "role": "small discovery context",
        "analysis_use": "initial state refinement",
        "data_extent": "9 samples",
    },
    "GSE189487": {
        "modality": "Visium spatial",
        "role": "small spatial context",
        "analysis_use": "pilot niche robustness",
        "data_extent": "6 sections",
    },
    "GSE164789": {
        "modality": "scRNA-seq",
        "role": "tumor-adjacent validation",
        "analysis_use": "tumor-adjacent immune-state contrast",
        "data_extent": "62 expression matrices",
    },
}

AXIS_SHORT_LABELS = {
    "mif_cd74_cxcr4": "MIF-CD74/CXCR4",
    "spp1_trem2_macrophage_epithelial": "SPP1/TREM2/PLA2G7",
    "cxcl9_cxcl10_cxcr3": "CXCL9/10-CXCR3",
    "c1q_apoe_trem2_lgals3": "C1Q/APOE/TREM2",
    "inflammatory_il1_tnf_cxcl8": "IL1B/TNF/CXCL8",
}

SIGNATURE_SHORT_LABELS = {
    "epithelial_progenitor_like": "Epi progenitor",
    "refined_epithelial_progenitor_like": "Refined epi progenitor",
    "spp1_macrophage": "SPP1 macrophage",
    "refined_spp1_macrophage": "Refined SPP1 macrophage",
    "c1q_macrophage": "C1Q macrophage",
    "refined_c1q_macrophage": "Refined C1Q macrophage",
    "inflammatory_macrophage": "Inflammatory macrophage",
    "refined_inflammatory_macrophage": "Refined inflammatory macrophage",
    "epithelial_progenitor_like_vs_other_epithelial": "Epi progenitor",
    "proliferating_epithelial_vs_other_epithelial": "Proliferating epi",
    "spp1_macrophage_vs_other_macrophage": "SPP1 macrophage",
    "c1q_macrophage_vs_other_macrophage": "C1Q macrophage",
    "inflammatory_macrophage_vs_other_macrophage": "Inflammatory macrophage",
    "resident_macrophage_vs_other_macrophage": "Resident macrophage",
}


def _ordered_existing(values: Iterable[str], preferred: Iterable[str]) -> list[str]:
    observed = [str(value) for value in values]
    preferred_present = [value for value in preferred if value in observed]
    extra = sorted(value for value in observed if value not in set(preferred_present))
    return preferred_present + extra


def stage_columns(table: pd.DataFrame) -> list[str]:
    """Return stage-count columns in manuscript display order."""
    return [stage for stage in STAGE_ORDER if stage in table.columns]


def summarize_dataset_composition(metadata: pd.DataFrame) -> pd.DataFrame:
    """Summarize per-dataset stage composition and analysis role annotations."""
    required = {"series_accession", "sample_accession", "interpreted_stage"}
    missing = required.difference(metadata.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"metadata is missing required columns: {missing_text}")

    working = metadata.copy()
    working["interpreted_stage"] = working["interpreted_stage"].fillna("Unknown").astype(str)

    counts = (
        working.groupby(["series_accession", "interpreted_stage"], dropna=False)
        .size()
        .unstack(fill_value=0)
    )
    ordered_stages = _ordered_existing(counts.columns, STAGE_ORDER)
    counts = counts.reindex(columns=ordered_stages, fill_value=0)
    counts["n_samples"] = counts.sum(axis=1)

    if "include_in_luad_progression" in working.columns:
        progression_counts = (
            working.assign(include_in_luad_progression=working["include_in_luad_progression"].fillna(False).astype(bool))
            .groupby("series_accession")["include_in_luad_progression"]
            .sum()
        )
        counts["n_progression_samples"] = progression_counts.reindex(counts.index, fill_value=0).astype(int)
    else:
        progression_cols = [column for column in PROGRESSION_STAGES if column in counts.columns]
        counts["n_progression_samples"] = counts[progression_cols].sum(axis=1).astype(int)

    annotations = pd.DataFrame.from_dict(DATASET_ANNOTATIONS, orient="index")
    annotations.index.name = "series_accession"
    summary = counts.join(annotations, how="left").reset_index()
    summary["modality"] = summary["modality"].fillna("other")
    summary["role"] = summary["role"].fillna("context")
    summary["analysis_use"] = summary["analysis_use"].fillna("supporting")
    summary["data_extent"] = summary["data_extent"].fillna(summary["n_samples"].astype(str) + " samples")

    dataset_order = {dataset: i for i, dataset in enumerate(DATASET_ORDER)}
    summary["_order"] = summary["series_accession"].map(dataset_order).fillna(len(dataset_order)).astype(int)
    summary = summary.sort_values(["_order", "series_accession"]).drop(columns="_order").reset_index(drop=True)
    return summary


def axis_short_label(axis_id: str) -> str:
    """Return a compact axis label for plotting."""
    return AXIS_SHORT_LABELS.get(axis_id, axis_id.replace("_", " "))


def prepare_axis_priority(evidence_matrix: pd.DataFrame) -> pd.DataFrame:
    """Prepare sorted axis priority rows for plotting."""
    required = {"axis_id", "priority_score", "evidence_grade"}
    missing = required.difference(evidence_matrix.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"evidence matrix is missing required columns: {missing_text}")

    table = evidence_matrix.copy()
    table["priority_score"] = pd.to_numeric(table["priority_score"], errors="coerce").fillna(0.0)
    table["axis_short_label"] = table["axis_id"].map(axis_short_label)
    table = table.sort_values("priority_score", ascending=False).reset_index(drop=True)
    table["display_rank"] = range(1, len(table) + 1)
    return table[["display_rank", "axis_id", "axis_short_label", "evidence_grade", "priority_score"]]


def _positive_scaled_mean(row: pd.Series, columns: list[str], scale: float) -> float:
    values = pd.to_numeric(row.reindex(columns), errors="coerce").dropna()
    if values.empty or scale <= 0:
        return 0.0
    score = max(float(values.mean()), 0.0) / scale
    return round(min(score, 1.0), 6)


def prepare_evidence_heatmap(evidence_matrix: pd.DataFrame) -> pd.DataFrame:
    """Build normalized 0-1 evidence components for a compact heatmap."""
    priority = prepare_axis_priority(evidence_matrix)
    working = evidence_matrix.set_index("axis_id")
    rows: list[dict[str, object]] = []
    for _, item in priority.iterrows():
        axis_id = item["axis_id"]
        row = working.loc[axis_id]
        rows.append(
            {
                "display_rank": item["display_rank"],
                "axis_id": axis_id,
                "axis_short_label": item["axis_short_label"],
                "Spatial niche": _positive_scaled_mean(
                    row,
                    ["source_spatial_delta", "target_spatial_delta"],
                    scale=0.20,
                ),
                "Specificity audit": _positive_scaled_mean(
                    row,
                    ["source_specificity_fraction", "target_specificity_fraction"],
                    scale=1.0,
                ),
                "Bulk trend": _positive_scaled_mean(row, ["bulk_delta_mean"], scale=40.0),
                "snRNA program": _positive_scaled_mean(row, ["snrna_program_delta"], scale=0.30),
                "Tumor-adjacent scRNA": _positive_scaled_mean(
                    row,
                    ["scrna_tumor_adjacent_delta"],
                    scale=0.10,
                ),
                "Target prioritization": _positive_scaled_mean(
                    row,
                    ["top_continuous_priority_score"],
                    scale=1.0,
                ),
            }
        )
    return pd.DataFrame(rows)


def simplify_evidence_type(value: str) -> str:
    """Shorten coupling evidence labels for figure annotations."""
    mapping = {
        "source_near_epithelial_progenitor": "source-side",
        "target_near_epithelial_progenitor": "receptor-side",
    }
    return mapping.get(str(value), str(value).replace("_", " "))


def prepare_top_perturbation_effects(
    perturbation_ranking: pd.DataFrame,
    top_n: int = 8,
    perturbation_type: str = "gene",
) -> pd.DataFrame:
    """Prepare strongest virtual-perturbation effects as positive coupling losses."""
    required = {
        "axis_id",
        "perturbation_type",
        "perturbed_genes",
        "perturbation_factor",
        "evidence_type",
        "n_samples",
        "coupling_relative_delta_mean",
        "continuous_priority_score",
    }
    missing = required.difference(perturbation_ranking.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"perturbation ranking is missing required columns: {missing_text}")

    table = perturbation_ranking.copy()
    table = table[table["perturbation_type"].eq(perturbation_type)].copy()
    table["coupling_relative_delta_mean"] = pd.to_numeric(
        table["coupling_relative_delta_mean"],
        errors="coerce",
    )
    table["continuous_priority_score"] = pd.to_numeric(
        table["continuous_priority_score"],
        errors="coerce",
    ).fillna(0.0)
    table["coupling_loss"] = -table["coupling_relative_delta_mean"]
    table = table[table["coupling_loss"].gt(0)].copy()
    table["axis_short_label"] = table["axis_id"].map(axis_short_label)
    table["evidence_short"] = table["evidence_type"].map(simplify_evidence_type)
    table["factor_label"] = table["perturbation_factor"].map(lambda value: f"x{float(value):g}")
    table["perturbation_label"] = (
        table["perturbed_genes"].astype(str) + " " + table["factor_label"] + " (" + table["evidence_short"] + ")"
    )
    table = table.sort_values(
        ["continuous_priority_score", "coupling_loss"],
        ascending=[False, False],
    ).head(top_n)
    return table[
        [
            "axis_id",
            "axis_short_label",
            "perturbed_genes",
            "perturbation_factor",
            "factor_label",
            "evidence_short",
            "n_samples",
            "coupling_relative_delta_mean",
            "coupling_loss",
            "continuous_priority_score",
            "perturbation_label",
        ]
    ].reset_index(drop=True)


def signature_short_label(signature: str) -> str:
    """Return compact display labels for signature-level figures."""
    return SIGNATURE_SHORT_LABELS.get(signature, signature.replace("_", " "))


def _class_label(value: str) -> str:
    mapping = {
        "epithelial": "Epithelial",
        "macrophage": "Macrophage",
    }
    text = str(value)
    return mapping.get(text, text.replace("_", " ").title())


def prepare_grn_target_ranking_source(target_ranking: pd.DataFrame) -> pd.DataFrame:
    """Prepare GSE308103 GRN-level target ranking rows for supplementary figures."""
    required = {
        "broad_class",
        "target_gene",
        "top_impacted_signature",
        "top_signature_mean_impact",
        "top_signature_max_impact",
        "top_signature_sum_impact",
        "top_impacted_genes",
    }
    missing = required.difference(target_ranking.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"GRN target ranking is missing required columns: {missing_text}")

    table = target_ranking.copy()
    for column in ("top_signature_mean_impact", "top_signature_max_impact", "top_signature_sum_impact"):
        table[column] = pd.to_numeric(table[column], errors="coerce").fillna(0.0)
    table["target_gene"] = table["target_gene"].astype(str)
    table["broad_class"] = table["broad_class"].astype(str)
    table["class_label"] = table["broad_class"].map(_class_label)
    table["top_impacted_signature"] = table["top_impacted_signature"].astype(str)
    table["signature_label"] = table["top_impacted_signature"].map(signature_short_label)
    table["top_impacted_genes"] = table["top_impacted_genes"].fillna("").astype(str)
    table = table.sort_values(
        ["top_signature_mean_impact", "top_signature_max_impact", "target_gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    table.insert(0, "display_rank", range(1, len(table) + 1))
    return table[
        [
            "display_rank",
            "broad_class",
            "class_label",
            "target_gene",
            "top_impacted_signature",
            "signature_label",
            "top_signature_mean_impact",
            "top_signature_max_impact",
            "top_signature_sum_impact",
            "top_impacted_genes",
        ]
    ]


def prepare_grn_top_gene_source(target_source: pd.DataFrame, top_n: int = 5) -> pd.DataFrame:
    """Split comma-delimited top impacted genes into long-form tile source data."""
    required = {
        "display_rank",
        "broad_class",
        "class_label",
        "target_gene",
        "top_impacted_signature",
        "signature_label",
        "top_impacted_genes",
    }
    missing = required.difference(target_source.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"GRN target source is missing required columns: {missing_text}")

    rows: list[dict[str, object]] = []
    for row in target_source.itertuples(index=False):
        genes = [gene.strip() for gene in str(row.top_impacted_genes).split(",") if gene.strip()]
        for position, gene in enumerate(genes[:top_n], start=1):
            rows.append(
                {
                    "display_rank": int(row.display_rank),
                    "broad_class": row.broad_class,
                    "class_label": row.class_label,
                    "target_gene": row.target_gene,
                    "top_impacted_signature": row.top_impacted_signature,
                    "signature_label": row.signature_label,
                    "gene_position": position,
                    "impacted_gene": gene,
                }
            )
    return pd.DataFrame(rows)


def prepare_signature_celltype_heatmap(
    celltype_summary: pd.DataFrame,
    signatures: Iterable[str],
    celltypes: Iterable[str] = CELLTYPE_ORDER,
) -> pd.DataFrame:
    """Prepare per-signature relative cell-type scores from GSE131907 summaries."""
    required = {"Cell_type.refined", "score", "mean_score"}
    missing = required.difference(celltype_summary.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"cell-type summary is missing required columns: {missing_text}")

    working = celltype_summary.copy()
    working["Cell_type.refined"] = working["Cell_type.refined"].fillna("Unassigned")
    working["mean_score"] = pd.to_numeric(working["mean_score"], errors="coerce").fillna(0.0)
    signatures = list(signatures)
    celltypes = list(celltypes)
    subset = working[working["score"].isin(signatures) & working["Cell_type.refined"].isin(celltypes)].copy()
    matrix = subset.pivot_table(
        index="score",
        columns="Cell_type.refined",
        values="mean_score",
        aggfunc="mean",
        fill_value=0.0,
    )
    matrix = matrix.reindex(index=signatures, columns=celltypes, fill_value=0.0)
    max_per_signature = matrix.max(axis=1).replace(0, 1.0)
    relative = matrix.div(max_per_signature, axis=0)
    relative.insert(0, "signature_label", [signature_short_label(item) for item in relative.index])
    relative.insert(0, "signature", relative.index)
    return relative.reset_index(drop=True)


def prepare_specificity_status_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Prepare expected/off-target/missing gene counts and fractions per signature."""
    required = {"signature", "specificity_status", "n_genes"}
    missing = required.difference(summary.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"specificity summary is missing required columns: {missing_text}")

    table = summary.copy()
    table["n_genes"] = pd.to_numeric(table["n_genes"], errors="coerce").fillna(0).astype(int)
    pivot = table.pivot_table(
        index="signature",
        columns="specificity_status",
        values="n_genes",
        aggfunc="sum",
        fill_value=0,
    )
    for column in ["expected", "off_target", "missing"]:
        if column not in pivot.columns:
            pivot[column] = 0
    ordered_signatures = [
        signature
        for signature in SIGNATURE_SHORT_LABELS
        if signature in pivot.index and signature.endswith("_vs_other_epithelial")
        or signature in pivot.index and signature.endswith("_vs_other_macrophage")
    ]
    remaining = [signature for signature in pivot.index if signature not in set(ordered_signatures)]
    pivot = pivot.reindex(ordered_signatures + sorted(remaining))
    pivot = pivot[["expected", "off_target", "missing"]].astype(int)
    pivot["total_genes"] = pivot.sum(axis=1)
    for column in ["expected", "off_target", "missing"]:
        pivot[f"{column}_fraction"] = (pivot[column] / pivot["total_genes"].replace(0, pd.NA)).fillna(0.0)
    pivot.insert(0, "signature_label", [signature_short_label(item) for item in pivot.index])
    pivot.insert(0, "signature", pivot.index)
    return pivot.reset_index(drop=True)


def prepare_signature_refinement_source(
    audit: pd.DataFrame,
    filtered_genes: Iterable[str],
    signature: str = "spp1_macrophage_vs_other_macrophage",
) -> pd.DataFrame:
    """Prepare gene-level source data for specificity-refinement supplementary figures."""
    required = {"signature", "rank", "gene", "expected_celltype", "top_celltype", "specificity_status"}
    missing = required.difference(audit.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"specificity audit is missing required columns: {missing_text}")

    filtered_set = {str(gene) for gene in filtered_genes}
    table = audit[audit["signature"].eq(signature)].copy()
    table["rank"] = pd.to_numeric(table["rank"], errors="coerce").fillna(0).astype(int)
    table["gene"] = table["gene"].astype(str)
    table["specificity_status"] = table["specificity_status"].fillna("missing").astype(str)
    table["expected_celltype"] = table["expected_celltype"].fillna("Unassigned").astype(str)
    table["top_celltype"] = table["top_celltype"].replace("", pd.NA).fillna("Missing in reference").astype(str)
    table["retained_after_audit"] = table["gene"].isin(filtered_set)
    table["retained_label"] = table["retained_after_audit"].map({True: "Retained", False: "Removed"})
    status_order = {"expected": 0, "off_target": 1, "missing": 2}
    table["status_order"] = table["specificity_status"].map(status_order).fillna(99).astype(int)
    table["specificity_label"] = table["specificity_status"].map(
        {
            "expected": "Expected myeloid",
            "off_target": "Off-target",
            "missing": "Missing",
        }
    ).fillna(table["specificity_status"].str.replace("_", " ").str.title())
    table = table.sort_values("rank").reset_index(drop=True)
    table["display_rank"] = range(1, len(table) + 1)
    return table[
        [
            "signature",
            "display_rank",
            "rank",
            "gene",
            "expected_celltype",
            "top_celltype",
            "specificity_status",
            "specificity_label",
            "status_order",
            "retained_after_audit",
            "retained_label",
        ]
    ]


def prepare_signature_refinement_status_summary(source: pd.DataFrame) -> pd.DataFrame:
    """Summarize original and retained specificity-status counts."""
    required = {"specificity_status", "retained_after_audit"}
    missing = required.difference(source.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"refinement source is missing required columns: {missing_text}")

    status_order = ["expected", "off_target", "missing"]
    table = source.copy()
    table["specificity_status"] = table["specificity_status"].fillna("missing").astype(str)
    table["retained_after_audit"] = table["retained_after_audit"].fillna(False).astype(bool)
    original = table.groupby("specificity_status").size().rename("original_n")
    retained = table[table["retained_after_audit"]].groupby("specificity_status").size().rename("retained_n")
    summary = pd.concat([original, retained], axis=1).fillna(0).astype(int)
    summary = summary.reindex(status_order + [item for item in summary.index if item not in status_order], fill_value=0)
    summary["removed_n"] = summary["original_n"] - summary["retained_n"]
    summary["original_fraction"] = (
        summary["original_n"] / summary["original_n"].sum() if int(summary["original_n"].sum()) else 0.0
    )
    summary["retained_fraction"] = (
        summary["retained_n"] / summary["retained_n"].sum() if int(summary["retained_n"].sum()) else 0.0
    )
    summary.insert(
        0,
        "specificity_label",
        [
            {
                "expected": "Expected myeloid",
                "off_target": "Off-target",
                "missing": "Missing",
            }.get(item, str(item).replace("_", " ").title())
            for item in summary.index
        ],
    )
    summary.insert(0, "specificity_status", summary.index)
    return summary.reset_index(drop=True)


def prepare_candidate_gene_specificity(
    gene_top_celltype: pd.DataFrame,
    genes: Iterable[str] = ("MIF", "CD74", "CXCR4", "SPP1", "TREM2", "PLA2G7", "IL1B", "TNF", "CXCL8"),
) -> pd.DataFrame:
    """Prepare top cell-type expression calls for candidate mechanism genes."""
    required = {"gene", "Cell_type.refined", "mean_expression"}
    missing = required.difference(gene_top_celltype.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"gene top-celltype table is missing required columns: {missing_text}")

    genes = list(genes)
    table = gene_top_celltype[gene_top_celltype["gene"].isin(genes)].copy()
    table["gene"] = pd.Categorical(table["gene"], categories=genes, ordered=True)
    table["Cell_type.refined"] = table["Cell_type.refined"].fillna("Unassigned")
    table["mean_expression"] = pd.to_numeric(table["mean_expression"], errors="coerce").fillna(0.0)
    table = table.sort_values("gene").reset_index(drop=True)
    return table[["gene", "Cell_type.refined", "mean_expression"]]


def prepare_axis_stage_heatmap(
    spatial_by_stage: pd.DataFrame,
    stages: Iterable[str] = SPATIAL_STAGE_ORDER,
) -> pd.DataFrame:
    """Prepare candidate-axis spatial deltas across ordered stages."""
    required = {"axis_id", "evidence_type", "stage", "enrichment_delta_mean", "n_samples"}
    missing = required.difference(spatial_by_stage.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"axis spatial table is missing required columns: {missing_text}")

    stages = list(stages)
    table = spatial_by_stage.copy()
    table = table[table["stage"].isin(stages)].copy()
    table["enrichment_delta_mean"] = pd.to_numeric(table["enrichment_delta_mean"], errors="coerce").fillna(0.0)
    table["axis_short_label"] = table["axis_id"].map(axis_short_label)
    table["evidence_short"] = table["evidence_type"].map(simplify_evidence_type)
    table["row_label"] = table["axis_short_label"] + " (" + table["evidence_short"] + ")"
    axis_order = [
        "mif_cd74_cxcr4",
        "spp1_trem2_macrophage_epithelial",
        "c1q_apoe_trem2_lgals3",
        "cxcl9_cxcl10_cxcr3",
        "inflammatory_il1_tnf_cxcl8",
    ]
    evidence_order = ["source-side", "receptor-side"]
    row_order = []
    for axis_id in axis_order:
        axis_label = axis_short_label(axis_id)
        for evidence in evidence_order:
            row_order.append(f"{axis_label} ({evidence})")
    matrix = table.pivot_table(
        index="row_label",
        columns="stage",
        values="enrichment_delta_mean",
        aggfunc="mean",
        fill_value=0.0,
    )
    observed_rows = [row for row in row_order if row in matrix.index]
    extra_rows = [row for row in matrix.index if row not in set(observed_rows)]
    matrix = matrix.reindex(index=observed_rows + sorted(extra_rows), columns=stages, fill_value=0.0)
    matrix.insert(0, "row_label", matrix.index)
    return matrix.reset_index(drop=True)


def _gene_order_map(genes: Iterable[str]) -> dict[str, int]:
    return {gene: i for i, gene in enumerate(genes)}


def _best_gene_rows(
    ranking: pd.DataFrame,
    genes: Iterable[str],
    perturbation_factor: float | None = None,
) -> pd.DataFrame:
    table = ranking.copy()
    table = table[table["perturbation_type"].eq("gene") & table["perturbed_genes"].isin(list(genes))].copy()
    if perturbation_factor is not None:
        table = table[pd.to_numeric(table["perturbation_factor"], errors="coerce").eq(float(perturbation_factor))].copy()
    table["continuous_priority_score"] = pd.to_numeric(
        table["continuous_priority_score"],
        errors="coerce",
    ).fillna(0.0)
    table = table.sort_values(
        ["perturbed_genes", "perturbation_factor", "continuous_priority_score"],
        ascending=[True, True, False],
    )
    if perturbation_factor is None:
        return table.groupby(["perturbed_genes", "perturbation_factor"], as_index=False).head(1).reset_index(drop=True)
    return table.groupby("perturbed_genes", as_index=False).head(1).reset_index(drop=True)


def prepare_perturbation_dose_response(
    continuous_ranking: pd.DataFrame,
    genes: Iterable[str] = ("MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7"),
) -> pd.DataFrame:
    """Prepare retained-factor dose-response rows from continuous perturbation summaries."""
    required = {
        "perturbation_type",
        "perturbed_genes",
        "perturbation_factor",
        "axis_id",
        "evidence_type",
        "n_samples",
        "coupling_relative_delta_mean",
        "continuous_priority_score",
    }
    missing = required.difference(continuous_ranking.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"continuous ranking is missing required columns: {missing_text}")

    genes = list(genes)
    order = _gene_order_map(genes)
    best = _best_gene_rows(continuous_ranking, genes)
    best["perturbation_factor"] = pd.to_numeric(best["perturbation_factor"], errors="coerce")
    best["coupling_relative_delta_mean"] = pd.to_numeric(
        best["coupling_relative_delta_mean"],
        errors="coerce",
    ).fillna(0.0)
    best["coupling_loss"] = (-best["coupling_relative_delta_mean"]).clip(lower=0.0)
    best["coupling_remaining"] = (1.0 + best["coupling_relative_delta_mean"]).clip(lower=0.0)
    best["evidence_short"] = best["evidence_type"].map(simplify_evidence_type)
    best["line_label"] = best["perturbed_genes"].astype(str) + " (" + best["evidence_short"] + ")"

    baseline_rows = []
    for gene in genes:
        gene_rows = best[best["perturbed_genes"].eq(gene)]
        if gene_rows.empty:
            continue
        template = gene_rows.sort_values("continuous_priority_score", ascending=False).iloc[0].to_dict()
        template["perturbation_id"] = f"{gene}_baseline_x1"
        template["perturbation_factor"] = 1.0
        template["coupling_relative_delta_mean"] = 0.0
        template["continuous_priority_score"] = 0.0
        template["coupling_loss"] = 0.0
        template["coupling_remaining"] = 1.0
        baseline_rows.append(template)

    table = pd.concat([pd.DataFrame(baseline_rows), best], ignore_index=True)
    table["_gene_order"] = table["perturbed_genes"].map(order)
    table = table.sort_values(["_gene_order", "perturbation_factor"], ascending=[True, False])
    return table[
        [
            "perturbation_id",
            "axis_id",
            "perturbed_genes",
            "evidence_type",
            "evidence_short",
            "line_label",
            "n_samples",
            "perturbation_factor",
            "coupling_relative_delta_mean",
            "coupling_loss",
            "coupling_remaining",
            "continuous_priority_score",
        ]
    ].reset_index(drop=True)


def prepare_perturbation_stage_loss(
    continuous_effects: pd.DataFrame,
    continuous_ranking: pd.DataFrame,
    genes: Iterable[str] = ("MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7"),
    stages: Iterable[str] = SPATIAL_STAGE_ORDER,
) -> pd.DataFrame:
    """Summarize full-dropout continuous coupling loss by stage for selected genes."""
    required_effects = {"perturbation_id", "stage", "coupling_relative_delta"}
    missing_effects = required_effects.difference(continuous_effects.columns)
    if missing_effects:
        missing_text = ", ".join(sorted(missing_effects))
        raise ValueError(f"continuous effects table is missing required columns: {missing_text}")

    genes = list(genes)
    stages = list(stages)
    best = _best_gene_rows(continuous_ranking, genes, perturbation_factor=0.0)
    best = best[["perturbation_id", "perturbed_genes", "axis_id", "evidence_type"]].copy()
    best["evidence_short"] = best["evidence_type"].map(simplify_evidence_type)
    best["row_label"] = best["perturbed_genes"].astype(str) + " (" + best["evidence_short"] + ")"

    effects = continuous_effects[continuous_effects["perturbation_id"].isin(best["perturbation_id"])].copy()
    effects = effects[effects["stage"].isin(stages)].copy()
    effects["coupling_relative_delta"] = pd.to_numeric(
        effects["coupling_relative_delta"],
        errors="coerce",
    )
    summary = (
        effects.groupby(["perturbation_id", "stage"], dropna=False)
        .agg(
            n_samples=("sample", "nunique") if "sample" in effects.columns else ("stage", "size"),
            coupling_loss=("coupling_relative_delta", lambda values: max(-float(values.mean()), 0.0)),
        )
        .reset_index()
        .merge(best, on="perturbation_id", how="left")
    )
    order = _gene_order_map(genes)
    summary["_gene_order"] = summary["perturbed_genes"].map(order)
    summary["stage"] = pd.Categorical(summary["stage"], categories=stages, ordered=True)
    summary = summary.sort_values(["_gene_order", "stage"]).drop(columns="_gene_order").reset_index(drop=True)
    return summary[
        [
            "perturbation_id",
            "perturbed_genes",
            "axis_id",
            "evidence_short",
            "row_label",
            "stage",
            "n_samples",
            "coupling_loss",
        ]
    ]


def prepare_perturbation_method_concordance(
    continuous_ranking: pd.DataFrame,
    virtual_ranking: pd.DataFrame,
    genes: Iterable[str] = ("MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7"),
    perturbation_factor: float = 0.0,
) -> pd.DataFrame:
    """Compare continuous-coupling loss with top-quantile dropout priority."""
    required_continuous = {
        "perturbation_id",
        "perturbation_type",
        "perturbed_genes",
        "perturbation_factor",
        "axis_id",
        "evidence_type",
        "coupling_relative_delta_mean",
        "continuous_priority_score",
    }
    required_virtual = {
        "perturbation_id",
        "axis_id",
        "evidence_type",
        "panel_relative_delta",
        "observed_fraction_delta",
        "dropout_priority_score",
    }
    missing_continuous = required_continuous.difference(continuous_ranking.columns)
    missing_virtual = required_virtual.difference(virtual_ranking.columns)
    if missing_continuous:
        missing_text = ", ".join(sorted(missing_continuous))
        raise ValueError(f"continuous ranking is missing required columns: {missing_text}")
    if missing_virtual:
        missing_text = ", ".join(sorted(missing_virtual))
        raise ValueError(f"virtual ranking is missing required columns: {missing_text}")

    genes = list(genes)
    best = _best_gene_rows(continuous_ranking, genes, perturbation_factor=perturbation_factor)
    merged = best.merge(
        virtual_ranking[
            [
                "perturbation_id",
                "axis_id",
                "evidence_type",
                "panel_relative_delta",
                "observed_fraction_delta",
                "dropout_priority_score",
            ]
        ],
        on=["perturbation_id", "axis_id", "evidence_type"],
        how="left",
    )
    merged["coupling_relative_delta_mean"] = pd.to_numeric(
        merged["coupling_relative_delta_mean"],
        errors="coerce",
    )
    merged["continuous_coupling_loss"] = (-merged["coupling_relative_delta_mean"]).clip(lower=0.0)
    merged["panel_relative_delta"] = pd.to_numeric(merged["panel_relative_delta"], errors="coerce")
    merged["top_quantile_panel_loss"] = (-merged["panel_relative_delta"]).clip(lower=0.0)
    merged["dropout_priority_score"] = pd.to_numeric(merged["dropout_priority_score"], errors="coerce").fillna(0.0)
    merged["evidence_short"] = merged["evidence_type"].map(simplify_evidence_type)
    merged["point_label"] = merged["perturbed_genes"].astype(str) + " (" + merged["evidence_short"] + ")"
    order = _gene_order_map(genes)
    merged["_gene_order"] = merged["perturbed_genes"].map(order)
    merged = merged.sort_values(["_gene_order", "continuous_coupling_loss"], ascending=[True, False])
    return merged[
        [
            "perturbation_id",
            "axis_id",
            "perturbed_genes",
            "evidence_type",
            "evidence_short",
            "point_label",
            "perturbation_factor",
            "continuous_coupling_loss",
            "top_quantile_panel_loss",
            "observed_fraction_delta",
            "dropout_priority_score",
            "continuous_priority_score",
        ]
    ].reset_index(drop=True)
