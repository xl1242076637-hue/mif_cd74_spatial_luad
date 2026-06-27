"""Dataset-specific metadata harmonization."""

from __future__ import annotations

import re


GSE282617_GROUP_MAP = {
    "ZCF": {
        "interpreted_stage": "Normal",
        "interpreted_condition": "normal_lung",
        "include_in_luad_progression": True,
    },
    "YWA": {
        "interpreted_stage": "AIS",
        "interpreted_condition": "adenocarcinoma_in_situ",
        "include_in_luad_progression": True,
    },
    "WJR": {
        "interpreted_stage": "MIA",
        "interpreted_condition": "minimally_invasive_adenocarcinoma",
        "include_in_luad_progression": True,
    },
    "JRX": {
        "interpreted_stage": "IAC",
        "interpreted_condition": "invasive_adenocarcinoma",
        "include_in_luad_progression": True,
    },
    "FLA": {
        "interpreted_stage": "LUSC",
        "interpreted_condition": "lung_squamous_cell_carcinoma_control",
        "include_in_luad_progression": False,
    },
}


def interpret_gse282617_group(group: str) -> dict:
    """Map GSE282617 group abbreviations to analysis-ready labels."""
    normalized = str(group).strip().upper()
    if normalized not in GSE282617_GROUP_MAP:
        raise ValueError(f"Unknown GSE282617 group code: {group!r}")
    return dict(GSE282617_GROUP_MAP[normalized])


def infer_lungpca_stage(title: str) -> str | None:
    """Infer Normal/AAH/AIS/MIA/LUAD labels from LungPCA GEO sample titles."""
    normalized = str(title).upper()
    for stage in ["NORMAL", "AAH", "AIS", "MIA", "LUAD"]:
        if re.search(rf"\b{stage}\b", normalized):
            return "Normal" if stage == "NORMAL" else stage
    return None


def interpret_gse164789_source(source_name: str) -> dict:
    """Map GSE164789 adjacent/tumor labels to coarse condition classes."""
    normalized = str(source_name).strip().lower()
    if "adjacent" in normalized:
        return {
            "interpreted_stage": "Adjacent",
            "interpreted_condition": "lung_neoplasm_adjacent",
            "include_in_luad_progression": False,
        }
    if "tumor" in normalized or "tumour" in normalized:
        return {
            "interpreted_stage": "Tumor",
            "interpreted_condition": "lung_neoplasm_tumor",
            "include_in_luad_progression": False,
        }
    return {}


def interpret_gse131907_source(source_name: str) -> dict:
    """Map GSE131907 source names to broad tissue-origin classes."""
    normalized = str(source_name).strip().lower()
    if "normal lung" in normalized:
        return {
            "interpreted_stage": "Normal",
            "interpreted_condition": "normal_lung_reference",
            "include_in_luad_progression": False,
        }
    if "tumour lung" in normalized or "tumor lung" in normalized:
        return {
            "interpreted_stage": "Primary tumor",
            "interpreted_condition": "primary_lung_tumor_reference",
            "include_in_luad_progression": False,
        }
    if "metastatic" in normalized or "pleural effusion" in normalized:
        return {
            "interpreted_stage": "Metastasis/effusion",
            "interpreted_condition": "metastatic_or_effusion_reference",
            "include_in_luad_progression": False,
        }
    if "normal lymph node" in normalized:
        return {
            "interpreted_stage": "Normal lymph node",
            "interpreted_condition": "normal_lymph_node_reference",
            "include_in_luad_progression": False,
        }
    return {}


def annotate_sample_record(record: dict) -> dict:
    """Add harmonized stage fields to one parsed GEO sample metadata record."""
    annotated = dict(record)
    series = annotated.get("series_accession")
    if series == "GSE282617" and annotated.get("group"):
        annotated.update(interpret_gse282617_group(annotated["group"]))
    elif series in {"GSE307534", "GSE308103"}:
        stage = infer_lungpca_stage(annotated.get("title", ""))
        if stage:
            annotated["interpreted_stage"] = stage
            annotated["interpreted_condition"] = "normal_precursor_luad_progression"
            annotated["include_in_luad_progression"] = True
    elif series == "GSE164789" and annotated.get("source_name"):
        annotated.update(interpret_gse164789_source(annotated["source_name"]))
    elif series == "GSE131907" and annotated.get("source_name"):
        annotated.update(interpret_gse131907_source(annotated["source_name"]))
    elif annotated.get("histological_type"):
        annotated["interpreted_stage"] = annotated["histological_type"]
        annotated["interpreted_condition"] = "lung_adenocarcinoma"
        annotated["include_in_luad_progression"] = True
    return annotated
