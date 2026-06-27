from luad_niche.metadata import (
    annotate_sample_record,
    infer_lungpca_stage,
    interpret_gse131907_source,
    interpret_gse164789_source,
    interpret_gse282617_group,
)


def test_interpret_gse282617_group_maps_geo_group_codes():
    assert interpret_gse282617_group("ZCF") == {
        "interpreted_stage": "Normal",
        "interpreted_condition": "normal_lung",
        "include_in_luad_progression": True,
    }
    assert interpret_gse282617_group("YWA")["interpreted_stage"] == "AIS"
    assert interpret_gse282617_group("WJR")["interpreted_stage"] == "MIA"
    assert interpret_gse282617_group("JRX")["interpreted_stage"] == "IAC"
    assert interpret_gse282617_group("FLA") == {
        "interpreted_stage": "LUSC",
        "interpreted_condition": "lung_squamous_cell_carcinoma_control",
        "include_in_luad_progression": False,
    }


def test_infer_lungpca_stage_from_geo_titles():
    assert infer_lungpca_stage("Normal lung tissue of patient 3 [snRNA-seq]") == "Normal"
    assert infer_lungpca_stage("Lung tissue with AAH of patient 1 [ST]") == "AAH"
    assert infer_lungpca_stage("Lung tissue with AIS of patient 3 [ST]") == "AIS"
    assert infer_lungpca_stage("Lung tissue with MIA of patient 10 [ST]") == "MIA"
    assert infer_lungpca_stage("Lung tissue with LUAD of patient 10 [ST]") == "LUAD"


def test_annotate_sample_record_marks_lungpca_progression_samples():
    annotated = annotate_sample_record(
        {
            "series_accession": "GSE307534",
            "title": "Lung tissue with MIA of patient 10 [ST]",
        }
    )

    assert annotated["interpreted_stage"] == "MIA"
    assert annotated["interpreted_condition"] == "normal_precursor_luad_progression"
    assert annotated["include_in_luad_progression"] is True


def test_interpret_gse164789_source_maps_adjacent_and_tumor():
    assert interpret_gse164789_source("lung neoplasm-adjacent") == {
        "interpreted_stage": "Adjacent",
        "interpreted_condition": "lung_neoplasm_adjacent",
        "include_in_luad_progression": False,
    }
    assert interpret_gse164789_source("lung neoplasm-tumor")["interpreted_stage"] == "Tumor"


def test_interpret_gse131907_source_maps_reference_tissues():
    assert interpret_gse131907_source("Normal lung")["interpreted_stage"] == "Normal"
    assert interpret_gse131907_source("Tumour lung")["interpreted_stage"] == "Primary tumor"
    assert interpret_gse131907_source("Metastatic brain")["interpreted_stage"] == "Metastasis/effusion"
    assert interpret_gse131907_source("Pleural effusion")["interpreted_stage"] == "Metastasis/effusion"
    assert interpret_gse131907_source("Normal lymph node")["interpreted_stage"] == "Normal lymph node"
