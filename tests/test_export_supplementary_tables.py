import importlib.util
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "export_supplementary_tables.py"
SPEC = importlib.util.spec_from_file_location("export_supplementary_tables", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_safe_filename_prefixes_stable_table_id():
    assert MODULE.safe_filename("ST03", "axis.csv") == "ST03_axis.csv"


def test_write_index_lists_packaged_tables(tmp_path):
    manifest = pd.DataFrame(
        [
            {
                "table_id": "ST01",
                "description": "Dataset composition",
                "source_path": "results/tables/source.csv",
                "output_path": "results/supplementary_tables/ST01_source.csv",
                "n_rows": 7,
                "n_columns": 4,
            }
        ]
    )
    output = tmp_path / "index.md"

    MODULE.write_index(output, manifest)

    text = output.read_text(encoding="utf-8")
    assert "`ST01`" in text
    assert "Dataset composition" in text
    assert "score-level in-silico target prioritization" in text
    assert "GRN-level virtual perturbation prioritization" in text


def test_table_specs_include_focused_orthogonal_validation_sources():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST18"] == "gse308103_snrna_candidate_gene_sample_summary.csv"
    assert table_specs["ST19"] == "gse308103_snrna_candidate_gene_stage_summary.csv"
    assert table_specs["ST20"] == "supplementary_figure_focused_orthogonal_validation_snrna_source.csv"
    assert table_specs["ST21"] == "supplementary_figure_focused_orthogonal_validation_bulk_source.csv"


def test_table_specs_include_mif_covariate_sensitivity_sources():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST22"] == "gse307534_mif_covariate_sensitivity_models.csv"
    assert table_specs["ST23"] == "gse307534_mif_covariate_paired_changes.csv"


def test_table_specs_include_spp1_signature_refinement_sources():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST24"] == "supplementary_figure_spp1_signature_refinement_gene_source.csv"
    assert table_specs["ST25"] == "supplementary_figure_spp1_signature_refinement_status_source.csv"
    assert table_specs["ST26"] == "supplementary_figure_spp1_signature_refinement_celltype_source.csv"


def test_table_specs_include_grn_virtual_perturbation_target_ranking():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST27"] == "gse308103_grn_virtual_perturbation_target_ranking.csv"


def test_table_specs_include_grn_virtual_perturbation_robustness_summary():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST28"] == "gse308103_grn_virtual_perturbation_robustness_summary.csv"


def test_table_specs_include_grn_cross_dataset_validation_tables():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST29"] == "grn_cross_dataset_signature_validation.csv"
    assert table_specs["ST30"] == "grn_cross_dataset_signature_validation_summary.csv"


def test_table_specs_include_original_sctenifoldknk_sensitivity_summary():
    table_specs = {table_id: source_name for table_id, _, source_name in MODULE.TABLE_SPECS}

    assert table_specs["ST31"] == "sctenifoldknk_original_expanded_interpretation.csv"
