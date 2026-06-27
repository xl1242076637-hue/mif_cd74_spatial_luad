import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "export_english_word_manuscript.py"
SPEC = importlib.util.spec_from_file_location("export_english_word_manuscript", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_english_figure_manifest_uses_current_boundary_wording():
    captions = [figure.caption for figure in MODULE.ENGLISH_FIGURES]

    assert any("receptor-side CD74" in caption for caption in captions)
    assert any("GRN-level virtual perturbation prioritization" in caption for caption in captions)
    assert all("ranks MIF and CD74" not in caption for caption in captions)


def test_english_table_manifest_has_five_review_tables():
    labels = [table.label for table in MODULE.ENGLISH_TABLES]

    assert labels == ["Table 1", "Table 2", "Table 3", "Table 4", "Table 5"]
    assert any("score-level in-silico target prioritization" in table.title for table in MODULE.ENGLISH_TABLES)
