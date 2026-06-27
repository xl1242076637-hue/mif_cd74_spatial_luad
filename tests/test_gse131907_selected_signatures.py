import importlib.util
import json
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "score_gse131907_selected_signatures.py"
SNRNA_SCRIPT_PATH = PROJECT_ROOT / "scripts" / "score_gse308103_snrna_states.py"
SCRNA_SCRIPT_PATH = PROJECT_ROOT / "scripts" / "score_gse164789_scrna_states.py"


def load_script_module(path=SCRIPT_PATH, module_name="score_gse131907_selected_signatures"):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_load_panels_includes_candidate_mechanism_genes(tmp_path):
    module = load_script_module()
    marker_path = tmp_path / "markers.yaml"
    refined_path = tmp_path / "refined.json"
    mechanism_path = tmp_path / "candidate_mechanisms.yaml"

    marker_path.write_text(
        yaml.safe_dump(
            {
                "broad_classes": {"epithelial": ["EPCAM"]},
                "state_panels": {"spp1_macrophage": ["SPP1"]},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    refined_path.write_text(json.dumps({}), encoding="utf-8")
    mechanism_path.write_text(
        yaml.safe_dump(
            {
                "axes": [
                    {
                        "source_genes": ["MIF", "SPP1"],
                        "target_genes": ["CD74", "CXCR4"],
                        "bulk_genes": ["ITGAV"],
                        "perturbation_genes": ["NLRP3", "MIF"],
                    }
                ]
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    panels = module.load_panels(marker_path, refined_path, mechanism_path)

    assert panels["candidate_mechanism_genes"] == ["MIF", "SPP1", "CD74", "CXCR4", "ITGAV", "NLRP3"]
    assert module.flatten_genes(panels) == ["EPCAM", "SPP1", "MIF", "CD74", "CXCR4", "ITGAV", "NLRP3"]


def test_single_cell_state_scripts_keep_c1q_scores_for_mechanism_ranking():
    snrna = load_script_module(SNRNA_SCRIPT_PATH, "score_gse308103_snrna_states")
    scrna = load_script_module(SCRNA_SCRIPT_PATH, "score_gse164789_scrna_states")

    for module in (snrna, scrna):
        assert "c1q_macrophage_score" in module.KEY_SCORE_COLUMNS
        assert "refined_c1q_macrophage_score" in module.KEY_SCORE_COLUMNS
