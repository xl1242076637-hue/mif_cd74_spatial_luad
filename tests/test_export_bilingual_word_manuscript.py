import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "export_bilingual_word_manuscript.py"
SPEC = importlib.util.spec_from_file_location("export_bilingual_word_manuscript", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_parse_markdown_blocks_preserves_headings_and_paragraphs():
    markdown = """# English Draft

## Results

### Score-level in-silico target prioritization ranks CD74

CD74 is prioritized as a receptor-side candidate.
"""

    blocks = MODULE.parse_markdown_blocks(markdown)

    assert blocks[0].kind == "h1"
    assert blocks[0].text == "English Draft"
    assert blocks[2].kind == "h3"
    assert "Score-level in-silico" in blocks[2].text
    assert blocks[3].kind == "p"
    assert "receptor-side" in blocks[3].text


def test_pair_bilingual_blocks_keeps_order_and_fills_missing_side():
    english = MODULE.parse_markdown_blocks("# Draft\n\n## Results\n\nEnglish paragraph.")
    chinese = MODULE.parse_markdown_blocks("# 草稿\n\n## 结果\n\n中文段落。\n\n额外中文段落。")

    pairs = MODULE.pair_bilingual_blocks(english, chinese)

    assert [(left.text if left else "", right.text if right else "") for left, right in pairs] == [
        ("Draft", "草稿"),
        ("Results", "结果"),
        ("English paragraph.", "中文段落。"),
        ("", "额外中文段落。"),
    ]


def test_figure_manifest_contains_grn_supplementary_figure():
    labels = [figure.label_zh for figure in MODULE.BILINGUAL_FIGURES]

    assert "补充图 4" in labels
    assert any("GRN-level virtual perturbation prioritization" in figure.caption_en for figure in MODULE.BILINGUAL_FIGURES)
