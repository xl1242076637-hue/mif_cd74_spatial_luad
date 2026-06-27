#!/usr/bin/env python
"""Package manuscripts, figures, tables, code and documentation for handoff."""

from __future__ import annotations

import argparse
import csv
import shutil
from datetime import datetime
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]


FILE_GROUPS = {
    "manuscripts": [
        "docs/manuscript_communications_biology_draft.md",
        "docs/manuscript_communications_biology_draft_zh.md",
        "docs/manuscript_communications_biology_draft_en_with_figures_tables.docx",
        "docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx",
        "docs/manuscript_communications_biology_draft_zh_zotero_fields.docx",
        "docs/manuscript_communications_biology_bilingual_with_figures_tables.docx",
    ],
    "documentation": [
        "README.md",
        "docs/reproducibility.md",
        "docs/analysis_log.md",
        "docs/decision_log.md",
        "docs/data_inventory.md",
        "docs/figure_legends.md",
        "docs/figure_plan.md",
        "docs/manuscript_figure_storyline_2026-06-03.md",
        "docs/supplementary_tables_index.md",
        "docs/grn_virtual_perturbation_2026-06-03.md",
        "docs/sctenifoldknk_original_smoke_2026-06-10.md",
        "docs/word_export_notes.md",
    ],
    "references": [
        "results/tables/communications_biology_reference_list.csv",
        "results/references/communications_biology_references.ris",
    ],
}

DIRECTORY_GROUPS = {
    "figures": [
        "results/figures",
    ],
    "tables": [
        "results/tables",
        "results/supplementary_tables",
    ],
    "code": [
        "scripts",
        "src",
        "tests",
        "config",
    ],
}

EXCLUDE_DIR_NAMES = {"__pycache__", ".pytest_cache"}
EXCLUDE_SUFFIXES = {".pyc"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=PROJECT_ROOT / "deliverables",
    )
    parser.add_argument(
        "--name",
        default="luad_epithelial_macrophage_niche_final_package_2026-06-16",
    )
    return parser.parse_args()


def should_skip(path: Path) -> bool:
    return any(part in EXCLUDE_DIR_NAMES for part in path.parts) or path.suffix in EXCLUDE_SUFFIXES


def copy_file(relative_path: str, destination_root: Path, manifest: list[dict[str, object]]) -> None:
    source = PROJECT_ROOT / relative_path
    destination = destination_root / relative_path
    if not source.exists():
        manifest.append(
            {
                "group": "missing",
                "source": relative_path,
                "destination": "",
                "bytes": 0,
                "status": "missing",
            }
        )
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    manifest.append(
        {
            "group": destination.parts[len(destination_root.parts)] if len(destination.parts) > len(destination_root.parts) else "",
            "source": relative_path,
            "destination": destination.relative_to(destination_root).as_posix(),
            "bytes": destination.stat().st_size,
            "status": "copied",
        }
    )


def copy_directory(relative_dir: str, destination_root: Path, manifest: list[dict[str, object]]) -> None:
    source_dir = PROJECT_ROOT / relative_dir
    if not source_dir.exists():
        manifest.append(
            {
                "group": "missing",
                "source": relative_dir,
                "destination": "",
                "bytes": 0,
                "status": "missing",
            }
        )
        return
    for source in sorted(path for path in source_dir.rglob("*") if path.is_file()):
        relative_path = source.relative_to(PROJECT_ROOT)
        if should_skip(relative_path):
            continue
        destination = destination_root / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        manifest.append(
            {
                "group": relative_dir,
                "source": relative_path.as_posix(),
                "destination": destination.relative_to(destination_root).as_posix(),
                "bytes": destination.stat().st_size,
                "status": "copied",
            }
        )


def write_manifest(destination_root: Path, manifest: list[dict[str, object]]) -> None:
    manifest_path = destination_root / "PACKAGE_MANIFEST.csv"
    with manifest_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["group", "source", "destination", "bytes", "status"])
        writer.writeheader()
        writer.writerows(manifest)


def write_readme(destination_root: Path, manifest: list[dict[str, object]]) -> None:
    copied = [row for row in manifest if row["status"] == "copied"]
    missing = [row for row in manifest if row["status"] == "missing"]
    lines = [
        "# LUAD epithelial-macrophage niche final package",
        "",
        f"Created: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "This folder collects the manuscript drafts, figures, source tables, supplementary tables, code, configuration, references and reproducibility notes for the public-data LUAD epithelial progenitor-macrophage niche project.",
        "",
        "## Main manuscript files",
        "",
        "- `docs/manuscript_communications_biology_draft_en_with_figures_tables.docx`: English Word manuscript with figures and tables.",
        "- `docs/manuscript_communications_biology_draft_zh_with_figures_tables.docx`: Chinese Word manuscript with figures and tables.",
        "- `docs/manuscript_communications_biology_draft.md`: English Markdown source.",
        "- `docs/manuscript_communications_biology_draft_zh.md`: Chinese Markdown source.",
        "",
        "## Key folders",
        "",
        "- `results/figures/`: main and supplementary figures.",
        "- `results/tables/`: analysis tables and figure source data.",
        "- `results/supplementary_tables/`: manuscript-facing supplementary CSV tables, including ST31 for the original scTenifoldKnk sensitivity summary.",
        "- `scripts/`, `src/`, `tests/`, `config/`: analysis code, package code, tests and configuration.",
        "- `docs/reproducibility.md`: command-level reproduction notes.",
        "",
        "## Interpretation boundary",
        "",
        "Score-level in-silico target prioritization and GRN-level virtual perturbation prioritization are computational target-ranking layers. They are not wet-lab perturbation experiments, causal knockout validation or treatment-response proof.",
        "",
        "## Manifest",
        "",
        f"`PACKAGE_MANIFEST.csv` lists {len(copied)} copied files.",
    ]
    if missing:
        lines.extend(
            [
                "",
                "## Missing expected files",
                "",
                *[f"- `{row['source']}`" for row in missing],
            ]
        )
    (destination_root / "README_PACKAGE.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    destination_root = args.output_root / args.name
    destination_root.mkdir(parents=True, exist_ok=True)

    manifest: list[dict[str, object]] = []
    for paths in FILE_GROUPS.values():
        for relative_path in paths:
            copy_file(relative_path, destination_root, manifest)
    for dirs in DIRECTORY_GROUPS.values():
        for relative_dir in dirs:
            copy_directory(relative_dir, destination_root, manifest)

    write_manifest(destination_root, manifest)
    write_readme(destination_root, manifest)
    copied = sum(1 for row in manifest if row["status"] == "copied")
    missing = sum(1 for row in manifest if row["status"] == "missing")
    print(f"Packaged {copied} files into {destination_root}")
    if missing:
        print(f"Missing expected files: {missing}")
    return 0 if missing == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
