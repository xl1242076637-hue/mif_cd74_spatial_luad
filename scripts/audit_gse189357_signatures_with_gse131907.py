#!/usr/bin/env python
"""Audit GSE189357 refined signatures using GSE131907 cell-type specificity."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.signature_specificity import audit_signature_specificity, filter_signatures_by_specificity  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--signatures",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_state_signature_genes.json",
    )
    parser.add_argument(
        "--top-celltypes",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse131907_selected_gene_top_celltype.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    signatures = json.loads(args.signatures.read_text(encoding="utf-8"))
    top_celltypes = pd.read_csv(args.top_celltypes)
    audit = audit_signature_specificity(signatures, top_celltypes)
    filtered = filter_signatures_by_specificity(signatures, audit)
    summary = (
        audit.groupby(["signature", "expected_celltype", "specificity_status"], dropna=False)
        .size()
        .reset_index(name="n_genes")
    )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    audit_output = args.table_dir / "gse189357_refined_signature_gse131907_specificity_audit.csv"
    summary_output = args.table_dir / "gse189357_refined_signature_gse131907_specificity_summary.csv"
    filtered_output = args.table_dir / "gse189357_refined_signature_genes_gse131907_specificity_filtered.json"
    audit.to_csv(audit_output, index=False, encoding="utf-8-sig")
    summary.to_csv(summary_output, index=False, encoding="utf-8-sig")
    filtered_output.write_text(json.dumps(filtered, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Wrote specificity audit: {audit_output}")
    print(f"Wrote specificity summary: {summary_output}")
    print(f"Wrote filtered signatures: {filtered_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
