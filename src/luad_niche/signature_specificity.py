"""Audit signature genes against external cell-type specificity evidence."""

from __future__ import annotations

import pandas as pd


def expected_celltype_for_signature(signature_name: str) -> str:
    """Return the expected broad cell type for a signature name."""
    if "macrophage" in signature_name:
        return "Myeloid cells"
    if "epithelial" in signature_name:
        return "Epithelial cells"
    return ""


def audit_signature_specificity(
    signatures: dict[str, list[str]],
    top_celltypes: pd.DataFrame,
) -> pd.DataFrame:
    """Audit whether each signature gene's top cell type matches expectation."""
    top_by_gene = top_celltypes.set_index("gene", drop=False)
    records = []
    for signature_name, genes in signatures.items():
        expected = expected_celltype_for_signature(signature_name)
        for rank, gene in enumerate(genes, start=1):
            if gene not in top_by_gene.index:
                records.append(
                    {
                        "signature": signature_name,
                        "rank": rank,
                        "gene": gene,
                        "expected_celltype": expected,
                        "top_celltype": "",
                        "top_mean_expression": None,
                        "specificity_status": "missing",
                    }
                )
                continue
            row = top_by_gene.loc[gene]
            top_celltype = row.get("Cell_type.refined", "")
            status = "expected" if expected and top_celltype == expected else "off_target"
            records.append(
                {
                    "signature": signature_name,
                    "rank": rank,
                    "gene": gene,
                    "expected_celltype": expected,
                    "top_celltype": top_celltype,
                    "top_mean_expression": row.get("mean_expression"),
                    "specificity_status": status,
                }
            )
    return pd.DataFrame(records)


def filter_signatures_by_specificity(
    signatures: dict[str, list[str]],
    audit: pd.DataFrame,
) -> dict[str, list[str]]:
    """Keep genes whose specificity audit status is expected."""
    keep = audit[audit["specificity_status"] == "expected"]
    filtered = {}
    for signature_name, genes in signatures.items():
        keep_genes = set(keep.loc[keep["signature"] == signature_name, "gene"])
        filtered[signature_name] = [gene for gene in genes if gene in keep_genes]
    return filtered
