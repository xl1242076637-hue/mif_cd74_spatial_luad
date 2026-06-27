import pandas as pd

from luad_niche.signature_specificity import audit_signature_specificity, filter_signatures_by_specificity


def test_audit_signature_specificity_marks_expected_top_celltype():
    signatures = {
        "spp1_macrophage_vs_other_macrophage": ["SPP1", "WFDC2", "MISSING"],
        "epithelial_progenitor_like_vs_other_epithelial": ["EPCAM", "LYZ"],
    }
    top_celltypes = pd.DataFrame(
        {
            "gene": ["SPP1", "WFDC2", "EPCAM", "LYZ"],
            "Cell_type.refined": ["Myeloid cells", "Epithelial cells", "Epithelial cells", "Myeloid cells"],
            "mean_expression": [1.0, 2.0, 3.0, 4.0],
        }
    )

    audit = audit_signature_specificity(signatures, top_celltypes)

    statuses = dict(zip(zip(audit["signature"], audit["gene"], strict=False), audit["specificity_status"], strict=False))
    assert statuses[("spp1_macrophage_vs_other_macrophage", "SPP1")] == "expected"
    assert statuses[("spp1_macrophage_vs_other_macrophage", "WFDC2")] == "off_target"
    assert statuses[("spp1_macrophage_vs_other_macrophage", "MISSING")] == "missing"
    assert statuses[("epithelial_progenitor_like_vs_other_epithelial", "EPCAM")] == "expected"
    assert statuses[("epithelial_progenitor_like_vs_other_epithelial", "LYZ")] == "off_target"


def test_filter_signatures_by_specificity_keeps_only_expected_genes():
    signatures = {
        "spp1_macrophage_vs_other_macrophage": ["SPP1", "WFDC2"],
        "epithelial_progenitor_like_vs_other_epithelial": ["EPCAM", "LYZ"],
    }
    audit = pd.DataFrame(
        {
            "signature": [
                "spp1_macrophage_vs_other_macrophage",
                "spp1_macrophage_vs_other_macrophage",
                "epithelial_progenitor_like_vs_other_epithelial",
                "epithelial_progenitor_like_vs_other_epithelial",
            ],
            "gene": ["SPP1", "WFDC2", "EPCAM", "LYZ"],
            "specificity_status": ["expected", "off_target", "expected", "off_target"],
        }
    )

    filtered = filter_signatures_by_specificity(signatures, audit)

    assert filtered == {
        "spp1_macrophage_vs_other_macrophage": ["SPP1"],
        "epithelial_progenitor_like_vs_other_epithelial": ["EPCAM"],
    }
