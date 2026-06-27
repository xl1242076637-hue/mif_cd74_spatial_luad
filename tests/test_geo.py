import pytest

from luad_niche.geo import geo_family_soft_url, geo_series_prefix, geo_supplementary_url


def test_geo_series_prefix_groups_accession_by_thousands():
    assert geo_series_prefix("GSE189357") == "GSE189nnn"
    assert geo_series_prefix("GSE282617") == "GSE282nnn"
    assert geo_series_prefix("GSE307534") == "GSE307nnn"


def test_geo_series_prefix_rejects_invalid_accessions():
    with pytest.raises(ValueError, match="GSE accession"):
        geo_series_prefix("189357")

    with pytest.raises(ValueError, match="GSE accession"):
        geo_series_prefix("GSM123")


def test_geo_supplementary_url_uses_ncbi_ftp_https_endpoint():
    assert (
        geo_supplementary_url("GSE189357")
        == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/"
    )


def test_geo_family_soft_url_uses_ncbi_ftp_https_endpoint():
    assert (
        geo_family_soft_url("GSE308103")
        == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308103/soft/GSE308103_family.soft.gz"
    )
