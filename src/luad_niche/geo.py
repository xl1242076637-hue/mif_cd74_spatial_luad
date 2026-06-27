"""GEO accession and supplementary-file URL helpers."""

from __future__ import annotations

import re


_GSE_RE = re.compile(r"^GSE(\d+)$", re.IGNORECASE)


def geo_series_prefix(accession: str) -> str:
    """Return the NCBI GEO series folder prefix for a GSE accession."""
    match = _GSE_RE.match(accession.strip())
    if not match:
        raise ValueError(f"Expected a GSE accession like 'GSE189357', got {accession!r}.")

    number = int(match.group(1))
    prefix = number // 1000
    return f"GSE{prefix}nnn"


def geo_supplementary_url(accession: str) -> str:
    """Return the HTTPS GEO supplementary-file directory URL."""
    normalized = accession.strip().upper()
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{geo_series_prefix(normalized)}/{normalized}/suppl/"
    )


def geo_family_soft_url(accession: str) -> str:
    """Return the HTTPS GEO family SOFT file URL."""
    normalized = accession.strip().upper()
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{geo_series_prefix(normalized)}/{normalized}/soft/{normalized}_family.soft.gz"
    )
