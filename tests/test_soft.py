from luad_niche.soft import normalize_characteristic_key, parse_soft_samples


def test_normalize_characteristic_key_fixes_common_geo_typo():
    assert normalize_characteristic_key("histolgical type") == "histological_type"
    assert normalize_characteristic_key("radiological type") == "radiological_type"


def test_parse_soft_samples_extracts_titles_and_characteristics():
    text = """
^SAMPLE = GSM1
!Sample_title = TD1 scRNA-seq
!Sample_geo_accession = GSM1
!Sample_source_name_ch1 = Lung adenocarcinoma
!Sample_characteristics_ch1 = histolgical type: IAC
!Sample_characteristics_ch1 = radiological type: SN
^SAMPLE = GSM2
!Sample_title = FLA1
!Sample_geo_accession = GSM2
!Sample_characteristics_ch1 = group: FLA
"""

    samples = parse_soft_samples(text, series_accession="GSETEST")

    assert samples == [
        {
            "series_accession": "GSETEST",
            "sample_accession": "GSM1",
            "title": "TD1 scRNA-seq",
            "source_name": "Lung adenocarcinoma",
            "histological_type": "IAC",
            "radiological_type": "SN",
        },
        {
            "series_accession": "GSETEST",
            "sample_accession": "GSM2",
            "title": "FLA1",
            "group": "FLA",
        },
    ]

