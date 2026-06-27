import gzip
from pathlib import Path

from luad_niche.sn_matrix import (
    parse_gse308103_filename,
    read_tabular_count_barcodes,
    read_tabular_selected_genes,
    read_tabular_selected_genes_with_totals,
)


def write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def test_parse_gse308103_filename_extracts_accession_and_sample_name():
    parsed = parse_gse308103_filename(Path("GSM9237926_P10_MIA.raw_counts.mtx.txt.gz"))

    assert parsed == {"sample_accession": "GSM9237926", "sample_name": "P10_MIA"}


def test_read_tabular_count_barcodes_reads_first_line(tmp_path):
    matrix_path = tmp_path / "sample.raw_counts.mtx.txt.gz"
    write_gzip(matrix_path, "cellA\tcellB\nGENE1\t1\t2\n")

    assert read_tabular_count_barcodes(matrix_path) == ["cellA", "cellB"]


def test_read_tabular_selected_genes_with_totals_streams_counts(tmp_path):
    matrix_path = tmp_path / "sample.raw_counts.mtx.txt.gz"
    write_gzip(
        matrix_path,
        "cellA\tcellB\tcellC\n"
        "GENE1\t1\t0\t2\n"
        "GENE2\t3\t4\t0\n"
        "GENE1\t2\t0\t1\n",
    )

    counts, totals = read_tabular_selected_genes_with_totals(matrix_path, ["GENE1", "MISSING"])

    assert counts.index.tolist() == ["cellA", "cellB", "cellC"]
    assert counts.columns.tolist() == ["GENE1"]
    assert counts["GENE1"].tolist() == [3.0, 0.0, 3.0]
    assert totals.tolist() == [6.0, 4.0, 3.0]


def test_read_tabular_selected_genes_skips_unrequested_count_rows(tmp_path):
    matrix_path = tmp_path / "sample.raw_counts.mtx.txt.gz"
    write_gzip(
        matrix_path,
        "cellA\tcellB\n"
        "GENE1\t1\t2\n"
        "UNREQUESTED\t100\t200\n"
        "GENE2\t3\t4\n",
    )

    counts = read_tabular_selected_genes(matrix_path, ["GENE2", "MISSING", "GENE1"])

    assert counts.index.tolist() == ["cellA", "cellB"]
    assert counts.columns.tolist() == ["GENE2", "GENE1"]
    assert counts["GENE2"].tolist() == [3.0, 4.0]
