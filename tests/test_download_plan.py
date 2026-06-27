from pathlib import Path

import pytest
import yaml

from luad_niche.download import (
    add_local_download_status,
    build_byte_ranges,
    build_dataset_targets,
    build_download_plan,
    filter_records,
    parse_geo_directory_listing,
)


def write_manifest(path: Path) -> None:
    manifest = {
        "datasets": [
            {
                "accession": "GSE189357",
                "supplementary_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/",
                "priority": 1,
            },
            {
                "accession": "GSE282617",
                "supplementary_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282617/suppl/",
                "priority": 1,
            },
        ]
    }
    path.write_text(yaml.safe_dump(manifest), encoding="utf-8")


def test_build_dataset_targets_uses_manifest_and_raw_accession_folders(tmp_path):
    manifest_path = tmp_path / "datasets.yaml"
    write_manifest(manifest_path)

    targets = build_dataset_targets(manifest_path, tmp_path, dataset="GSE189357")

    assert len(targets) == 1
    assert targets[0]["accession"] == "GSE189357"
    assert targets[0]["supplementary_url"].endswith("/GSE189357/suppl/")
    assert targets[0]["raw_dir"] == tmp_path / "data" / "raw" / "GSE189357"


def test_build_dataset_targets_rejects_unknown_dataset(tmp_path):
    manifest_path = tmp_path / "datasets.yaml"
    write_manifest(manifest_path)

    with pytest.raises(ValueError, match="not found"):
        build_dataset_targets(manifest_path, tmp_path, dataset="GSE000001")


def test_parse_geo_directory_listing_extracts_file_links():
    html = """
    <html><body>
    <a href="../">Parent Directory</a>
    <a href="GSE189357_counts.csv.gz">GSE189357_counts.csv.gz</a>
    <a href="GSE189357_metadata.tsv.gz">GSE189357_metadata.tsv.gz</a>
    <a href="?C=N;O=D">Name</a>
    </body></html>
    """

    files = parse_geo_directory_listing(
        html,
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/",
    )

    assert files == [
        {
            "filename": "GSE189357_counts.csv.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_counts.csv.gz",
        },
        {
            "filename": "GSE189357_metadata.tsv.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/GSE189357_metadata.tsv.gz",
        },
    ]


def test_parse_geo_directory_listing_ignores_external_footer_links():
    html = """
    <html><body>
    <a href="GSE282617_processed_data.csv.gz">GSE282617_processed_data.csv.gz</a>
    <a href="https://www.hhs.gov/vulnerability-disclosure-policy/index.html">HHS disclosure</a>
    </body></html>
    """

    files = parse_geo_directory_listing(
        html,
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282617/suppl/",
    )

    assert files == [
        {
            "filename": "GSE282617_processed_data.csv.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282617/suppl/GSE282617_processed_data.csv.gz",
        }
    ]


def test_build_download_plan_maps_remote_files_to_local_paths(tmp_path):
    target = {
        "accession": "GSE189357",
        "raw_dir": tmp_path / "data" / "raw" / "GSE189357",
    }
    remote_files = [
        {
            "filename": "GSE189357_counts.csv.gz",
            "url": "https://example.org/GSE189357_counts.csv.gz",
        }
    ]

    plan = build_download_plan(target, remote_files)

    assert plan == [
        {
            "accession": "GSE189357",
            "filename": "GSE189357_counts.csv.gz",
            "url": "https://example.org/GSE189357_counts.csv.gz",
            "local_path": str(tmp_path / "data" / "raw" / "GSE189357" / "GSE189357_counts.csv.gz"),
        }
    ]


def test_filter_records_keeps_matching_filenames():
    records = [
        {"filename": "GSE189357_RAW.tar"},
        {"filename": "filelist.txt"},
    ]

    assert filter_records(records, include_regex=r"filelist\.txt$") == [
        {"filename": "filelist.txt"}
    ]


def test_add_local_download_status_marks_missing_partial_and_complete(tmp_path):
    missing = add_local_download_status(
        {
            "filename": "missing.tar",
            "local_path": str(tmp_path / "missing.tar"),
            "remote_bytes": 3,
        }
    )
    assert missing["local_bytes"] == 0
    assert missing["download_status"] == "missing"
    assert missing["download_complete"] is False

    partial_path = tmp_path / "partial.tar"
    partial_path.write_bytes(b"ab")
    partial = add_local_download_status(
        {
            "filename": "partial.tar",
            "local_path": str(partial_path),
            "remote_bytes": 3,
        }
    )
    assert partial["local_bytes"] == 2
    assert partial["download_status"] == "partial"
    assert partial["download_complete"] is False

    complete_path = tmp_path / "complete.tar"
    complete_path.write_bytes(b"abc")
    complete = add_local_download_status(
        {
            "filename": "complete.tar",
            "local_path": str(complete_path),
            "remote_bytes": 3,
        }
    )
    assert complete["local_bytes"] == 3
    assert complete["download_status"] == "complete"
    assert complete["download_complete"] is True


def test_build_byte_ranges_returns_inclusive_ranges():
    assert build_byte_ranges(total_bytes=10, range_chunk_size=4) == [
        (0, 0, 3),
        (1, 4, 7),
        (2, 8, 9),
    ]
