from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_third_party_files_are_not_committed():
    excluded_paths = [
        ROOT / "CRspectra.py",
        ROOT / "CRdata",
        ROOT / "crs-and-exoplanets-main",
    ]

    assert not any(path.exists() for path in excluded_paths)


def test_reproducibility_docs_exist():
    required = [
        ROOT / "README.md",
        ROOT / "requirements.txt",
        ROOT / "CITATION.cff",
        ROOT / "docs" / "third_party_data.md",
        ROOT / "docs" / "reproducibility.md",
    ]

    assert all(path.exists() for path in required)
