# Reproducibility Guide

This project has two reproducibility levels:

1. Public smoke-test reproducibility, using only files in this repository.
2. Full scientific workflow reproducibility, using the external CR spectra code/data that are not redistributed here.

## 1. Public Smoke Test

Create an environment and install the package with development dependencies:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -e ".[dev]"
```

Run tests:

```bash
pytest
```

This verifies:

- habitable-zone boundary calculations
- simulation-output filename parsing
- simulation-output grid reading
- third-party CR spectra files are not accidentally committed
- core documentation is present

## 2. Full Scientific Workflow

Install runtime dependencies:

```bash
pip install -r requirements.txt
```

Add the external CR spectra dependency locally:

```text
CRspectra.py
CRdata/
  CR-modulation_LIS
  CR-modulation_Earth
```

These files are intentionally ignored by Git. See `docs/third_party_data.md`.

Run the legacy Solar reference model from the original public code:

```bash
python Fullcode.py Sol 5780 6.8e-14 1 0 400 1 22 false
```

Expected result:

- a generated output file is created under the model output directory
- generated output contains three whitespace-separated columns:
  - radius in cm
  - energy in erg
  - modulation value

## Notes

The legacy research scripts still contain exploratory plotting sections and publication-specific plotting choices. The tested package code in `stellar_bubbles/` captures reusable pieces that can be safely imported and extended.
