# Third-Party CR Spectra Inputs

The full cosmic-ray modulation workflow imports `CRspectra.py`:

```python
import CRspectra as cr
```

The local working copy also uses files under `CRdata/`.

These files are not included in this public repository because they are not my code/data to redistribute. They appear to come from a separate CR spectra project, and the local copy does not include a license that would make public redistribution safe.

## Required Location

If you have permission to use the external CR spectra code and data, place them here:

```text
scripts/
  CRspectra.py
  CRdata/
    CR-modulation_LIS
    CR-modulation_Earth
```

Then run the model from `scripts/`.

## Optional External Summary Table

The script `Modulation_data_ploting.py` expects:

```text
scripts/crs-and-exoplanets-main/modulation_results.csv
```

That table is also treated as an external/local dependency and is not redistributed here.

## Git Ignore Policy

The repository ignores:

```text
CRspectra.py
CRdata/
crs-and-exoplanets-main/
scripts/CRspectra.py
scripts/CRdata/
scripts/crs-and-exoplanets-main/
data/external/
```

This prevents accidental publication of third-party material while still documenting exactly where authorised users should place the files for full reproduction.

## Before Vendoring

Before adding `CRspectra.py`, `CRdata/`, or `crs-and-exoplanets-main/` to this repository, confirm:

- original author/source
- license
- redistribution permission
- citation requirements
