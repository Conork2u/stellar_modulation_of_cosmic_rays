# Third-Party CR Spectra Inputs

`Fullcode.py` and the newer local model workflow import `CRspectra.py` and expect that module to read local files from `CRdata/`.

Those files are not included in this public repository because the current local copies appear to come from a separate `crs-and-exoplanets` codebase and do not include an explicit license in the working folder.

## Expected Local Files

If you have permission to use the external CR spectra code/data, place the files in the repository root during local runs:

```text
CRspectra.py
CRdata/
  CR-modulation_LIS
  CR-modulation_Earth
```

Some plotting workflows also expect:

```text
external/modulation_results.csv
```

`external/`, `CRspectra.py`, `CRdata/`, and `crs-and-exoplanets-main/` are ignored by Git so they are not accidentally redistributed.

## Attribution Note

The local `CRspectra.py` copy references external cosmic-ray spectra, including the Eta Carinae non-thermal process constraints from White et al. (2020), DOI `10.1051/0004-6361/201937031`.

Before publishing or vendoring this dependency, confirm the original source, authorship, and license.
