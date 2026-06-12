# Reproducibility Guide

This project has two reproducibility levels:

1. Public repository reproducibility: install the project, inspect the scripts, and run the import-safe tests.
2. Full scientific workflow reproducibility: regenerate the model outputs using external CR spectra files that are not redistributed here.

## 1. Clone And Install

```bash
git clone https://github.com/Conork2u/stellar_modulation_of_cosmic_rays.git
cd stellar_modulation_of_cosmic_rays
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

On macOS/Linux:

```bash
source .venv/bin/activate
```

For tests:

```bash
pip install -e ".[dev]"
pytest
```

## 2. Add Third-Party CR Inputs

Full runs require:

```text
scripts/
  CRspectra.py
  CRdata/
    CR-modulation_LIS
    CR-modulation_Earth
```

These files are intentionally ignored by Git because they come from a separate CR spectra code/data source. If you have access to them, copy them into `scripts/` before running the model. See [third_party_data.md](third_party_data.md).

## 3. Run A Single Model

Run from `scripts/`:

```bash
cd scripts
python Bubble_model.py Sol 5780 6.8e-14 1 0 400 1 22
```

Argument order:

```text
Spectral_Class Teff_K Mass_Loss_Msun_per_year Stellar_Radius_Rsun log10_Luminosity_Lsun Wind_Speed_km_s Magnetic_Field_G Stellar_Velocity_km_s
```

The model writes generated output to `scripts/Data/`.

## 4. Run The Batch Grid

From `scripts/`:

```bash
python Simulation_code.py
```

The batch runner calls `Bubble_model.py` for:

- Solar reference parameters
- representative F5V, G5V, K5V, and M5V stellar classes
- magnetic-field and stellar-velocity sweeps
- TRAPPIST-1, Proxima, and Teegarden's Star cases

## 5. Read And Plot Results

After `Data/` has been generated:

```bash
python Simulation_results_reader.py
```

This reads each generated model file, parses the parameter values from the filename, reshapes the radius-energy grid, and produces the comparison plots.

## 6. Modulation-Ratio Summary Plots

`modulation_data_plotting.py` requires the external table:

```text
scripts/crs-and-exoplanets-main/modulation_results.csv
```

Run it only when that local third-party workflow is available:

```bash
python modulation_data_plotting.py
```

## Output Policy

Generated outputs are not committed because they can be regenerated and can quickly make the repository noisy:

```text
scripts/Data/
Data/
```

The repository commits source code, lightweight reference files, documentation, tests, and citation metadata.
