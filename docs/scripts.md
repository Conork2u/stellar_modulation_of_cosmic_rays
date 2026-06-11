# Script Guide

This repository keeps the current research workflow in `scripts/`. The scripts are intentionally close to the working research code used to regenerate figures, while the repository documentation explains how to run them safely and reproducibly.

## `Bubble_model.py`

Main model entry point.

Responsibilities:

- reads one stellar/system parameter set from command-line arguments
- defines constants and stellar/ISM parameters in cgs units using `astropy`
- calculates habitable-zone boundaries using Kopparapu et al. coefficients
- computes wind ram pressure, ISM pressure, wind termination shock radius, and contact discontinuity radius
- builds radial velocity and magnetic-field profiles
- defines cosmic-ray modulation and diffusion functions across the stellar bubble
- writes generated modulation grids to `Data/` when plotting is disabled
- produces diagnostic plots when plotting sections are enabled

Key functions:

- `compute_hz_boundaries(teff, L)`: returns HZ boundary distances in AU.
- `pRAM(rho_ISM, v_star)`: ram-pressure contribution from stellar motion through the ISM.
- `pthermal(n, k_B, T)`: thermal-pressure contribution from the ISM.
- `pISM(pthermal, pRAM)`: combined ISM pressure term used in the shock-radius calculation.
- `R1(lw, vw, pISM)`: termination shock radius.
- `v(radius)`: wind velocity profile.
- `B(radius)`: magnetic-field profile.
- `N_B(radius, E)`, `N_A(radius, E)`, `N(radius, E)`: cosmic-ray modulation profiles in the different regions.
- `D_c_B(radius, E)`, `D_c_A(radius, E)`, `D_c(radius, E)`: diffusion coefficient profiles.

Full runs require `CRspectra.py` and `CRdata/` beside the script.

## `Simulation_code.py`

Batch runner for the main model.

Responsibilities:

- creates `Data/` if it does not already exist
- defines the F5V, G5V, K5V, and M5V parameter grid
- defines magnetic-field and stellar-velocity sweeps
- defines selected nearby systems: TRAPPIST-1, Proxima, and Teegarden's Star
- calls `Bubble_model.py` with each parameter set
- captures model output for inspection

The runner now calls the current model file:

```python
command = [sys.executable, "Bubble_model.py", *map(str, parameters)]
```

## `Simulation_results_reader.py`

Reader and plotting utility for generated model outputs.

Responsibilities:

- parses filenames in `Data/` to recover spectral class, magnetic field, and stellar velocity
- identifies the available parameter grid
- reads whitespace-separated model outputs
- reshapes generated radius, energy, and modulation arrays
- plots modulation ranges by spectral class
- plots fixed and varying parameter comparisons
- plots extreme modulation cases and Solar reference comparisons

Expected generated file format:

```text
radius_cm energy_erg modulation_value
```

## `Modulation_data_ploting.py`

Summary plotting script for modulation-ratio tables.

Responsibilities:

- reads `crs-and-exoplanets-main/modulation_results.csv`
- groups rows by spectral class, magnetic field, velocity, and spectral index
- compares ratios above 30 GeV and across the wider energy range
- creates publication-style error-bar plots by spectral class

This script depends on the external `crs-and-exoplanets-main` output table and should only be run when that third-party/local workflow is available.

## Data Products

`Bubble_model.py` writes generated model outputs to:

```text
scripts/Data/
```

Those files can be regenerated and are not committed to the repository.
