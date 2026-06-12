# Stellar Modulation of Cosmic Rays

Scientific Python code for modelling how stellar winds modulate cosmic-ray spectra around main-sequence stars and nearby exoplanet systems. The project estimates stellar wind, magnetic field, diffusion, and habitable-zone conditions, then uses those profiles to calculate how the incident cosmic-ray spectrum changes across the circumstellar environment.

This repository is the public project home for the Stellar Bubbles / stellar modulation work associated with My thesis and the IAU proceedings article on cosmic-ray induced radiation in exoplanetary environments.

## Research Context

This code supports:

- Conor Kennedy, *Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments*, thesis supervised by Dr. Robert Brose, Dublin City University, 2024.
- Conor Kennedy and Robert Brose, "Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments", *Proceedings of the International Astronomical Union*, Symposium 387, 2026. DOI: [10.1017/S1743921324003296](https://doi.org/10.1017/S1743921324003296).

The public IAU article is available through Cambridge Core: [Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C6A093A430EEA36A668661F8E1F747D3/S1743921324003296a.pdf/investigating_cosmicray_induced_radiation_in_exoplanetary_environments.pdf).

## What The Model Does

The core model combines stellar and interstellar-medium parameters to:

- estimate habitable-zone boundaries using the Kopparapu et al. prescriptions
- calculate wind ram pressure, thermal pressure, wind termination shock radius, and contact discontinuity radius
- build radial wind-velocity and magnetic-field profiles
- define energy-dependent cosmic-ray modulation and diffusion profiles
- run parameter sweeps for F5V, G5V, K5V, M5V, Solar, and selected nearby systems
- read generated simulation grids and recreate comparison plots

## Repository Layout

```text
scripts/
  Bubble_model.py                # Main stellar-wind and cosmic-ray modulation model
  Simulation_code.py             # Batch runner for stellar classes and specific systems
  Simulation_results_reader.py   # Reads generated Data/ files and plots modulation ranges
  modulation_data_plotting.py    # Plotting script for modulation-ratio summary tables

stellar_bubble_characterisation/
  Standalone earlier scripts for calculating and plotting one stellar bubble

data/
  reference/                     # Small reference tables used alongside the model
  external/                      # Local-only third-party CR inputs go here; not committed

docs/
  reproducibility.md             # Step-by-step reproduction workflow
  scripts.md                     # What each script does and how it fits into the workflow
  third_party_data.md            # CRspectra.py / CRdata requirements and licensing boundary

stellar_bubbles/
  Import-safe helper modules and tests for reusable calculations.
```

Generated simulation outputs are written to `scripts/Data/` and are intentionally ignored by Git.

## Standalone Stellar-Bubble Characterisation

The folder `stellar_bubble_characterisation/` contains the earlier standalone scripts:

- `stellar_parameters.py`
- `hz_calculations.py`
- `function_definitions.py`
- `plotting.py`

Together these scripts characterise the basic properties of one stellar bubble: stellar inputs, habitable-zone boundaries, wind pressure, termination shock radius, contact discontinuity, velocity profile, magnetic-field profile, and diagnostic plots. This is a smaller single-system workflow. It does not run the full multi-system simulation grid; that is handled by `scripts/Bubble_model.py` and `scripts/Simulation_code.py`.

## Installation

Clone the repository:

```bash
git clone https://github.com/Conork2u/stellar_modulation_of_cosmic_rays.git
cd stellar_modulation_of_cosmic_rays
```

Create and activate a Python environment:

```bash
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

On macOS/Linux, activate with:

```bash
source .venv/bin/activate
```

For development checks, install the package with test dependencies:

```bash
pip install -e ".[dev]"
pytest
```

## Third-Party CR Spectra Requirement

Full scientific reproduction requires two external inputs that are **not included in this public repository** because they are not my files to redistribute:

```text
CRspectra.py
CRdata/
  CR-modulation_LIS
  CR-modulation_Earth
```

If you have permission to use those files, place them beside the model scripts:

```text
scripts/
  Bubble_model.py
  CRspectra.py
  CRdata/
    CR-modulation_LIS
    CR-modulation_Earth
```

The model imports `CRspectra.py` directly:

```python
import CRspectra as cr
```

Without those files, the repository can still be inspected, tested, and understood, but the full cosmic-ray spectrum calculations cannot be reproduced. See [docs/third_party_data.md](docs/third_party_data.md) for details.

## Running The Model

Run commands from `scripts/` so the model can find its local reference files and write outputs to `scripts/Data/`:

```bash
cd scripts
python Bubble_model.py Sol 5780 6.8e-14 1 0 400 1 22
```

Argument order:

```text
Spectral_Class Teff_K Mass_Loss_Msun_per_year Stellar_Radius_Rsun log10_Luminosity_Lsun Wind_Speed_km_s Magnetic_Field_G Stellar_Velocity_km_s
```

Example batch run:

```bash
python Simulation_code.py
```

`Simulation_code.py` calls the current model file, `Bubble_model.py`, for Solar, F/G/K/M grid cases, and selected nearby systems.

Read and plot generated model outputs:

```bash
python Simulation_results_reader.py
```

## Reproducing The Work

1. Clone the repository and install dependencies.
2. Add the external `CRspectra.py` and `CRdata/` files under `scripts/` if you have access to them.
3. Run `python Simulation_code.py` from `scripts/` to regenerate the model grid.
4. Run `python Simulation_results_reader.py` from `scripts/` to inspect the generated `Data/` files and recreate modulation plots.
5. Use `modulation_data_plotting.py` only when the external `crs-and-exoplanets-main/modulation_results.csv` table is available.

For a more detailed workflow, see [docs/reproducibility.md](docs/reproducibility.md). For a script-by-script map, see [docs/scripts.md](docs/scripts.md).

## Notes On Included And Excluded Files

Included:

- current working simulation scripts
- standalone stellar-bubble characterisation scripts grouped under `stellar_bubble_characterisation/`
- lightweight reference tables
- documentation, citation metadata, tests, and install configuration

Excluded:

- `OldCode` / `Old Code` / `OldCode Folders`
- generated simulation output grids
- literature PDFs and personal files
- `CRspectra.py`, `CRdata/`, and `crs-and-exoplanets-main/` until their authorship and license are confirmed

## References

- Kopparapu et al. (2013), "Habitable Zones Around Main-Sequence Stars: New Estimates", *The Astrophysical Journal*, 765, 131.
- Kopparapu et al. (2014), "Habitable Zones Around Main-Sequence Stars: Dependence on Planetary Mass", *The Astrophysical Journal Letters*, 787, L29.
- Kennedy and Brose (2026), "Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments", *Proceedings of the International Astronomical Union*, Symposium 387.
