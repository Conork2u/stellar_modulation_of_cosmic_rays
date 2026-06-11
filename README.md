# Stellar Modulation of Cosmic Rays

Scientific Python research code for modelling how stellar-wind environments modulate cosmic-ray spectra around F, G, K, and M stars, with a focus on exoplanets in habitable zones.

The repository provides the public, citable home for the Stellar Bubbles / stellar modulation project associated with thesis and IAU proceedings work on cosmic-ray induced radiation in exoplanetary environments.

## Research Context

This code supports the research project:

- Conor Kennedy, *Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments*, thesis supervised by Dr. Robert Brose, Dublin City University, 2024.
- Conor Kennedy and Robert Brose, "Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments", *Proceedings of the International Astronomical Union*, Symposium 387, 2026. DOI: [10.1017/S1743921324003296](https://doi.org/10.1017/S1743921324003296).

The public IAU article is available through Cambridge Core: [Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C6A093A430EEA36A668661F8E1F747D3/S1743921324003296a.pdf/investigating_cosmicray_induced_radiation_in_exoplanetary_environments.pdf).

## Highlights

- Models stellar-wind modulation for representative F, G, K, M, Solar, and selected nearby systems.
- Computes habitable-zone boundaries using the Kopparapu et al. prescriptions.
- Uses a tested Python package layer for reusable calculations and output parsing.
- Documents external cosmic-ray spectra inputs without redistributing unclear-license files.
- Includes CI and tests for the public, import-safe code.

## Repository Contents

- `stellar_bubbles/`: reusable Python package code.
- `tests/`: pytest coverage for habitable-zone calculations, simulation-output parsing, and repository checks.
- `docs/reproducibility.md`: setup and reproduction notes.
- `docs/third_party_data.md`: notes on the external `CRspectra.py` / `CRdata` dependency.
- `CITATION.cff`: citation metadata for the associated IAU proceedings paper.

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -e ".[dev]"
```

On macOS/Linux, use `source .venv/bin/activate` instead of the Windows activation command.

## Tests

```bash
pytest
```

The test suite covers reusable calculations and repository checks. Full model runs require the external CR spectra files described in [docs/third_party_data.md](docs/third_party_data.md).

## Full Workflow

The full model workflow depends on external CR spectra files that are not redistributed here. If you have permission to use those files, place them locally as described in [docs/third_party_data.md](docs/third_party_data.md), then follow [docs/reproducibility.md](docs/reproducibility.md).

## References

- Kopparapu et al. (2013), "Habitable Zones Around Main-Sequence Stars: New Estimates", *The Astrophysical Journal*, 765, 131.
- Kopparapu et al. (2014), "Habitable Zones Around Main-Sequence Stars: Dependence on Planetary Mass", *The Astrophysical Journal Letters*, 787, L29.
- Kennedy and Brose (2026), "Investigating Cosmic-Ray Induced Radiation in Exoplanetary Environments", *Proceedings of the International Astronomical Union*, Symposium 387.
