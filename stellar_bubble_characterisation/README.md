# Stellar Bubble Characterisation

This folder contains the earlier standalone scripts for characterising the basic structure of an individual stellar bubble. They are kept together because they are useful for single-system exploration without running the full simulation grid.

## Purpose

The characterisation workflow computes and plots the physical environment around one selected star:

- stellar input parameters
- habitable-zone limits
- wind and ISM pressure terms
- wind termination shock radius
- contact discontinuity radius
- radial velocity profile
- magnetic-field profile
- simple cosmic-ray modulation profile

This is distinct from the full simulation pipeline in `scripts/`, where `Bubble_model.py` is run repeatedly across stellar classes, magnetic-field values, stellar velocities, and nearby systems.

## Files

- `stellar_parameters.py`: defines constants, stellar choices, and converted cgs quantities.
- `hz_calculations.py`: calculates Kopparapu et al. habitable-zone boundaries for the selected star.
- `function_definitions.py`: defines pressure, radius, wind, magnetic-field, and modulation helper functions.
- `plotting.py`: generates diagnostic plots for the selected stellar bubble.

## Notes

These scripts reflect the earlier single-system research workflow. They are kept for transparency and for smaller exploratory calculations, while the current reproducible multi-system workflow is documented from `scripts/Bubble_model.py`.
