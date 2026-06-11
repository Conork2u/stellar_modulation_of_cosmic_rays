# Repository Consolidation Plan

Target public repository:

`Conork2u/stellar_modulation_of_cosmic_rays`

Reason:

- It is already public.
- It has existing stars and history.
- Its name is descriptive and aligned with the thesis/proceedings paper.
- It can become the canonical home for the current Stellar Bubbles code.

## Repositories/Folders To Consolidate

- `stellar_modulation_of_cosmic_rays`: original public GitHub repository.
- `Stellar_Bubbles_Sim`: newest local working code from `C:\Stellar Book\Stellar Bubbles\Stellar_Bubbles_Sim`.
- `Stellar-Bubble-SIM`: private GitHub repository; appears empty after connector access was granted.
- `Stellar_Bubbles_Full`: private GitHub repository; appears to contain only a placeholder README.
- `Stellar-Profiles`: private GitHub repository; early 2023 ancestor containing `function_definitions.py`, `hz_calculations.py`, `plotting.py`, and `stellar_parameters.py`.

## Proposed Structure

```text
stellar_modulation_of_cosmic_rays/
  stellar_bubbles/
    habitable_zone.py
    simulation_data.py
  tests/
  docs/
    reproducibility.md
    third_party_data.md
    consolidation_plan.md
  Fullcode.py
  README.md
  CITATION.cff
  pyproject.toml
  requirements.txt
```

## What To Do With Old Repos

After the canonical repo is published and checked:

1. Keep `stellar_modulation_of_cosmic_rays` as the main repo.
2. For duplicate public/private repos, replace their READMEs with a short redirect to the canonical repo.
3. Archive duplicates on GitHub once the redirect is in place.
4. Keep private/local archives only if they contain unpublished, private, or third-party data.

## Rule For Third-Party Material

Do not merge `CRspectra.py`, `CRdata`, or `crs-and-exoplanets-main` until the original source and license are confirmed.
