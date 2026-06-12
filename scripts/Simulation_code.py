# -*- coding: utf-8 -*-
"""
Batch runner for the stellar bubble model.

Runs the current Bubble_model.py script for the Solar reference case,
selected nearby systems, and the F/G/K/M parameter grid.
"""

from __future__ import annotations

import itertools
import os
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)

DATA_DIR = SCRIPT_DIR / "Data"
DATA_DIR.mkdir(exist_ok=True)


spectral_classes = [
    ("F5V", 6550, 5e-14, 1.473, 0.56, 500),
    ("G5V", 5660, 2.4e-14, 0.977, -0.05, 400),
    ("K5V", 4440, 2.99e-15, 0.701, -0.76, 300),
    ("M5V", 3060, 1.92e-16, 0.196, -2.52, 200),
]

magnetic_fields = {
    "F5V": (0.8, 2.2, 3.3),
    "G5V": (0.2, 1, 47.9),
    "K5V": (1.50, 36.5, 78),
    "M5V": (16, 51, 406),
}

stellar_velocities = (10, 20, 30, 60, 120)

specific_systems = [
    ("Trappist1", 2570, 6.9e-17, 0.114, -3.28, 200, 600, 80.81),
    ("Proxima", 2930, 3.5e-17, 0.156, -2.79, 400, 600, 25.0),
    ("Teegarden", 2680, 6.8e-16, 0.120, -3.19, 200, 1200, 92.0),
]

solar_reference = ("Sol", 5780, 6.8e-14, 1, 0, 400, 1, 22)


def grid_parameter_sets():
    """Yield the model parameters for the full stellar-class grid."""

    for stellar_params in spectral_classes:
        spectral_class = stellar_params[0]
        for b_field, v_star in itertools.product(
            magnetic_fields[spectral_class],
            stellar_velocities,
        ):
            yield stellar_params + (b_field, v_star)


def run_model(parameters):
    """Run Bubble_model.py for one parameter set."""

    command = [sys.executable, "Bubble_model.py", *map(str, parameters)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)

    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    result.check_returncode()


def main():
    """Run the full documented simulation set."""

    runs = [*specific_systems, solar_reference, *grid_parameter_sets()]
    for parameters in runs:
        print(f"Running: {parameters}")
        run_model(parameters)


if __name__ == "__main__":
    main()
