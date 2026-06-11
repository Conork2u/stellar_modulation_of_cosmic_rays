"""Habitable-zone boundary calculations from Kopparapu et al."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HabitableZoneBoundaries:
    """Habitable-zone boundaries in astronomical units."""

    recent_venus: float
    runaway_greenhouse: float
    max_greenhouse: float
    early_mars: float
    runaway_greenhouse_5_me: float
    runaway_greenhouse_01_me: float


S_EFF_SUN = (1.776, 1.107, 0.356, 0.320, 1.188, 0.99)
A = (2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4)
B = (2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8)
C = (-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12)
D = (-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15)


def _effective_flux(teff: float, index: int) -> float:
    t_star = teff - 5780.0
    return (
        S_EFF_SUN[index]
        + A[index] * t_star
        + B[index] * t_star**2
        + C[index] * t_star**3
        + D[index] * t_star**4
    )


def compute_hz_boundaries(teff: float, luminosity_solar: float) -> HabitableZoneBoundaries:
    """Return Kopparapu habitable-zone boundaries in AU.

    Args:
        teff: Stellar effective temperature in Kelvin.
        luminosity_solar: Stellar luminosity in units of solar luminosity.
    """

    if luminosity_solar <= 0:
        raise ValueError("luminosity_solar must be positive")

    distances = [
        (luminosity_solar / _effective_flux(teff, index)) ** 0.5
        for index in range(len(S_EFF_SUN))
    ]
    return HabitableZoneBoundaries(*distances)
