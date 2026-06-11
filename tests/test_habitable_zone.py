import pytest

from stellar_bubbles import compute_hz_boundaries


def test_solar_habitable_zone_boundaries_are_in_expected_order():
    hz = compute_hz_boundaries(teff=5780, luminosity_solar=1.0)

    assert hz.recent_venus < hz.runaway_greenhouse
    assert hz.runaway_greenhouse < hz.max_greenhouse
    assert hz.max_greenhouse < hz.early_mars


def test_solar_boundaries_match_kopparapu_scale():
    hz = compute_hz_boundaries(teff=5780, luminosity_solar=1.0)

    assert hz.runaway_greenhouse == pytest.approx(0.95, rel=0.02)
    assert hz.max_greenhouse == pytest.approx(1.68, rel=0.02)


def test_habitable_zone_rejects_non_positive_luminosity():
    with pytest.raises(ValueError):
        compute_hz_boundaries(teff=5780, luminosity_solar=0)
