import numpy as np
import pytest

from stellar_bubbles import parse_output_filename, read_simulation_output


def test_parse_output_filename():
    assert parse_output_filename("G5V_B1.0_V30.0") == ("G5V", 1.0, 30.0)


def test_parse_output_filename_rejects_unknown_shape():
    with pytest.raises(ValueError):
        parse_output_filename("not-a-simulation-output")


def test_read_simulation_output_handles_small_fixture(tmp_path):
    fixture = tmp_path / "Sol_B1.0_V22.0"
    rows = np.array([
        [1.0, 10.0, 0.8],
        [2.0, 10.0, 0.7],
        [3.0, 10.0, 0.6],
        [1.0, 20.0, 0.5],
        [2.0, 20.0, 0.4],
        [3.0, 20.0, 0.3],
    ])
    np.savetxt(fixture, rows)

    output = read_simulation_output(fixture, energy_bins=2)

    assert output.radius_cm.shape == (2, 3)
    assert output.energy_erg.shape == (2, 3)
    assert output.modulation.shape == (2, 3)
    assert output.modulation.min() >= 0
