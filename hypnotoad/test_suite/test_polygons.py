from hypnotoad.utils import polygons
import numpy as np


def test_area():
    assert np.isclose(polygons.area([(0, 0), (0, 1), (1, 1), (1, 0)]), 1.0)
    assert np.isclose(polygons.area([(0, 0), (0, 1), (1, 0)]), 0.5)


def test_clockwise():
    assert not polygons.clockwise([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert polygons.clockwise([(0, 0), (0, 1), (1, 1), (1, 0)])
