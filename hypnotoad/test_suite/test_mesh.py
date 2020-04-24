import numpy
import pytest
from hypnotoad.core import mesh
from .utils_for_tests import tight_approx


class TestMultiLocationArray:
    nx = 4
    ny = 5

    @pytest.fixture
    def MLArray(self):
        MLArray = mesh.MultiLocationArray(self.nx, self.ny)
        return MLArray

    def test_zero(self, MLArray):
        a = MLArray.zero()

        # don't use @property getters for the tests, as these will set the arrays to zero
        # when called, and we want to check that they are already zero.
        assert a._centre_array == tight_approx(numpy.zeros([self.nx, self.ny]))
        assert a._xlow_array == tight_approx(numpy.zeros([self.nx + 1, self.ny]))
        assert a._ylow_array == tight_approx(numpy.zeros([self.nx, self.ny + 1]))
        assert a._corners_array == tight_approx(numpy.zeros([self.nx + 1, self.ny + 1]))
