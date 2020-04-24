# Copyright 2019 J.T. Omotani
#
# Contact John Omotani john.omotani@ukaea.uk
#
# This file is part of Hypnotoad 2.
#
# Hypnotoad 2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Hypnotoad 2 is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# Hypnotoad 2.  If not, see <http://www.gnu.org/licenses/>.

import pytest

import numpy

from .utils_for_tests import tight_approx
from hypnotoad.utils.dct_interpolation import DCT_2D


def test_DCT_2D():
    def f(R, Z):
        return (R - 0.5) ** 2 - (Z - 0.1) ** 2

    nR = 10
    nZ = 15
    Rmin = 0.2
    Rmax = 1.2
    Zmin = -0.3
    Zmax = 0.4

    R_array = numpy.linspace(Rmin, Rmax, nR)
    Z_array = numpy.linspace(Zmin, Zmax, nZ)

    f_array = f(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])

    f_dct = DCT_2D(R_array, Z_array, f_array)

    # check the input array values are correctly reproduced
    f_reconstructed = f_dct(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])
    assert f_reconstructed == tight_approx(f_array)

    # check a few random points
    R = numpy.array([0.2784230, 0.357892, 0.578237, 0.732580, 1.1326794])[
        numpy.newaxis, :
    ]
    Z = numpy.array([-0.232123, -0.178594, -0.053789, 0.172530, 0.375072])[
        :, numpy.newaxis
    ]
    # can't use tight tolerance because interpolation does not reproduce exactly the
    # values away from the grid points
    assert f_dct(R, Z) == pytest.approx(f(R, Z), abs=1.0e-2)


def test_DCT_2D_ddR():
    # check R-derivative
    def f(R, Z):
        return (R - 0.5) ** 2 - (Z - 0.1) ** 2

    def dfdR(R, Z):
        return 2.0 * (R - 0.5) + 0.0 * Z

    nR = 60
    nZ = 23
    Rmin = 0.2
    Rmax = 1.2
    Zmin = -0.3
    Zmax = 0.4

    R_array = numpy.linspace(Rmin, Rmax, nR)
    Z_array = numpy.linspace(Zmin, Zmax, nZ)

    f_array = f(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])
    dfdR_array = dfdR(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])

    f_dct = DCT_2D(R_array, Z_array, f_array)

    # check on the input array
    dfdR_reconstructed = f_dct.ddR(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])
    # exclude edge points because gradient reconstruction is poor there
    assert dfdR_reconstructed[:, nR // 10 : -nR // 10] == pytest.approx(
        dfdR_array[:, nR // 10 : -nR // 10], abs=1.0e-2
    )

    # check a few random points
    R = numpy.array([0.2784230, 0.357892, 0.578237, 0.732580, 1.0326794])[
        numpy.newaxis, :
    ]
    Z = numpy.array([-0.232123, -0.178594, -0.053789, 0.172530, 0.375072])[
        :, numpy.newaxis
    ]
    # can't use tight tolerance because interpolation does not reproduce exactly the
    # values away from the grid points
    assert f_dct.ddR(R, Z) == pytest.approx(dfdR(R, Z), abs=1.0e-2)


def test_DCT_2D_ddZ():
    # check Z-derivative
    def f(R, Z):
        return (R - 0.5) ** 2 - (Z - 0.1) ** 2

    def dfdZ(R, Z):
        return 0.0 * R - 2.0 * (Z - 0.1)

    nR = 11
    nZ = 40
    Rmin = 0.2
    Rmax = 1.2
    Zmin = -0.3
    Zmax = 0.4

    R_array = numpy.linspace(Rmin, Rmax, nR)
    Z_array = numpy.linspace(Zmin, Zmax, nZ)

    f_array = f(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])
    dfdZ_array = dfdZ(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])

    f_dct = DCT_2D(R_array, Z_array, f_array)

    # check on the input array
    dfdZ_reconstructed = f_dct.ddZ(R_array[numpy.newaxis, :], Z_array[:, numpy.newaxis])
    # exclude edge points because gradient reconstruction is poor there
    assert dfdZ_reconstructed[nZ // 10 : -nZ // 10] == pytest.approx(
        dfdZ_array[nZ // 10 : -nZ // 10], abs=1.0e-2
    )

    # check a few random points
    R = numpy.array([0.2784230, 0.357892, 0.578237, 0.732580, 1.1326794])[
        numpy.newaxis, :
    ]
    Z = numpy.array([-0.222123, -0.178594, -0.053789, 0.172530, 0.275072])[
        :, numpy.newaxis
    ]
    # can't use tight tolerance because interpolation does not reproduce exactly the
    # values away from the grid points
    assert f_dct.ddZ(R, Z) == pytest.approx(dfdZ(R, Z), abs=1.0e-2)
