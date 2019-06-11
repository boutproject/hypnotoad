import pytest
from ..equilibrium import Equilibrium
from ..torpex import *
from .utils_for_tests import *

class TestTORPEX:

    @pytest.fixture
    def equilib(self):
        equilib = Equilibrium()
        equilib.psi = lambda R, Z: (R - .98)**2 - .7*(Z - .04)**2
        equilib.psi_sep = [0.]
        equilib.fpol = lambda psi: 1.
        equilib.TORPEX_wall = lambda theta: Point2D( .65 + .85*numpy.cos(theta),
                .8 + .85*numpy.sin(theta))

        return equilib

    def test_Eqdsk(self, equilib):
        testfile = 'test_torpex_gridgen.g'

        nR = 20
        Rmin = -1.
        Rmax = 1.4
        nZ = 40
        Zmin = -.4
        Zmax = 2.

        createEqdsk(equilib, nR=nR, Rmin=Rmin, Rmax=Rmax, nZ=nZ, Zmin=Zmin, Zmax=Zmax,
                filename=testfile)

        grid_equilib = TORPEXMagneticField({'gfile':testfile}, {})

        R = numpy.linspace(Rmin, Rmax, nR)[numpy.newaxis, :]
        Z = numpy.linspace(Zmin, Zmax, nZ)[:, numpy.newaxis]

        assert equilib.psi(R, Z) == pytest.approx(grid_equilib.psi(R,Z), rel=1.e-10)

        import os
        os.remove(testfile)
