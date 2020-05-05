import pytest

import numpy

from hypnotoad.core.equilibrium import Equilibrium, Point2D
from hypnotoad.cases import torpex


class ThisEquilibrium(Equilibrium):
    def __init__(self):
        self.user_options = self.user_options_factory.create({})
        super().__init__({})


class TestTORPEX:
    @pytest.fixture
    def equilib(self):
        equilib = ThisEquilibrium()
        equilib.psi = lambda R, Z: (R - 0.98) ** 2 - 0.7 * (Z - 0.04) ** 2
        equilib.psi_sep = [0.0]
        equilib.fpol = lambda psi: 1.0
        equilib.TORPEX_wall = lambda theta: Point2D(
            0.65 + 0.85 * numpy.cos(theta), 0.8 + 0.85 * numpy.sin(theta)
        )

        return equilib

    def test_Eqdsk(self, equilib):
        testfile = "test_torpex_gridgen.g"

        nR = 20
        Rmin = -1.0
        Rmax = 1.4
        nZ = 40
        Zmin = -0.4
        Zmax = 2.0

        torpex.createEqdsk(
            equilib,
            nR=nR,
            Rmin=Rmin,
            Rmax=Rmax,
            nZ=nZ,
            Zmin=Zmin,
            Zmax=Zmax,
            filename=testfile,
        )

        grid_equilib = torpex.TORPEXMagneticField(
            {"gfile": testfile}, {"psi_core": 0.0, "psi_sol": 1.0}
        )

        R = numpy.linspace(Rmin, Rmax, nR)[numpy.newaxis, :]
        Z = numpy.linspace(Zmin, Zmax, nZ)[:, numpy.newaxis]

        assert equilib.psi(R, Z) == pytest.approx(grid_equilib.psi(R, Z), rel=1.0e-9)

        import os

        os.remove(testfile)
