import numpy as np
import numpy.testing as npt
import pytest

from hypnotoad.cases.circular import CircularEquilibrium
from hypnotoad.core.mesh import BoutMesh


class TestCircular:
    R0 = 2.3
    B0 = 3.2
    r_inner = 0.1
    r_outer = 1.4
    q = 4.1
    test_settings = {
        "R0": R0,
        "B0": B0,
        "poloidal_spacing_method": "linear",
        "q_coefficients": [q],
        "r_inner": r_inner,
        "r_outer": r_outer,
        "number_of_processors": 20,
    }
    rtol = 2.0e-5
    atol = 1.0e-6
    # Fine, Cartesian grid of R and Z used for testing Equilibrium methods
    R_cartesian = (
        np.linspace(R0 - r_outer, R0 + r_outer, 3200)[:, np.newaxis]
        * np.ones(3210)[np.newaxis, :]
    )
    dR = R_cartesian[1, 0] - R_cartesian[0, 0]
    Z_cartesian = (
        np.linspace(-r_outer, +r_outer, 3210)[np.newaxis, :]
        * np.ones(3200)[:, np.newaxis]
    )
    dZ = Z_cartesian[0, 1] - Z_cartesian[0, 0]
    DDR_slice = (slice(1, -1), slice(None))
    DDZ_slice = (slice(None), slice(1, -1))
    D2DRDZ_slice = (slice(1, -1), slice(1, -1))
    r_cartesian = np.sqrt((R_cartesian - R0) ** 2 + Z_cartesian ** 2)
    theta_cartesian = np.arctan2(Z_cartesian, R_cartesian - R0)
    epsilon_cartesian = r_cartesian / R0
    qbar_cartesian = q * np.sqrt(1.0 - epsilon_cartesian ** 2)
    Bp_cartesian = B0 * r_cartesian / (qbar_cartesian * R_cartesian)
    Bt_cartesian = B0 * R0 / R_cartesian
    BR_cartesian = Bp_cartesian * np.sin(theta_cartesian)
    BZ_cartesian = -Bp_cartesian * np.cos(theta_cartesian)
    B2_cartesian = Bt_cartesian ** 2 + Bp_cartesian ** 2
    B_cartesian = np.sqrt(B2_cartesian)

    @pytest.fixture
    def equilib(self):
        """
        CircularEquilibrium object used to test methods of Equilibrium
        """
        return CircularEquilibrium(self.test_settings)

    def get_mesh(self, settings):
        """
        Create BoutMesh to test geometric quantities created by BoutMesh
        """
        equilib = CircularEquilibrium(settings)

        mesh = BoutMesh(equilib, settings)
        mesh.geometry()

        return mesh

    def DDR(self, f):
        return (f[2:, :] - f[:-2, :]) / (2.0 * self.dR)

    def D2DR2(self, f):
        return (f[2:, :] - 2.0 * f[1:-1, :] + f[:-2, :]) / self.dR ** 2

    def DDZ(self, f):
        return (f[:, 2:] - f[:, :-2]) / (2.0 * self.dZ)

    def D2DZ2(self, f):
        return (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / self.dZ ** 2

    def D2DRDZ(self, f):
        return self.DDR(self.DDZ(f))

    def test_Equilibrium_Bp_R(self, equilib):
        """
        Major-radial component of the magnetic field
        """
        npt.assert_allclose(
            equilib.Bp_R(self.R_cartesian, self.Z_cartesian),
            self.BR_cartesian,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_Bp_Z(self, equilib):
        """
        Vertical component of the magnetic field
        """
        npt.assert_allclose(
            equilib.Bp_Z(self.R_cartesian, self.Z_cartesian),
            self.BZ_cartesian,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_d2psidR2(self, equilib):
        """
        d2psi/dR2
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        d2psidR2 = self.D2DR2(equilib.psi(R, Z))
        npt.assert_allclose(
            equilib.d2psidR2(R[self.DDR_slice], Z[self.DDR_slice]),
            d2psidR2,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_d2psidZ2(self, equilib):
        """
        d2psi/dZ2
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        d2psidZ2 = self.D2DZ2(equilib.psi(R, Z))
        npt.assert_allclose(
            equilib.d2psidZ2(R[self.DDZ_slice], Z[self.DDZ_slice]),
            d2psidZ2,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_d2psidRdZ(self, equilib):
        """
        d2psi/dRdZ
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        d2psidRdZ = self.D2DRDZ(equilib.psi(R, Z))
        npt.assert_allclose(
            equilib.d2psidRdZ(R[self.D2DRDZ_slice], Z[self.D2DRDZ_slice]),
            d2psidRdZ,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_Bzeta(self, equilib):
        """
        Toroidal magnetic field
        """
        npt.assert_allclose(
            equilib.Bzeta(self.R_cartesian, self.Z_cartesian),
            self.Bt_cartesian,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_B2(self, equilib):
        """
        B^2
        """
        B2 = self.Bt_cartesian ** 2 + self.Bp_cartesian ** 2
        npt.assert_allclose(
            equilib.B2(self.R_cartesian, self.Z_cartesian),
            B2,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBzetadR(self, equilib):
        """
        d(Bzeta)/dR
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBzetadR = self.DDR(self.Bt_cartesian)
        npt.assert_allclose(
            equilib.dBzetadR(R, Z)[self.DDR_slice],
            dBzetadR,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBzetadZ(self, equilib):
        """
        d(Bzeta)/dZ
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBzetadZ = self.DDZ(self.Bt_cartesian)
        npt.assert_allclose(
            equilib.dBzetadZ(R, Z)[self.DDZ_slice],
            dBzetadZ,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBRdR(self, equilib):
        """
        d(BR)/dR
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBRdR = self.DDR(self.BR_cartesian)
        npt.assert_allclose(
            equilib.dBRdR(R, Z)[self.DDR_slice],
            dBRdR,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBRdZ(self, equilib):
        """
        d(BR)/dZ
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBRdZ = self.DDZ(self.BR_cartesian)
        npt.assert_allclose(
            equilib.dBRdZ(R, Z)[self.DDZ_slice],
            dBRdZ,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBZdR(self, equilib):
        """
        d(BZ)/dR
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBZdR = self.DDR(self.BZ_cartesian)
        npt.assert_allclose(
            equilib.dBZdR(R, Z)[self.DDR_slice],
            dBZdR,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBZdZ(self, equilib):
        """
        d(BZ)/dZ
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBZdZ = self.DDZ(self.BZ_cartesian)
        npt.assert_allclose(
            equilib.dBZdZ(R, Z)[self.DDZ_slice],
            dBZdZ,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dB2dR(self, equilib):
        """
        d(B2)/dR
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dB2dR = self.DDR(self.B2_cartesian)
        npt.assert_allclose(
            equilib.dB2dR(R, Z)[self.DDR_slice],
            dB2dR,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dB2dZ(self, equilib):
        """
        d(B2)/dZ
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dB2dZ = self.DDZ(self.B2_cartesian)
        npt.assert_allclose(
            equilib.dB2dZ(R, Z)[self.DDZ_slice],
            dB2dZ,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBdR(self, equilib):
        """
        dB/dR
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBdR = self.DDR(self.B_cartesian)
        npt.assert_allclose(
            equilib.dBdR(R, Z)[self.DDR_slice],
            dBdR,
            rtol=self.rtol,
            atol=self.atol,
        )

    def test_Equilibrium_dBdZ(self, equilib):
        """
        dB/dZ
        """
        R = self.R_cartesian
        Z = self.Z_cartesian
        dBdZ = self.DDZ(self.B_cartesian)
        npt.assert_allclose(
            equilib.dBdZ(R, Z)[self.DDZ_slice],
            dBdZ,
            rtol=self.rtol,
            atol=self.atol,
        )
