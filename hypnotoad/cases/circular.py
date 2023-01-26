from collections import OrderedDict
from collections.abc import Sequence
import numpy as np
from optionsfactory import WithMeta
from optionsfactory.checks import is_positive

from ..core.equilibrium import Equilibrium, EquilibriumRegion, Point2D
from ..utils.utils import with_default


class CircularEquilibrium(Equilibrium):
    """
    Magnetic configuration with circular, concentric flux surfaces

    Geometry based on the description of large but finite aspect ratio in
    [S.  Jolliet, F. D. Halpern, J. Loizu, A. Mosetto and P. Ricci, Phys.
    Plasmas 21, 022303 (2014)].

    Can be either SOL-only with a limiter at the inboard midplane or core-only.
    """

    # Circular-geometry-specific options and default values
    user_options_factory = Equilibrium.user_options_factory.add(
        B0=WithMeta(
            1.0,
            doc="Toroidal magnetic field at the magnetic axis",
            value_type=float,
            check_all=is_positive,
        ),
        limiter=WithMeta(
            False,
            doc=(
                "Core only, no limiter (False, default) or SOL only, with a limiter "
                "(True)"
            ),
            value_type=bool,
        ),
        nx=WithMeta(
            6,
            doc="Number of radial points (including boundary points)",
            value_type=int,
            check_all=is_positive,
        ),
        ny=WithMeta(
            8,
            doc="Number of poloidal points (excluding boundary points)",
            value_type=int,
            check_all=is_positive,
        ),
        q_coefficients=WithMeta(
            [3.4567890123456789],
            doc=(
                "Coefficients determining the safety factor q. Value is a float "
                "(treated as a length-0 sequence) or a sequence of floats. The values "
                "a[i] are used as the coefficients of a polynomial series "
                "q(r) = a[0] + a[1]*r**2 + a[2]*r**4 + ..."
            ),
            value_type=(float, Sequence),
            check_all=lambda x: all(is_positive(y) for y in tuple(x)),
        ),
        R0=WithMeta(
            1.0,
            doc="Major radius of the magnetic axis",
            value_type=float,
            check_all=is_positive,
        ),
        r_inner=WithMeta(
            0.1,
            doc="Minor radius of the inner edge of the grid.",
            value_type=float,
            check_all=is_positive,
        ),
        r_outer=WithMeta(
            0.3,
            doc="Minor radius of the outer edge of the grid.",
            value_type=float,
            check_all=is_positive,
        ),
    )

    def __init__(self, settings=None, nonorthogonal_settings=None):
        """
        Create a circular flux surface geometry

        Parameters
        ----------
        settings : dict
            A dict that will be used to set non-default values of options
            (self.user_options)
        """
        if settings is None:
            settings = {}
        if nonorthogonal_settings is None:
            nonorthogonal_settings = {}

        self.user_options = self.user_options_factory.create(settings)

        self.makeWall()

        super().__init__(nonorthogonal_settings)

        print(self.user_options.as_table(), flush=True)
        if not self.user_options.orthogonal:
            print(self.nonorthogonal_options.as_table(), flush=True)

        self.poloidal_spacing_delta_psi = with_default(
            self.user_options.poloidal_spacing_delta_psi,
            np.abs(
                (
                    self.psi_r(self.user_options.r_inner)
                    - self.psi_r(self.user_options.r_outer)
                )
                / 20.0
            ),
        )

        self.regions = OrderedDict(circular=self.makeRegion())

        if self.user_options.limiter:
            self.psi_sep = [self.psi_r(self.user_options.r_inner)]
        else:
            self.psi_sep = [self.psi_r(self.user_options.r_outer)]

            # Core-only geometry: add connection so domain is periodic in y
            self.makeConnection("circular", 0, "circular", 0)

    def makeRegion(self):
        """
        Create the EquilibriumRegion to build the grid around
        """
        r0 = 0.5 * (self.user_options.r_inner + self.user_options.r_outer)

        if self.user_options.limiter:
            kind = "wall.wall"
        else:
            # No actual 'X-point' in this circular geometry, but this 'kind' will create
            # a y-periodic region
            kind = "X.X"

        n_points = 100
        thetas = np.linspace(np.pi, -np.pi, n_points)
        R0 = self.user_options.R0
        points = [
            Point2D(R0 + r0 * np.cos(theta), r0 * np.sin(theta)) for theta in thetas
        ]

        eq_reg = EquilibriumRegion(
            equilibrium=self,
            name="circular",
            nSegments=1,  # The number of radial regions
            nx=[self.user_options.nx],
            ny=self.user_options.ny,
            ny_total=self.user_options.ny,
            kind=kind,
            # The following arguments are passed through to PsiContour
            points=points,  # list of Point2D objects on the line
            psival=self.psi_r(r0),
            Rrange=(-float("inf"), float("inf")),
            Zrange=(-float("inf"), float("inf")),
        )

        eq_reg.psi_vals = [
            np.linspace(
                self.psi_r(self.user_options.r_inner),
                self.psi_r(self.user_options.r_outer),
                2 * self.user_options.nx + 1,
            )
        ]

        if self.user_options.limiter:
            eq_reg.separatrix_radial_index = 0
        else:
            eq_reg.separatrix_radial_index = -1

        return eq_reg

    def makeWall(self):
        """
        Create a wall that includes the limiter (if present)

        Wall points go around anti-clockwise.
        """
        R0 = self.user_options.R0
        rmax = 1.1 * self.user_options.r_outer
        self.Rmin = R0 - rmax
        self.Rmax = R0 + rmax
        self.Zmin = -rmax
        self.Zmax = rmax
        if self.user_options.limiter:
            self.wall = [
                Point2D(self.Rmin, self.Zmin),
                Point2D(self.Rmax, self.Zmin),
                Point2D(self.Rmax, self.Zmax),
                Point2D(self.Rmin, self.Zmax),
                Point2D(self.Rmin, 0.0),
                Point2D(R0, 0.0),
                Point2D(self.Rmin, 0.0),
            ]
        else:
            self.wall = [
                Point2D(self.Rmin, self.Zmin),
                Point2D(self.Rmax, self.Zmin),
                Point2D(self.Rmax, self.Zmax),
                Point2D(self.Rmin, self.Zmax),
            ]

        return self.wall

    def r(self, R, Z):
        """
        r = sqrt((R - R0)**2 + Z**2)
        """
        R0 = self.user_options.R0
        return np.sqrt((R - R0) ** 2 + Z**2)

    def theta(self, R, Z):
        """
        theta = arctan(Z/(R - R0))
        """
        R0 = self.user_options.R0
        return np.arctan2(Z, R - R0)

    def q(self, r):
        """
        Evaluate the safety factor as a function of minor radius
        """
        if not hasattr(self, "_q"):
            coef_list = self.user_options.q_coefficients
            exponent_list = range(0, 2 * len(coef_list), 2)

            def func(x):
                return sum(c * x**e for c, e in zip(coef_list, exponent_list))

            self._q = func

        return self._q(r)

    def dqdr(self, r):
        """
        Evaluate dq/dr
        """
        if not hasattr(self, "_dqdr"):
            coef_list = self.user_options.q_coefficients
            if len(coef_list) == 1:

                def func(x):
                    return 0.0

                self._dqdr = func
            else:
                coef_list = coef_list[1:]
                exponent_list = range(1, 2 * len(coef_list) + 1, 2)

                def func(x):
                    return sum(
                        c * e * x ** (e - 1) for c, e in zip(coef_list, exponent_list)
                    )

                self._dqdr = func

        return self._dqdr(r)

    def dpsidr_r(self, r):
        """
        dpsi/dr as a function of r
        """
        if not hasattr(self, "_dpsidr_r"):
            B0 = self.user_options.B0
            R0 = self.user_options.R0

            def func(x):
                return B0 * x / (np.sqrt(1.0 - x**2 / R0**2) * self.q(x))

            self._dpsidr_r = func

        return self._dpsidr_r(r)

    def d2psidr2_r(self, r):
        """
        d2psi/dr2 as a function of r

        ::

            d2psi/dr2 = d/dr(B0 r / (sqrt(1 - r**2 / R0**2) q))
                      = B0 / (sqrt(1 - r**2 / R0**2) q)
                        + B0 r**2 / (R0**2 (1 - r**2 / R0**2)**1.5 q)
                        - B0 r dq/dr / (sqrt(1 - r**2 / R0**2) q**2)
        """
        if not hasattr(self, "_d2psidr2_r"):
            R0 = self.user_options.R0
            B0 = self.user_options.B0

            def func(x):
                return (
                    B0 / (np.sqrt(1.0 - x**2 / R0**2) * self.q(x))
                    + B0
                    * x**2
                    / (R0**2 * (1.0 - x**2 / R0**2) ** 1.5 * self.q(x))
                    - B0
                    * x
                    * self.dqdr(x)
                    / (np.sqrt(1.0 - x**2 / R0**2) * self.q(x) ** 2)
                )

            self._d2psidr2_r = func

        return self._d2psidr2_r(r)

    def psi_r(self, r):
        """
        Evaluate the poloidal flux function psi as a function of minor radius

        Calculated by integrating psi'(r) = B0*r/(sqrt(1-r**2/R0**2)*q(r)), as given
        above eq (24) in Jolliet et al.

        Integrals done with Wolfram Alpha for the first few values of the maximum order
        in the polynomial q(r) = a0 + a1 * r**2 + a2 * r**4 + ...
        where only even powers of r are used so that q(r) is even at r=0 so there is no
        discontinuity in gradient at the magnetic axis.
        """
        if not hasattr(self, "_psi_r"):
            coef_array = np.array(self.user_options.q_coefficients)
            R0 = self.user_options.R0
            B0 = self.user_options.B0
            if len(coef_array) == 1:

                def func(x):
                    return (
                        B0
                        * R0**2
                        / coef_array[0]
                        * (1.0 - np.sqrt(1.0 - x**2 / R0**2))
                    )

                self._psi_r = func
            elif len(coef_array) == 2:
                a0, a1 = coef_array

                def func(x):
                    return (
                        B0
                        * R0
                        * np.log(
                            (
                                (1.0 + R0 * np.sqrt(a1 / (a0 + a1 * R0**2)))
                                * (
                                    -1.0
                                    + np.sqrt(
                                        (a1 * (-(x**2) + R0**2))
                                        / (a0 + a1 * R0**2)
                                    )
                                )
                            )
                            / (
                                (-1.0 + R0 * np.sqrt(a1 / (a0 + a1 * R0**2)))
                                * (
                                    1.0
                                    + np.sqrt(
                                        (a1 * (-(x**2) + R0**2))
                                        / (a0 + a1 * R0**2)
                                    )
                                )
                            )
                        )
                        / (2.0 * np.sqrt(a1 * (a0 + a1 * R0**2)))
                    )

                self._psi_r = func
            # Note, don't seem to be able to get closed-form expression for integral
            # with three q-coefficients (quartic q) from Mathematica without it being
            # written with complex numbers, even though the answer is real.
            else:
                raise ValueError(
                    f"This number of coefficients in the polynomial expansion of q not "
                    f"supported yet, got {self.user_options.q_coefficients}"
                )

        return self._psi_r(r)

    def drdR(self, R, Z):
        """
        r = sqrt((R - R0)**2 + Z**2)
        dr/dR = (R - R0) / sqrt((R - R0)**2 + Z**2)
        """
        R0 = self.user_options.R0
        return (R - R0) / np.sqrt((R - R0) ** 2 + Z**2)

    def drdZ(self, R, Z):
        """
        r = sqrt((R - R0)**2 + Z**2)
        dr/dZ = Z / sqrt((R - R0)**2 + Z**2)
        """
        R0 = self.user_options.R0
        return Z / np.sqrt((R - R0) ** 2 + Z**2)

    def dRdr(self, R, Z):
        """
        R = R0 + r cos(theta)
        dRdr = cos(theta)
        """
        return np.cos(self.theta(R, Z))

    def dZdr(self, R, Z):
        """
        Z = r sin(theta)
        dZdr = sin(theta)
        """
        return np.sin(self.theta(R, Z))

    def psi(self, R, Z):
        """
        Evaluate the poloidal flux function psi as a function of R and Z
        """
        return self.psi_r(self.r(R, Z))

    def f_R(self, R, Z):
        """
        R component of the vector :math:`\\nabla\\psi/|\\nabla\\psi|^2`.
        This is in the minor radius direction for concentric, circular flux surface
        geometry.
        """
        R0 = self.user_options.R0
        r = self.r(R, Z)
        return (R - R0) / r / self.dpsidr_r(r)

    def f_Z(self, R, Z):
        """
        Z component of the vector :math:`\\nabla\\psi/|\\nabla\\psi|^2`.
        This is in the minor radius direction for concentric, circular flux surface
        geometry.
        """
        r = self.r(R, Z)
        return Z / r / self.dpsidr_r(r)

    def Bp_R(self, R, Z):
        """
        Bp_R = dpsi/dZ / R
             = dpsi/dr dr/dZ / R
        """
        return self.dpsidr_r(self.r(R, Z)) * self.drdZ(R, Z) / R

    def Bp_Z(self, R, Z):
        """
        Bp_Z = -dpsi/dR / R
             = -dpsi/dr dr/dR / R
        """
        return -self.dpsidr_r(self.r(R, Z)) * self.drdR(R, Z) / R

    def d2psidR2(self, R, Z):
        """
        d2psi/dR2
        """
        # d2psi/dR2 = d/dR(dpsi/dR)
        #           = d/dR(dpsi/dr dr/dR)
        #           = d/dR(dpsi/dr (R - R0)/sqrt((R - R0)**2 + Z**2))
        #           = d/dR(dpsi/dr) (R - R0)/sqrt((R - R0)**2 + Z**2)
        #             + dpsi/dr d/dR((R - R0)/sqrt((R - R0)**2 + Z**2))
        #           = d2psi/dr2 dr/dR (R - R0)/sqrt((R - R0)**2 + Z**2)
        #             + dpsi/dr (1/sqrt((R - R0)**2 + Z**2)
        #                        - (R - R0)**2/((R - R0)**2 + Z**2)**1.5)
        #           = d2psi/dr2 (R - R0)**2/((R - R0)**2 + Z**2)
        #             + dpsi/dr (1/sqrt((R - R0)**2 + Z**2)
        #                        - (R - R0)**2/((R - R0)**2 + Z**2)**1.5)
        #           = d2psi/dr2 (R - R0)**2/((R - R0)**2 + Z**2)
        #             + dpsi/dr ((R - R0)**2 + Z**2 - (R - R0)**2)
        #               / ((R - R0)**2 + Z**2)**1.5
        #           = d2psi/dr2 (R - R0)**2/((R - R0)**2 + Z**2)
        #             + dpsi/dr Z**2/((R - R0)**2 + Z**2)**1.5
        R0 = self.user_options.R0
        r = self.r(R, Z)
        return (
            self.d2psidr2_r(r) * (R - R0) ** 2 / ((R - R0) ** 2 + Z**2)
            + self.dpsidr_r(r) * Z**2 / ((R - R0) ** 2 + Z**2) ** 1.5
        )

    def d2psidZ2(self, R, Z):
        """
        d2psi/dZ2
        """
        # d2psi/dZ2 = d/dZ(dpsi/dZ)
        #           = d/dZ(dpsi/dr dr/dZ)
        #           = d/dZ(dpsi/dr Z/sqrt((R - R0)**2 + Z**2))
        #           = d/dZ(dpsi/dr) Z/sqrt((R - R0)**2 + Z**2)
        #             + dpsi/dr d/dZ(Z/sqrt((R - R0)**2 + Z**2))
        #           = d2psi/dr2 dr/dZ Z/sqrt((R - R0)**2 + Z**2)
        #             + dpsi/dr (1/sqrt((R - R0)**2 + Z**2)
        #                        - Z**2/((R - R0)**2 + Z**2)**1.5)
        #           = d2psi/dr2 Z**2/((R - R0)**2 + Z**2)
        #             + dpsi/dr (1/sqrt((R - R0)**2 + Z**2)
        #                        - Z**2/((R - R0)**2 + Z**2)**1.5)
        #           = d2psi/dr2 Z**2/((R - R0)**2 + Z**2)
        #             + dpsi/dr ((R - R0)**2 + Z**2 - Z**2)/((R - R0)**2 + Z**2)**1.5
        #           = d2psi/dr2 Z**2/((R - R0)**2 + Z**2)
        #             + dpsi/dr (R - R0)**2/((R - R0)**2 + Z**2)**1.5
        R0 = self.user_options.R0
        r = self.r(R, Z)
        return (
            self.d2psidr2_r(r) * Z**2 / ((R - R0) ** 2 + Z**2)
            + self.dpsidr_r(r) * (R - R0) ** 2 / ((R - R0) ** 2 + Z**2) ** 1.5
        )

    def d2psidRdZ(self, R, Z):
        """
        d2psi/dR2
        """
        # d2psi/dRdZ = d/dR(dpsi/dZ)
        #            = d/dR(dpsi/dr dr/dZ)
        #            = d/dR(dpsi/dr Z/sqrt((R - R0)**2 + Z**2))
        #            = d/dR(dpsi/dr) Z/sqrt((R - R0)**2 + Z**2)
        #              + dpsi/dr d/dR(Z/sqrt((R - R0)**2 + Z**2))
        #            = d2psi/dr2 dr/dR Z/sqrt((R - R0)**2 + Z**2)
        #              - dpsi/dr (R - R0) Z / ((R - R0)**2 + Z**2)**1.5
        #            = d2psi/dr2 (R - R0) Z / ((R - R0)**2 + Z**2)
        #              - dpsi/dr (R - R0) Z / ((R - R0)**2 + Z**2)**1.5
        #            = (d2psi/dr2 - dpsi/dr / sqrt((R - R0)**2 + Z**2))
        #              * (R - R0) Z / ((R - R0)**2 + Z**2)
        R0 = self.user_options.R0
        r = self.r(R, Z)
        return (
            (self.d2psidr2_r(r) - self.dpsidr_r(r) / np.sqrt((R - R0) ** 2 + Z**2))
            * (R - R0)
            * Z
            / ((R - R0) ** 2 + Z**2)
        )

    def fpol(self, psi):
        """poloidal current function, returns fpol such that B_toroidal = fpol/R"""
        return self.user_options.R0 * self.user_options.B0

    def fpolprime(self, psi):
        """psi-derivative of fpol"""
        return 0.0

    @property
    def Bt_axis(self):
        """Calculate toroidal field on axis"""
        return self.user_options.B0
