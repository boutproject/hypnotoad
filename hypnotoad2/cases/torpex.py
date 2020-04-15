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

from collections import OrderedDict
import warnings

import numpy

from ..core.mesh import BoutMesh
from ..core.equilibrium import (
    setDefault,
    Equilibrium,
    PsiContour,
    Point2D,
    EquilibriumRegion,
    SolutionError,
)
from ..utils.hypnotoad_options import (
    HypnotoadOptions,
    HypnotoadInternalOptions,
    optionsTableString,
)

# type for manipulating information about magnetic field coils
from collections import namedtuple

Coil = namedtuple("Coil", "R, Z, I")


class TORPEXMagneticField(Equilibrium):
    """
    Magnetic configuration defined by coil positions and currents for the TORPEX device
    """

    # TORPEX wall is a circle radius 0.2 m around (1 m, 0 m)
    awall = 0.2
    Rcentre = 1.0
    Zcentre = 0.0

    # Bounding box
    Rmin = 0.8
    Rmax = 1.2
    Zmin = -0.2
    Zmax = 0.2

    # Add TORPEX-specific options and default values
    user_options = HypnotoadOptions.add(
        nx_core=None,
        nx_sol=None,
        ny_inner_lower_divertor=None,
        ny_inner_upper_divertor=None,
        ny_outer_upper_divertor=None,
        ny_outer_lower_divertor=None,
        psi_core=None,
        psi_sol=None,
        psi_sol_inner=None,
        psi_pf=None,
        psi_pf_lower=None,
        psi_pf_upper=None,
        saddle_point_p1=[0.85, -0.15],
        saddle_point_p2=[0.85, 0.15],
    )

    def __init__(self, equilibOptions, meshOptions, **kwargs):

        self.equilibOptions = equilibOptions

        # Set up options read from user input
        self.user_options = TORPEXMagneticField.user_options

        # Set sensible defaults for options
        self.user_options.set(
            xpoint_poloidal_spacing_length=5.0e-2,
            nonorthogonal_xpoint_poloidal_spacing_length=5.0e-2,
            follow_perpendicular_rtol=2.0e-8,
            follow_perpendicular_atol=1.0e-8,
            refine_width=1.0e-5,
            refine_atol=2.0e-8,
        )

        default_options = self.user_options.copy()
        self.user_options.set(**meshOptions)
        self.user_options = self.user_options.push(kwargs)

        setDefault(self.user_options, "psi_pf", self.user_options.psi_core)
        setDefault(self.user_options, "psi_pf_lower", self.user_options.psi_pf)
        setDefault(self.user_options, "psi_pf_upper", self.user_options.psi_pf)

        setDefault(self.user_options, "psi_sol_inner", self.user_options.psi_sol)

        setDefault(
            self.user_options,
            "poloidal_spacing_delta_psi",
            numpy.abs((self.user_options.psi_core - self.user_options.psi_sol) / 20.0),
        )

        print(optionsTableString(self.user_options, default_options))

        # Call Equilibrium constructor after adding stuff to options
        super().__init__(**kwargs)

        if "Coils" in self.equilibOptions:
            self.coils = [Coil(**c) for c in self.equilibOptions["Coils"]]

            self.magneticFunctionsFromCoils()

            self.Bt_axis = self.equilibOptions["Bt_axis"]
        elif "gfile" in self.equilibOptions:
            from ..utils.dct_interpolation import DCT_2D

            # load a g-file
            try:
                from pyEquilibrium.geqdsk import Geqdsk

                gfile = Geqdsk(self.equilibOptions["gfile"])
            except AttributeError:
                from boututils.geqdsk import Geqdsk

                gfile = Geqdsk()
                gfile.openFile(self.equilibOptions["gfile"])

            R = numpy.linspace(
                gfile["rleft"], gfile["rleft"] + gfile["rdim"], gfile["nw"]
            )
            Z = numpy.linspace(
                gfile["zmid"] - 0.5 * gfile["zdim"],
                gfile["zmid"] + 0.5 * gfile["zdim"],
                gfile["nh"],
            )
            psirz = gfile["psirz"]

            # check sign of psirz is consistent with signs of psi_axis, psi_bndry and
            # plasma current
            try:
                # called siaxis by boututils.geqdsk.Geqdsk
                psi_axis = gfile["siaxis"]
            except KeyError:
                # called simag by pyEquilibrium.geqdsk.Geqdsk
                psi_axis = gfile["simag"]
            R_axis = gfile["rmaxis"]
            Z_axis = gfile["zmaxis"]
            psi_bndry = gfile["sibry"]
            Ip = gfile["current"]
            if psi_axis < psi_bndry:
                # psi increases outward radially, so Bp is clockwise in the poloidal plane
                # (outward in major radius at the top of the torus, inward at the bottom).
                # This corresponds to plasma current in the anti-clockwise toroidal
                # direction looking from above, so current should be positive
                assert (
                    Ip >= 0.0
                ), "direction of plasma current should be anti-clockwise to be consistent with sign of grad(psi)"
            else:
                # psi decreases outward radially, so current should be in the opposite
                # direction
                assert (
                    Ip <= 0.0
                ), "direction of plasma current should be clockwise to be consistent with sign of grad(psi)"
            # index of a point close to the magnetic axis
            i_axis = numpy.searchsorted(R, R_axis)
            j_axis = numpy.searchsorted(Z, Z_axis)

            # Approximate value of psi at the magnetic axis from the psi array
            grid_psi_axis = psirz[j_axis, i_axis]
            if numpy.abs(psi_axis) > 1.0e-10:
                if numpy.abs(psi_axis + grid_psi_axis) / numpy.abs(psi_axis) < 1.0e-2:
                    # In some EFIT files, psirz might be defined with the 'wrong' sign, so
                    # that it increases radially, in which case we would have
                    # psi_axis ~ -grid_psi_axis and we need to flip the sign of psirz
                    psirz = -psirz
                else:
                    if (
                        numpy.abs(psi_axis - grid_psi_axis) / numpy.abs(psi_axis)
                        < 1.0e-2
                    ):
                        warnings.warn(
                            "psi_axis is not consistent with the value of psirz at (rmaxis, zmaxis)"
                        )
            else:
                # psi_axis is zero
                if numpy.abs(grid_psi_axis) < 1.0e-10:
                    # [i_axis, j_axis] must have been right on the magnetic axis, so go
                    # one point away to get a non-zero value
                    test_psi = psirz[j_axis + 1, i_axis + 1]
                else:
                    test_psi = grid_psi_axis

                if (test_psi > 0.0 and psi_bndry < 0.0) or (
                    test_psi < 0.0 and psi_bndry > 0.0
                ):
                    # psi is changing (from 0) in the 'wrong' direction away from the
                    # axis. Should be going towards psi_bndry. So flip the sign
                    psirz = -psirz

            self.magneticFunctionsFromGrid(R, Z, psirz)

            self.Bt_axis = gfile["bcentr"]
        elif "matfile" in self.equilibOptions:
            # Loading directly from the TORPEX-provided matlab file should be slightly
            # more accurate than going via a g-file because g-files don't save full
            # double-precision
            from ..utils.dct_interpolation import DCT_2D

            from scipy.io import loadmat

            eqfile = loadmat(self.equilibOptions["matfile"])["eq"]

            R = eqfile["R"][0, 0]
            Z = eqfile["Z"][0, 0]
            # TORPEX psi uses different sign convention from us
            psi = -eqfile["psi"][0, 0]

            extra = 0.04
            Rinds = (R[0, :] >= self.Rmin - extra) * (R[0, :] <= self.Rmax + extra)
            Zinds = (Z[:, 0] >= self.Zmin - extra) * (Z[:, 0] <= self.Zmax + extra)

            R = R[:, Rinds]
            R = R[Zinds, :]
            Z = Z[:, Rinds]
            Z = Z[Zinds, :]
            psi = psi[:, Rinds]
            psi = psi[Zinds, :]

            self.magneticFunctionsFromGrid(R[0, :], Z[:, 0], psi)

            Bt = eqfile["Bphi"][0, 0]
            RindMid = Bt.shape[1] // 2
            ZindMid = Bt.shape[0] // 2
            assert eqfile["R"][0, 0][ZindMid, RindMid] == 1.0
            assert eqfile["Z"][0, 0][ZindMid, RindMid] == 0.0
            self.Bt_axis = Bt[ZindMid, RindMid]
        else:
            raise ValueError("Failed to initialise psi function from inputs")

        # TORPEX plasma pressure so low fpol is constant
        self.fpol = lambda psi: self.Bt_axis / self.Rcentre
        self.fpolprime = lambda psi: 0.0

        # Make a set of points representing the wall
        self.wall = [
            self.TORPEX_wall(theta)
            for theta in numpy.linspace(0.0, 2.0 * numpy.pi, 100, endpoint=False)
        ]

        try:
            self.x_points = [
                self.findSaddlePoint(
                    Point2D(*(self.user_options.saddle_point_p1)),
                    Point2D(*(self.user_options.saddle_point_p2)),
                )
            ]
            self.psi_sep = [self.psi(*self.x_points[0])]
        except SolutionError:
            warnings.warn(
                "Warning: failed to find X-point. Equilibrium generation will fail"
            )

    def TORPEX_wall(self, theta):
        """
        Return the location of the TORPEX wall parameterized by the angle theta
        anticlockwise around the centre of the vacuum vessel
        """
        return Point2D(
            self.Rcentre + self.awall * numpy.cos(theta),
            self.Zcentre + self.awall * numpy.sin(theta),
        )

    def addWallToPlot(self, npoints=None):
        if npoints is not None:
            theta = numpy.linspace(0.0, 2.0 * numpy.pi, npoints + 1)
            pyplot.plot(*self.TORPEX_wall(theta))
        else:
            R = [p.R for p in self.wall]
            Z = [p.Z for p in self.wall]
            R.append(R[0])
            Z.append(Z[0])
            pyplot.plot(R, Z)

    def magneticFunctionsFromCoils(self):
        """
        Calculate the poloidal magnetic flux function psi = -R*A_phi, where A_phi is the
        toroidal (anti-clockwise) component of magnetic vector potential due to coils.
        See for example http://physics.usask.ca/~hirose/p812/notes/Ch3.pdf

        The currents in the coils are taken to be positive in the anti-clockwise direction
        here.

        Note e_R x e_phi = e_Z

        A radially increasing psi results in Bp going clockwise in the poloidal plane.
        """
        import sympy
        from sympy.functions.special.elliptic_integrals import elliptic_k, elliptic_e
        import scipy.special

        R, Z = sympy.symbols("R Z")
        mu0 = 4.0e-7 * sympy.pi

        A_phi = 0 * R

        for coil in self.coils:
            # little-r is the vector position from the centre of the coil to (R,Z)
            # sinTheta is the angle between r and the axis through the centre of the coil
            rSquared = R ** 2 + (Z - coil.Z) ** 2
            r = sympy.sqrt(rSquared)
            sinTheta = R / r
            kSquared = (
                4
                * coil.R
                * r
                * sinTheta
                / (rSquared + coil.R ** 2 + 2 * coil.R * r * sinTheta)
            )
            A_phi += (
                coil.I
                * coil.R
                / sympy.sqrt(r ** 2 + coil.R ** 2 + 2 * coil.R * r * sinTheta)
                / kSquared
                * ((2 - kSquared) * elliptic_k(kSquared) - 2 * elliptic_e(kSquared))
            )

        # multiply by costant pre-factor
        A_phi *= mu0 / sympy.pi

        psi = -R * A_phi
        dpsidR = sympy.diff(psi, R)
        dpsidZ = sympy.diff(psi, Z)
        modGradpsiSquared = dpsidR ** 2 + dpsidZ ** 2
        B_R = dpsidZ / R
        B_Z = -dpsidR / R
        d2psidR2 = sympy.diff(psi, R, R)
        d2psidZ2 = sympy.diff(psi, Z, Z)
        d2psidRdZ = sympy.diff(psi, R, Z)

        self.psi = sympy.lambdify(
            [R, Z],
            psi,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.f_R = sympy.lambdify(
            [R, Z],
            dpsidR / modGradpsiSquared,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.f_Z = sympy.lambdify(
            [R, Z],
            dpsidZ / modGradpsiSquared,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.Bp_R = sympy.lambdify(
            [R, Z],
            B_R,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.Bp_Z = sympy.lambdify(
            [R, Z],
            B_Z,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.d2psidR2 = sympy.lambdify(
            [R, Z],
            d2psidR2,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.d2psidZ2 = sympy.lambdify(
            [R, Z],
            d2psidZ2,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )
        self.d2psidRdZ = sympy.lambdify(
            [R, Z],
            d2psidRdZ,
            modules=[
                "numpy",
                {
                    "elliptic_k": scipy.special.ellipk,
                    "elliptic_e": scipy.special.ellipe,
                },
            ],
        )

    def makeRegions(self, npoints=100):
        """
        Find the separatrix and create the regions to grid.

        For TORPEX, follow 4 legs away from the x-point, starting with a rough guess and
        then refining to the separatrix value of A_toroidal.
        """
        assert len(self.x_points) == 1, "should be one X-point for TORPEX configuration"
        xpoint = self.x_points[0]

        boundary = self.findRoots_1d(
            lambda s: self.psi(*self.wallPosition(s)) - self.psi_sep[0], 4, 0.0, 1.0
        )

        # put lower left leg first in list, go clockwise
        boundary = boundary[2::-1] + [boundary[3]]

        legnames = [
            "inner_lower_divertor",
            "inner_upper_divertor",
            "outer_upper_divertor",
            "outer_lower_divertor",
        ]
        kinds = ["wall.X", "X.wall", "wall.X", "X.wall"]

        # create input options for EquilibriumRegions
        legoptions = {}
        for i, name in enumerate(legnames):
            options = {}
            options["nx"] = [self.user_options.nx_core, self.user_options.nx_sol]
            options["ny"] = self.user_options["ny_" + name]
            options["kind"] = kinds[i]
            legoptions[name] = options

        # set hard-wired poloidal grid spacing options
        ny_total = sum(opt["ny"] for opt in legoptions.values())

        setDefault(self.options, "N_norm", ny_total)
        self.regions = OrderedDict()
        wall_vectors = OrderedDict()
        s = numpy.linspace(10.0 * PsiContour.options.refine_atol, 1.0, npoints)
        for i, boundary_position in enumerate(boundary):
            name = legnames[i]
            boundaryPoint = self.wallPosition(boundary_position)
            legR = xpoint.R + s * (boundaryPoint.R - xpoint.R)
            legZ = xpoint.Z + s * (boundaryPoint.Z - xpoint.Z)
            leg = EquilibriumRegion(
                self,
                legnames[i],
                2,
                self.user_options,
                self.options.push(legoptions[name]),
                [Point2D(R, Z) for R, Z in zip(legR, legZ)],
                self.psi,
                self.psi_sep[0],
            )
            self.regions[name] = leg.getRefined()
            wall_vectors[name] = self.wallVector(boundary_position)

        # Make the SeparatrixContours go around clockwise
        # Record the x-point position
        # Record the psi-values of segment boundaries
        # Record the desired radial grid spacing dpsidi at internal boundaries

        dpsidi_sep_inner = (
            self.user_options.psi_sol_inner - self.psi_sep[0]
        ) / self.user_options.nx_sol
        dpsidi_sep_outer = (
            self.user_options.psi_sol - self.psi_sep[0]
        ) / self.user_options.nx_sol
        dpsidi_sep_lower = (
            self.psi_sep[0] - self.user_options.psi_pf_lower
        ) / self.user_options.nx_core
        dpsidi_sep_upper = (
            self.psi_sep[0] - self.user_options.psi_pf_upper
        ) / self.user_options.nx_core
        if self.user_options.psi_pf_lower < self.user_options.psi_sol:
            dpsidi_sep = min(
                dpsidi_sep_inner, dpsidi_sep_outer, dpsidi_sep_lower, dpsidi_sep_upper
            )
        else:
            dpsidi_sep = max(
                dpsidi_sep_inner, dpsidi_sep_outer, dpsidi_sep_lower, dpsidi_sep_upper
            )

        # decrease (assuming the factor is <1) the spacing around the separatrix by the
        # factor psi_spacing_separatrix_multiplier
        if self.user_options.psi_spacing_separatrix_multiplier is not None:
            dpsidi_sep = (
                self.user_options.psi_spacing_separatrix_multiplier * dpsidi_sep
            )

        # lower PF
        lower_psi_func = self.getPolynomialGridFunc(
            self.user_options.nx_core,
            self.user_options.psi_pf_lower,
            self.psi_sep[0],
            grad_upper=dpsidi_sep,
        )
        lower_psi_vals = self.make1dGrid(self.user_options.nx_core, lower_psi_func)

        # upper PF
        upper_psi_func = self.getPolynomialGridFunc(
            self.user_options.nx_core,
            self.user_options.psi_pf_upper,
            self.psi_sep[0],
            grad_upper=dpsidi_sep,
        )
        upper_psi_vals = self.make1dGrid(self.user_options.nx_core, upper_psi_func)

        # inner SOL
        inner_psi_func = self.getPolynomialGridFunc(
            self.user_options.nx_sol,
            self.psi_sep[0],
            self.user_options.psi_sol_inner,
            grad_lower=dpsidi_sep,
        )
        inner_psi_vals = self.make1dGrid(self.user_options.nx_sol, inner_psi_func)

        # outer SOL
        outer_psi_func = self.getPolynomialGridFunc(
            self.user_options.nx_sol,
            self.psi_sep[0],
            self.user_options.psi_sol,
            grad_lower=dpsidi_sep,
        )
        outer_psi_vals = self.make1dGrid(self.user_options.nx_sol, outer_psi_func)

        def setupRegion(name, psi_vals1, psi_vals2, reverse):
            r = self.regions[name]
            r.psi_vals = [psi_vals1, psi_vals2]
            r.separatrix_radial_index = 1
            if reverse:
                r.reverse()
                r.xPointsAtEnd[1] = xpoint
                r.wallSurfaceAtStart = wall_vectors[name]
            else:
                r.xPointsAtStart[1] = xpoint
                r.wallSurfaceAtEnd = wall_vectors[name]

        setupRegion("inner_lower_divertor", lower_psi_vals, inner_psi_vals, True)
        setupRegion("inner_upper_divertor", upper_psi_vals, inner_psi_vals, False)
        setupRegion("outer_upper_divertor", upper_psi_vals, outer_psi_vals, True)
        setupRegion("outer_lower_divertor", lower_psi_vals, outer_psi_vals, False)

        # inner lower PF -> outer lower PF
        self.makeConnection("inner_lower_divertor", 0, "outer_lower_divertor", 0)

        # inner lower SOL -> inner upper SOL
        self.makeConnection("inner_lower_divertor", 1, "inner_upper_divertor", 1)

        # outer upper PF -> inner upper PF
        self.makeConnection("outer_upper_divertor", 0, "inner_upper_divertor", 0)

        # outer upper SOL -> outer lower SOL
        self.makeConnection("outer_upper_divertor", 1, "outer_lower_divertor", 1)


def parseInput(filename):
    import yaml

    with open(filename, "r") as inputfile:
        inputs = yaml.safe_load(inputfile)
    mesh_inputs = inputs["Mesh"]
    del inputs["Mesh"]
    equilib_inputs = inputs
    if "Coils" in equilib_inputs:
        print("Coils:", equilib_inputs["Coils"])
    elif "gfile" in equilib_inputs:
        print("gfile:", equilib_inputs["gfile"])

    return equilib_inputs, mesh_inputs


def createMesh(filename, **kwargs):
    # parse input file
    equilibOptions, meshOptions = parseInput(filename)

    equilibrium = TORPEXMagneticField(equilibOptions, meshOptions, **kwargs)

    print("X-point", equilibrium.x_points[0], "with psi=" + str(equilibrium.psi_sep[0]))

    equilibrium.makeRegions()

    return BoutMesh(equilibrium)


def createEqdsk(
    equilib,
    *,
    nR=None,
    Rmin=None,
    Rmax=None,
    nZ=None,
    Zmin=None,
    Zmax=None,
    filename="torpex_test.g"
):
    from pyEquilibrium.geqdsk import Geqdsk

    R = numpy.linspace(Rmin, Rmax, nR)[numpy.newaxis, :]
    Z = numpy.linspace(Zmin, Zmax, nZ)[:, numpy.newaxis]

    gout = Geqdsk()
    gout.set("nw", nR)
    gout.set("nh", nZ)
    gout.set("rdim", Rmax - Rmin)
    gout.set("zdim", Zmax - Zmin)
    gout.set("rcentr", 0.5 * (Rmax - Rmin))
    gout.set("rleft", Rmin)
    gout.set("zmid", 0.5 * (Zmax + Zmin))
    gout.set("rmaxis", 1.0)
    gout.set("zmaxis", 0.0)
    # these values very arbitrary as don't have a magnetic axis
    gout.set("simag", equilib.psi(1.0, Zmax))
    gout.set("sibry", equilib.psi_sep[0])
    gout.set("bcentr", equilib.fpol(0.0) / 1.0)
    gout.set("current", 0.0)
    gout.set("xdum", 0.0)

    gout.set(
        "fpol", equilib.fpol(0.0) * numpy.ones(nR)
    )  # works for TORPEX because we assume fpol is constant - plasma response neglected
    gout.set("pres", numpy.zeros(nR))
    gout.set("ffprime", numpy.zeros(nR))
    gout.set("pprime", numpy.zeros(nR))
    gout.set("psirz", equilib.psi(R, Z))

    gout.set("rbbbs", [Rmin, Rmax, Rmax, Rmin])
    gout.set("zbbbs", [Zmin, Zmin, Zmax, Zmax])

    theta = numpy.linspace(0.0, 2.0 * numpy.pi, 100, endpoint=False)
    gout.set("rlim", [equilib.TORPEX_wall(t).R for t in theta])
    gout.set("zlim", [equilib.TORPEX_wall(t).Z for t in theta])

    gout.dump(filename)
