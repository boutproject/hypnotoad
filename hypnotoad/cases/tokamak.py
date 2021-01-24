# Generate grids for tokamak configurations
#

import numpy as np
from optionsfactory import WithMeta
from optionsfactory.checks import (
    is_non_negative,
    is_positive,
    NoneType,
)
from scipy import interpolate
from scipy.integrate import solve_ivp
import warnings
from collections import OrderedDict
import functools

from ..core.equilibrium import Equilibrium, EquilibriumRegion, Point2D
from ..core.mesh import MultiLocationArray

from ..utils import critical, polygons
from ..utils.utils import with_default


class TokamakEquilibrium(Equilibrium):
    """
    Represents an axisymmetric tokamak equilibrium

    Data members
    - x_points: list of Point2D objects giving the position of the X-points ordered
                from primary X-point (nearest the core) outward
    - o_point: Point2D object for the magnetic O-point
    - psi_sep: values of psi on the separatrices ordered the same as self.x_points
    - Rmin, Rmax, Zmin, Zmax: positions of the corners of a bounding
                              box for the gridding
    - regions: OrderedDict of EquilibriumRegion objects that specify this equilibrium
    - wall: list of Point2D giving vertices of polygon representing the wall, in
            anti-clockwise order; assumed to be closed so last element and first are
            taken to be connected
    """

    # Tokamak-specific options and default values
    user_options_factory = Equilibrium.user_options_factory.add(
        nx_core=WithMeta(
            5,
            doc="Number of radial points in the core",
            value_type=int,
            check_all=is_positive,
        ),
        nx_pf=WithMeta(
            "nx_core",
            doc=(
                "Number of radial points in the PF region Note: Currently can't be "
                "varied due to BOUT++ limitations"
            ),
            value_type=int,
            check_all=is_positive,
        ),
        nx_inter_sep=WithMeta(
            0,
            doc="Number of radial points in the inter-separatrix region",
            value_type=int,
            check_all=is_non_negative,
        ),
        nx_sol=WithMeta(
            5,
            doc="Number of radial points in the SOL",
            value_type=int,
            check_all=is_positive,
        ),
        nx_sol_inner=WithMeta(
            "nx_sol",
            doc=(
                "Number of radial points in the outer SOL. Note: Currently can't be "
                "varied due to BOUT++ limitations"
            ),
            value_type=int,
            check_all=is_positive,
        ),
        nx_sol_outer=WithMeta(
            "nx_sol",
            doc=(
                "Number of radial points in the inner SOL. Note: Currently can't be "
                "varied due to BOUT++ limitations"
            ),
            value_type=int,
            check_all=is_positive,
        ),
        ny_inner_divertor=WithMeta(
            4,
            doc="Number of poloidal points in the inner divertor(s)",
            value_type=int,
            check_all=is_positive,
        ),
        ny_inner_lower_divertor=WithMeta(
            "ny_inner_divertor",
            doc="Number of poloidal points in the inner, lower divertor",
            value_type=int,
            check_all=is_positive,
        ),
        ny_inner_upper_divertor=WithMeta(
            "ny_inner_divertor",
            doc="Number of poloidal points in the inner, upper divertor",
            value_type=int,
            check_all=is_positive,
        ),
        ny_outer_divertor=WithMeta(
            4,
            doc="Number of poloidal points in the outer divertor(s)",
            value_type=int,
            check_all=is_positive,
        ),
        ny_outer_lower_divertor=WithMeta(
            "ny_outer_divertor",
            doc="Number of poloidal points in the outer, lower divertor",
            value_type=int,
            check_all=is_positive,
        ),
        ny_outer_upper_divertor=WithMeta(
            "ny_outer_divertor",
            doc="Number of poloidal points in the outer, upper divertor",
            value_type=int,
            check_all=is_positive,
        ),
        ny_sol=WithMeta(
            8,
            doc="Number of poloidal points in the SOL upstream of the X-point(s)",
            value_type=int,
            check_all=is_positive,
        ),
        ny_inner_sol=WithMeta(
            lambda options: options.ny_sol // 2,
            doc="Number of poloidal points in the inner SOL upstream of the X-point(s)",
            value_type=int,
            check_all=is_positive,
        ),
        ny_outer_sol=WithMeta(
            lambda options: options.ny_sol - options.ny_inner_sol,
            doc="Number of poloidal points in the outer SOL upstream of the X-point(s)",
            value_type=int,
            check_all=is_positive,
        ),
        psinorm_core=WithMeta(
            0.9,
            doc="Normalised psi of the inner radial (core) boundary",
            value_type=[float, int],
        ),
        psinorm_sol=WithMeta(
            1.1,
            doc="Normalised psi of the outer radial (SOL) boundary",
            value_type=[float, int],
        ),
        psinorm_sol_inner=WithMeta(
            "psinorm_sol",
            doc="Normalised psi of the outer radial boundary in the inner SOL",
            value_type=[float, int],
        ),
        psinorm_pf=WithMeta(
            "psinorm_core",
            doc="Normalised psi of the inner radial boundary in the PFR",
            value_type=[float, int],
        ),
        psinorm_pf_lower=WithMeta(
            "psinorm_pf",
            doc="Normalised psi of the inner radial boundary in the lower PFR",
            value_type=[float, int],
        ),
        psinorm_pf_upper=WithMeta(
            "psinorm_pf",
            doc="Normalised psi of the inner radial boundary in the upper PFR",
            value_type=[float, int],
        ),
        # Poloidal flux ranges.
        # These are the values which are used in the mesh generation
        # The default values comes from the psinorm values, but the user
        # can override these defaults.
        psi_core=WithMeta(
            None,
            doc=(
                "Unnormalised poloidal flux value at the core boundary, used in mesh "
                "generation. Overrides psinorm_core if this value is given, if the "
                "option is none, then calculated from psinorm_core"
            ),
            value_type=[float, int, NoneType],
        ),
        psi_sol=WithMeta(
            None,
            doc=(
                "Unnormalised poloidal flux value at the SOL boundary, used in mesh "
                "generation. Overrides psinorm_sol if this value is given, if the "
                "option is none, then calculated from psinorm_sol"
            ),
            value_type=[float, int, NoneType],
        ),
        psi_sol_inner=WithMeta(
            None,
            doc=(
                "Unnormalised poloidal flux value at the inner SOL boundary, used in "
                "mesh generation. Overrides psinorm_sol_inner if this value is given, "
                "if the option is none, then calculated from psinorm_sol_inner"
            ),
            value_type=[float, int, NoneType],
        ),
        psi_pf_lower=WithMeta(
            None,
            doc=(
                "Unnormalised poloidal flux value at the lower PFR boundary, used in "
                "mesh generation. Overrides psinorm_pf_lower if this value is given, "
                "if the option is none, then calculated from psinorm_pf_lower"
            ),
            value_type=[float, int, NoneType],
        ),
        psi_pf_upper=WithMeta(
            None,
            doc=(
                "Unnormalised poloidal flux value at the upper PFR boundary, used in "
                "mesh generation. Overrides psinorm_pf_upper if this value is given, "
                "if the option is none, then calculated from psinorm_pf_upper"
            ),
            value_type=[float, int, NoneType],
        ),
        start_at_upper_outer=WithMeta(
            False,
            doc=(
                "Start gridding double-null at upper-outer divertor instead of "
                "lower-inner. Warning: this option was added to enable backward "
                "compatibility with restart files from simulations in "
                "upper-disconnected-double-null configuration using grid files from "
                "the IDL hypnotoad; it is not well tested and not recommended to use."
            ),
            value_type=bool,
        ),
        # Tolerance for positioning points that should be at X-point, but need to be
        # slightly displaced from the null so code can follow Grad(psi).
        # Number between 0. and 1.
        xpoint_offset=WithMeta(
            0.1,
            doc=(
                "Tolerance for placing intial positions for tracing perpendiculars "
                "that should start exactly at an X-point, but initial positions to be "
                "slightly displaced from the null so code can follow Grad(psi).  This "
                "is a numerical fudge factor that may need to be increased for "
                "low-resolution input equilibria."
            ),
            value_type=float,
            check_all=[is_positive, lambda x: x < 1.0],
        ),
    )

    def __init__(
        self,
        R1D,
        Z1D,
        psi2D,
        psi1D,
        fpol1D,
        pressure=None,
        wall=None,
        psi_axis=None,
        dct=False,
        make_regions=True,
        settings=None,
        nonorthogonal_settings=None,
    ):
        """
        Create a Tokamak equilibrium.

        Inputs
        ------

        R1D[nx]       1D array of major radius [m]
        Z1D[ny]       1D array of height [m]
        psi2D[nx,ny]  2D array of poloidal flux [Wb]
        psi1D[nf]     1D array of poloidal flux [Wb]
        fpol1D[nf]    1D array of f=R*Bt [mT]

        Keywords
        --------

        pressure[nf] = 1D array of pressure as a function of psi1D [Pa]

        wall = [(R0,Z0), (R1, Z1), ...]
               A list of coordinate pairs, defining the vessel wall.
               The wall is closed, so the last point connects to the first.
        psi_axis = float
               The value of poloidal flux on the magnetic axis. If not
               given, the value found at the axis will be used.
        dct = bool
               EXPERIMENTAL: If true, use a DCT to interpolate and differentiate
               poloidal flux. By default a cubic spline is used.
        make_regions = bool
               Generate the regions to be meshed. The default (True)
               means that the object is complete after initialisation.
               If set to False, e.g. for testing, self.makeRegions() should
               be called to generate the regions.
        settings = A dict that will be used to set non-default values of options
               (self.user_options)
        nonorthogonal_settings = A dict that will be used to set non-default values of
               options (self.nonorthogonal_options)

        """
        if dct:
            # Create an interpolation
            # This sets the functions
            #   self.psi
            #   self.f_R
            #   self.f_Z
            #   self.Bp_R
            #   self.Bp_Z
            #   self.d2psidR2
            #   self.d2psidZ2
            #   self.d2psidRdZ
            self.magneticFunctionsFromGrid(R1D, Z1D, psi2D)
        else:
            self.psi_func = interpolate.RectBivariateSpline(R1D, Z1D, psi2D)

        self.f_psi_sign = 1.0
        if len(fpol1D) > 0:
            # Spline for interpolation of f = R*Bt

            # Note: psi1D must be increasing
            if psi1D[-1] < psi1D[0]:
                self.f_psi_sign = -1.0

            self.f_spl = interpolate.InterpolatedUnivariateSpline(
                psi1D * self.f_psi_sign, fpol1D, ext=3
            )
            # ext=3 specifies that boundary values are used outside range

            # Spline representing the derivative of f
            self.fprime_spl = self.f_spl.derivative()
        else:
            self.f_spl = lambda psi: 0.0
            self.fprime_spl = lambda psi: 0.0

        # Optional pressure profile
        if pressure is not None:
            self.p_spl = interpolate.InterpolatedUnivariateSpline(
                psi1D * self.f_psi_sign, pressure, ext=3
            )
        else:
            # If no pressure, then not output to grid file
            self.p_spl = None

        # Find critical points (O- and X-points)
        R2D, Z2D = np.meshgrid(R1D, Z1D, indexing="ij")
        opoints, xpoints = critical.find_critical(R2D, Z2D, psi2D)
        if len(opoints) == 0:
            warnings.warn("No O-points found in TokamakEquilibrium input")
        else:
            if psi_axis is None:
                psi_axis = opoints[0][2]  # Psi on magnetic axis
            self.o_point = Point2D(opoints[0][0], opoints[0][1])
        self.psi_axis = psi_axis

        if len(xpoints) == 0:
            warnings.warn("No X-points found in TokamakEquilibrium input")

        self.x_points = [Point2D(r, z) for r, z, psi in xpoints]
        self.psi_sep = [psi for r, z, psi in xpoints]

        # Bounding box for domain
        self.Rmin = min(R1D)
        self.Rmax = max(R1D)
        self.Zmin = min(Z1D)
        self.Zmax = max(Z1D)

        # Wall geometry. Note: should be anti-clockwise
        if wall is None:
            # No wall given, so add one which is just inside the domain edge
            offset = 1e-2  # in m
            wall = [
                (self.Rmin + offset, self.Zmin + offset),
                (self.Rmax - offset, self.Zmin + offset),
                (self.Rmax - offset, self.Zmax - offset),
                (self.Rmin + offset, self.Zmax - offset),
            ]
        elif len(wall) < 3:
            raise ValueError(
                f"Wall must be a polygon, so should have at least 3 points. Got "
                f"wall={wall}"
            )

        if polygons.clockwise(wall):
            wall = wall[
                ::-1
            ]  # Reverse, without modifying input list (which .reverse() would)
        self.wall = [Point2D(r, z) for r, z in wall]

        # Take the default settings, then the options keyword, then
        # any additional keyword arguments
        self.user_options = self.user_options_factory.create(settings)

        self.equilibOptions = {}

        super().__init__(nonorthogonal_settings)

        # Print the table of options
        print(self.user_options.as_table(), flush=True)

        if make_regions:
            # Create self.regions
            self.makeRegions()

    def findLegs(self, xpoint, radius=0.01, step=0.01):
        """Find the divertor legs coming from a given X-point

        xpoint   Point2D object giving position
        radius   Search distance from X-point, in meters
        step     Integration step size, in meters
        """

        psi_sep = self.psi(xpoint.R, xpoint.Z)  # Value of psi

        # Draw a half-circle around this X-point
        theta0 = np.arctan2(xpoint.Z - self.o_point.Z, xpoint.R - self.o_point.R)
        angles = np.linspace(theta0 - np.pi / 2.0, theta0 + np.pi / 2.0)

        rvals = xpoint.R + radius * np.cos(angles)
        zvals = xpoint.Z + radius * np.sin(angles)
        psivals = self.psi(rvals, zvals)
        # Note: If psivals crosses psi_sep, the value becomes negative
        inds = np.nonzero((psivals[1:] - psi_sep) * (psivals[:-1] - psi_sep) < 0.0)[0]

        # Currently only handle standard X-points (no snowflakes)
        assert len(inds) == 2

        # Divide-and-conquer to get a points on the leg
        # This goes into a list leg_points = [(r,z),..]
        # These are arranged anticlockwise
        leg_points = []
        for ind in inds:
            # Define the line along which the flux surface lies
            # e.g r = r0 + dr * s
            r0 = rvals[ind]
            z0 = zvals[ind]
            dr = rvals[ind + 1] - r0
            dz = zvals[ind + 1] - z0

            s1 = 0.0
            s2 = 1.0
            psi1 = psivals[ind]  # s = 0

            while s2 - s1 > 1e-5:
                smid = 0.5 * (s1 + s2)
                psi_mid = self.psi(r0 + smid * dr, z0 + smid * dz)

                if (psi_mid - psi_sep) * (psi1 - psi_sep) < 0.0:
                    # Between psi_mid and psi1
                    s2 = smid
                else:
                    psi1 = psi_mid
                    s1 = smid
            smid = 0.5 * (s1 + s2)
            r = r0 + smid * dr
            z = z0 + smid * dz
            leg_points.append((r, z))

        # For each leg, follow the magnetic field away from the X-point
        # until the line intersects the wall
        leg_lines = []
        for leg in leg_points:
            line = [xpoint]  # Start with the X-point

            Br = self.Bp_R(*leg)
            Bz = self.Bp_Z(*leg)
            # Dot product vector from X-point to leg with Bp
            # The sign is used to tell which way to integrate
            sign = np.sign((leg[0] - xpoint.R) * Br + (leg[1] - xpoint.Z) * Bz)

            # Integrate in this direction until the wall is intersected
            # This is affected by sign, which determines which way to integrate
            def dpos_dl(distance, pos):
                r = pos[0]
                z = pos[1]
                Br = self.Bp_R(r, z)
                Bz = self.Bp_Z(r, z)
                B = np.sqrt(Br ** 2 + Bz ** 2)
                return [sign * Br / B, sign * Bz / B]

            pos = leg  # Starting position
            while True:
                # Integrate a distance "step" along the leg
                solve_result = solve_ivp(dpos_dl, (0.0, step), pos)
                newpos = (solve_result.y[0][1], solve_result.y[1][1])

                # Check if we have crossed the boundary
                # somewhere between pos and newpos

                intersect = self.wallIntersection(Point2D(*pos), Point2D(*newpos))
                if intersect is not None:
                    line.append(intersect)  # Put the intersection in the line
                    break
                pos = newpos
                line.append(Point2D(*pos))

            # Should now have an intersect with the wall
            # which is a Point2D object
            leg_lines.append(line)

        # Now have a list of 2 legs. Check which one is the inner leg
        # by comparing the major radius of the strike points
        if leg_lines[0][-1].R > leg_lines[1][-1].R:
            leg_lines = leg_lines[::-1]
        return {"inner": leg_lines[0], "outer": leg_lines[1]}

    def makeRegions(self):
        """Main region generation function. Regions are logically
        rectangular ranges in poloidal angle; segments are
        ranges of poloidal flux (radial coordinate). Poloidal ranges,
        radial segments, and the connections between them describe
        the topology and geometry of the mesh.

        This function is called by __init__ to generate regions unless
        make_regions is set to False.

        The main steps in doing this are:
        1. Set defaults if not already set by user
        2. Identify whether single or double null
        3. Describe the leg and core regions, depending on the topology
        4. Follow flux surfaces based on core region descriptions
           (self.coreRegionToRegion)
        5. Process all regions into a set of EquilibriumObjects
           (self.createRegionObjects)
        6. Sort the EquilibriumRegion objects for BoutMesh output
           (self.createRegionObjects).
        7. Connect regions together

        Inputs
        ------

        Keywords set options in user_options. Note that not all options
        will have any effect, because they were already used in __init__.

        Modifies
        --------

        - self.user_options    Sets default values if not set by user
        - self.regions         OrderedDict of EquilibriumRegion objects

        """
        assert self.psi_axis is not None

        # psi values
        def psinorm_to_psi(psinorm):
            if psinorm is None:
                return None
            return self.psi_axis + psinorm * (self.psi_sep[0] - self.psi_axis)

        def psi_to_psinorm(psi):
            if psi is None:
                return None
            return (psi - self.psi_axis) / (self.psi_sep[0] - self.psi_axis)

        self.psi_core = with_default(
            self.user_options.psi_core, psinorm_to_psi(self.user_options.psinorm_core)
        )
        self.psi_sol = with_default(
            self.user_options.psi_sol, psinorm_to_psi(self.user_options.psinorm_sol)
        )
        self.psi_sol_inner = with_default(
            self.user_options.psi_sol_inner,
            psinorm_to_psi(self.user_options.psinorm_sol_inner),
        )
        self.psi_pf_lower = with_default(
            self.user_options.psi_pf_lower,
            psinorm_to_psi(self.user_options.psinorm_pf_lower),
        )
        self.psi_pf_upper = with_default(
            self.user_options.psi_pf_upper,
            psinorm_to_psi(self.user_options.psinorm_pf_upper),
        )

        self.poloidal_spacing_delta_psi = with_default(
            self.user_options.poloidal_spacing_delta_psi,
            np.abs((self.psi_core - self.psi_sol) / 20.0),
        )

        # Filter out the X-points not in range.
        # Keep only those with normalised psi < psinorm_sol
        self.psi_sep, self.x_points = zip(
            *(
                (psi, xpoint)
                for psi, xpoint in zip(self.psi_sep, self.x_points)
                if psi_to_psinorm(psi) < self.user_options.psinorm_sol
            )
        )

        # Check that there are only one or two left
        assert 0 < len(self.x_points) <= 2

        if len(self.x_points) == 1:
            # Generate the specifications for a lower or upper single null
            leg_regions, core_regions, segments, connections = self.describeSingleNull()
        else:
            # Specifications for a double null (connected or disconnected)
            leg_regions, core_regions, segments, connections = self.describeDoubleNull()

        # Create a new dictionary, which will contain all regions
        # including core and legs
        all_regions = leg_regions.copy()
        all_regions.update(self.coreRegionToRegion(core_regions))

        # Create the regions in an OrderedDict, assign to self.regions
        self.regions = self.createRegionObjects(all_regions, segments)

        # Make the connections between regions
        for connection in connections:
            self.makeConnection(*connection)

    def describeSingleNull(self):
        """
        Create the specifications for a single null configuration

        Returns
        -------

        leg_regions    Dictionary describing poloidal regions in legs
        core_regions   Dictionary describing poloidal regions between X-points
        segments       Dictionary describing radial segments
                        nx          Number of radial (x) cells
                        psi_vals    1D array of poloidal flux values. Length 2*nx+1
        connections    List of connections between regions
        """
        assert len(self.x_points) == 1
        # Single null. Could be lower or upper

        # Find lines along the legs from X-point to target
        legs = self.findLegs(self.x_points[0])

        # Move the first point of each leg slightly away from the X-point
        diff = self.user_options.xpoint_offset
        for leg in legs.values():
            leg[0] = diff * leg[1] + (1.0 - diff) * leg[0]

        psi_sep = self.psi_sep[0]
        xpoint = self.x_points[0]

        # Average radial grid spacing in each region
        dpsidi_sol = (self.psi_sol - self.psi_sep[0]) / self.user_options.nx_sol
        dpsidi_core = (self.psi_sep[0] - self.psi_core) / self.user_options.nx_core
        if self.x_points[0].Z < self.o_point.Z:
            # Lower single null
            dpsidi_pf = (self.psi_sep[0] - self.psi_pf_lower) / self.user_options.nx_pf
        else:
            # Upper single null
            dpsidi_pf = (self.psi_sep[0] - self.psi_pf_upper) / self.user_options.nx_pf

        # Get the smallest absolute grid spacing for the separatrix
        dpsidi_sep = min([dpsidi_sol, dpsidi_core, dpsidi_pf], key=abs)

        # decrease (assuming the factor is <1) the spacing around the separatrix by the
        # factor psi_spacing_separatrix_multiplier
        if self.user_options.psi_spacing_separatrix_multiplier is not None:
            dpsidi_sep *= self.user_options.psi_spacing_separatrix_multiplier

        # Radial segments common to both lower and upper single null geometries
        segments = {
            "core": {
                "nx": self.user_options.nx_core,
                "psi_start": self.psi_core,
                "psi_end": psi_sep,
                "grad_end": dpsidi_sep,
            },
            "sol": {
                "nx": self.user_options.nx_sol,
                "psi_start": psi_sep,
                "psi_end": self.psi_sol,
                "grad_start": dpsidi_sep,
            },
        }

        core_regions = {
            "core": {
                "segments": ["core", "sol"],
                "ny": self.user_options.ny_inner_sol + self.user_options.ny_outer_sol,
                "kind": "X.X",
                "xpoints_at_start": [None, xpoint, None],
                "xpoints_at_end": [None, xpoint, None],
                "psi_at_start": psi_sep,
                "psi_at_end": psi_sep,
            }
        }

        if self.x_points[0].Z < self.o_point.Z:
            print("Generating a lower single null", flush=True)

            segments["lower_pf"] = {
                "nx": self.user_options.nx_pf,
                "psi_start": self.psi_pf_lower,
                "psi_end": psi_sep,
                "grad_end": dpsidi_sep,
            }

            leg_regions = {
                "inner_lower_divertor": {
                    "segments": ["lower_pf", "sol"],
                    "ny": self.user_options.ny_inner_lower_divertor,
                    "kind": "wall.X",
                    "points": legs["inner"][
                        ::-1
                    ],  # Reversed because leg is from X-point to target
                    "psi": psi_sep,
                    "wall_at_start": [0, 0],
                    "xpoints_at_end": [None, xpoint, None],
                },
                "outer_lower_divertor": {
                    "segments": ["lower_pf", "sol"],
                    "ny": self.user_options.ny_outer_lower_divertor,
                    "kind": "X.wall",
                    "points": legs["outer"],
                    "psi": psi_sep,
                    "xpoints_at_start": [None, xpoint, None],
                    "wall_at_end": [0, 0],
                },
            }
            # Connections between regions
            connections = [
                (
                    "inner_lower_divertor",
                    0,
                    "outer_lower_divertor",
                    0,
                ),  # inner lower PF -> outer lower PF
                ("inner_lower_divertor", 1, "core", 1),  # inner lower PF -> Core
                ("core", 1, "outer_lower_divertor", 1),  # Core -> outer lower PF
                ("core", 0, "core", 0),  # Core -> core
            ]
        else:
            print("Upper Single Null", flush=True)

            # The mesh is still arranged clockwise, meaning that y indexing starts
            # at the outer upper divertor, and ends at the inner upper divertor.
            #
            # The points in the legs are ordered from X-point to target,
            # so here the outer leg needs to be reversed rather than the inner leg

            segments["upper_pf"] = {
                "nx": self.user_options.nx_pf,
                "psi_start": self.psi_pf_upper,
                "psi_end": psi_sep,
                "grad_end": dpsidi_sep,
            }

            leg_regions = {
                "inner_upper_divertor": {
                    "segments": ["upper_pf", "sol"],
                    "ny": self.user_options.ny_inner_upper_divertor,
                    "kind": "X.wall",
                    "points": legs["inner"],
                    "psi": psi_sep,
                    "xpoints_at_start": [None, xpoint, None],
                    "wall_at_end": [0, 0],
                },
                "outer_upper_divertor": {
                    "segments": ["upper_pf", "sol"],
                    "ny": self.user_options.ny_outer_upper_divertor,
                    "kind": "wall.X",
                    "points": legs["outer"][::-1],
                    "psi": psi_sep,
                    "wall_at_start": [0, 0],
                    "xpoints_at_end": [None, xpoint, None],
                },
            }
            # Connections between regions
            connections = [
                (
                    "outer_upper_divertor",
                    0,
                    "inner_upper_divertor",
                    0,
                ),  # outer lower PF -> inner lower PF
                ("outer_upper_divertor", 1, "core", 1),  # outer lower PF -> Core
                ("core", 1, "inner_upper_divertor", 1),  # Core -> inner upper PF
                ("core", 0, "core", 0),  # Core -> core
            ]

        return (
            leg_regions,
            core_regions,
            self.segmentsWithPsivals(segments),
            connections,
        )

    def describeDoubleNull(self):
        """
        Create the specifications for a double null configuration,
        either connected or disconnected.

        Returns
        -------

        leg_regions    Dictionary describing poloidal regions in legs
        core_regions   Dictionary describing poloidal regions between X-points
        segments       Dictionary describing radial segments
                        nx          Number of radial (x) cells
                        psi_vals    1D array of poloidal flux values. Length 2*nx+1
        connections    List of connections between regions
        """

        assert len(self.x_points) == 2

        if self.x_points[0].Z < self.o_point.Z:
            # Lower double null
            lower_xpt_ind, upper_xpt_ind = 0, 1
        else:
            # Upper double null
            lower_xpt_ind, upper_xpt_ind = 1, 0

        lower_x_point = self.x_points[lower_xpt_ind]
        upper_x_point = self.x_points[upper_xpt_ind]
        lower_psi = self.psi_sep[lower_xpt_ind]
        upper_psi = self.psi_sep[upper_xpt_ind]

        assert np.isclose(lower_psi, self.psi(*lower_x_point))
        assert np.isclose(upper_psi, self.psi(*upper_x_point))

        # Find lines along the legs from X-point to target
        lower_legs = self.findLegs(lower_x_point)
        upper_legs = self.findLegs(upper_x_point)

        # Move the first point of each leg slightly away from the X-point
        diff = self.user_options.xpoint_offset
        for legs in [lower_legs, upper_legs]:
            for leg in legs.values():
                leg[0] = diff * leg[1] + (1.0 - diff) * leg[0]

        # Average radial grid spacing in each region
        dpsidi_sol_inner = (
            self.psi_sol_inner - self.psi_sep[0]
        ) / self.user_options.nx_sol_inner
        dpsidi_sol_outer = (
            self.psi_sol - self.psi_sep[0]
        ) / self.user_options.nx_sol_outer
        dpsidi_core = (self.psi_sep[0] - self.psi_core) / self.user_options.nx_core
        dpsidi_pf_upper = (upper_psi - self.psi_pf_upper) / self.user_options.nx_pf
        dpsidi_pf_lower = (lower_psi - self.psi_pf_lower) / self.user_options.nx_pf

        # Get the smallest absolute grid spacing for the separatrix
        dpsidi_sep = min(
            [
                dpsidi_sol_inner,
                dpsidi_sol_outer,
                dpsidi_core,
                dpsidi_pf_upper,
                dpsidi_pf_lower,
            ],
            key=abs,
        )

        # decrease (assuming the factor is <1) the spacing around the separatrix by the
        # factor psi_spacing_separatrix_multiplier
        if self.user_options.psi_spacing_separatrix_multiplier is not None:
            dpsidi_sep *= self.user_options.psi_spacing_separatrix_multiplier

        # Number of points in the inter-separatrix region
        nx_inter_sep = self.user_options.nx_inter_sep

        # Adjust the number of points in upper and lower PF regions,
        # to keep nx constant between regions. This is because some regions
        # will have a near SOL but not others
        if self.x_points[0] == lower_x_point:
            nx_pf_lower = self.user_options.nx_pf
            nx_pf_upper = self.user_options.nx_pf + nx_inter_sep
        else:
            nx_pf_lower = self.user_options.nx_pf + nx_inter_sep
            nx_pf_upper = self.user_options.nx_pf

        # Radial segments i.e. gridded ranges of poloidal flux
        # These are common to both connected and disconnected double null
        segments = {
            "core": {
                "nx": self.user_options.nx_core,
                "psi_start": self.psi_core,
                "psi_end": self.psi_sep[0],
                "grad_end": dpsidi_sep,
            },
            "upper_pf": {
                "nx": nx_pf_upper,
                "psi_start": self.psi_pf_upper,
                "psi_end": upper_psi,
                "grad_end": dpsidi_sep,
            },
            "lower_pf": {
                "nx": nx_pf_lower,
                "psi_start": self.psi_pf_lower,
                "psi_end": lower_psi,
                "grad_end": dpsidi_sep,
            },
            "inner_sol": {
                "nx": self.user_options.nx_sol_inner,
                "psi_start": self.psi_sep[-1],
                "psi_end": self.psi_sol_inner,
                "grad_start": dpsidi_sep,
            },
            "outer_sol": {
                "nx": self.user_options.nx_sol_outer,
                "psi_start": self.psi_sep[-1],
                "psi_end": self.psi_sol,
                "grad_start": dpsidi_sep,
            },
        }

        if nx_inter_sep == 0:
            print("Generating a connected double null", flush=True)

            # Only use psi of inner separatrix, not the outer one (if it is slightly
            # different)
            segments["upper_pf"]["psi_end"] = self.psi_sep[0]
            segments["lower_pf"]["psi_end"] = self.psi_sep[0]
            segments["inner_sol"]["psi_start"] = self.psi_sep[0]
            segments["outer_sol"]["psi_start"] = self.psi_sep[0]

            # Description of each poloidal region
            leg_regions = {
                "inner_lower_divertor": {
                    "segments": ["lower_pf", "inner_sol"],
                    "ny": self.user_options.ny_inner_lower_divertor,
                    "kind": "wall.X",
                    "points": lower_legs["inner"][::-1],
                    "psi": lower_psi,
                    "wall_at_start": [0, 0],
                    "xpoints_at_end": [None, lower_x_point, None],
                },
                "outer_lower_divertor": {
                    "segments": ["lower_pf", "outer_sol"],
                    "ny": self.user_options.ny_outer_lower_divertor,
                    "kind": "X.wall",
                    "points": lower_legs["outer"],
                    "psi": lower_psi,
                    "xpoints_at_start": [None, lower_x_point, None],
                    "wall_at_end": [0, 0],
                },
                "inner_upper_divertor": {
                    "segments": ["upper_pf", "inner_sol"],
                    "ny": self.user_options.ny_inner_upper_divertor,
                    "kind": "X.wall",
                    "points": upper_legs["inner"],
                    "psi": upper_psi,
                    "xpoints_at_start": [None, upper_x_point, None],
                    "wall_at_end": [0, 0],
                },
                "outer_upper_divertor": {
                    "segments": ["upper_pf", "outer_sol"],
                    "ny": self.user_options.ny_outer_upper_divertor,
                    "kind": "wall.X",
                    "points": upper_legs["outer"][::-1],
                    "psi": upper_psi,
                    "wall_at_start": [0, 0],
                    "xpoints_at_end": [None, upper_x_point, None],
                },
            }

            core_regions = {
                "inner_core": {
                    "segments": ["core", "inner_sol"],
                    "ny": self.user_options.ny_inner_sol,
                    "kind": "X.X",
                    "xpoints_at_start": [None, lower_x_point, None],
                    "xpoints_at_end": [None, upper_x_point, None],
                    "psi_at_start": lower_psi,
                    "psi_at_end": upper_psi,
                },
                "outer_core": {
                    "segments": ["core", "outer_sol"],
                    "ny": self.user_options.ny_outer_sol,
                    "kind": "X.X",
                    "xpoints_at_start": [None, upper_x_point, None],
                    "xpoints_at_end": [None, lower_x_point, None],
                    "psi_at_start": upper_psi,
                    "psi_at_end": lower_psi,
                },
            }

            connections = [
                ("inner_lower_divertor", 0, "outer_lower_divertor", 0),
                ("outer_upper_divertor", 0, "inner_upper_divertor", 0),
                ("inner_core", 0, "outer_core", 0),
                ("outer_core", 0, "inner_core", 0),
                ("inner_lower_divertor", 1, "inner_core", 1),
                ("inner_core", 1, "inner_upper_divertor", 1),
                ("outer_upper_divertor", 1, "outer_core", 1),
                ("outer_core", 1, "outer_lower_divertor", 1),
            ]

            return (
                leg_regions,
                core_regions,
                self.segmentsWithPsivals(segments),
                connections,
            )

        else:
            print("Generating a disconnected double null", flush=True)

            # Disconnected double null -> Additional radial segment
            segments["near_sol"] = {
                "nx": nx_inter_sep,
                "psi_start": self.psi_sep[0],
                "psi_end": self.psi_sep[1],
                "grad_start": dpsidi_sep,
                "grad_end": dpsidi_sep,
            }

            if self.x_points[0] == lower_x_point:
                print("Lower double null", flush=True)

                inner_lower_segments = ["lower_pf", "near_sol", "inner_sol"]
                outer_lower_segments = ["lower_pf", "near_sol", "outer_sol"]

                inner_upper_segments = ["upper_pf", "upper_pf2", "inner_sol"]
                outer_upper_segments = ["upper_pf", "upper_pf2", "outer_sol"]

                # Lower X-point between 1st region (PF or core) and near SOL
                lower_core_xpoints = [None, lower_x_point, None, None]
                lower_pf_xpoints = [None, lower_x_point, None, None]

                # Upper X-point between near SOL and inner/outer SOL
                upper_core_xpoints = [None, None, upper_x_point, None]
                upper_pf_xpoints = [None, None, upper_x_point, None]

                connections = [
                    # PFR
                    ("inner_lower_divertor", 0, "outer_lower_divertor", 0),
                    ("outer_upper_divertor", 0, "inner_upper_divertor", 0),  # upper_pf
                    ("outer_upper_divertor", 1, "inner_upper_divertor", 1),  # upper_pf2
                    # core
                    ("inner_core", 0, "outer_core", 0),
                    ("outer_core", 0, "inner_core", 0),
                    # near SOL
                    ("inner_lower_divertor", 1, "inner_core", 1),
                    ("inner_core", 1, "outer_core", 1),
                    ("outer_core", 1, "outer_lower_divertor", 1),
                    # inner SOL
                    ("inner_lower_divertor", 2, "inner_core", 2),
                    ("inner_core", 2, "inner_upper_divertor", 2),
                    # outer SOL
                    ("outer_upper_divertor", 2, "outer_core", 2),
                    ("outer_core", 2, "outer_lower_divertor", 2),
                ]

                # Perform radial gridding
                segments = self.segmentsWithPsivals(segments)

                # Split the secondary PFR segment into two segments.
                # This is so that all regions have the same number of segments
                # (currently needed by BoutMesh)

                segments["upper_pf2"] = {
                    "nx": nx_inter_sep,
                    "psi_vals": segments["upper_pf"]["psi_vals"][
                        -(2 * nx_inter_sep + 1) : None
                    ],
                }

                segments["upper_pf"]["psi_vals"] = segments["upper_pf"]["psi_vals"][
                    : (-2 * nx_inter_sep)
                ]
                segments["upper_pf"]["nx"] -= nx_inter_sep

            else:
                print("Upper double null", flush=True)
                inner_lower_segments = ["lower_pf", "lower_pf2", "inner_sol"]
                outer_lower_segments = ["lower_pf", "lower_pf2", "outer_sol"]

                inner_upper_segments = ["upper_pf", "near_sol", "inner_sol"]
                outer_upper_segments = ["upper_pf", "near_sol", "outer_sol"]

                upper_core_xpoints = [None, upper_x_point, None, None]
                upper_pf_xpoints = [None, upper_x_point, None, None]

                lower_core_xpoints = [None, None, lower_x_point, None]
                lower_pf_xpoints = [None, None, lower_x_point, None]

                connections = [
                    # PFR
                    ("inner_lower_divertor", 0, "outer_lower_divertor", 0),  # lower_pf
                    ("inner_lower_divertor", 1, "outer_lower_divertor", 1),  # lower_pf2
                    ("outer_upper_divertor", 0, "inner_upper_divertor", 0),
                    # core
                    ("inner_core", 0, "outer_core", 0),
                    ("outer_core", 0, "inner_core", 0),
                    # near SOL
                    ("outer_upper_divertor", 1, "outer_core", 1),
                    ("outer_core", 1, "inner_core", 1),
                    ("inner_core", 1, "inner_upper_divertor", 1),
                    # inner SOL
                    ("inner_lower_divertor", 2, "inner_core", 2),
                    ("inner_core", 2, "inner_upper_divertor", 2),
                    # outer SOL
                    ("outer_upper_divertor", 2, "outer_core", 2),
                    ("outer_core", 2, "outer_lower_divertor", 2),
                ]

                # Perform radial gridding
                segments = self.segmentsWithPsivals(segments)

                # Split the secondary PFR segment into two segments.
                # This is so that all regions have the same number of segments
                # (currently needed by BoutMesh)
                segments["lower_pf2"] = {
                    "nx": nx_inter_sep,
                    "psi_vals": segments["lower_pf"]["psi_vals"][
                        -(2 * nx_inter_sep + 1) : None
                    ],
                }

                segments["lower_pf"]["psi_vals"] = segments["lower_pf"]["psi_vals"][
                    : (-2 * nx_inter_sep)
                ]
                segments["lower_pf"]["nx"] -= nx_inter_sep

            # Here these are for both lower and upper disconnected double null
            core_regions = {
                "inner_core": {
                    "segments": ["core", "near_sol", "inner_sol"],
                    "ny": self.user_options.ny_inner_sol,
                    "kind": "X.X",
                    "xpoints_at_start": lower_core_xpoints,
                    "xpoints_at_end": upper_core_xpoints,
                    "psi_at_start": lower_psi,
                    "psi_at_end": upper_psi,
                },
                "outer_core": {
                    "segments": ["core", "near_sol", "outer_sol"],
                    "ny": self.user_options.ny_outer_sol,
                    "kind": "X.X",
                    "xpoints_at_start": upper_core_xpoints,
                    "xpoints_at_end": lower_core_xpoints,
                    "psi_at_start": upper_psi,
                    "psi_at_end": lower_psi,
                },
            }

            leg_regions = {
                "inner_lower_divertor": {
                    "segments": inner_lower_segments,
                    "ny": self.user_options.ny_inner_lower_divertor,
                    "kind": "wall.X",
                    "points": lower_legs["inner"][::-1],
                    "psi": lower_psi,
                    "wall_at_start": [0, 0],
                    "xpoints_at_end": lower_pf_xpoints,
                },
                "outer_lower_divertor": {
                    "segments": outer_lower_segments,
                    "ny": self.user_options.ny_outer_lower_divertor,
                    "kind": "X.wall",
                    "points": lower_legs["outer"],
                    "psi": lower_psi,
                    "xpoints_at_start": lower_pf_xpoints,
                    "wall_at_end": [0, 0],
                },
                "inner_upper_divertor": {
                    "segments": inner_upper_segments,
                    "ny": self.user_options.ny_inner_upper_divertor,
                    "kind": "X.wall",
                    "points": upper_legs["inner"],
                    "psi": upper_psi,
                    "xpoints_at_start": upper_pf_xpoints,
                    "wall_at_end": [0, 0],
                },
                "outer_upper_divertor": {
                    "segments": outer_upper_segments,
                    "ny": self.user_options.ny_outer_upper_divertor,
                    "kind": "wall.X",
                    "points": upper_legs["outer"][::-1],
                    "psi": upper_psi,
                    "wall_at_start": [0, 0],
                    "xpoints_at_end": upper_pf_xpoints,
                },
            }

            return leg_regions, core_regions, segments, connections

    def coreRegionToRegion(self, core_regions, npoints=100):
        """
        For each poloidal arc along a separatrix between two X-points
        (core region), find a set of points between the X-points.
        The result is returned as a dict of regions (like leg regions)

        Inputs
        ------

        core_regions  A dictionary containing definitions of core regions.
                      Keys are:
                        segments   A list of segment names
                        ny         Number of poloidal (y) points
                        kind       A string e.g. "wall.X"
                        xpoints_at_start    A list of Point2D objects or None
                        xpoints_at_end      A list of Point2D objects or None
                        psi_at_start    Poloidal flux at the start of the line
                        psi_at_end      Poloidal flux at the end of the line

        npoints   number of points in each core region

        Returns
        -------
        A dictionary of region definitions, compatible with processing code
        for leg regions.
        """

        # Loop through core regions, calculate points along the lines,
        # and put into the result dictionary
        result = {}
        for name, region in core_regions.items():
            region = region.copy()  # So we don't modify the input

            def any_value(values):
                "Return the first non-None value in list, or None"
                return next((val for val in values if val is not None), None)  # Default

            start_x = any_value(region["xpoints_at_start"])
            end_x = any_value(region["xpoints_at_end"])
            start_psi = region["psi_at_start"]
            end_psi = region["psi_at_end"]

            # Range of angles. Note: This angle goes anticlockwise
            # so core regions need to be reversed
            start_angle = np.arctan2(
                start_x.Z - self.o_point.Z, start_x.R - self.o_point.R
            )
            end_angle = np.arctan2(end_x.Z - self.o_point.Z, end_x.R - self.o_point.R)
            if end_angle >= start_angle:
                end_angle -= 2 * np.pi

            # Angle offset from the X-point. This is to reduce the chances
            # of missing the X-point, passing through to the other side.
            dtheta = 0.5 * (end_angle - start_angle) / npoints
            r0, z0 = self.o_point.R, self.o_point.Z  # Location of O-point

            # If the start and end X-point are different, interpolate
            # from one to the other. This helps ensure that constructed
            # coordinate lines don't go the wrong side of X-points.
            def psival(angle):
                "Interpolation in psi with angle"
                norm = (angle - start_angle) / (end_angle - start_angle)

                # Smoother step function (Ken Perlin)
                # https://en.wikipedia.org/wiki/Smoothstep
                norm = 6.0 * norm ** 5 - 15.0 * norm ** 4 + 10.0 * norm ** 3

                return norm * end_psi + (1.0 - norm) * start_psi

            # Iterate in angle from start to end
            points = [
                Point2D(
                    *critical.find_psisurface(
                        self,
                        r0,
                        z0,
                        r0 + 8.0 * np.cos(angle),
                        z0 + 8.0 * np.sin(angle),
                        psival=psival(angle),
                    )
                )
                for angle in np.linspace(
                    start_angle + dtheta, end_angle - dtheta, npoints
                )
            ]

            # Add points to the beginning and end near (but not at) the X-points
            diff = self.user_options.xpoint_offset
            if diff < 0.0 or diff > 1.0:
                raise ValueError(f"xpoint_offset={diff} should be between 0 and 1.")

            region["points"] = (
                [(1.0 - diff) * start_x + diff * points[0]]
                + points
                + [(1.0 - diff) * end_x + diff * points[-1]]
            )

            region["psi"] = None  # Not all points on the same flux surface

            result[name] = region
        return result

    def segmentsWithPsivals(self, segments):
        """
        Grids radial segments

        Input
        -----

        segments   A dict of segments, each of which is a dictionary containing
                      nx    Number of points in psi (x)
                      psi_start   The poloidal flux at the start of the segment
                      psi_end     The poloidal flux at the end of the segment
                      grad_start [optional]  Cell spacing at the start
                      grad_end   [optional] Cell spacing at the end

        The input is not modified

        Returns
        -------

        A dictionary of segments, with an additional key "psi_vals"
        """
        result = {}
        for name, segment in segments.items():
            segment_with_psival = segment.copy()

            psi_func = self.getPolynomialGridFunc(
                segment["nx"],
                segment["psi_start"],
                segment["psi_end"],
                grad_lower=segment.get("grad_start", None),
                grad_upper=segment.get("grad_end", None),
            )

            segment_with_psival["psi_vals"] = self.make1dGrid(segment["nx"], psi_func)
            result[name] = segment_with_psival
        return result

    def createRegionObjects(self, all_regions, segments):
        """
        Create an OrderedDict of EquilibriumRegion objects,
        using specifications for the regions and segments
        in the all_regions and segments dictionaries.

        These regions need to be sorted so that BoutMesh
        can generate branch cut indices. To do this ordering,
        a limited set of region names should be used:
        - 'inner_lower_divertor'
        - 'core' (for single null)
        - 'inner_core' and 'outer_core' (for double null)
        - 'inner_upper_divertor'
        - 'outer_upper_divertor'
        - 'outer_lower_divertor'

        Inputs
        ------

        all_regions   Dictionary containing specification for each region
        segments      Dictionary of radial segment definitions
                        nx            Number of radial cells
                        psi_vals      1D array of psi values, length 2*nx+1

        Returns
        -------

        None. Modifies self.regions
        """

        # Set the total number of grid cells in y
        self.ny_total = sum([region["ny"] for region in all_regions.values()])

        # Loop through all regions. For each one create an EquilibriumRegion
        region_objects = {}
        for name, region in all_regions.items():
            eqreg = EquilibriumRegion(
                equilibrium=self,
                name=name,
                nSegments=len(region["segments"]),  # The number of radial regions
                nx=[segments[seg_name]["nx"] for seg_name in region["segments"]],
                ny=region["ny"],
                ny_total=self.ny_total,
                kind=region["kind"],
                # The following arguments are passed through to PsiContour
                points=region["points"],  # list of Point2D objects on the line
                psival=region["psi"],
            )

            # Grids of psi values in each segment
            eqreg.psi_vals = [segments[name]["psi_vals"] for name in region["segments"]]

            eqreg.separatrix_radial_index = 1

            if "xpoints_at_start" in region:
                eqreg.xPointsAtStart = region["xpoints_at_start"]

            if "xpoints_at_end" in region:
                eqreg.xPointsAtEnd = region["xpoints_at_end"]

            if "wall_at_start" in region:
                eqreg.wallSurfaceAtStart = region["wall_at_start"]

            if "wall_at_end" in region:
                eqreg.wallSurfaceAtEnd = region["wall_at_end"]

            # Pressure profiles
            if self.p_spl is not None:
                if "wall" in region["kind"]:
                    # A leg region. Reflect the pressure in poloidal flux
                    # so that the pressure in the private flux region falls
                    # away from the separatrix

                    # Determine if poloidal flux is increasing or decreasing with radius
                    sign = np.sign(self.psi_sep[0] - self.psi_axis)

                    assert region["psi"] is not None
                    leg_psi = region["psi"]
                    eqreg.pressure = lambda psi: self.pressure(
                        leg_psi + sign * abs(psi - leg_psi)
                    )
                else:
                    # Core region, so use the core pressure
                    eqreg.pressure = self.pressure

            region_objects[name] = eqreg
        # The region objects need to be sorted, so that the
        # BoutMesh generator can use jyseps indices to introduce branch cuts

        if "inner_lower_divertor" in region_objects:
            if not self.user_options.start_at_upper_outer:
                ordering = [
                    "inner_lower_divertor",
                    # For single null; in double null this will be ignored
                    "core",
                    # For double null; in single null these will be ignored
                    "inner_core",
                    "inner_upper_divertor",
                    "outer_upper_divertor",
                    "outer_core",
                    #
                    "outer_lower_divertor",
                ]
            else:
                # Special case intended for backward compatibility with simulations
                # using upper-disconnected-double-null IDL-hypnotoad grid files
                ordering = [
                    "outer_upper_divertor",
                    "outer_core",
                    "outer_lower_divertor",
                    "inner_lower_divertor",
                    "core",
                    "inner_core",
                    "inner_upper_divertor",
                ]
        else:
            # Upper single null special case
            ordering = ["outer_upper_divertor", "core", "inner_upper_divertor"]

        # Check that all regions are in the ordering
        for key in region_objects:
            if key not in ordering:
                raise ValueError("Region '" + key + "' is not in ordering")

        # Sort the region objects, ignoring regions which are not present
        return OrderedDict(
            [(key, region_objects[key]) for key in ordering if key in region_objects]
        )

    def handleMultiLocationArray(getResult):
        @functools.wraps(getResult)
        # Define a function which handles MultiLocationArray arguments
        def handler(self, *args):
            if isinstance(args[0], MultiLocationArray):
                for arg in args[1:]:
                    assert isinstance(arg, MultiLocationArray), (
                        "if first arg is a MultiLocationArray, then others must be as "
                        "well"
                    )
                result = MultiLocationArray(args[0].nx, args[0].ny)

                if all(arg.centre is not None for arg in args):
                    result.centre = getResult(self, *(arg.centre for arg in args))

                if all(arg.xlow is not None for arg in args):
                    result.xlow = getResult(self, *(arg.xlow for arg in args))

                if all(arg.ylow is not None for arg in args):
                    result.ylow = getResult(self, *(arg.ylow for arg in args))

                if all(arg.corners is not None for arg in args):
                    result.corners = getResult(self, *(arg.corners for arg in args))
            else:
                result = getResult(self, *args)
            return result

        return handler

    @handleMultiLocationArray
    def psi(self, R, Z):
        "Return the poloidal flux at the given (R,Z) location"
        return self.psi_func(R, Z, grid=False)

    @handleMultiLocationArray
    def f_R(self, R, Z):
        """returns the R component of the vector Grad(psi)/|Grad(psi)|**2."""
        dpsidR = self.psi_func(R, Z, dx=1, grid=False)
        dpsidZ = self.psi_func(R, Z, dy=1, grid=False)
        return dpsidR / (dpsidR ** 2 + dpsidZ ** 2)

    @handleMultiLocationArray
    def f_Z(self, R, Z):
        """returns the Z component of the vector Grad(psi)/|Grad(psi)|**2."""
        dpsidR = self.psi_func(R, Z, dx=1, grid=False)
        dpsidZ = self.psi_func(R, Z, dy=1, grid=False)
        return dpsidZ / (dpsidR ** 2 + dpsidZ ** 2)

    @handleMultiLocationArray
    def Bp_R(self, R, Z):
        """returns the R component of the poloidal magnetic field."""
        return self.psi_func(R, Z, dy=1, grid=False) / R

    @handleMultiLocationArray
    def Bp_Z(self, R, Z):
        """returns the Z component of the poloidal magnetic field."""
        return -self.psi_func(R, Z, dx=1, grid=False) / R

    @handleMultiLocationArray
    def fpol(self, psi):
        """poloidal current function,
        returns fpol such that B_toroidal = fpol/R"""
        return self.f_spl(psi * self.f_psi_sign)

    @handleMultiLocationArray
    def fpolprime(self, psi):
        """psi-derivative of fpol"""
        return self.fprime_spl(psi * self.f_psi_sign)

    @handleMultiLocationArray
    def pressure(self, psi):
        """Plasma pressure in Pascals"""
        if self.p_spl is None:
            return None
        return self.p_spl(psi * self.f_psi_sign)

    @property
    def Bt_axis(self):
        """Calculate toroidal field on axis"""
        return self.fpol(self.psi_axis) / self.o_point.R


def read_geqdsk(
    filehandle, settings=None, nonorthogonal_settings=None, make_regions=True
):
    """
    Read geqdsk formatted data from a file object, returning
    a TokamakEquilibrium object

    Inputs
    ------
    filehandle   A file handle to read
    settings     dict passed to TokamakEquilibrium
    nonorthogonal_settings  dict passed to TokamakEquilibrium

    Options
    -------
    reverse_current = bool  Changes the sign of poloidal flux psi
    extrapolate_profiles = bool   Extrapolate pressure using exponential
    """

    if settings is None:
        settings = {}

    from ..geqdsk._geqdsk import read as geq_read

    data = geq_read(filehandle)

    # Range of psi normalises psi derivatives
    psi_boundary = data["sibdry"]
    psi_axis = data["simagx"]

    # 1D grid on which fpol is defined. Goes from normalised psi 0 to 1
    psi1D = np.linspace(psi_axis, psi_boundary, data["nx"], endpoint=True)

    R1D = np.linspace(
        data["rleft"], data["rleft"] + data["rdim"], data["nx"], endpoint=True
    )

    Z1D = np.linspace(
        data["zmid"] - 0.5 * data["zdim"],
        data["zmid"] + 0.5 * data["zdim"],
        data["ny"],
        endpoint=True,
    )

    psi2D = data["psi"]

    if settings.get("reverse_current", False):
        warnings.warn("Reversing the sign of the poloidal field")
        psi2D *= -1.0
        psi1D *= -1.0

    # Get the wall
    if "rlim" in data and "zlim" in data:
        wall = list(zip(data["rlim"], data["zlim"]))
    else:
        wall = None

    pressure = data["pres"]
    fpol = data["fpol"]

    if settings.get("extrapolate_profiles", False):
        # Use an exponential decay for the pressure, based on
        # the value and gradient at the plasma edge
        dpdpsi = (pressure[-1] - pressure[-2]) / (psi1D[-1] - psi1D[-2])
        p0 = pressure[-1]
        # Extend the array out to normalised psi of 1.2
        # Exclude first point since duplicates last point in core
        psiSOL = np.linspace(0.0, 0.2 * (psi1D[-1] - psi1D[0]), 50)[1:]

        psi1D = np.concatenate([psi1D, psiSOL])

        # p = p0 * exp( (psi - psi0) * dpdpsi / p0)
        pressure = np.concatenate([pressure, p0 * np.exp(psiSOL * dpdpsi / p0)])

        # fpol constant in SOL
        fpol = np.concatenate([fpol, np.full(psiSOL.shape, fpol[-1])])

    result = TokamakEquilibrium(
        R1D,
        Z1D,
        psi2D,
        psi1D,
        fpol,
        pressure=pressure,
        wall=wall,
        make_regions=make_regions,
        settings=settings,
        nonorthogonal_settings=nonorthogonal_settings,
    )

    # Store geqdsk input as a string in the TokamakEquilibrium object so we can save it
    # in BoutMesh.writeGridFile
    # reset to beginning of file
    filehandle.seek(0)
    # read file as a single string and store in result
    result.geqdsk_input = filehandle.read()
    # also save filename, if it exists
    if hasattr(filehandle, "name"):
        result.geqdsk_filename = filehandle.name

    return result
