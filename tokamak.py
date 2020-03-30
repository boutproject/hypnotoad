# Generate grids for tokamak configurations
#
# 

import numpy as np
from scipy import interpolate
from scipy.integrate import solve_ivp
import warnings
from collections import OrderedDict
import functools

from .equilibrium import Equilibrium, EquilibriumRegion, Point2D, setDefault
from .hypnotoad_options import HypnotoadOptions
from .mesh import MultiLocationArray

from . import critical
from . import polygons

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
    default_options = HypnotoadOptions.add(

        # Radial grid cell numbers
        nx_core = 5,   # Number of radial points in the core
        nx_pf = None,  # Number of radial points in the PF region
        nx_pf_lower = None,
        nx_pf_upper = None,
        nx_sol = 5, # Number of radial points in the SOL
        
        # Poloidal grid cell numbers
        ny_inner_lower_divertor = 4,
        ny_inner_sol = 4,
        ny_inner_upper_divertor = 4,
        ny_outer_upper_divertor = 4,
        ny_outer_sol = 4,
        ny_outer_lower_divertor = 4,

        # Normalised poloidal flux ranges.
        psinorm_core = 0.9,
        psinorm_sol = 1.1,
        psinorm_sol_inner = None,    # Default: psinorm_sol
        psinorm_pf = None,           # Default: psinorm_core
        psinorm_pf_lower = None,     # Default: psinorm_pf
        psinorm_pf_upper = None,     # Default: psinorm_pf

        # Poloidal flux ranges.
        # These are the values which are used in the mesh generation
        # The default values comes from the psinorm values, but the user
        # can override these defaults.
        psi_core = None,         # Default: psinorm_core
        psi_sol = None,          # Default: psnorm_sol
        psi_sol_inner = None,    # Default: psinorm_sol_inner
        psi_pf_lower = None,     # Default: psinorm_pf_lower
        psi_pf_upper = None)     # Default: psinorm_pf_upper
    
    def __init__(self, R1D, Z1D, psi2D, psi1D, fpol1D,
                 wall=None, psi_axis=None, dct=False, **kwargs):
        
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

        if len(fpol1D) > 0:
            # Spline for interpolation of f = R*Bt
            self.f_spl = interpolate.InterpolatedUnivariateSpline(psi1D, fpol1D, ext=3)
            # ext=3 specifies that boundary values are used outside range
            
            # Spline representing the derivative of f
            self.fprime_spl = self.f_spl.derivative()
        else:
            self.f_spl = lambda psi: 0.0
            self.fprime_spl = lambda psi: 0.0

        # Find critical points (O- and X-points)
        R2D, Z2D = np.meshgrid(R1D, Z1D, indexing='ij')
        opoints, xpoints = critical.find_critical(R2D, Z2D, psi2D)
        if len(opoints) == 0:
            warnings.warn("No O-points found in TokamakEquilibrium input")
        else:
            if psi_axis is None:
                psi_axis = opoints[0][2] # Psi on magnetic axis
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
            offset = 1e-2 # in m
            wall = [(self.Rmin + offset, self.Zmin + offset),
                    (self.Rmax - offset, self.Zmin + offset),
                    (self.Rmax - offset, self.Zmax - offset),
                    (self.Rmin + offset, self.Zmax - offset)]
            
        if polygons.clockwise(wall):
            wall = wall[::-1] # Reverse, without modifying input list (which .reverse() would)
        self.wall = [Point2D(r,z) for r,z in wall]

        self.user_options = TokamakEquilibrium.default_options
        # Set sensible defaults for options
        self.user_options.set(
            xpoint_poloidal_spacing_length = 5.e-2,
            nonorthogonal_xpoint_poloidal_spacing_length = 5.e-2,
            follow_perpendicular_rtol = 2.e-8,
            follow_perpendicular_atol = 1.e-8,
            refine_width = 1.e-2,
            refine_atol = 2.e-8,
            refine_methods = "integrate+newton",
            finecontour_Nfine = 100,
            finecontour_diagnose = False,
            finecontour_maxits = 200,
            curvature_type = 'curl(b/B) with x-y derivatives')
        
        super().__init__(**kwargs)
        
    def optionsTableStr(self):
        """Return a string containing a table of options set"""
        formatstring = '{:<50}|  {:<30}\n'

        # Header
        result = ('\nOptions\n=======\n'
                  + formatstring.format('Name', 'Value')
                  + '-'*80
                  + "\n")

        # Row for each value
        for name, value in sorted(self.user_options.items()):
            valuestring = str(value)
            if value == self.default_options[name]:
                valuestring += '\t(default)'
            result += formatstring.format(name, valuestring)
        return result
    
    def findLegs(self, xpoint, radius=0.01, step=0.01):
        """Find the divertor legs coming from a given X-point
        
        xpoint   Point2D object giving position
        radius   Search distance from X-point, in meters
        step     Integration step size, in meters
        """

        psi_sep = self.psi(xpoint.R, xpoint.Z) # Value of psi

        # Draw a half-circle around this X-point
        if xpoint.Z < self.o_point.Z:
            # If the X-point is below the O-point, choose only angles
            # between pi and 2pi
            angles = np.linspace(np.pi, 2.*np.pi, 50, endpoint=True)
        else:
            # X-point above O-point
            angles = np.linspace(0.0, np.pi, 50, endpoint=True)
            
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
            dr = rvals[ind+1] - r0
            dz = zvals[ind+1] - z0

            s1 = 0.0
            s2 = 1.0
            psi1 = psivals[ind]  # s = 0
            psi2 = psivals[ind+1] # s = 1
            
            while s2 - s1 > 1e-5:
                smid = 0.5 * (s1 + s2)
                psi_mid = self.psi(r0 + smid * dr,
                                   z0 + smid * dz)

                if (psi_mid - psi_sep) * (psi1 - psi_sep) < 0.0:
                    # Between psi_mid and psi1
                    psi2 = psi_mid
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
            line = [xpoint] # Start with the X-point
            
            Br = self.Bp_R(*leg)
            Bz = self.Bp_Z(*leg)
            # Dot product vector from X-point to leg with Bp
            # The sign is used to tell which way to integrate
            sign = np.sign((leg[0] - xpoint.R) * Br +
                           (leg[1] - xpoint.Z) * Bz)

            # Integrate in this direction until the wall is intersected
            # This is affected by sign, which determines which way to integrate
            def dpos_dl(distance, pos):
                r = pos[0]
                z = pos[1]
                Br = self.Bp_R(r, z)
                Bz = self.Bp_Z(r, z)
                B = np.sqrt(Br**2 + Bz**2)
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
                    line.append(intersect) # Put the intersection in the line
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
        return {"inner": leg_lines[0],
                "outer": leg_lines[1]}
        
    def makeRegions(self, **kwargs):
        assert 0 < len(self.x_points) <= 2
        assert self.psi_axis is not None
        
        self.user_options = self.default_options.copy()
        self.user_options = self.user_options.push(kwargs)

        # Radial grid cells in PF regions
        setDefault(self.user_options, 'nx_pf', self.user_options.nx_core)
        setDefault(self.user_options, 'nx_pf_lower', self.user_options.nx_pf)
        setDefault(self.user_options, 'nx_pf_upper', self.user_options.nx_pf)

        # Normalised psi values
        setDefault(self.user_options, 'psinorm_pf', self.user_options.psinorm_core)
        setDefault(self.user_options, 'psinorm_pf_lower', self.user_options.psinorm_pf)
        setDefault(self.user_options, 'psinorm_pf_upper', self.user_options.psinorm_pf)
        
        setDefault(self.user_options, 'psinorm_sol_inner', self.user_options.psinorm_sol)
        
        # psi values
        def psinorm_to_psi(psinorm):
            if psinorm is None:
                return None
            return self.psi_axis + psinorm * (self.psi_sep[0] - self.psi_axis)
        
        setDefault(self.user_options, 'psi_core',
                   psinorm_to_psi(self.user_options.psinorm_core))
        setDefault(self.user_options, 'psi_sol',
                   psinorm_to_psi(self.user_options.psinorm_sol))
        setDefault(self.user_options, 'psi_sol_inner',
                   psinorm_to_psi(self.user_options.psinorm_sol_inner))
        setDefault(self.user_options, 'psi_pf_lower',
                   psinorm_to_psi(self.user_options.psinorm_pf_lower))
        setDefault(self.user_options, 'psi_pf_upper',
                   psinorm_to_psi(self.user_options.psinorm_pf_upper))

        setDefault(self.user_options, 'poloidal_spacing_delta_psi',
                np.abs((self.user_options.psi_core - self.user_options.psi_sol)/20.))
        
        # Print the table of options
        print(self.optionsTableStr())
        
        self.regions = OrderedDict()
        
        if len(self.x_points) == 1:
            # Single null. Could be lower or upper

            # Lower Single Null
            
            # Default options for equilibrium regions
            # Here the total number of poloidal grid cells
            setDefault(self.options, 'N_norm', (self.user_options.ny_inner_lower_divertor +
                                                self.user_options.ny_inner_sol +
                                                self.user_options.ny_outer_sol +
                                                self.user_options.ny_outer_lower_divertor))

            # Find lines along the legs from X-point to target
            legs = self.findLegs(self.x_points[0])

            # Move the first point of each leg slightly away from the X-point
            diff = 0.1 #100.*self.user_options.refine_atol
            for leg in legs.values():
                leg[0] = diff * leg[1] + (1.0 - diff) * leg[0] 

            # wall_vector = {"inner": self.wallVector(legs["inner"][-1]),
            #               "outer": self.wallVector(legs["outer"][-1])}
            
            self.regions['inner_lower_divertor'] = EquilibriumRegion(
                self,
                'inner_lower_divertor',     # Name
                2,   # nSegments, the number of radial regions
                self.user_options,
                self.options.push(
                    {"nx": [self.user_options.nx_pf_lower, self.user_options.nx_sol],
                     "ny": self.user_options.ny_inner_lower_divertor,
                     "kind": "wall.X"}),
                # The following arguments are passed through to PsiContour
                legs["inner"],  # list of Point2D objects on the line
                self.psi,   # Function to calculate the poloidal flux
                self.psi_sep[0])
            
            # Average radial grid spacing in each region
            dpsidi_sol = (self.user_options.psi_sol - self.psi_sep[0]) / self.user_options.nx_sol
            dpsidi_core = (self.psi_sep[0] - self.user_options.psi_core) / self.user_options.nx_core
            dpsidi_pf = (self.psi_sep[0] - self.user_options.psi_pf_lower) / self.user_options.nx_pf_lower

            # Get the smallest absolute grid spacing for the separatrix
            dpsidi_sep = min([dpsidi_sol, dpsidi_core, dpsidi_pf], key=abs)

            # decrease (assuming the factor is <1) the spacing around the separatrix by the
            # factor psi_spacing_separatrix_multiplier
            if self.user_options.psi_spacing_separatrix_multiplier is not None:
                dpsidi_sep = self.user_options.psi_spacing_separatrix_multiplier * dpsidi_sep

            # Private Flux
            pf_psi_func = self.getPolynomialGridFunc(
                self.user_options.nx_pf_lower,
                self.user_options.psi_pf_lower, self.psi_sep[0], grad_upper=dpsidi_sep)

            pf_psi_vals = self.make1dGrid(self.user_options.nx_pf_lower, pf_psi_func)
            
            # Core
            core_psi_func = self.getPolynomialGridFunc(
                self.user_options.nx_core,
                self.user_options.psi_core, self.psi_sep[0], grad_upper=dpsidi_sep)

            core_psi_vals = self.make1dGrid(self.user_options.nx_core, core_psi_func)

            # SOL
            sol_psi_func = self.getPolynomialGridFunc(
                self.user_options.nx_sol,
                self.psi_sep[0], self.user_options.psi_sol, grad_upper=dpsidi_sep)

            sol_psi_vals = self.make1dGrid(self.user_options.nx_sol, sol_psi_func)
            
            def setupRegion(name, psi_vals1, psi_vals2, reverse):
                r = self.regions[name]
                r.psi_vals = [psi_vals1, psi_vals2]
                r.separatrix_radial_index = 1
                if reverse:
                    r.reverse()
                    r.xPointsAtEnd[1] = self.x_points[0]
                    r.wallSurfaceAtStart = wall_vectors[name]
                else:
                    r.xPointsAtStart[1] = self.x_points[0]
                    r.wallSurfaceAtEnd = wall_vectors[name]

            r = self.regions['inner_lower_divertor']
            r.psi_vals = [pf_psi_vals, sol_psi_vals]
            r.separatrix_radial_index = 1
            r.reverse()
            r.xPointsAtEnd[1] = self.x_points[0]
            r.wallSurfaceAtStart = [0., 0.] #wall_vector["inner"]
            
            
        else:
            # Double null
            
            # Regions ordered clockwise, starting lower inner divertor
            legnames = ['inner_lower_divertor', 'inner_sol', 'inner_upper_divertor',
                        'outer_upper_divertor', 'outer_sol', 'outer_lower_divertor']
            kinds = ['wall.X', 'X.X', 'X.wall',
                     'wall.X', 'X.X', 'X.wall']
            
    def handleMultiLocationArray(getResult):
        @functools.wraps(getResult)
        # Define a function which handles MultiLocationArray arguments
        def handler(self, R, Z):
            if isinstance(R, MultiLocationArray):
                assert isinstance(Z, MultiLocationArray), 'if R is a MultiLocationArray, then Z must be as well'
                result = MultiLocationArray(R.nx, R.ny)
                if R.centre is not None and Z.centre is not None:
                    result.centre = getResult(self, R.centre, Z.centre)
                    
                if R.xlow is not None and Z.xlow is not None:
                    result.xlow = getResult(self, R.xlow, Z.xlow)
                        
                if R.ylow is not None and Z.ylow is not None:
                    result.ylow = getResult(self, R.ylow, Z.ylow)

                if R.corners is not None and Z.corners is not None:
                    result.corners = getResult(self, R.corners, Z.corners)
            else:
                result = getResult(self, R, Z)
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
        return dpsidR / np.sqrt(dpsidR**2 + dpsidZ**2)

    @handleMultiLocationArray
    def f_Z(self, R, Z):
        """returns the Z component of the vector Grad(psi)/|Grad(psi)|**2."""
        dpsidR = self.psi_func(R, Z, dx=1, grid=False)
        dpsidZ = self.psi_func(R, Z, dy=1, grid=False)
        return dpsidZ / np.sqrt(dpsidR**2 + dpsidZ**2)

    @handleMultiLocationArray
    def Bp_R(self, R, Z):
        """returns the R component of the poloidal magnetic field."""
        return -self.psi_func(R, Z, dy=1, grid=False) / R

    @handleMultiLocationArray
    def Bp_Z(self, R, Z):
        """returns the Z component of the poloidal magnetic field."""
        return self.psi_func(R, Z, dx=1, grid=False) / R

    def fpol(self, psi):
        """poloidal current function, 
        returns fpol such that B_toroidal = fpol/R"""
        return self.f_spl(psi)

    def fpolprime(self, psi):
        """psi-derivative of fpol
        """
        return self.fprime_spl(psi)
    

def read_geqdsk(filehandle):
    """
    Read geqdsk formatted data from a file object, returning
    a TokamakEquilibrium object
    """

    from ._geqdsk import read as geq_read
    
    data = geq_read(filehandle)
    
    # Range of psi normalises psi derivatives
    psi_boundary = data["sibdry"]
    psi_axis = data["simagx"]

    # 1D grid on which fpol is defined. Goes from normalised psi 0 to 1
    psi1D = np.linspace(psi_axis, psi_boundary,
                        data["nx"], endpoint=True)

    R1D = np.linspace(data["rleft"], data["rleft"] + data["rdim"],
                      data["nx"], endpoint=True)
                            
    Z1D = np.linspace(data["zmid"] - 0.5*data["zdim"],
                      data["zmid"] + 0.5*data["zdim"],
                      data["ny"], endpoint=True)

    # Get the wall 
    if "rlim" in data and "zlim" in data:
        wall = list(zip(data["rlim"], data["zlim"]))
    else:
        wall = []
    
    return TokamakEquilibrium(R1D,
                              Z1D,
                              data["psi"],  # psi2D
                              psi1D, 
                              data["fpol"], # fpol1D
                              wall=wall)

def example():
    nx = 65
    ny = 65
    
    r1d = np.linspace(1.2, 1.8, nx)
    z1d = np.linspace(-0.5, 0.5, ny)
    r2d, z2d = np.meshgrid(r1d, z1d, indexing='ij')

    r0 = 1.5
    z0 = -0.3

    # This has two O-points, and one x-point at (r0, z0)
    def psi_func(R,Z):
        return np.exp(-((R - r0)**2 + (Z - z0 - 0.3)**2)/0.3**2) + np.exp(-((R - r0)**2 + (Z - z0 + 0.3)**2)/0.3**2)
    
    eq = TokamakEquilibrium(r1d, z1d, psi_func(r2d, z2d),
                            [], []) # psi1d, fpol

    eq.makeRegions(psinorm_pf=0.9, psinorm_sol=1.1)
    
    from .mesh import BoutMesh
    mesh = BoutMesh(eq)
    mesh.geometry()
    
    import matplotlib.pyplot as plt

    eq.plotPotential()
    #eq.addWallToPlot()
    plt.plot(*eq.x_points[0], 'rx')
    
    mesh.plotPoints(xlow=True, ylow=True, corners=True)
    
    plt.show()
