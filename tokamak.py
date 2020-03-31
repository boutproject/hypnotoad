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
        nx_pf_lower = None,   # Default nx_pf
        nx_pf_upper = None,   # Default nx_pf
        nx_sol = 5, # Number of radial points in the SOL
        nx_sol_inner = None,  # Default nx_sol
        nx_sol_outer = None,  # Default nx_sol
        
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
            refine_methods = "integrate+newton, integrate, none",
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

        # Radial grid cells in SOL regions (for double null)
        setDefault(self.user_options, 'nx_sol_inner', self.user_options.nx_sol)
        setDefault(self.user_options, 'nx_sol_outer', self.user_options.nx_sol)
        
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

            if self.x_points[0].Z < self.o_point.Z:
                # Lower Single Null

                mesh_type = "lower single null"
                inner_leg_name = "inner_lower_divertor"
                outer_leg_name = "outer_lower_divertor"
                ny_inner_leg = self.user_options.ny_inner_lower_divertor
                ny_outer_leg = self.user_options.ny_outer_lower_divertor
                nx_pf = self.user_options.nx_pf_lower
                psi_pf = self.user_options.psi_pf_lower
            else:
                # Upper Single Null

                mesh_type = "upper single null"
                inner_leg_name = 'inner_upper_divertor'
                outer_leg_name = 'outer_upper_divertor'
                ny_inner_leg = self.user_options.ny_inner_upper_divertor
                ny_outer_leg = self.user_options.ny_outer_upper_divertor
                nx_pf = self.user_options.nx_pf_upper
                psi_pf = self.user_options.psi_pf_upper

            print("Generating tokamak mesh type: " + mesh_type)
                
            # Default options for equilibrium regions
            # Here the total number of poloidal grid cells
            setDefault(self.options, 'N_norm', (ny_inner_leg +
                                                self.user_options.ny_inner_sol +
                                                self.user_options.ny_outer_sol +
                                                ny_outer_leg))

            # Find lines along the legs from X-point to target
            legs = self.findLegs(self.x_points[0])

            # Move the first point of each leg slightly away from the X-point
            diff = 0.1 #100.*self.user_options.refine_atol
            for leg in legs.values():
                leg[0] = diff * leg[1] + (1.0 - diff) * leg[0] 

            # wall_vector = {"inner": self.wallVector(legs["inner"][-1]),
            #               "outer": self.wallVector(legs["outer"][-1])}
            
            self.regions[inner_leg_name] = EquilibriumRegion(
                self,
                inner_leg_name,
                2,   # nSegments, the number of radial regions
                self.user_options,
                self.options.push(
                    {"nx": [nx_pf, self.user_options.nx_sol],
                     "ny": ny_inner_leg,
                     "kind": "wall.X"}),
                # The following arguments are passed through to PsiContour
                legs["inner"],  # list of Point2D objects on the line
                self.psi,   # Function to calculate the poloidal flux
                self.psi_sep[0])

            self.regions[outer_leg_name] = EquilibriumRegion(
                self,
                outer_leg_name,
                2,   # nSegments, the number of radial regions
                self.user_options,
                self.options.push(
                    {"nx": [nx_pf, self.user_options.nx_sol],
                     "ny": ny_outer_leg,
                     "kind": "X.wall"}),
                # The following arguments are passed through to PsiContour
                legs["outer"],  # list of Point2D objects on the line
                self.psi,   # Function to calculate the poloidal flux
                self.psi_sep[0])
            
            # Core region
            ax = None
            if False:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                self.plotPotential(axis=ax)
            
            r0, z0 = self.o_point.R, self.o_point.Z # Location of O-point
            # Find the angle of the X-point from the O-point
            angle_xpt = np.arctan2(self.x_points[0].Z -self.o_point.Z,
                                   self.x_points[0].R -self.o_point.R)
            
            npoints = 100
            dtheta = np.pi / npoints

            points = [Point2D(*critical.find_psisurface(self,
                                                        r0, z0,
                                                        r0 + 8.*np.cos(angle),
                                                        z0 + 8.*np.sin(angle),
                                                        psival=self.psi_sep[0],
                                                        axis=ax))
                      for angle in np.linspace(angle_xpt + dtheta,
                                               angle_xpt + 2*np.pi - dtheta,
                                               npoints)]
            if ax is not None:
                plt.show()

            # Add points to the beginning and end near (but not at) the X-point
            diff = 0.1
            points = ([(1.0 - diff) * self.x_points[0] + diff * points[0]] +
                      points +
                      [(1.0 - diff) * self.x_points[0] + diff * points[-1]])

            self.regions['core'] = EquilibriumRegion(
                self,
                'core',     # Name
                2,   # nSegments, the number of radial regions
                self.user_options,
                self.options.push(
                    {"nx": [self.user_options.nx_core, self.user_options.nx_sol],
                     "ny": self.user_options.ny_inner_sol + self.user_options.ny_outer_sol,
                     "kind": "X.X"}),
                # The following arguments are passed through to PsiContour
                points,  # list of Point2D objects on the line
                self.psi,   # Function to calculate the poloidal flux
                self.psi_sep[0])
            
            # Average radial grid spacing in each region
            dpsidi_sol = (self.user_options.psi_sol - self.psi_sep[0]) / self.user_options.nx_sol
            dpsidi_core = (self.psi_sep[0] - self.user_options.psi_core) / self.user_options.nx_core
            dpsidi_pf = (self.psi_sep[0] - psi_pf) / nx_pf

            # Get the smallest absolute grid spacing for the separatrix
            dpsidi_sep = min([dpsidi_sol, dpsidi_core, dpsidi_pf], key=abs)

            # decrease (assuming the factor is <1) the spacing around the separatrix by the
            # factor psi_spacing_separatrix_multiplier
            if self.user_options.psi_spacing_separatrix_multiplier is not None:
                dpsidi_sep = self.user_options.psi_spacing_separatrix_multiplier * dpsidi_sep

            # Private Flux
            pf_psi_func = self.getPolynomialGridFunc(
                nx_pf,
                psi_pf, self.psi_sep[0], grad_upper=dpsidi_sep)

            pf_psi_vals = self.make1dGrid(nx_pf, pf_psi_func)
            
            # Core
            core_psi_func = self.getPolynomialGridFunc(
                self.user_options.nx_core,
                self.user_options.psi_core, self.psi_sep[0], grad_upper=dpsidi_sep)

            core_psi_vals = self.make1dGrid(self.user_options.nx_core, core_psi_func)

            # SOL
            sol_psi_func = self.getPolynomialGridFunc(
                self.user_options.nx_sol,
                self.psi_sep[0], self.user_options.psi_sol, grad_lower=dpsidi_sep)

            sol_psi_vals = self.make1dGrid(self.user_options.nx_sol, sol_psi_func)


            if mesh_type == "lower single null":
                # Note: The grid is arranged so that the y index goes clockwise
                # 
                r = self.regions['inner_lower_divertor']
                r.psi_vals = [pf_psi_vals, sol_psi_vals]
                r.separatrix_radial_index = 1
                r.reverse()  # Leg points are ordered from X-point -> target, so reverse
                r.xPointsAtEnd[1] = self.x_points[0]
                r.wallSurfaceAtStart = [0., 0.] #wall_vector["inner"]
                
                r = self.regions['outer_lower_divertor']
                r.psi_vals = [pf_psi_vals, sol_psi_vals]
                r.separatrix_radial_index = 1
                r.xPointsAtStart[1] = self.x_points[0]
                r.wallSurfaceAtEnd = [0., 0.] #wall_vector["outer"]
                
                r = self.regions['core']
                r.psi_vals = [core_psi_vals, sol_psi_vals]
                r.separatrix_radial_index = 1
                r.reverse()  # Core points are created anti-clockwise, so reverse
                r.xPointsAtEnd[1] = self.x_points[0]
                r.xPointsAtStart[1] = self.x_points[0]
                
                # Connections between regions
                # inner lower PF -> outer lower PF
                self.makeConnection('inner_lower_divertor', 0, 'outer_lower_divertor', 0)
                
                # inner lower PF -> Core
                self.makeConnection('inner_lower_divertor', 1, 'core', 1)
                
                # Core -> outer lower PF
                self.makeConnection('core', 1, 'outer_lower_divertor', 1)
                
                # Core -> core
                self.makeConnection('core', 0, 'core', 0)
            else:
                # Upper single null
                # The mesh is still arranged clockwise, meaning that y indexing starts
                # at the outer upper divertor, and ends at the inner upper divertor.
                #
                # The points in the legs are ordered from X-point to target,
                # so here the outer leg needs to be reversed rather than the inner leg
                # 
                r = self.regions['inner_upper_divertor']
                r.psi_vals = [pf_psi_vals, sol_psi_vals]
                r.separatrix_radial_index = 1
                r.xPointsAtEnd[1] = self.x_points[0]
                r.wallSurfaceAtStart = [0., 0.] #wall_vector["inner"]
                
                r = self.regions['outer_upper_divertor']
                r.psi_vals = [pf_psi_vals, sol_psi_vals]
                r.separatrix_radial_index = 1
                r.reverse()
                r.xPointsAtStart[1] = self.x_points[0]
                r.wallSurfaceAtEnd = [0., 0.] #wall_vector["outer"]
                
                r = self.regions['core']
                r.psi_vals = [core_psi_vals, sol_psi_vals]
                r.separatrix_radial_index = 1
                r.reverse()
                r.xPointsAtEnd[1] = self.x_points[0]
                r.xPointsAtStart[1] = self.x_points[0]
                
                # Connections between regions
                # outer lower PF -> inner lower PF
                self.makeConnection('outer_upper_divertor', 0, 'inner_upper_divertor', 0)
                
                # outer lower PF -> Core
                self.makeConnection('outer_upper_divertor', 1, 'core', 1)
                
                # Core -> inner upper PF
                self.makeConnection('core', 1, 'inner_upper_divertor', 1)
                
                # Core -> core
                self.makeConnection('core', 0, 'core', 0)
            
        else:
            # Double null, connected or disconnected    

            if self.x_points[0].Z < self.o_point.Z:
                lower_xpt_ind, upper_xpt_ind = 0, 1
            else:
                lower_xpt_ind, upper_xpt_ind = 1, 0
            
            lower_x_point = self.x_points[lower_xpt_ind]
            upper_x_point = self.x_points[upper_xpt_ind]
            lower_psi = self.psi_sep[lower_xpt_ind]
            upper_psi = self.psi_sep[upper_xpt_ind]
            
            # Find lines along the legs from X-point to target
            lower_legs = self.findLegs(lower_x_point)
            upper_legs = self.findLegs(upper_x_point)
            
            # Average radial grid spacing in each region
            dpsidi_sol_inner = (self.user_options.psi_sol_inner - self.psi_sep[0]) / self.user_options.nx_sol_inner
            dpsidi_sol_outer = (self.user_options.psi_sol - self.psi_sep[0]) / self.user_options.nx_sol_outer
            dpsidi_core = (self.psi_sep[0] - self.user_options.psi_core) / self.user_options.nx_core
            dpsidi_pf_upper = (upper_psi - self.user_options.psi_pf_upper) / self.user_options.nx_pf_upper
            dpsidi_pf_lower = (lower_psi - self.user_options.psi_pf_lower) / self.user_options.nx_pf_lower
            
            # Get the smallest absolute grid spacing for the separatrix
            dpsidi_sep = min([dpsidi_sol_inner, dpsidi_sol_outer,
                              dpsidi_core,
                              dpsidi_pf_upper, dpsidi_pf_lower], key=abs)

            # Number of points in the inter-separatrix region
            nx_inter_sep = np.rint(abs((upper_psi - lower_psi) / dpsidi_sep))
            if nx_inter_sep == 0:
                print("Generating a connected double null")

                # Radial segments i.e. gridded ranges of poloidal flux
                segments = {"core": {'nx': self.user_options.nx_core,
                                     'psi_start': self.user_options.psi_core,
                                     'psi_end': self.psi_sep[0],
                                     'grad_end': dpsidi_sep},
                            "upper_pf": {'nx': self.user_options.nx_pf_upper,
                                         'psi_start':self.user_options.psi_pf_upper,
                                         'psi_end': upper_psi,
                                         'grad_end': dpsidi_sep},
                            "lower_pf": {'nx': self.user_options.nx_pf_lower,
                                         'psi_start': self.user_options.psi_pf_lower,
                                         'psi_end': lower_psi,
                                         'grad_end': dpsidi_sep},
                            "inner_sol": {'nx': self.user_options.nx_sol_inner,
                                          'psi_start': self.psi_sep[0],
                                          'psi_end': self.user_options.psi_sol_inner,
                                          'grad_start': dpsidi_sep},
                            "outer_sol": {'nx': self.user_options.nx_sol_outer,
                                          'psi_start': self.psi_sep[0],
                                          'psi_end': self.user_options.psi_sol,
                                          'grad_start': dpsidi_sep}}
                
                # Description of each poloidal region
                leg_regions = {'inner_lower_divertor': {'segments': ["lower_pf", "inner_sol"],
                                                        'ny': self.user_options.ny_inner_lower_divertor,
                                                        'kind': "wall.X",
                                                        'points': lower_legs["inner"][::-1],
                                                        'psi': lower_psi,
                                                        'wall_at_start': [0,0],
                                                        'xpoint_at_end': lower_x_point},
                               
                               'outer_lower_divertor': {'segments': ["lower_pf", "outer_sol"],
                                                        'ny': self.user_options.ny_outer_lower_divertor,
                                                        'kind': "X.wall",
                                                        'points': lower_legs["outer"],
                                                        'psi': lower_psi,
                                                        'xpoint_at_start': lower_x_point,
                                                        'wall_at_end': [0,0]},
                               
                               'inner_upper_divertor': {'segments':["upper_pf", "inner_sol"],
                                                        'ny': self.user_options.ny_inner_upper_divertor,
                                                        'kind': "X.wall",
                                                        'points': upper_legs["inner"],
                                                        'psi': upper_psi,
                                                        'xpoint_at_start': upper_x_point,
                                                        'wall_at_end': [0,0]},
                               
                               'outer_upper_divertor': {'segments': ["upper_pf", "outer_sol"],
                                                        'ny': self.user_options.ny_outer_upper_divertor,
                                                        'kind': "wall.X",
                                                        'points': upper_legs["outer"][::-1],
                                                        'psi': upper_psi,
                                                        'wall_at_start': [0,0],
                                                        'xpoint_at_end': upper_x_point}}
                
                core_regions = {'inner_core': {'segments': ["core", "inner_sol"],
                                               'ny': self.user_options.ny_inner_sol,
                                               'kind': "X.X",
                                               'xpoint_at_start': lower_x_point,
                                               'xpoint_at_end': upper_x_point,
                                               'psi_at_start': lower_psi,
                                               'psi_at_end': upper_psi},
                                               
                               'outer_core': {'segments': ["core", "outer_sol"],
                                              'ny': self.user_options.ny_outer_sol,
                                              'kind': "X.X",
                                              'xpoint_at_start': upper_x_point,
                                              'xpoint_at_end': lower_x_point,
                                              'psi_at_start': upper_psi,
                                              'psi_at_end': lower_psi}}
                
            else:
                print("Generating a disconnected double null")
                
                # Radial segments i.e. gridded ranges of poloidal flux
                segments = {"core": {'nx': self.user_options.nx_core,
                                     'psi_start': self.user_options.psi_core,
                                     'psi_end': self.psi_sep[0],
                                     'grad_end': dpsidi_sep},
                            "upper_pf": {'nx': self.user_options.nx_pf_upper,
                                         'psi_start':self.user_options.psi_pf_upper,
                                         'psi_end': upper_psi,
                                         'grad_end': dpsidi_sep},
                            "lower_pf": {'nx': self.user_options.nx_pf_lower,
                                         'psi_start': self.user_options.psi_pf_lower,
                                         'psi_end': lower_psi,
                                         'grad_end': dpsidi_sep},
                            "near_sol": {'nx': nx_inter_sep,
                                         'psi_start': self.psi_sep[0],
                                         'psi_end': self.psi_sep[1],
                                         'grad_start': dpsidi_sep,
                                         'grad_end': dpsidi_sep},
                            "inner_sol": {'nx': self.user_options.nx_sol_inner - nx_inter_sep,
                                          'psi_start': self.psi_sep[1],
                                          'psi_end': self.user_options.psi_sol_inner,
                                          'grad_start': dpsidi_sep},
                            "outer_sol": {'nx': self.user_options.nx_sol_outer - nx_inter_sep,
                                          'psi_start': self.psi_sep[1],
                                          'psi_end': self.user_options.psi_sol,
                                          'grad_start': dpsidi_sep}}
                
                core_regions = {'inner_core': {'segments': ["core", "near_sol", "inner_sol"],
                                               'ny': self.user_options.ny_inner_sol,
                                               'kind': "X.X",
                                               'xpoint_at_start': lower_x_point,
                                               'xpoint_at_end': upper_x_point,
                                               'psi_at_start': lower_psi,
                                               'psi_at_end': upper_psi},
                                               
                               'outer_core': {'segments': ["core", "near_sol", "outer_sol"],
                                              'ny': self.user_options.ny_outer_sol,
                                              'kind': "X.X",
                                              'xpoint_at_start': upper_x_point,
                                              'xpoint_at_end': lower_x_point,
                                              'psi_at_start': upper_psi,
                                              'psi_at_end': lower_psi}} 
                
                if self.x_points[0] == lower_x_point:
                    # Lower double null
                    inner_lower_segments = ["lower_pf", "near_sol", "inner_sol"]
                    outer_lower_segments = ["lower_pf", "near_sol", "outer_sol"]
                    
                    inner_upper_segments = ["upper_pf", "inner_sol"]
                    outer_upper_segments = ["upper_pf", "outer_sol"]
                    
                else:
                    # Upper double null
                    inner_lower_segments = ["lower_pf", "inner_sol"]
                    outer_lower_segments = ["lower_pf", "outer_sol"]
                    
                    inner_upper_segments = ["upper_pf", "near_sol", "inner_sol"]
                    outer_upper_segments = ["upper_pf", "near_sol", "outer_sol"]
                    
                    
                leg_regions = {'inner_lower_divertor': {'segments': inner_lower_segments,
                                                        'ny': self.user_options.ny_inner_lower_divertor,
                                                        'kind': "wall.X",
                                                        'points': lower_legs["inner"][::-1],
                                                        'psi': lower_psi,
                                                        'wall_at_start': [0,0],
                                                        'xpoint_at_end': lower_x_point},
                               
                               'outer_lower_divertor': {'segments': outer_lower_segments,
                                                        'ny': self.user_options.ny_outer_lower_divertor,
                                                        'kind': "X.wall",
                                                        'points': lower_legs["outer"],
                                                        'psi': lower_psi,
                                                        'xpoint_at_start': lower_x_point,
                                                        'wall_at_end': [0,0]},
                               
                               'inner_upper_divertor': {'segments': inner_upper_segments,
                                                        'ny': self.user_options.ny_inner_upper_divertor,
                                                        'kind': "X.wall",
                                                        'points': upper_legs["inner"],
                                                        'psi': upper_psi,
                                                        'xpoint_at_start': upper_x_point,
                                                        'wall_at_end': [0,0]},
                               
                               'outer_upper_divertor': {'segments': outer_upper_segments,
                                                        'ny': self.user_options.ny_outer_upper_divertor,
                                                        'kind': "wall.X",
                                                        'points': upper_legs["outer"][::-1],
                                                        'psi': upper_psi,
                                                        'wall_at_start': [0,0],
                                                        'xpoint_at_end': upper_x_point}}
                                                        

                    
            
            # Create a new dictionary, which will contain all regions
            # including core and legs
            all_regions = leg_regions.copy()
            
            # number of points in each core region
            npoints = 100
            
            # Loop through core regions, calculate points along the lines,
            # and put into the all_regions dictionary
            for name, region in core_regions.items(): 

                start_x = region['xpoint_at_start']
                end_x = region['xpoint_at_end']
                start_psi = region['psi_at_start']
                end_psi = region['psi_at_end']
                
                # Range of angles. Note: This angle goes anticlockwise
                # so core regions need to be reversed
                start_angle = np.arctan2(start_x.Z - self.o_point.Z,
                                         start_x.R - self.o_point.R)
                end_angle = np.arctan2(end_x.Z - self.o_point.Z,
                                       end_x.R - self.o_point.R)
                if end_angle > start_angle:
                    end_angle -= 2*np.pi

                # Angle offset from the X-point
                dtheta = 0.5 * (end_angle - start_angle) / npoints
                r0, z0 = self.o_point.R, self.o_point.Z  # Location of O-point

                def psival(angle):
                    "Linear interpolation in psi with angle"
                    norm = (angle - start_angle) / (end_angle - start_angle)
                    return norm * end_psi + (1. - norm) * start_psi

                # Iterate in angle from start to end
                points = [Point2D(*critical.find_psisurface(self,
                                                            r0, z0,
                                                            r0 + 8.*np.cos(angle),
                                                            z0 + 8.*np.sin(angle),
                                                            psival = psival(angle)))
                          for angle in np.linspace(start_angle + dtheta,
                                                   end_angle - dtheta,
                                                   npoints)]
                
                # Add points to the beginning and end near (but not at) the X-points
                diff = 0.1
                
                region["points"] = ([(1.0 - diff) * start_x + diff * points[0]] +
                                    points +
                                    [(1.0 - diff) * end_x + diff * points[-1]])

                region["psi"] = 0.5 * (start_psi + end_psi) # Average psi along line
                
                all_regions[name] = region

            
            # Grid each radial segment, put result in psi_vals
            psi_vals = {}
            for name, segment in segments.items():
                psi_func = self.getPolynomialGridFunc(segment["nx"],
                                                      segment["psi_start"],
                                                      segment["psi_end"],
                                                      grad_lower=segment.get("grad_start", None),
                                                      grad_upper=segment.get("grad_end", None))
                
                psi_vals[name] = self.make1dGrid(segment["nx"], psi_func)

            # Normalisation of the grid cells number
            # Set N_norm to the total number of grid cells in y
            setDefault(self.options, 'N_norm', sum([region["ny"]
                                                    for region in all_regions.values()]))

            # Loop through all regions. For each one create a 
            for name, region in all_regions.items():
                eqreg = EquilibriumRegion(self,
                                          name,
                                          len(region["segments"]), # The number of radial regions
                                          self.user_options,
                                          self.options.push(
                                              {"nx": [segments[seg_name]["nx"]
                                                      for seg_name in region["segments"]],
                                               "ny": region["ny"],
                                               "kind": region["kind"]}),
                                          # The following arguments are passed through to PsiContour
                                          region["points"],  # list of Point2D objects on the line
                                          self.psi,   # Function to calculate the poloidal flux
                                          region["psi"])

                # Grids of psi values in each segment
                eqreg.psi_vals = [psi_vals[segment]
                                 for segment in region["segments"]]
                
                eqreg.separatrix_radial_index = 1

                if 'xpoint_at_start' in region:
                    eqreg.xPointsAtStart[1] = region['xpoint_at_start']
                
                if 'xpoint_at_end' in region:
                    eqreg.xPointsAtEnd[1] = region['xpoint_at_end']

                if 'wall_at_start' in region:
                    eqreg.wallSurfaceAtStart = region['wall_at_start']

                if 'wall_at_end' in region:
                    eqreg.wallSurfaceAtEnd = region['wall_at_end']

                self.regions[name] = eqreg
            
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

    if False:
        r0 = 1.5
        z0 = 0.3
        
        # This has two O-points, and one x-point at (r0, z0)
        def psi_func(R,Z):
            return np.exp(-((R - r0)**2 + (Z - z0 - 0.3)**2)/0.3**2) + np.exp(-((R - r0)**2 + (Z - z0 + 0.3)**2)/0.3**2)
    else:
        r0 = 1.5
        z0 = 0.3
        
        # This has two X-points
        def psi_func(R,Z):
            return np.exp(-((R - r0)**2 + Z**2)/0.3**2) + np.exp(-((R - r0)**2 + (Z + 2*z0)**2)/0.3**2) + np.exp(-((R - r0)**2 + (Z - 2*z0)**2)/0.3**2)

        
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
    
    #for region in eq.regions.values():
    #    plt.plot([p.R for p in region.points], [p.Z for p in region.points], '-o')
    
    mesh.plotPoints(xlow=True, ylow=True, corners=True)
    
    plt.show()
