# Generate grids for tokamak configurations
#
# 

import numpy as np
from scipy import interpolate
import warnings
from collections import OrderedDict

from .equilibrium import Equilibrium  # Base class 
from .equilibrium import Point2D 
from .equilibrium import setDefault # For options 
from .hypnotoad_options import HypnotoadOptions

from . import critical
from . import polygons

class TokamakEquilibrium(Equilibrium):
    """
    Represents an axisymmetric tokamak equilibrium

    
    Data members
    - x_points: list of Point2D objects giving the position of the X-points ordered
                from primary X-point (nearest the core) outward
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
    
    def __init__(self, R1D, Z1D, psi2D, psi1D, fpol1D, wall=[], psi_axis=None):
        # Create a 2D spline interpolation for psi
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
        if polygons.clockwise(wall):
            wall = wall[::-1] # Reverse, without modifying input list (which .reverse() would)
        self.wall = [Point2D(r,z) for r,z in wall]
        

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
        
        # Print the table of options
        print(self.optionsTableStr())
        
        self.regions = OrderedDict()
        
        if len(self.x_points) == 1:
            # Single null. Could be lower or upper

            # Lower Single Null
            # Regions ordered clockwise, starting lower inner divertor
            legnames = ['inner_lower_divertor', 'sol', 'outer_lower_divertor']
            kinds = ['wall.X', 'X.X', 'X.wall']

            # Average grid spacing
            dpsidi_sol = (self.user_options.psi_sol - self.psi_sep[0]) / self.user_options.nx_sol
            dpsidi_core = (self.psi_sep[0] - self.user_options.psi_core) / self.user_options.nx_core
            dpsidi_pf = (self.psi_sep[0] - self.user_options.psi_pf_lower) / self.user_options.nx_pf_lower

            
        else:
            # Double null
            
            # Regions ordered clockwise, starting lower inner divertor
            legnames = ['inner_lower_divertor', 'inner_sol', 'inner_upper_divertor',
                        'outer_upper_divertor', 'outer_sol', 'outer_lower_divertor']
            kinds = ['wall.X', 'X.X', 'X.wall',
                     'wall.X', 'X.X', 'X.wall']
            
        

        
    def psi(self, R, Z):
        "Return the poloidal flux at the given (R,Z) location"
        return self.psi_func(R, Z, grid=False)

    def f_R(self, R, Z):
        """returns the R component of the vector Grad(psi)/|Grad(psi)|**2."""
        dpsidR = self.psi_func(R, Z, dx=1, grid=False)
        dpsidZ = self.psi_func(R, Z, dy=1, grid=False)
        return dpsidR / np.sqrt(dpsidR**2 + dpsidZ**2)

    def f_Z(self, R, Z):
        """returns the Z component of the vector Grad(psi)/|Grad(psi)|**2."""
        dpsidR = self.psi_func(R, Z, dx=1, grid=False)
        dpsidZ = self.psi_func(R, Z, dy=1, grid=False)
        return dpsidZ / np.sqrt(dpsidR**2 + dpsidZ**2)

    def Bp_R(self, R, Z):
        """returns the R component of the poloidal magnetic field."""
        return -self.psi_func(R, Z, dy=1, grid=False) / R

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

    import matplotlib.pyplot as plt

    eq.plotPotential()
    plt.plot(*eq.x_points[0], 'rx')
    plt.show()
    
    eq.makeRegions()
    
