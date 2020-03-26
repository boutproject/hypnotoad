# Generate grids for tokamak configurations
#
# 

import numpy as np
from scipy import interpolate
import warnings

from .equilibrium import Equilibrium  # Base class 
from .equilibrium import Point2D 

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

    def __init__(self, R1D, Z1D, psi2D, psi1D, fpol1D, wall=[]):
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

    
