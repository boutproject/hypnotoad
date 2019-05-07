#!/usr/bin/env python3
"""
Create a BOUT++ grid for TORPEX from an input file giving coil currents and positions

Input file should contain coil parameters, for each coil:
    R: major radius in metres
    Z: major radius in metres
    I: clockwise current in Amps

Note: positions of cell corners are generated first, grid points are then put in the
centre of the cell.
"""

plotStuff = True

# Torpex parameters
Rmin = 0.8
Rmax = 1.2
Zmin = -.2
Zmax = .2

import numpy
from mesh import Mesh, Point2D
from analyse_equilibrium import findSaddlePoint, findSeparatrix
if plotStuff:
    from matplotlib import pyplot

def TORPEX_wall(theta):
    """
    Return the location of the TORPEX wall parameterized by the angle theta
    anticlockwise around the centre of the vacuum vessel
    """
    # TORPEX wall is a circle radius 0.2 m around (1 m, 0 m)
    awall = 0.2
    Rcentre = 1.
    Zcentre = 0.
    return Point2D(Rcentre + awall*numpy.cos(theta), Zcentre + awall*numpy.sin(theta))

def addWallToPlot(npoints=100):
    theta = numpy.linspace(0., 2.*numpy.pi, npoints+1)
    pyplot.plot(*TORPEX_wall(theta))

def parseInput(filename):
    import yaml
    from collections import namedtuple

    with open(filename, 'r') as inputfile:
        coil_inputs, mesh_inputs = yaml.safe_load_all(inputfile)
    print('Coils:',coil_inputs['Coils'])
    
    Coil = namedtuple('Coil', 'R, Z, I')
    return [Coil(**c) for c in coil_inputs['Coils']], coil_inputs['Bt_axis'], mesh_inputs['Mesh']

def magneticFunctions(coils):
    """
    Calculate toroidal (anticlockwise) component of magnetic vector potential due to coils
    See for example http://physics.usask.ca/~hirose/p812/notes/Ch3.pdf

    Note e_R x e_phi = e_Z
    """
    import sympy
    from sympy.functions.special.elliptic_integrals import elliptic_k, elliptic_e
    import scipy.special

    R,Z = sympy.symbols('R Z')
    mu0 = 4.e-7*sympy.pi

    potential = 0*R

    for coil in coils:
        # little-r is the vector position from the centre of the coil to (R,Z)
        # sinTheta is the angle between r and the axis through the centre of the coil
        rSquared = R**2 + (Z - coil.Z)**2
        r = sympy.sqrt(rSquared)
        sinTheta = R / r
        kSquared = 4*coil.R*r*sinTheta / (rSquared + coil.R**2 + 2*coil.R*r*sinTheta)
        potential += -(
          coil.I*coil.R / sympy.sqrt(r**2 + coil.R**2 + 2*coil.R*r*sinTheta) / kSquared
          * ( (2-kSquared)*elliptic_k(kSquared) - 2*elliptic_e(kSquared) )
          )

    # multiply by costant pre-factor
    potential *= mu0/sympy.pi

    dAdR = sympy.diff(potential, R)
    dAdZ = sympy.diff(potential, Z)
    modGradASquared = dAdR**2 + dAdZ**2

    A_func = sympy.lambdify([R,Z], potential, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    f_R_func = sympy.lambdify([R,Z], dAdR/modGradASquared, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    f_Z_func = sympy.lambdify([R,Z], dAdZ/modGradASquared, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    Bp_R_func = sympy.lambdify([R,Z], dAdZ, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    Bp_Z_func = sympy.lambdify([R,Z], -1/R * sympy.diff(R*potential, R), modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    return (A_func, f_R_func, f_Z_func, Bp_R_func, Bp_Z_func)

def createMesh(filename):
    # parse input file
    coils, Bt_axis, meshOptions = parseInput(filename)

    A_toroidal, f_R, f_Z, Bp_R, Bp_Z = magneticFunctions(coils)

    xpoint = findSaddlePoint(A_toroidal, Rmin, Rmax, Zmin, Zmax)
    A_xpoint = A_toroidal(*xpoint)
    print('X-point',xpoint,'with A_toroidal='+str(A_xpoint))

    # note legs are ordered in theta
    separatrixLegs = findSeparatrix(A_toroidal, lambda s: TORPEX_wall(s*2.*numpy.pi), xpoint, A_xpoint)

    separatrix = {'inner_lower_divertor': separatrixLegs[0],
                  'inner_upper_divertor': separatrixLegs[1],
                  'outer_upper_divertor': separatrixLegs[2],
                  'outer_lower_divertor': separatrixLegs[3]}

    fpol = Bt_axis / 1. # Major radius of TORPEX axis is 1m

    return (Mesh(meshOptions, A_toroidal, f_R, f_Z, Bp_R, Bp_Z, fpol, A_xpoint,
                 separatrix),
            {'xpoint':xpoint})

if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]
    gridname = 'torpex.grd.nc'

    mesh, info = createMesh(filename)

    mesh.geometry()

    if plotStuff:
        mesh.plotPotential(Rmin, Rmax, Zmin, Zmax)
        addWallToPlot()
        pyplot.plot(*info['xpoint'], 'rx')
        mesh.plotPoints()
        pyplot.show()

    mesh.writeGridfile(gridname)

    exit(0)
