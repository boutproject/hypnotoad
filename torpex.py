#!/usr/bin/env python3
"""
Create a BOUT++ grid for TORPEX from an input file giving coil currents and positions

Input file should contain coil parameters, for each coil:
    R: major radius in metres
    Z: major radius in metres
    I: anti-clockwise current in Amps
"""

import numpy

plotStuff = True

# Torpex parameters
Rmin = 0.8
Rmax = 1.2
Zmin = -.2
Zmax = .2

from scipy.optimize import minimize_scalar
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
    return (Rcentre + awall*numpy.cos(theta), Zcentre + awall*numpy.sin(theta))

def addWallToPlot(npoints=100):
    theta = numpy.linspace(0., 2.*numpy.pi, npoints+1, endpoint=True)
    pyplot.plot(*TORPEX_wall(theta))

def parseInput(filename):
    import yaml
    from collections import namedtuple

    with open(filename, 'r') as inputfile:
        inputs = yaml.safe_load(inputfile)
    print(inputs['Coils'])
    
    Coil = namedtuple('Coil', 'R, Z, I')
    return [Coil(**c) for c in inputs['Coils']]

def potentialFunction(coils):
    """
    Calculate toroidal (anticlockwise) component of magnetic vector potential due to coils
    See for example http://physics.usask.ca/~hirose/p812/notes/Ch3.pdf
    """
    import sympy
    from sympy.functions.special.elliptic_integrals import elliptic_k, elliptic_e

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
        potential += (
          coil.I*coil.R / sympy.sqrt(r**2 + coil.R**2 + 2*coil.R*r*sinTheta) / kSquared
          * ( (2-kSquared)*elliptic_k(kSquared) - 2*elliptic_e(kSquared) )
          )

    # multiply by costant pre-factor
    potential *= mu0/sympy.pi
    print(potential)

    return numpy.vectorize(sympy.lambdify([R,Z], potential, 'mpmath'))

def plotPotential(potential, npoints=100, ncontours=40):
    pyplot.figure()
    R = numpy.linspace(Rmin, Rmax, npoints)
    Z = numpy.linspace(Zmin, Zmax, npoints)
    contours = pyplot.contour(
            R, Z, potential(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T, ncontours)
    pyplot.clabel(contours, inline=False, fmt='%1.3g')

def distance(p1, p2):
    return numpy.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

def findMinimum_1d(pos1, pos2, f, atol=1.e-14):
    coords = lambda s: (pos1[0] + s*(pos2[0]-pos1[0]), (pos1[1] + s*(pos2[1]-pos1[1])))
    result = minimize_scalar(lambda s: f(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
    if result.success:
        return coords(result.x)
    else:
        print('findMinimum_1d failed?')
        return coords(result.x)

def findMaximum_1d(pos1, pos2, f, atol=1.e-14):
    coords = lambda s: (pos1[0] + s*(pos2[0]-pos1[0]), (pos1[1] + s*(pos2[1]-pos1[1])))
    # minimize -f to find maximum
    result = minimize_scalar(lambda s: -f(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
    if result.success:
        return coords(result.x)
    else:
        print('findMaximum_1d failed?')
        return coords(result.x)

def findExtremum_1d(pos1, pos2, f, rtol=1.e-5, atol=1.e-14):
    smallDistance = 10.*rtol*distance(pos1, pos2)

    minpos = findMinimum_1d(pos1, pos2, f, atol)
    if distance(pos1,minpos) > smallDistance and distance(pos2,minpos) > smallDistance:
        # minimum is not at either end of the interval
        return minpos, True

    maxpos = findMaximum_1d(pos1, pos2, f, atol)
    if distance(pos1,maxpos) > smallDistance and distance(pos2,maxpos) > smallDistance:
        return maxpos, False

    raise ValueError("Neither minimum nor maximum found in interval")

def findSaddlePoint(f, atol=5.e-8):
    posTop, minTop = findExtremum_1d((Rmin, Zmax), (Rmax, Zmax), f)
    posBottom, minBottom = findExtremum_1d((Rmin, Zmin), (Rmax, Zmin), f)
    posLeft, minLeft = findExtremum_1d((Rmin, Zmin), (Rmin, Zmax), f)
    posRight, minRight = findExtremum_1d((Rmax, Zmin), (Rmax, Zmax), f)

    assert minTop == minBottom
    assert minLeft == minRight
    assert minTop != minLeft

    if minTop:
        vertSearch = findMaximum_1d
    else:
        vertSearch = findMinimum_1d

    if minLeft:
        horizSearch = findMaximum_1d
    else:
        horizSearch = findinximum_1d

    extremumVert = (Rmin, Zmin)
    extremumHoriz = (Rmax, Zmax)

    count = 0
    while distance(extremumVert, extremumHoriz) > atol:
        print(count, extremumVert, extremumHoriz, distance(extremumVert, extremumHoriz))
        count = count+1

        extremumVert = vertSearch(posBottom, posTop, f, 0.5*atol)
        posLeft = (posLeft[0], extremumVert[1])
        posRight = (posRight[0], extremumVert[1])

        extremumHoriz = horizSearch(posLeft, posRight, f, 0.5*atol)
        posBottom = (extremumHoriz[0], posBottom[1])
        posTop = (extremumHoriz[0], posTop[1])

    return ((extremumVert[0]+extremumHoriz[0])/2., (extremumVert[1]+extremumHoriz[1])/2.)

if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]

    # parse input file
    coils = parseInput(filename)

    A_toroidal = potentialFunction(coils)

    R = numpy.linspace(Rmin, Rmax, 50)
    Z = numpy.linspace(Zmin, Zmax, 50)
    AtorArray = A_toroidal(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T

    xpoint = findSaddlePoint(A_toroidal)

    if plotStuff:
        #plotPotential(A_toroidal)
        plotPotential(lambda R,Z: A_toroidal(R,Z)-A_xpoint)
        pyplot.axes().set_aspect('equal')
        addWallToPlot()
        pyplot.plot(*xpoint, 'rx')
        pyplot.show()

    exit(0)
