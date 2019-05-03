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
Rmin = 0.9
Rmax = 1.1
Zmin = -.1
Zmax = .1

# Golden ratio
oneOverPhi = 2./(1. + numpy.sqrt(5.))

from scipy.optimize import minimize_scalar
if plotStuff:
    from matplotlib import pyplot

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
    minpos = findMinimum_1d(pos1, pos2, f, atol)

    smallDistance = 10.*rtol*distance(pos1, pos2)
    if distance(pos1,minpos) > smallDistance and distance(pos2,minpos) > smallDistance:
        # minimum is not at either end of the interval
        return minpos, True

    return findMaximum_1d(pos1, pos2, f, atol), False

def findSaddlePoint(f, atol=2.e-14):
    posTop, minTop = findExtremum_1d((Rmin, Zmax), (Rmax, Zmax), f)
    posBottom, minBottom = findExtremum_1d((Rmin, Zmin), (Rmax, Zmin), f)
    posLeft, minLeft = findExtremum_1d((Rmin, Zmin), (Rmin, Zmax), f)
    posRight, minRight = findExtremum_1d((Rmax, Zmin), (Rmax, Zmax), f)

    assert minTop == minBottom
    assert minLeft == minRight
    assert minTop != minLeft

    if minTop:
        vertSearch = findMinimum_1d
    else:
        vertSearch = findMaximum_1d

    if minLeft:
        horizSearch = findMinimum_1d
    else:
        horizSearch = findMaximum_1d

    def isAbove(p, pLeft, pRight):
        # is p above the line connecting pLeft and pRight?

        # Z at R=p[0]
        Z = pLeft[1] + (p[0]-pLeft[0])/(pRight[0] - pLeft[0]) * (pRight[1] - pLeft[1])

        return p[1] > Z

    def isToRight(p, pBottom, pTop):
        # is p to the right of the line connecting pBottom and pTop?

        # R at Z=p[1]
        R = pBottom[0] + (p[1]-pBottom[1])/(pTop[1] - pBottom[1]) * (pTop[0] - pBottom[0])

        return p[0] > R

    def checkLims(midpoint, p1, p2, atol):
        # if the extremum is at p1 or p2, push that point away
        if distance(midpoint, p1) < atol:
            print('found end 1')
            p1 = (2.*p1[0]-p2[0], 2.*p1[1]-p2[1])
        if distance(midpoint, p2) < atol:
            print('found end 2')
            p2 = (2.*p2[0]-p1[0], 2.*p2[1]-p1[1])
        return p1, p2

    def updateFurther(midpoint, p1, p2):
        # advance the further away of p1, p2 towards midpoint
        if distance(midpoint, p1) > distance(midpoint, p2):
            p1 = ( midpoint[0] - (midpoint[0] - p1[0])*oneOverPhi,
                   midpoint[1] - (midpoint[1] - p1[1])*oneOverPhi )
        else:
            p2 = ( midpoint[0] - (midpoint[0] - p2[0])*oneOverPhi,
                   midpoint[1] - (midpoint[1] - p2[1])*oneOverPhi )
        return p1, p2

    count = 0
    while distance(posBottom, posTop) > atol or distance(posLeft, posRight) > atol:
        print(count, distance(posBottom, posTop), distance(posLeft, posRight))
        count = count+1

        midpoint = vertSearch(posBottom, posTop, f, 0.5*atol)
        posBottom, posTop = checkLims(midpoint, posBottom, posTop, atol)
        posLeft, posRight = updateFurther(midpoint, posLeft, posRight)

        midpoint = horizSearch(posLeft, posRight, f, 0.5*atol)
        posLeft, posRight = checkLims(midpoint, posLeft, posRight, atol)
        posBottom, posTop = updateFurther(midpoint, posBottom, posTop)

    return ((posLeft[0]+posRight[0])/2., (posBottom[1]+posTop[1])/2.)

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
        plotPotential(A_toroidal)
        pyplot.plot(*xpoint, 'rx')
        pyplot.show()

    exit(0)
