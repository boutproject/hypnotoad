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
Rmin = 0.9
Rmax = 1.1
Zmin = -.1
Zmax = .1

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

    return numpy.vectorize(sympy.lambdify([R,Z], potential, 'sympy'))

def plotPotential(potential, npoints=100, ncontours=40):
    pyplot.figure()
    R = numpy.linspace(Rmin, Rmax, npoints)
    Z = numpy.linspace(Zmin, Zmax, npoints)
    contours = pyplot.contour(
            R, Z, potential(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T, ncontours)
    pyplot.clabel(contours, inline=False, fmt='%1.3g')

def findMinimum_1d(pos1, pos2, f, rtol=1.e-5, atol=1.e-14):
    # Golden-section search: https://en.wikipedia.org/wiki/Golden-section_search

    # Golden ratio
    oneOverPhi = 2./(1. + numpy.sqrt(5.))

    # coordinates of line between pos1 and pos2, parameterized by 0<=s<=1
    coords = lambda s: (pos1[0] + s*(pos2[0]-pos1[0]), (pos1[1] + s*(pos2[1]-pos1[1])))
    def realdist(s1, s2):
        R1,Z1 = coords(s1)
        R2,Z2 = coords(s2)
        return numpy.sqrt((R2-R1)**2 + (Z2-Z1)**2)

    # extremum is in an interval rtol times smaller than the distance between pos1 and
    # pos2
    sLeft = 0.
    sRight = 1.
    fLeft = f(*coords(sLeft))
    fRight = f(*coords(sRight))
    newsLeft = sRight - (sRight - sLeft)*oneOverPhi
    newsRight = sLeft + (sRight - sLeft)*oneOverPhi
    newfLeft = f(*coords(newsLeft))
    newfRight = f(*coords(newsRight))
    while sRight - sLeft > rtol and realdist(sLeft, sRight) > atol:
        print('findmin',newsLeft,newfLeft,newsRight,newfRight)
        if newfLeft < newfRight:
            sRight = newsRight
            fRight = newfRight
            newsRight = newsLeft
            newfRight = newfLeft
            newsLeft = sRight - (sRight - sLeft)*oneOverPhi
            newfLeft = f(*coords(newsLeft))
        else:
            sLeft = newsLeft
            fLeft = newsLeft
            newsLeft = newsRight
            newfLeft = newfRight
            newsRight = sLeft + (sRight - sLeft)*oneOverPhi
            newfRight = f(*coords(newsRight))

    return coords((sRight+sLeft)/2.)

def findMaximum_1d(pos1, pos2, f, rtol=1.e-5, atol=1.e-14):
    # Golden-section search: https://en.wikipedia.org/wiki/Golden-section_search

    # Golden ratio
    oneOverPhi = 2./(1. + numpy.sqrt(5.))

    # coordinates of line between pos1 and pos2, parameterized by 0<=s<=1
    coords = lambda s: (pos1[0] + s*(pos2[0]-pos1[0]), (pos1[1] + s*(pos2[1]-pos1[1])))
    def realdist(s1, s2):
        R1,Z1 = coords(s1)
        R2,Z2 = coords(s2)
        return numpy.sqrt((R2-R1)**2 + (Z2-Z1)**2)

    # extremum is in an interval rtol times smaller than the distance between pos1 and
    # pos2
    sLeft = 0.
    sRight = 1.
    fLeft = f(*coords(sLeft))
    fRight = f(*coords(sRight))
    newsLeft = sRight - (sRight - sLeft)*oneOverPhi
    newsRight = sLeft + (sRight - sLeft)*oneOverPhi
    newfLeft = f(*coords(newsLeft))
    newfRight = f(*coords(newsRight))
    while sRight - sLeft > rtol and realdist(sLeft, sRight) > atol:
        print('findmax',newsLeft,newfLeft,newsRight,newfRight)
        if newfLeft > newfRight:
            sRight = newsRight
            fRight = newfRight
            newsRight = newsLeft
            newfRight = newfLeft
            newsLeft = sRight - (sRight - sLeft)*oneOverPhi
            newfLeft = f(*coords(newsLeft))
        else:
            sLeft = newsLeft
            fLeft = newsLeft
            newsLeft = newsRight
            newfLeft = newfRight
            newsRight = sLeft + (sRight - sLeft)*oneOverPhi
            newfRight = f(*coords(newsRight))

    return coords((sRight+sLeft)/2.)

if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]

    # parse input file
    coils = parseInput(filename)

    Ator = potentialFunction(coils)

    # test extremum?
    pos1 = (.9,-0.1)
    pos2 = (.95,.1)
    R = numpy.linspace(pos1[0],pos2[0],50)
    Z = numpy.linspace(pos1[1],pos2[1],50)
    minpoint = findMinimum_1d(pos1, pos2, Ator)
    maxpoint = findMaximum_1d(pos1, pos2, Ator)
    pyplot.plot(R,Ator(R,Z))
    pyplot.xlabel('R')
    pyplot.axvline(minpoint[0], c='r')
    pyplot.axvline(maxpoint[0], c='k')
    pyplot.figure()
    pyplot.plot(Z,Ator(R,Z))
    pyplot.xlabel('Z')
    pyplot.axvline(minpoint[1], c='r')
    pyplot.axvline(maxpoint[1], c='k')
    pyplot.show()
    exit(0)

    if plotStuff:
        plotPotential(Ator)
        pyplot.show()

    exit(0)
