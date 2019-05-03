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

# Could try to calculate this from range of A_toroidal in the domain or something, but
# easier just to eyeball it from contour plots.
# Only used to calculate some relative tolerances
typical_A = 1.e-7

from scipy.optimize import minimize_scalar, brentq, root
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
    theta = numpy.linspace(0., 2.*numpy.pi, npoints+1, endpoint=True)
    pyplot.plot(*TORPEX_wall(theta))

class Point2D:
    """
    A point in 2d space.
    Can be added, subtracted, multiplied by scalar
    """
    def __init__(self, R, Z):
        self.R = R
        self.Z = Z

    def __add__(self, other):
        return Point2D(self.R+other.R, self.Z+other.Z)

    def __sub__(self, other):
        return Point2D(self.R-other.R, self.Z-other.Z)

    def __mul__(self, other):
        return Point2D(self.R*other, self.Z*other)

    def __rmul__(self, other):
        return Point2D(self.R*other, self.Z*other)

    def __truediv__(self, other):
        return Point2D(self.R/other, self.Z/other)

    def __iter__(self):
        """
        Along with __next__() allows Point2D class to be treated like a tuple, e.g.
        p = Point2D(1., 0.)
        val = f(*p)
        where f is a function that takes two arguments
        """
        self.iterStep = 0
        return self

    def __next__(self):
        if self.iterStep == 0:
            self.iterStep = 1
            return self.R
        elif self.iterStep == 1:
            self.iterStep = 2
            return self.Z
        else:
            raise StopIteration

    def __repr__(self):
        """
        Allow Point2D to be printed
        """
        return 'Point2D('+str(self.R)+','+str(self.Z)+')'

class MeshContour:
    """
    Represents a contour as a collection of points.
    Includes methods for interpolation.
    Mostly behaves like a list
    """
    def __init__(self, points):
        self.points = points

    def __iter__(self):
        return self.points.__iter__()

    def append(self, point):
        self.points.append(point)

    def getRefined(self, A_target, width=.2, atol=2.e-8):
        f = lambda R,Z: A_toroidal(R, Z) - A_target

        def perpLine(p, tangent, w):
            # p - point through which to draw perpLine
            # tangent - vector tangent to original curve, result will be perpendicular to this
            # w - width on either side of p to draw the perpLine to
            modTangent = numpy.sqrt(tangent.R**2 + tangent.Z**2)
            perpIdentityVector = Point2D(tangent.Z/modTangent, -tangent.R/modTangent)
            return lambda s: p + 2.*(s-0.5)*w*perpIdentityVector

        def refinePoint(p, tangent):
            converged = False
            w = width
            sp = []
            ep = []
            while not converged:
                try:
                    pline = perpLine(p, tangent, w)
                    sp.append(pline(0.))
                    ep.append(pline(1.))
                    snew, info = brentq(lambda s: f(*pline(s)), 0., 1., xtol=atol, full_output=True)
                    converged = info.converged
                except ValueError:
                    pass
                w /= 2.
                if w < atol:
                    if numpy.abs(f(*pline(0.))) < atol*typical_A and numpy.abs(f(*pline(1.))) < atol*typical_A:
                        # f is already so close to 0 that the point does not need refining
                        snew = 0.5
                        pass
                    else:
                        raise ValueError("Could not find interval to refine point")

            return pline(snew)

        newpoints = []
        newpoints.append(refinePoint(self.points[0], self.points[1] - self.points[0]))
        for i,p in enumerate(self.points[1:-1]):
            newpoints.append(refinePoint(p, self.points[i+1] - self.points[i-1]))
        newpoints.append(refinePoint(self.points[-1], self.points[-1] - self.points[-2]))

        return MeshContour(newpoints)

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

    return numpy.vectorize(sympy.lambdify([R,Z], potential, 'mpmath'))

def plotPotential(potential, npoints=100, ncontours=40):
    pyplot.figure()
    R = numpy.linspace(Rmin, Rmax, npoints)
    Z = numpy.linspace(Zmin, Zmax, npoints)
    contours = pyplot.contour(
            R, Z, potential(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T, ncontours)
    pyplot.clabel(contours, inline=False, fmt='%1.3g')

def distance(p1, p2):
    d = p2 - p1
    return numpy.sqrt(d.R**2 + d.Z**2)

def findMinimum_1d(pos1, pos2, f, atol=1.e-14):
    coords = lambda s: pos1 + s*(pos2-pos1)
    result = minimize_scalar(lambda s: f(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
    if result.success:
        return coords(result.x)
    else:
        print('findMinimum_1d failed?')
        return coords(result.x)

def findMaximum_1d(pos1, pos2, f, atol=1.e-14):
    coords = lambda s: pos1 + s*(pos2-pos1)
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

def findSaddlePoint(f, atol=2.e-8):
    posTop, minTop = findExtremum_1d(Point2D(Rmin, Zmax), Point2D(Rmax, Zmax), f)
    posBottom, minBottom = findExtremum_1d(Point2D(Rmin, Zmin), Point2D(Rmax, Zmin), f)
    posLeft, minLeft = findExtremum_1d(Point2D(Rmin, Zmin), Point2D(Rmin, Zmax), f)
    posRight, minRight = findExtremum_1d(Point2D(Rmax, Zmin), Point2D(Rmax, Zmax), f)

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

    extremumVert = Point2D(Rmin, Zmin)
    extremumHoriz = Point2D(Rmax, Zmax)

    count = 0
    while distance(extremumVert, extremumHoriz) > atol:
        count = count+1

        extremumVert = vertSearch(posBottom, posTop, f, 0.5*atol)
        posLeft.Z = extremumVert.Z
        posRight.Z = extremumVert.Z

        extremumHoriz = horizSearch(posLeft, posRight, f, 0.5*atol)
        posBottom.R = extremumHoriz.R
        posTop.R = extremumHoriz.R

    print('findSaddlePoint took',count,'iterations to converge')

    return (extremumVert+extremumHoriz)/2.

def findRoots_1d(f, n, xmin, xmax, atol = 2.e-8, rtol = 1.e-5, maxintervals=1024):
    """
    Find n roots of a scalar function f(x) in the range xmin<=x<=xmax
    Assume they're not too close to each other - exclude a small region around each found
    root when searching for more.
    """
    smallDistance = rtol * (xmax - xmin)
    foundRoots = 0
    roots = []
    n_intervals = n
    while True:
        interval_points = numpy.linspace(xmin, xmax, n_intervals+1, endpoint=True)
        interval_f = f(interval_points)
        lucky_roots = numpy.where(interval_f == 0.)
        if len(lucky_roots[0]) > 0:
            raise NotImplementedError("Don't handle interval points that happen to land "
                    "on a root yet!")
        intervals_with_roots = numpy.where(numpy.sign(interval_f[:-1]) !=
                                           numpy.sign(interval_f[1:]))[0]
        if len(intervals_with_roots) >= n:
            break
        n_intervals *= 2
        if n_intervals > maxintervals:
            raise ValueError("Could not find", n, "roots when checking", maxintervals,
                             "intervals")

    # find roots in the intervals
    for i in intervals_with_roots:
        root, info = brentq(f, interval_points[i], interval_points[i+1], xtol=atol,
                full_output=True)
        if not info.converged:
            raise ValueError("Root finding failed in {" + str(interval_points[i]) + "," +
                    str(interval_points[i+1]) + "} with end values {" + str(interval_f[i])
                    + "," + str(interval_f[i+1]))
        roots.append(root)
        foundRoots += 1

    return roots

def findSeparatrix(xpoint, A_x, atol = 2.e-8, npoints=100):
    """
    Follow 4 legs away from the x-point, starting with a rough guess and then refining to
    the separatrix value of A_toroidal.
    """
    boundaryThetas = findRoots_1d(lambda theta: A_toroidal(*TORPEX_wall(theta)) - A_x, 4,
            0., 2.*numpy.pi)
    boundaryPoints = tuple(TORPEX_wall(theta) for theta in boundaryThetas)

    legs = []
    s = numpy.linspace(10.*atol, 1., npoints, endpoint=True)
    for point in boundaryPoints:
        legR = xpoint.R + s*(point.R - xpoint.R)
        legZ = xpoint.Z + s*(point.Z - xpoint.Z)
        leg = MeshContour([Point2D(R,Z) for R,Z in zip(legR, legZ)])
        leg = leg.getRefined(A_x, atol=atol)
        legs.append(leg)

    return legs

if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]

    # parse input file
    coils = parseInput(filename)

    A_toroidal = potentialFunction(coils)

    #R = numpy.linspace(Rmin, Rmax, 50)
    #Z = numpy.linspace(Zmin, Zmax, 50)

    xpoint = findSaddlePoint(A_toroidal)
    A_xpoint = A_toroidal(*xpoint)
    print('X-point',xpoint,'with A_toroidal='+str(A_xpoint))

    # note legs are ordered in theta
    separatrixLegs = findSeparatrix(xpoint, A_xpoint)

    if plotStuff:
        #plotPotential(A_toroidal)
        plotPotential(lambda R,Z: A_toroidal(R,Z)-A_xpoint)
        pyplot.axes().set_aspect('equal')
        addWallToPlot()
        pyplot.plot(*xpoint, 'rx')
        for l in separatrixLegs:
            pyplot.plot([x.R for x in l], [x.Z for x in l])
        pyplot.show()

    exit(0)
