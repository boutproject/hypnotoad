"""
Functions for analysing an equilibrium for which an interpolating function is given for
the potential.
"""

import numpy
from scipy.optimize import minimize_scalar, brentq, root
from scipy.interpolate import interp1d
from collections import OrderedDict
from copy import deepcopy
import warnings

class SolutionError(Exception):
    """
    Solution was not found
    """
    pass

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

def calc_distance(p1, p2):
    d = p2 - p1
    return numpy.sqrt(d.R**2 + d.Z**2)

class PsiContour:
    """
    Represents a contour as a collection of points.
    Includes methods for interpolation.
    Mostly behaves like a list
    """

    # Number of points for high-resolution intermediate contours used in re-gridding
    Nfine = 1000

    def __init__(self, points, psi, psival):
        self.points = points

        self.startInd = 0
        self.endInd = len(points) - 1

        self.recalculateDistance()

        # Function that evaluates the vector potential at R,Z
        self.psi = psi

        # Value of vector potential on this contour
        self.psival = psival

    def __iter__(self):
        return self.points.__iter__()

    def __str__(self):
        return self.points.__str__()

    def __getitem__(self, key):
        return self.points.__getitem__(key)

    def __len__(self):
        return self.points.__len__()

    def append(self, point):
        self.points.append(point)
        self.distance.append(
                self.distance[-1] + calc_distance(self.points[-2], self.points[-1]))

    def recalculateDistance(self):
        self.distance = [0.]
        for i in range(1, len(self.points)):
            self.distance.append(
                    self.distance[-1] + calc_distance(self.points[i-1], self.points[i]))

    def reverse(self):
        self.points.reverse()
        old_start = self.startInd
        self.startInd = len(self) - 1 - self.endInd
        self.endInd = len(self) - 1 - old_start
        self.distance.reverse()
        self.distance = [self.distance[0] - d for d in self.distance]

    def refine(self, *args, **kwargs):
        new = self.getRefined(*args, **kwargs)
        self.points = new.points
        self.distance = new.distance

    def getRefined(self, width=1.e-5, atol=2.e-8):
        f = lambda R,Z: self.psi(R, Z) - self.psival

        def perpLine(p, tangent, w):
            # p - point through which to draw perpLine
            # tangent - vector tangent to original curve, result will be perpendicular to this
            # w - width on either side of p to draw the perpLine to
            modTangent = numpy.sqrt(tangent.R**2 + tangent.Z**2)
            perpIdentityVector = Point2D(tangent.Z/modTangent, -tangent.R/modTangent)
            return lambda s: p + 2.*(s-0.5)*w*perpIdentityVector

        def refinePoint(p, tangent):
            if numpy.abs(f(*p)) < atol*numpy.abs(self.psival):
                # don't need to refine
                return p
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
                    from matplotlib import pyplot
                    pline0 = perpLine(p, tangent, width)
                    Rbox = numpy.linspace(p.R-.05,p.R+.05,100)[numpy.newaxis,:]
                    Zbox = numpy.linspace(p.Z-.05,p.Z+.05,100)[:,numpy.newaxis]
                    svals = numpy.linspace(0., 1., 40)
                    pyplot.figure()
                    pyplot.contour(Rbox+0.*Zbox,Zbox+0.*Rbox,self.psi(Rbox,Zbox))
                    pyplot.plot([pline0(s).R for s in svals], [pline0(s).Z for s in svals], 'x')
                    pyplot.figure()
                    pyplot.plot([f(*pline0(s)) for s in svals])
                    pyplot.show()
                    raise SolutionError("Could not find interval to refine point at "+str(p))

            return pline(snew)

        newpoints = []
        newpoints.append(refinePoint(self.points[0], self.points[1] - self.points[0]))
        for i,p in enumerate(self.points[1:-1]):
            newpoints.append(refinePoint(p, self.points[i+1] - self.points[i-1]))
        newpoints.append(refinePoint(self.points[-1], self.points[-1] - self.points[-2]))

        return PsiContour(newpoints, self.psi, self.psival)

    def interpFunction(self):
        distance = numpy.array(numpy.float64(self.distance)) - self.distance[self.startInd]
        R = numpy.array(numpy.float64([p.R for p in self.points]))
        Z = numpy.array(numpy.float64([p.Z for p in self.points]))
        interpR = interp1d(distance, R, kind='cubic',
                           assume_sorted=True, fill_value='extrapolate')
        interpZ = interp1d(distance, Z, kind='cubic',
                           assume_sorted=True, fill_value='extrapolate')
        return lambda s: Point2D(interpR(s), interpZ(s))

    def getRegridded(self, npoints, *, width=1.e-5, atol=2.e-8, sfunc=None,
            extend_lower=0, extend_upper=0):
        """
        Interpolate onto set of npoints points, then refine positions.
        By default points are uniformly spaced, this can be changed by passing 'sfunc'
        which replaces the uniform interval 's' with 's=sfunc(s)'.
        'extend_lower' and 'extend_upper' extend the contour past its existing ends by a
        number of points.
        Returns a new PsiContour.

        Note: '*,' in the arguments list forces the following arguments to be passed as
        keyword, not positional, arguments
        """
        # To make the new points accurate, first regrid onto a high-resolution contour,
        # then interpolate.
        # Extend further than will be needed in the final contour, because extrapolation
        # past the end of the fine contour is very bad.
        extend_lower_fine = 2*(extend_lower * self.Nfine) // npoints
        extend_upper_fine = 2*(extend_upper * self.Nfine) // npoints

        indices_fine = numpy.linspace(-extend_lower_fine,
                (self.Nfine - 1 + extend_upper_fine),
                self.Nfine + extend_lower_fine + extend_upper_fine)

        if sfunc is not None:
            sfine = sfunc(indices_fine*(npoints - 1)/(self.Nfine - 1))
        else:
            sfine = (self.distance[self.endInd] - self.distance[self.startInd]) / (self.Nfine - 1) * indices_fine

        interp_self = self.interpFunction()

        fine_contour = PsiContour([interp_self(x) for x in sfine], self.psi, self.psival)
        fine_contour.startInd = extend_lower_fine
        fine_contour.endInd = len(fine_contour) - 1 - extend_upper_fine
        fine_contour.refine(width, atol)

        # fine_contour does not have exactly the same total distance as the original,
        # because it integrates on a finer grid. But sfunc was defined with the total
        # distance of the original, so scale the distance of fine_contour to make the
        # totals the same
        scale_factor = (self.distance[self.endInd] - self.distance[self.startInd]) / (fine_contour.distance[fine_contour.endInd] - fine_contour.distance[fine_contour.startInd])
        fine_contour.distance = [scale_factor*d for d in fine_contour.distance]

        indices = numpy.linspace(-extend_lower, (npoints - 1 + extend_upper),
                npoints + extend_lower + extend_upper)
        if sfunc is not None:
            s = sfunc(indices)
        else:
            s = (self.distance[self.endInd] - self.distance[self.startInd]) / (npoints - 1) * indices

        interp_fine = fine_contour.interpFunction()

        new_contour = PsiContour([interp_fine(x) for x in s], self.psi, self.psival)
        new_contour.startInd = extend_lower
        new_contour.endInd = len(new_contour) - 1 - extend_upper
        # new_contour was interpolated from a high-resolution contour, so should not need
        # a large width for refinement - use 1.e-7 instead of 'width'
        return new_contour.getRefined(1.e-7, atol)

    def plot(self, *args, **kwargs):
        from matplotlib import pyplot
        pyplot.plot([x.R for x in self], [x.Z for x in self], *args, **kwargs)

class EquilibriumRegion(PsiContour):
    """
    Specialization of PsiContour for representing an equilibrium segment, which is a
    poloidal segment based around a contour (normally a segment of a separatrix). Includes
    members giving the connections to other regions and to list the X-points at the
    boundaries where the contour starts or ends.
    """
    def __init__(self, name, nSegments, options, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name
        self.nSegments = nSegments
        self.options = options
        self.nx = self.options['nx']
        self.ny_noguards = self.options['ny']
        self.y_boundary_guards = self.options['y_boundary_guards']

        self.xPointsAtStart = []
        self.xPointsAtEnd = []
        self.connections = []
        self.psi_vals = []
        self.sfunc = None

        # xPointsAtStart and xPointsAtEnd should have an entry at the lower and upper side
        # of each segment, so they both have length=nSegments+1
        self.xPointsAtStart.append(None)
        self.xPointsAtEnd.append(None)
        for i in range(nSegments):
            c = {'inner':None, 'outer':None, 'lower':None, 'upper':None}
            if i > 0:
                c['inner'] = (self.name, i - 1)
            if i < nSegments - 1:
                c['outer'] = (self.name, i + 1)
            self.connections.append(c)
            self.xPointsAtStart.append(None)
            self.xPointsAtEnd.append(None)

    def fromPsiContour(self, contour):
        result = EquilibriumRegion(self.name, self.nSegments, self.options,
                contour.points, contour.psi, contour.psival)
        result.xPointsAtStart = deepcopy(self.xPointsAtStart)
        result.xPointsAtEnd = deepcopy(self.xPointsAtEnd)
        result.connections = deepcopy(self.connections)
        result.psi_vals = deepcopy(self.psi_vals)
        return result

    def ny(self, radialIndex):
        # Get ny for a segment of this EquilibriumRegion, including any y-boundary guard
        # cells
        result = self.ny_noguards
        if self.connections[radialIndex]['lower'] is None:
            result += self.y_boundary_guards
        if self.connections[radialIndex]['upper'] is None:
            result += self.y_boundary_guards
        return result

    def getRefined(self, *args, **kwargs):
        return self.fromPsiContour(super().getRefined(*args, **kwargs))

    def getRegridded(self, *, radialIndex, **kwargs):
        for wrong_argument in ['npoints', 'extend_lower', 'extend_upper', 'sfunc']:
            # these are valid arguments to PsiContour.getRegridded, but not to
            # EquilibriumRegion.getRegridded. EquilibriumRegion.getRegridded knows its own
            # ny and connections, so must use these
            if wrong_argument in kwargs:
                raise ValueError("'"+wrong_argument+"' should not be given as an "
                        "argument to EquilibriumRegion.getRegridded")
        if self.connections[radialIndex]['lower'] is None:
            extend_lower = 2*self.y_boundary_guards
        else:
            extend_lower = 0
        if self.connections[radialIndex]['upper'] is None:
            extend_upper = 2*self.y_boundary_guards
        else:
            extend_upper = 0
        return self.fromPsiContour(super().getRegridded(2*self.ny_noguards + 1,
            extend_lower=extend_lower, extend_upper=extend_upper, sfunc=self.sfunc, **kwargs))

class Equilibrium:
    """
    Base class to provide an interface to an interpolating function for the flux function psi that defines
    the magnetic equilibrium, along with some useful methods.

    psi is the magnetic flux function.

    f_R and f_Z are the components of a vector Grad(psi)/|Grad(psi)|**2. This vector
    points along a path perpendicular to psi-contours, and its value is ds/dpsi where s is
    the coordinate along the path, so we can follow the path by integrating this vector:
    R(psi) = \\int_0^\\psi f_R
    and
    Z(psi) = \\int_0^\\psi f_Z

    Derived classes must provide:
      - self.psi: function which takes two arguments, {R,Z}, and returns the value of psi
        at that position.
      - self.f_R: function which takes two arguments, {R,Z}, and returns the R
        component of the vector Grad(psi)/|Grad(psi)|**2.
      - self.f_Z: function which takes two arguments, {R,Z}, and returns the Z
        component of the vector Grad(psi)/|Grad(psi)|**2.
      - self.Bp_R: function which takes two arguments, {R,Z}, and returns the R
        component of the poloidal magnetic field.
      - self.Bp_Z: function which takes two arguments, {R,Z}, and returns the Z
        component of the poloidal magnetic field.
      - self.x_points: list of Point2D objects giving the position of the X-points ordered
        from primary X-point (nearest the core) outward
      - self.psi_sep: values of psi on the separatrices ordered the same as self.x_points
      - self.fpol: poloidal current function, takes one argument, psi, and returns fpol
        (function such that B_toroidal = fpol/R)
      - self.Rmin, self.Rmax, self.Zmin, self.Zmax: positions of the corners of a bounding
        box for the gridding
      - self.regions: OrderedDict of EquilibriumRegion objects that specify this equilibrium
    """
    def makeConnection(self, lowerRegion, lowerSegment, upperRegion, upperSegment):
        """
        Make a connection between the upper edge of a certain segment of lowerRegion and
        the lower edge of a certain segment of upperRegion.
        """
        # Needs to be OrderedDict so that Mesh can iterate through it in consistent order
        assert type(self.regions) == OrderedDict, 'self.regions should be OrderedDict'

        lRegion = self.regions[lowerRegion]
        uRegion = self.regions[upperRegion]

        assert lRegion.connections[lowerSegment]['upper'] is None, 'lRegion.connections[\'upper\'] should not have been set already'
        assert uRegion.connections[upperSegment]['lower'] is None, 'uRegion.connections[\'lower\'] should not have been set already'

        # Check nx of both segments is the same - otherwise the connection must be between
        # some wrong regions
        assert lRegion.nx[lowerSegment] == uRegion.nx[upperSegment], 'nx should match across connection'

        lRegion.connections[lowerSegment]['upper'] = (upperRegion, upperSegment)
        uRegion.connections[upperSegment]['lower'] = (lowerRegion, lowerSegment)

    def magneticFunctionsFromGrid(self, R, Z, psiRZ):
        from dct_interpolation import DCT_2D

        self._dct = DCT_2D(R, Z, psiRZ)

        self.psi = lambda R, Z: self._dct(R, Z)
        modGradpsiSquared = lambda R, Z: self._dct.ddR(R, Z)**2 + self._dct.ddZ(R, Z)**2
        self.f_R = lambda R, Z: self._dct.ddR(R, Z) / modGradpsiSquared(R, Z)
        self.f_Z = lambda R, Z: self._dct.ddZ(R, Z) / modGradpsiSquared(R, Z)
        self.Bp_R = lambda R, Z: self._dct.ddZ(R, Z) / R
        self.Bp_Z = lambda R, Z: -self._dct.ddR(R, Z) / R

    def findMinimum_1d(self, pos1, pos2, atol=1.e-14):
        coords = lambda s: pos1 + s*(pos2-pos1)
        result = minimize_scalar(lambda s: self.psi(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
        if result.success:
            return coords(result.x)
        else:
            raise SolutionError('findMinimum_1d failed')

    def findMaximum_1d(self, pos1, pos2, atol=1.e-14):
        coords = lambda s: pos1 + s*(pos2-pos1)
        # minimize -f to find maximum
        result = minimize_scalar(lambda s: -self.psi(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
        if result.success:
            return coords(result.x)
        else:
            raise SolutionError('findMaximum_1d failed')

    def findExtremum_1d(self, pos1, pos2, rtol=1.e-5, atol=1.e-14):
        smallDistance = 10.*rtol*calc_distance(pos1, pos2)

        minpos = self.findMinimum_1d(pos1, pos2, atol)
        if calc_distance(pos1,minpos) > smallDistance and calc_distance(pos2,minpos) > smallDistance:
            # minimum is not at either end of the interval
            return minpos, True

        maxpos = self.findMaximum_1d(pos1, pos2, atol)
        if calc_distance(pos1,maxpos) > smallDistance and calc_distance(pos2,maxpos) > smallDistance:
            return maxpos, False

        raise SolutionError("Neither minimum nor maximum found in interval")

    def findSaddlePoint(self, Rmin, Rmax, Zmin, Zmax, atol=2.e-8):
        """
        Find a saddle point in the function self.psi atol is the tolerance on the position
        of the saddle point {Rmin,Rmax}, {Zmin,Zmax} are the bounding values of the box to
        search for a saddle point in.
        """

        posTop, minTop = self.findExtremum_1d(Point2D(Rmin, Zmax), Point2D(Rmax, Zmax))
        posBottom, minBottom = self.findExtremum_1d(Point2D(Rmin, Zmin), Point2D(Rmax, Zmin))
        posLeft, minLeft = self.findExtremum_1d(Point2D(Rmin, Zmin), Point2D(Rmin, Zmax))
        posRight, minRight = self.findExtremum_1d(Point2D(Rmax, Zmin), Point2D(Rmax, Zmax))

        assert minTop == minBottom, 'if minumum is found at top, should also be found at bottom'
        assert minLeft == minRight, 'if minumum is found at left, should also be found at right'
        assert minTop != minLeft, 'if minimum is found at top, maximum should be found at left'

        if minTop:
            vertSearch = self.findMaximum_1d
        else:
            vertSearch = self.findMinimum_1d

        if minLeft:
            horizSearch = self.findMaximum_1d
        else:
            horizSearch = self.findMinimum_1d

        extremumVert = Point2D(Rmin, Zmin)
        extremumHoriz = Point2D(Rmax, Zmax)

        count = 0
        while calc_distance(extremumVert, extremumHoriz) > atol:
            count = count+1

            extremumVert = vertSearch(posBottom, posTop, 0.5*atol)
            posLeft.Z = extremumVert.Z
            posRight.Z = extremumVert.Z

            extremumHoriz = horizSearch(posLeft, posRight, 0.5*atol)
            posBottom.R = extremumHoriz.R
            posTop.R = extremumHoriz.R

        print('findSaddlePoint took',count,'iterations to converge')

        return (extremumVert+extremumHoriz)/2.

    def findRoots_1d(self, f, n, xmin, xmax, atol = 2.e-8, rtol = 1.e-5, maxintervals=1024):
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
            interval_points = numpy.linspace(xmin, xmax, n_intervals+1)
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
                raise SolutionError("Could not find", n, "roots when checking", maxintervals,
                                 "intervals")

        # find roots in the intervals
        for i in intervals_with_roots:
            root, info = brentq(f, interval_points[i], interval_points[i+1], xtol=atol,
                    full_output=True)
            if not info.converged:
                raise SolutionError("Root finding failed in {" + str(interval_points[i]) + "," +
                        str(interval_points[i+1]) + "} with end values {" + str(interval_f[i])
                        + "," + str(interval_f[i+1]))
            roots.append(root)
            foundRoots += 1

        if foundRoots > n:
            warnings.warn('Warning: found',foundRoots,'roots but expected only',n)

        return roots

    def make1dGrid(self, n, spacingFunc):
        """
        Make a 1d grid:
        - Start by generating grid of cell-face values, with values from spacingFunc.
        - Place cell-centre values half-way between cell-faces.

        spacingFunc should take an index between 0 and n, and returns the desired
        coordinate value.
        """
        face_vals = [spacingFunc(i) for i in range(n + 1)]

        result = numpy.zeros(2*n + 1)
        result[::2] = face_vals
        result[1::2] = 0.5*(result[:-1:2] + result[2::2])

        return result

    def getPolynomialGridFunc(self, n, lower, upper, *, grad_lower=None, grad_upper=None):
        """
        A polynomial function with value 'lower' at 0 and 'upper' at n, used to
        non-uniformly place grid point values in index space. Optionally matches the
        gradient grad_lower at the lower end and grad_upper at the upper end. Linear if
        neither gradient given, quadritic if one given, cubic if both given.
        """
        if grad_lower is None and grad_upper is None:
            return lambda i: lower + (upper - lower)*i/n
        elif grad_lower is None:
            # psi(i) = a*i^2 + b*i + c
            # psi(0) = lower = c
            # psi(n) = upper = a*n**2 + b*n + c
            # dpsidi(n) = grad_upper = 2*a*n + b
            # n*a - c/n = grad_upper - upper/n
            c = lower
            a = (grad_upper - upper/n + c/n) / n
            b = grad_upper - 2.*a*n
            return lambda i: a*i**2 + b*i + c
        elif grad_upper is None:
            # psi(i) = a*i^2 + b*i + c
            # psi(0) = lower = c
            # psi(n) = upper = a*n**2 + b*n + c
            # dpsidi(0) = grad_lower = b
            c = lower
            b = grad_lower
            a = (upper - b*n - c) / n**2
            return lambda i: a*i**2 + b*i + c
        else:
            # psi(i) = a*i^3 + b*i^2 + c*i + d
            # psi(0) = lower = d
            # dpsidi(0) = grad_lower = c
            # dpsidi(n) = grad_upper = 3*a*n^2 + 2*b*n + c
            # psi(n) = upper = a*n^3 + b*n^2 + c*n + d
            # grad_upper - 2*upper/n = a*n^2 - c - d/n
            d = lower
            c = grad_lower
            a = (grad_upper - 2*upper/n + c + d/n) / n**2
            b = (grad_upper - 3*a*n**2 - c) / (2*n)
            return lambda i: a*i**3 + b*i**2 + c*i + d

    def getSqrtPoloidalDistanceFunc(self, length, N, N_norm, *, d_lower=None, d_sqrt_lower=None,
            d_upper=None, d_sqrt_upper=None):
        """
        Return a function s(i) giving poloidal distance as a function of index-number.
        Construct s(i)=sN(iN) as a function of the normalized iN = i/N_norm so that it has the
        same form when resolution is changed. The total Ny in the grid might be a good
        choice for N_norm.
        sN(0) = 0
        sN(N/N_norm) = L
        ds/diN(0) ~ d_sqrt_lower/sqrt(iN)+d_lower at iN=0 (if d_lower not None, else no
                                                           sqrt(iN) term)
        ds/diN(N/N_norm) ~ d_sqrt_upper/sqrt(N/N_norm-iN)+d_upper at iN=N_norm (if d_upper
                                                                              is not None)

        By default d_sqrt_lower=d_lower and d_sqrt_upper=d_upper, unless both are
        specified explicitly
        """
        if d_lower is None and d_upper is None:
            assert d_sqrt_lower is None, 'cannot set d_sqrt_lower unless d_lower is set'
            assert d_sqrt_upper is None, 'cannot set d_sqrt_upper unless d_upper is set'
            # always monotonic
            return lambda i: i*length/N
        elif d_lower is None:
            assert d_sqrt_lower is None, 'cannot set d_sqrt_lower unless d_lower is set'
            if d_sqrt_upper is None:
                d_sqrt_upper = d_upper
            # s(iN) = -b*sqrt(N/N_norm-iN) + c + d*iN + e*(iN)^2
            # s(0) = 0 = -b*sqrt(N/N_norm) + c
            # ds/diN(N/N_norm) = b/(2*sqrt(N/N_norm-iN))+d+2*e*N/N_norm ~ d_sqrt_upper/sqrt(N/N_norm-iN)+d_upper
            # b = 2*d_sqrt_upper
            # d + 2*e*N/N_norm = d_upper
            # d = d_upper - 2*e*N/N_norm
            # s(N/N_norm) = L = c + d*N/N_norm + e*(N/N_norm)^2
            # L = c + d_upper*N/N_norm - 2*e*(N/N_norm)^2 + e*(N/N_norm)^2
            # e = (c + d_upper*N/N_norm - L) / (N/N_norm)^2
            b = 2.*d_sqrt_upper
            c = b*numpy.sqrt(N/N_norm)
            e = (c + d_upper*N/N_norm - length) / (N/N_norm)**2
            d = d_upper - 2*e*N/N_norm

            # check function is monotonic: gradients at beginning and end should both be
            # positive.
            # lower boundary:
            assert b/(2.*numpy.sqrt(N/N_norm)) + d > 0., 'gradient at start should be positive'
            # upper boundary:
            assert b > 0., 'sqrt part of function should be positive at end'
            assert d + 2.*e*N/N_norm >= 0., 'gradient of polynomial part should be positive at end'

            return lambda i: -b*numpy.sqrt((N-i)/N_norm) + c + d*i/N_norm + e*(i/N_norm)**2
        elif d_upper is None:
            if d_sqrt_lower is None:
                d_sqrt_lower = d_lower
            assert d_sqrt_upper is None
            # s(iN) = a*sqrt(iN) + c + d*iN + e*iN^2
            # s(0) = 0 = c
            # ds/di(0) = a/(2*sqrt(iN))+d ~ d_sqrt_lower/sqrt(iN)+d_lower
            # a = 2*d_sqrt_lower
            # d = d_lower
            # s(N/N_norm) = L = a*sqrt(N/N_norm) + c + d*N/N_norm + e*(N/N_norm)^2
            a = 2.*d_sqrt_lower
            d = d_lower
            e = (length - a*numpy.sqrt(N/N_norm) - d*N/N_norm)/(N/N_norm)**2

            # check function is monotonic: gradients at beginning and end should both be
            # positive.
            # lower boundary:
            assert a > 0., 'sqrt part of function should be positive at start'
            assert d >= 0., 'gradient of polynomial part should be positive at start'
            # upper boundary:
            assert a/(2.*numpy.sqrt(N/N_norm)) + d + 2.*e*N/N_norm > 0., 'gradient at end should be positive'

            return lambda i: a*numpy.sqrt(i/N_norm) + d*i/N_norm + e*(i/N_norm)**2
        else:
            if d_sqrt_lower is None:
                d_sqrt_lower = d_lower
            if d_sqrt_upper is None:
                d_sqrt_upper = d_upper
            # s(iN) = a*sqrt(iN) - b*sqrt(N/N_norm-iN) + c + d*iN + e*iN^2 + f*iN^3
            # s(0) = 0 = -b*sqrt(N/N_norm) + c
            # c = b*sqrt(N/N_norm)
            # ds/diN(0) = a/(2*sqrt(iN))+b/(2*sqrt(N/N_norm))+d ~ d_sqrt_lower/sqrt(iN)+d_lower
            # a = 2*d_sqrt_lower
            # b/(2*sqrt(N/N_norm)) + d = d_lower
            # d = d_lower - b/(2*sqrt(N/N_norm)
            # ds/di(N) =
            # b/(2*sqrt(N/N_norm-iN))+a/(2*sqrt(N))+d+2*e*N/N_norm+3*f*(N/N_norm)^2 ~ d_sqrt_upper/sqrt(N/N_norm-i)+d_upper
            # b = 2*d_sqrt_upper
            # a/(2*sqrt(N/N_norm) + d + 2*e*N/N_norm + 3*f*(N/N_norm)^2 = d_upper
            # e = (d_upper - a/(2*sqrt(N/N_norm)) - d)/(2*N/N_norm) - 3/2*f*N/N_norm
            # s(N/N_norm) = L = a*sqrt(N/N_norm) + c + d*N/N_norm + e*(N/N_norm)^2 + f*(N/N_norm)^3
            # L = a*sqrt(N/N_norm) + c + d*N/N_norm + (d_upper - a/(2*sqrt(N/N_norm)) - d)*N/(2*N_norm) - 3/2*f*(N/N_norm)^3 + f*(N/N_norm)^3
            # f = 2*(a*sqrt(N/N_norm) + c + d*N/(2*N_norm) + d_upper*N/(2*N_norm) - a*sqrt(N/N_norm)/4 - L)*N_norm^3/N^3
            a = 2.*d_sqrt_lower
            b = 2.*d_sqrt_upper
            c = b*numpy.sqrt(N/N_norm)
            d = d_lower - b/2./numpy.sqrt(N/N_norm)
            f = 2.*(a*numpy.sqrt(N/N_norm) + c + d*N/N_norm/2. + d_upper*N/N_norm/2. - a*numpy.sqrt(N/N_norm)/4. - length)*N_norm**3/N**3
            e = (d_upper - a/2./numpy.sqrt(N/N_norm) - d)*N_norm/2./N - 1.5*f*N/N_norm

            # check function is monotonic: gradients at beginning and end should both be
            # positive. Only check the boundaries here, should really add a check that
            # gradient does not reverse in the middle somewhere...
            # lower boundary:
            assert a > 0., 'sqrt part of function should be positive at start'
            assert b/(2.*numpy.sqrt(N/N_norm)) + d >= 0., 'gradient of non-singular part should be positive at start'
            # upper boundary:
            assert b > 0., 'sqrt part of function should be positive at end'
            assert a/(2.*numpy.sqrt(N/N_norm)) + d + 2.*e*N/N_norm + 3.*f*(N/N_norm)**2 >= 0., 'gradient of non-singular part should be positive at end'

            return lambda i: a*numpy.sqrt(i/N_norm) - b*numpy.sqrt((N-i)/N_norm) + c + d*i/N_norm + e*(i/N_norm)**2 + f*(i/N_norm)**3

    def readOption(self, name, default=None):
        print('reading option', name, end='')
        try:
            value = self.options[name]
            print(':', value)
            return value
        except KeyError:
            print(' - not set - setting default:', default)
            return default

    def plotPotential(self, Rmin=None, Rmax=None, Zmin=None, Zmax=None, npoints=100,
            ncontours=40):
        from matplotlib import pyplot

        if Rmin is None: Rmin = self.Rmin
        if Rmax is None: Rmax = self.Rmax
        if Zmin is None: Zmin = self.Zmin
        if Zmax is None: Zmax = self.Zmax

        R = numpy.linspace(Rmin, Rmax, npoints)
        Z = numpy.linspace(Zmin, Zmax, npoints)
        ax = pyplot.axes(aspect='equal')
        contours = ax.contour(
                R, Z, self.psi(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T, ncontours)
        pyplot.clabel(contours, inline=False, fmt='%1.3g')

    def plotSeparatrix(self):
        from matplotlib import pyplot
        for region in self.regions.values():
            R = [p.R for p in region]
            Z = [p.Z for p in region]
            pyplot.scatter(R, Z, marker='x', label=region.name)

class DoubleNull(Equilibrium):
    """
    Analyse tokamak equilibrium - single-null, connected double-null or disconnected
    double-null - with regions lined up according to BOUT++ requirements.
    """

    def __init__(self, options):
        self.options = options

        self.nx_core = self.readOption('nx_core')
        self.nx_between = self.readOption('nx_between')
        self.nx_sol = self.readOption('nx_sol')
        self.ny_inner_lower_divertor = self.readOption('ny_inner_lower_divertor')
        self.ny_inner_core = self.readOption('ny_inner_core')
        self.ny_inner_upper_divertor = self.readOption('ny_inner_upper_divertor')
        self.ny_outer_upper_divertor = self.readOption('ny_outer_upper_divertor')
        self.ny_outer_core = self.readOption('ny_outer_core')
        self.ny_outer_lower_divertor = self.readOption('ny_outer_lower_divertor')

        self.psi_core = self.readOption('psi_core')
        self.psi_lower_pf = self.readOption('psi_lower_pf', self.psi_core)
        self.psi_upper_pf = self.readOption('psi_upper_pf', self.psi_core)

        self.psi_sol = self.readOption('psi_sol')
        self.psi_inner_sol = self.readOption('psi_inner_sol', self.psi_sol)
        # this option can only be set different from psi_sol in a double-null
        # configuration (i.e. if there are upper divertor legs)
        if self.psi_sol != self.psi_inner_sol:
            assert self.ny_inner_upper_divertor > 0 or self.ny_outer_upper_divertor > 0, 'inner and outer SOL should be separated by an upper divertor, i.e. topology should be double-null'

        self.psi_spacing_separatrix_multiplier = self.readOption('psi_spacing_separatrix_multiplier', None)
        if self.psi_spacing_separatrix_multiplier is not None:
            if self.nx_between > 0:
                raise ValueError("Cannot use psi_spacing_separatrix_multiplier when "
                                 "there are points between two separatrices - to get "
                                 "the same effect, increase the number of points in the "
                                 "between-separatrix region.")

        self.y_boundary_guards = self.readOption('y_boundary_guards')

        if ( (self.ny_inner_upper_divertor > 0 or self.ny_outer_upper_divertor > 0)
                and (self.ny_inner_lower_divertor > 0 or self.ny_inner_core > 0
                     or self.ny_outer_core > 0 or self.ny_outer_lower_divertor > 0) ):
            # there are two distinct divertor/limiter targets
            #  - the upper divertor legs are not empty, and also the other regions are not
            #    all empty
            self.upper_target_y_boundary_guards = self.y_boundary_guards
