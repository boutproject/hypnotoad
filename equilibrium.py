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

# Dictionary of global parameters for contours
ContourParameters = {'Nfine':1000}

# tolerance used to try and avoid missed intersections between lines
# also if two sets of lines appear to intersect twice, only count it once if the
# distance between the intersections is less than this
intersect_tolerance = 1.e-14

class PoloidalSpacingParameters:
    def __init__(self):
        # Method to use for poloidal spacing function:
        #  - 'sqrt' for getSqrtPoloidalSpacingFunction
        #  - 'polynomial' for getPolynomialPoloidalSpacingFunction
        self.method = 'sqrt'

        ## Parameters for sqrt spacing function
        # Distance for polynomial part of spacing function at lower end
        self.sqrt_b_lower = None

        # Distance for polynomial part of spacing function at upper end
        self.sqrt_b_upper = None

        # Distance for sqrt part of spacing function (if used) at lower end
        self.sqrt_a_lower = None

        # Distance for sqrt part of spacing function (if used) at upper end
        self.sqrt_a_upper = None

        ## Parameters for polynomial spacing function
        # Distance for spacing function at lower end
        self.polynomial_d_lower = None

        # Distance for spacing function at upper end
        self.polynomial_d_upper = None

        # Normalization factor for number of points in contours, used to scale grid
        # spacing with total number of points, to keep functions consistent when
        # resolution is changed
        self.N_norm = None

        # Method to use for spacing of points in non-orthogonal grid
        self.nonorthogonal_method = 'combined'

        # Distance for transition between fixed-poloidal-spacing grid and orthogonal grid
        # at the lower end. If 'None' then the value of polynomial_d_lower will be used instead.
        self.nonorthogonal_range_lower = None

        # Distance for transition between fixed-poloidal-spacing grid and orthogonal grid
        # at the upper end. If 'None' then the value of polynomial_d_upper will be used instead.
        self.nonorthogonal_range_upper = None

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

    def as_ndarray(self):
        return numpy.array((self.R, self.Z))

def calc_distance(p1, p2):
    d = p2 - p1
    return numpy.sqrt(d.R**2 + d.Z**2)

def swap_points(p1, p2):
    tempR = p1.R
    tempZ = p1.Z
    p1.R = p2.R
    p1.Z = p2.Z
    p2.R = tempR
    p2.Z = tempZ

def find_intersection(l1start, l1end, l2start, l2end):
    """
    Find the intersection (if there is one) between the lines 'l1' and 'l2'
    """
    # Copy so we don't change the inputs
    l1start = deepcopy(l1start)
    l1end = deepcopy(l1end)
    l2start = deepcopy(l2start)
    l2end = deepcopy(l2end)

    R2 = l2start.R
    Z2 = l2start.Z
    dR2 = l2end.R - l2start.R
    dZ2 = l2end.Z - l2start.Z

    if numpy.abs(l1end.R - l1start.R) > numpy.abs(l1end.Z - l1start.Z):
        # if l1 is sensible, dR1 shouldn't be too small as it's bigger than dZ1
        # l1 is Z = Z1 + dZ1/dR1 * (R - R1)

        # sort l1 points in R
        if l1start.R > l1end.R:
            swap_points(l1start, l1end)
        R1 = l1start.R
        Z1 = l1start.Z
        dR1 = l1end.R - l1start.R
        dZ1 = l1end.Z - l1start.Z

        if numpy.abs(l2end.R - l2start.R) > numpy.abs(l2end.Z - l2start.Z):
            # if l2 is sensible, dR2 shouldn't be too small as it's bigger than dZ2
            # l2 is Z = Z2 + dZ2/dR2 * (R - R2)

            # sort l2 points in R
            if l2start.R > l2end.R:
                swap_points(l2start, l2end)
            R2 = l2start.R
            Z2 = l2start.Z
            dR2 = l2end.R - l2start.R
            dZ2 = l2end.Z - l2start.Z

            if numpy.abs(dZ1/dR1 - dZ2/dR2) < 1.e-15:
                # lines are parallel so cannot intersect
                return None

            # intersection where
            # Z1 + dZ1/dR1 * (R - R1) = Z2 + dZ2/dR2 * (R - R2)
            # (dZ1/dR1 - dZ2/dR2)*R = Z2 - Z1 + dZ1/dR1*R1 - dZ2/dR2*R2
            R = (Z2 - Z1 + dZ1/dR1*R1 - dZ2/dR2*R2) / (dZ1/dR1 - dZ2/dR2)
            if (R >= R1 - intersect_tolerance and R <= l1end.R + intersect_tolerance and
                    R >= R2 - intersect_tolerance and R <= l2end.R + intersect_tolerance):
                Z = Z1 + dZ1/dR1 * (R - R1)
                return Point2D(R, Z)
            else:
                return None
        else:
            # if l2 is sensible, dZ2 shouldn't be too small as it's bigger than dR2
            # l2 is R = R2 + dR2/dZ2 * (Z - Z2)

            # sort l2 points in Z
            if l2start.Z > l2end.Z:
                swap_points(l2start, l2end)
            R2 = l2start.R
            Z2 = l2start.Z
            dR2 = l2end.R - l2start.R
            dZ2 = l2end.Z - l2start.Z

            # intersection where
            # Z = Z1 + dZ1/dR1 * (R2 + dR2/dZ2 * (Z - Z2) - R1)
            # (1 - dZ1*dR2/dR1/dZ2) * Z = Z1 + dZ1/dR1 * (R2 - dR2/dZ2*Z2 - R1)
            Z = (Z1 + dZ1/dR1 * (R2 - dR2/dZ2*Z2 - R1)) / (1. - dZ1*dR2/(dR1*dZ2))
            R = R2 + dR2/dZ2 * (Z - Z2)
            if (R >= R1 - intersect_tolerance and R <= l1end.R + intersect_tolerance and
                    Z >= Z2 - intersect_tolerance and Z <= l2end.Z + intersect_tolerance):
                return Point2D(R, Z)
            else:
                return None
    else:
        # if l1 is sensible, dZ1 shouldn't be too small as it's bigger than dR1
        # l1 is R = R1 + dR1/dZ1 * (Z - Z1)

        # sort l1 points in Z
        if l1start.Z > l1end.Z:
            swap_points(l1start, l1end)
        R1 = l1start.R
        Z1 = l1start.Z
        dR1 = l1end.R - l1start.R
        dZ1 = l1end.Z - l1start.Z

        if numpy.abs(dR2) > numpy.abs(dZ2):
            # if l2 is sensible, dR2 shouldn't be too small as it's bigger than dZ2
            # l2 is Z = Z2 + dZ2/dR2 * (R - R2)

            # sort l2 points in R
            if l2start.R > l2end.R:
                swap_points(l2start, l2end)
            R2 = l2start.R
            Z2 = l2start.Z
            dR2 = l2end.R - l2start.R
            dZ2 = l2end.Z - l2start.Z

            # intersection where
            # R = R1 + dR1/dZ1 * (Z2 + dZ2/dR2 * (R - R2) - Z1)
            # (1 - dR1/dZ1*dZ2/dR2) * R = R1 + dR1/dZ1 * (Z2 - dZ2/dR2*R2 - Z1)
            R = (R1 + dR1/dZ1 * (Z2 - dZ2/dR2*R2 - Z1)) / (1. - dR1/dZ1 * dZ2/dR2)
            Z = Z2 + dZ2/dR2 * (R - R2)
            if (Z >= Z1 - intersect_tolerance and Z <= l1end.Z + intersect_tolerance and
                    R >= R2 - intersect_tolerance and R <= l2end.R + intersect_tolerance):
                return Point2D(R, Z)
            else:
                return None
        else:
            # if l2 is sensible, dZ2 shouldn't be too small as it's bigger than dR2
            # l2 is R = R2 + dR2/dZ2 * (Z - Z2)

            # sort l2 points in Z
            if l2start.Z > l2end.Z:
                swap_points(l2start, l2end)
            R2 = l2start.R
            Z2 = l2start.Z
            dR2 = l2end.R - l2start.R
            dZ2 = l2end.Z - l2start.Z

            if numpy.abs(dR2/dZ2 - dR1/dZ1) < 1.e-15:
                # lines are parallel so cannot intersect
                return None

            # intersection where
            # R2 + dR2/dZ2 * (Z - Z2) = R1 + dR1/dZ1 * (Z - Z1)
            # (dR2/dZ2 - dR1*dZ1) * Z = R1 - R2 + dR2/dZ2*Z2 - dR1/dZ1*Z1
            Z = (R1 - R2 + dR2/dZ2*Z2 - dR1/dZ1*Z1) / (dR2/dZ2 - dR1/dZ1)
            if (Z >= Z1 - intersect_tolerance and Z <= l1end.Z + intersect_tolerance and
                    Z >= Z2 - intersect_tolerance and Z <= l2end.Z + intersect_tolerance):
                R = R2 + dR2/dZ2 * (Z - Z2)
                return Point2D(R, Z)
            else:
                return None

class PsiContour:
    """
    Represents a contour as a collection of points.
    Includes methods for interpolation.
    Mostly behaves like a list
    """
    def __init__(self, points, psi, psival):
        self.points = points

        self.startInd = 0
        self.endInd = len(points) - 1

        self.recalculateDistance()

        # Function that evaluates the vector potential at R,Z
        self.psi = psi

        # Value of vector potential on this contour
        self.psival = psival

        # Number of boundary guard cells at either end
        # This may be set even if the contour has not been extended yet, to specify how
        # many guard cells should be added when it is - this is extra information to
        # startInd and endInd.
        self.extend_lower = 0
        self.extend_upper = 0

    def __iter__(self):
        return self.points.__iter__()

    def __str__(self):
        return self.points.__str__()

    def __getitem__(self, key):
        return self.points.__getitem__(key)

    def __len__(self):
        return self.points.__len__()

    def setSelfToContour(self, contour):
        """
        Copy the state of this object from contour
        """
        self.points = deepcopy(contour.points)
        self.startInd = contour.startInd
        self.endInd = contour.endInd
        self.distance = contour.distance
        self.psi = contour.psi
        self.psival = contour.psival
        self.extend_lower = contour.extend_lower
        self.extend_upper = contour.extend_upper

    def newContourFromSelf(self, *, points=None, psival=None):
        if points is None:
            points = deepcopy(self.points)
        if psival is None:
            psival = self.psival
        new_contour = PsiContour(points, self.psi, psival)

        new_contour.startInd = self.startInd
        new_contour.endInd = self.endInd
        new_contour.extend_lower = self.extend_lower
        new_contour.extend_upper = self.extend_upper

        return new_contour

    def append(self, point):
        self.points.append(point)
        if len(self.distance) is 0:
            self.distance = [0.]
        else:
            self.distance.append(
                    self.distance[-1] + calc_distance(self.points[-2], self.points[-1]))

    def prepend(self, point):
        self.points.insert(0, point)
        self.recalculateDistance()

    def insert(self, index, point):
        # Make sure index is positive, following behaviour of list.insert()
        if index < 0:
            index += len(self) + 1
            if index < 0:
                index = 0

        self.points.insert(index, point)
        if index <= self.startInd:
            self.startInd += 1
        if index <= self.endInd:
            self.endInd += 1
        if self.endInd < 0 and index > len(self) + self.endInd:
            self.endInd -= 1

        self.recalculateDistance()

    def recalculateDistance(self):
        if len(self.points) is 0:
            self.distance = []
        else:
            self.distance = [0.]
            for i in range(1, len(self.points)):
                self.distance.append(
                        self.distance[-1] + calc_distance(self.points[i-1], self.points[i]))

    def totalDistance(self):
        return self.distance[self.endInd] - self.distance[self.startInd]

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
                    print('width =',width)
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

        return self.newContourFromSelf(points=newpoints)

    def interpFunction(self):
        distance = numpy.array(numpy.float64(self.distance)) - self.distance[self.startInd]
        R = numpy.array(numpy.float64([p.R for p in self.points]))
        Z = numpy.array(numpy.float64([p.Z for p in self.points]))
        interpR = interp1d(distance, R, kind='cubic',
                           assume_sorted=True, fill_value='extrapolate')
        interpZ = interp1d(distance, Z, kind='cubic',
                           assume_sorted=True, fill_value='extrapolate')
        return lambda s: Point2D(interpR(s), interpZ(s))

    def contourSfunc(self):
        """
        Function interpolating distance as a function of index for the current state of
        this contour. When outside [startInd, endInd], set to constant so the results
        aren't affected by extrapolation errors.
        """
        interpS = interp1d(numpy.arange(len(self), dtype=float), self.distance, kind='cubic', assume_sorted=True,
                fill_value='extrapolate')
        thisStartInd = self.startInd
        thisEndInd = self.endInd
        if thisEndInd < 0:
            # endInd might be negative, which would mean relative to the end of the list,
            # but we need the actual index below
            thisEndInd += len(self)
        startDistance = self.distance[thisStartInd]
        endDistance = self.distance[thisEndInd]
        return lambda i: numpy.piecewise(i, [i <= 0., i >= thisEndInd - thisStartInd],
                [0., endDistance - startDistance,
                 lambda i: interpS(i + thisStartInd) - startDistance])

    def interpSSperp(self, vec):
        """
        Returns:
        1. a function s(s_perp) for interpolating the poloidal distance along the contour
           from the distance perpendicular to vec.
           's_perp' is modified to be a monotonically increasing function along the
           contour.
        2. the total perpendicular distance between startInd and endInd of the contour.
        """

        # vec_perp is a vector in the direction of either increasing or decreasing sperp
        vec_perp = numpy.zeros(2)
        vec_perp[0] = -vec[1]
        vec_perp[1] = vec[0]

        # make vec_perp a unit vector
        vec_perp = vec_perp / numpy.sqrt(numpy.sum(vec_perp**2))
        Rstart = self[self.startInd].R
        Zstart = self[self.startInd].Z

        # s_perp = (vec_perp).(r) where r is the displacement vector of each point from self[self.startInd]
        s_perp = numpy.array([vec_perp[0]*(p.R - Rstart) + vec_perp[1]*(p.Z - Zstart) for p in self])

        # s_perp might not be monotonic in which case s(s_perp) is not well defined.
        # To get around this, if d(s_perp) between two points is negative, flip its sign
        # to make a fake 's_perp' that is always increasing.
        # Note we only need s_perp to be good near one of the ends, the function using it
        # will be multiplied by a weight that goes to zero far from the end.
        # This correction means s_perp is always increasing, regardless of sign of
        # vec_perp, so don't need to check sign of vec_perp when creating it.
        for i in range(self.startInd + 1, len(s_perp)):
            ds = s_perp[i] - s_perp[i-1]
            if ds < 0.:
                s_perp[i:] = 2.*s_perp[i-1] - s_perp[i:]
        for i in range(self.startInd - 1, -1, -1):
            ds = s_perp[i+1] - s_perp[i]
            if ds < 0.:
                s_perp[:i+1] = 2.*s_perp[i+1] - s_perp[:i+1]

        s_perp_total = s_perp[self.endInd] - s_perp[self.startInd]

        distance = numpy.array(self.distance) - self.distance[self.startInd]
        s_of_sperp = interp1d(s_perp, distance, kind='cubic', assume_sorted=True,
                fill_value = 'extrapolate')

        return s_of_sperp, s_perp_total

    def regrid(self, *args, **kwargs):
        """
        Regrid this contour, modifying the object
        """
        self.setSelfToContour(self.getRegridded(*args, **kwargs))
        return self

    def getRegridded(self, npoints, *, width=1.e-5, atol=2.e-8, sfunc=None,
            extend_lower=None, extend_upper=None):
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
        if extend_lower is not None:
            self.extend_lower = extend_lower
        if extend_upper is not None:
            self.extend_upper = extend_upper

        Nfine = ContourParameters['Nfine']
        # To make the new points accurate, first regrid onto a high-resolution contour,
        # then interpolate.
        # Extend further than will be needed in the final contour, because extrapolation
        # past the end of the fine contour is very bad.
        extend_lower_fine = 2*(self.extend_lower * Nfine) // npoints
        extend_upper_fine = 2*(self.extend_upper * Nfine) // npoints

        indices_fine = numpy.linspace(-extend_lower_fine,
                (Nfine - 1 + extend_upper_fine),
                Nfine + extend_lower_fine + extend_upper_fine)

        if sfunc is not None:
            sfine = sfunc(indices_fine*(npoints - 1)/(Nfine - 1))
        else:
            sfine = (self.distance[self.endInd] - self.distance[self.startInd]) / (Nfine - 1) * indices_fine

        interp_self = self.interpFunction()

        fine_contour = PsiContour([interp_self(x) for x in sfine], self.psi, self.psival)
        fine_contour.startInd = extend_lower_fine
        fine_contour.endInd = len(fine_contour) - 1 - extend_upper_fine
        fine_contour.refine(width, atol)

        # fine_contour does not have exactly the same total distance as the original,
        # because it integrates on a finer grid. But sfunc was defined with the total
        # distance of the original, so scale the distance of fine_contour to make the
        # totals the same
        if sfunc is not None:
            scale_factor = (sfunc(npoints - 1.) - sfunc(0.)) / (fine_contour.distance[fine_contour.endInd] - fine_contour.distance[fine_contour.startInd])
        else:
            scale_factor = (self.distance[self.endInd] - self.distance[self.startInd]) / (fine_contour.distance[fine_contour.endInd] - fine_contour.distance[fine_contour.startInd])
        fine_contour.distance = [scale_factor*d for d in fine_contour.distance]

        indices = numpy.linspace(-self.extend_lower, (npoints - 1 + self.extend_upper),
                npoints + self.extend_lower + self.extend_upper)
        if sfunc is not None:
            s = sfunc(indices)

            # offset fine_contour.interpFunction in case sfunc(0.)!=0.
            sbegin = sfunc(0.)
        else:
            s = (self.distance[self.endInd] - self.distance[self.startInd]) / (npoints - 1) * indices
            sbegin = 0.

        interp_fine_unadjusted = fine_contour.interpFunction()
        interp_fine = lambda s: interp_fine_unadjusted(s - sbegin)

        new_contour = self.newContourFromSelf(points=[interp_fine(x) for x in s])
        new_contour.startInd = self.extend_lower
        new_contour.endInd = len(new_contour) - 1 - self.extend_upper
        # new_contour was interpolated from a high-resolution contour, so should not need
        # a large width for refinement - use 1.e-5 instead of 'width'
        return new_contour.getRefined(1.e-5, atol)

    def temporaryExtend(self, *, extend_lower=0, extend_upper=0):
        """
        Add temporary guard-cell points to the beginning and/or end of a contour
        """
        if extend_lower > 0:
            ds = self.distance[1] - self.distance[0]
            for i in range(extend_lower):
                new_point = self.interpFunction()(self.distance[0] - ds)
                self.points.insert(0, new_point)
                self.startInd += 1
                self.endInd += 1
                self.refine() # also re-calculates self.distance
        if extend_upper > 0:
            ds = self.distance[-1] - self.distance[-2]
            for i in range(extend_upper):
                new_point = self.interpFunction()(self.distance[-1] + ds)
                self.points.append(new_point)
                self.refine()

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
    def __init__(self, equilibrium, name, nSegments, options, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.equilibrium = equilibrium
        self.name = name
        self.nSegments = nSegments
        self.options = options
        self.nx = self.options['nx']
        self.ny_noguards = self.options['ny']
        self.y_boundary_guards = self.options['y_boundary_guards']

        self.xPointsAtStart = []
        self.xPointsAtEnd = []

        # Set if this segment starts on a wall, with value of vector along wall
        self.surfaceAtStart = None

        # Set if this segment ends on a wall, with value of vector along wall
        self.surfaceAtEnd = None

        self.connections = []
        self.psi_vals = []
        self.separatrix_radial_index = 0

        # parameters for poloidal distance functions
        self.poloidalSpacingParameters = PoloidalSpacingParameters()

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

    def newRegionFromPsiContour(self, contour):
        result = EquilibriumRegion(self.equilibrium, self.name, self.nSegments, self.options,
                contour.points, contour.psi, contour.psival)
        result.xPointsAtStart = deepcopy(self.xPointsAtStart)
        result.xPointsAtEnd = deepcopy(self.xPointsAtEnd)
        result.connections = deepcopy(self.connections)
        result.psi_vals = deepcopy(self.psi_vals)
        result.poloidalSpacingParameters = self.poloidalSpacingParameters
        result.separatrix_radial_index = self.separatrix_radial_index
        result.startInd = contour.startInd
        result.endInd = contour.endInd
        result.extend_lower = contour.extend_lower
        result.extend_upper = contour.extend_upper
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
        return self.newRegionFromPsiContour(super().getRefined(*args, **kwargs))

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
        sfunc = self.getSfuncFixedSpacing(2*self.ny_noguards+1,
                self.distance[self.endInd] - self.distance[self.startInd])
        return self.newRegionFromPsiContour(super().getRegridded(2*self.ny_noguards + 1,
            extend_lower=extend_lower, extend_upper=extend_upper, sfunc=sfunc, **kwargs))

    def getSfuncFixedSpacing(self, npoints, distance, *, method=None):
        if method is None:
            if self.equilibrium.orthogonal:
                method = self.poloidalSpacingParameters.method
            else:
                method = 'nonorthogonal'

        if method == 'sqrt':
            return self.getSqrtPoloidalDistanceFunc(distance, npoints-1,
                    self.poloidalSpacingParameters.N_norm,
                    b_lower=self.poloidalSpacingParameters.sqrt_b_lower,
                    a_lower=self.poloidalSpacingParameters.sqrt_a_lower,
                    b_upper=self.poloidalSpacingParameters.sqrt_b_upper,
                    a_upper=self.poloidalSpacingParameters.sqrt_a_upper)
        elif method == 'polynomial':
            return self.getPolynomialPoloidalDistanceFunc(distance, npoints-1,
                    self.poloidalSpacingParameters.N_norm,
                    d_lower=self.poloidalSpacingParameters.polynomial_d_lower,
                    d_upper=self.poloidalSpacingParameters.polynomial_d_upper)
        elif method == 'nonorthogonal':
            nonorth_method = self.poloidalSpacingParameters.nonorthogonal_method
            if nonorth_method == 'poloidal_orthogonal_combined':
                raise ValueError('nonorth_method='+nonorth_method+' not implemented yet')
            elif nonorth_method == 'combined':
                if self.surfaceAtStart is not None:
                    # surface is a wall
                    lower_surface = self.surfaceAtStart
                else:
                    # No wall, surface is orthogonal to flux surfaces
                    lower_surface = [self.equilibrium.f_R(*self[self.startInd]),
                                     self.equilibrium.f_Z(*self[self.startInd])]

                if self.surfaceAtEnd is not None:
                    # surface is a wall
                    upper_surface = self.surfaceAtEnd
                else:
                    # No wall, surface is orthogonal to flux surfaces
                    upper_surface = [self.equilibrium.f_R(*self[self.endInd]),
                                     self.equilibrium.f_Z(*self[self.endInd])]

                return self.combineSfuncsPerpSpacing(self, lower_surface, upper_surface,
                        None)
            elif nonorth_method == 'orthogonal':
                orth_method = self.poloidalSpacingParameters.method
                return self.getSfuncFixedSpacing(npoints, distance, method=orth_method)
            else:
                return self.getSfuncFixedSpacing(npoints, distance, method=nonorth_method)
        else:
            raise ValueError('Unrecognized option '
                             +str(self.poloidalSpacingParameters.method)
                             +' for poloidal spacing method')

    def combineSfuncsPoloidalSpacing(self, sfunc_orthogonal, total_distance):

        # this sfunc gives a fixed poloidal spacing at beginning and end of contours
        sfunc_fixed_spacing = self.equilibriumRegion.getSfuncFixedSpacing(
                2*self.ny_noguards + 1, total_distance, method='polynomial')

        if self.poloidalSpacingParameters.nonorthogonal_range_lower is not None:
            range_lower = self.poloidalSpacingParameters.nonorthogonal_range_lower
        else:
            range_lower = self.poloidalSpacingParameters.polynomial_d_lower

        if self.poloidalSpacingParameters.nonorthogonal_range_upper is not None:
            range_upper = self.poloidalSpacingParameters.nonorthogonal_range_upper
        else:
            range_upper = self.poloidalSpacingParameters.polynomial_d_upper

        if range_lower is not None and range_upper is not None:
            def new_sfunc(i):
                sfixed = sfunc_fixed_spacing(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_lower so it is 1. at the lower boundary and 0. at the
                # upper boundary and the gradient is zero at both boundaries
                weight_lower = numpy.piecewise(sfixed,
                        [sfixed < 0., sfixed > total_distance],
                        [1., 0., lambda s: numpy.exp(-(s/range_lower)**2)
                            * (0.5 + 0.5*numpy.cos(numpy.pi*s/total_distance))])

                # define weight_upper so it is 1. at the upper boundary and 0. at the
                # lower boundary and the gradient is zero at both boundaries
                weight_upper = numpy.piecewise(sfixed,
                    [sfixed < 0., sfixed > total_distance],
                    [0., 1., lambda s: numpy.exp(-((total_distance - s)/range_upper)**2)
                        * (0.5 - 0.5*numpy.cos(numpy.pi*s/total_distance))])

                # make sure weight_lower + weight_upper <= 1
                weight = weight_lower + weight_upper
                weight_over_slice = weight[weight > 1.]
                weight_lower[weight > 1.] /= weight_over_slice
                weight_upper[weight > 1.] /= weight_over_slice

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncsPoloidalSpacing again
                    # for this contour with sfunc_orthogonal from self.contourSfunc() then
                    # we get the same spacing again. This is used to make the contours
                    # along the separatrix keep the same values when pushing the other
                    # contours towards orthogonal positions
                    # s = weight_lower*sfixed + weight_upper*sfixed + (1. - weight_lower - weight_upper)*s
                    sorth = (weight_lower*sfixed + weight_upper*sfixed) / (weight_lower + weight_upper)

                return weight_lower*sfixed + weight_upper*sfixed + (1. - weight_lower - weight_upper)*sorth
        elif range_lower is not None:
            def new_sfunc(i):
                sfixed = sfunc_fixed_spacing(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_lower so it is 1. at the lower boundary and 0. at the
                # upper boundary
                weight_lower = numpy.piecewise(sfixed,
                        [sfixed < 0., sfixed > total_distance],
                        [1., 0., lambda s: numpy.exp(-s/range_lower)
                                           * (1. - s/total_distance)])

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncsPoloidalSpacing again
                    # for this contour with sfunc_orthogonal from self.contourSfunc() then
                    # we get the same spacing again. This is used to make the contours
                    # along the separatrix keep the same values when pushing the other
                    # contours towards orthogonal positions
                    # s = weight_lower*sfixed + (1. - weight_lower)*s
                    sorth = sfixed

                return (weight_lower)*sfixed + (1. - weight_lower) * sorth
        elif range_upper is not None:
            def new_sfunc(i):
                sfixed = sfunc_fixed_spacing(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_upper so it is 1. at the upper boundary and 0. at the
                # lower boundary
                weight_upper = numpy.piecewise(sfixed,
                    [sfixed < 0., sfixed > total_distance],
                    [0., 1., lambda s: numpy.exp(-(total_distance - s)/range_upper)
                                       * s/total_distance])

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncsPoloidalSpacing again
                    # for this contour with sfunc_orthogonal from self.contourSfunc() then
                    # we get the same spacing again. This is used to make the contours
                    # along the separatrix keep the same values when pushing the other
                    # contours towards orthogonal positions
                    # s = weight_upper*sfixed + (1. - weight_upper)*s
                    sorth = sfixed

                return (weight_upper)*sfixed + (1. - weight_upper) * sorth
        else:
            assert sfunc_orthogonal is not None, 'Without range_lower or range_upper, cannot use with sfunc_orthogonal=None'
            def new_sfunc(i):
                return sfunc_orthogonal(i)

        # Check new_sfunc is monotonically increasing
        indices = numpy.arange(-self.extend_lower, len(self)-self.extend_lower, dtype=float)
        scheck = new_sfunc(indices)
        if numpy.any(scheck[1:] < scheck[:-1]):
            from matplotlib import pyplot
            print('check ranges',range_lower, range_upper)
            pyplot.figure()
            pyplot.plot(scheck)
            pyplot.axhline(0.)
            pyplot.axhline(total_distance)
            pyplot.show()
            decreasing = numpy.where(scheck[1:] < scheck[:-1])[0] + 1
            raise ValueError('In region ' + self.name + 'combined spacing function is '
                    + 'decreasing at indices ' + str(decreasing) + ' on contour of '
                    + 'length ' + str(len(self)) + '. It may help to increase '
                    + 'nonorthogonal_range_target or nonorthogonal_range_xpoint.')

        return new_sfunc

    def combineSfuncsPerpSpacing(self, contour, vec_lower, vec_upper, sfunc_orthogonal):

        sfunc_fixed_lower, sperp_func_lower = self.getSfuncFixedPerpSpacing(
                2*self.ny_noguards + 1, contour, vec_lower, True)
        total_perp_distance_lower = sperp_func_lower(2.*self.ny_noguards + 1.)

        sfunc_fixed_upper, sperp_func_upper = self.getSfuncFixedPerpSpacing(
                2*self.ny_noguards + 1, contour, vec_upper, False)
        total_perp_distance_upper = sperp_func_upper(2.*self.ny_noguards + 1.)

        if self.poloidalSpacingParameters.nonorthogonal_range_lower is not None:
            range_lower = self.poloidalSpacingParameters.nonorthogonal_range_lower
        else:
            range_lower = self.poloidalSpacingParameters.polynomial_d_lower

        if self.poloidalSpacingParameters.nonorthogonal_range_upper is not None:
            range_upper = self.poloidalSpacingParameters.nonorthogonal_range_upper
        else:
            range_upper = self.poloidalSpacingParameters.polynomial_d_upper

        if range_lower is not None and range_upper is not None:
            def new_sfunc(i):
                sfixed_lower = sfunc_fixed_lower(i)
                sperp_lower = sperp_func_lower(i)

                sfixed_upper = sfunc_fixed_upper(i)
                sperp_upper = sperp_func_upper(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_lower so it is 1. at the lower boundary and 0. at the
                # upper boundary and the gradient is zero at both boundaries
                weight_lower = numpy.piecewise(sperp_lower,
                        [sperp_lower < 0., sperp_lower > total_perp_distance_lower],
                        [1., 0., lambda s: numpy.exp(-(s/range_lower)**2)
                            * (0.5 + 0.5*numpy.cos(numpy.pi*s/total_perp_distance_lower))])

                # define weight_upper so it is 1. at the upper boundary and 0. at the
                # lower boundary and the gradient is zero at both boundaries
                weight_upper = numpy.piecewise(sperp_upper,
                    [sperp_upper < 0., sperp_upper > total_perp_distance_upper],
                    [0., 1., lambda s: numpy.exp(-((total_perp_distance_upper - s)/range_upper)**2)
                        * (0.5 - 0.5*numpy.cos(numpy.pi*s/total_perp_distance_upper))])

                # make sure weight_lower + weight_upper <= 1
                weight = weight_lower + weight_upper
                weight_over_slice = weight[weight > 1.]
                weight_lower[weight > 1.] /= weight_over_slice
                weight_upper[weight > 1.] /= weight_over_slice

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncs again for this contour
                    # with sfunc_orthogonal from self.contourSfunc() then we get the same
                    # spacing again. This is used to make the contours along the
                    # separatrix keep the same values when pushing the other contours
                    # towards orthogonal positions
                    # s = weight_lower*sfixed_lower + weight_upper*sfixed_upper + (1. - weight_lower - weight_upper)*s
                    sorth = (weight_lower*sfixed_lower + weight_upper*sfixed_upper) / (weight_lower + weight_upper)

                return weight_lower*sfixed_lower + weight_upper*sfixed_upper + (1. - weight_lower - weight_upper)*sorth
        elif range_lower is not None:
            def new_sfunc(i):
                sfixed_lower = sfunc_fixed_lower(i)
                sperp_lower = sperp_func_lower(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_lower so it is 1. at the lower boundary and 0. at the
                # upper boundary
                weight_lower = numpy.piecewise(sperp_lower,
                        [sperp_lower < 0., sperp_lower > total_perp_distance_lower],
                        [1., 0., lambda s: numpy.exp(-s/range_lower)
                                           * (1. - s/total_perp_distance_lower)])

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncs again for this contour
                    # with sfunc_orthogonal from self.contourSfunc() then we get the same
                    # spacing again. This is used to make the contours along the
                    # separatrix keep the same values when pushing the other contours
                    # towards orthogonal positions
                    # s = weight_lower*sperp_lower + (1. - weight_lower)*s
                    sorth = sperp_lower

                return (weight_lower)*sfixed_lower + (1. - weight_lower) * sorth
        elif range_upper is not None:
            def new_sfunc(i):
                sfixed_upper = sfunc_fixed_upper(i)
                sperp_upper = sperp_func_upper(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_upper so it is 1. at the upper boundary and 0. at the
                # lower boundary
                weight_upper = numpy.piecewise(sperp_upper,
                    [sperp_upper < 0., sperp_upper > total_perp_distance_upper],
                    [0., 1., lambda s: numpy.exp(-(total_perp_distance_upper - s)/range_upper)
                                       * s/total_perp_distance_upper])

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncs again for this contour
                    # with sfunc_orthogonal from self.contourSfunc() then we get the same
                    # spacing again. This is used to make the contours along the
                    # separatrix keep the same values when pushing the other contours
                    # towards orthogonal positions
                    # s = weight_upper*sperp_upper + (1. - weight_upper)*s
                    sorth = sperp_upper

                return (weight_upper)*sfixed_upper + (1. - weight_upper) * sorth
        else:
            assert sfunc_orthogonal is not None, 'Without range_lower or range_upper, cannot use with sfunc_orthogonal=None'
            def new_sfunc(i):
                return sfunc_orthogonal(i)

        # Check new_sfunc is monotonically increasing
        indices = numpy.arange(-self.extend_lower, len(self)-self.extend_lower, dtype=float)
        scheck = new_sfunc(indices)
        if numpy.any(scheck[1:] < scheck[:-1]):
            from matplotlib import pyplot
            print('check ranges',range_lower, range_upper)
            pyplot.figure()
            pyplot.plot(scheck)
            pyplot.axhline(0.)
            pyplot.axhline(contour.totalDistance())
            pyplot.plot(sperp_func_lower(indices))
            pyplot.plot(sperp_func_upper(indices))
            pyplot.show()
            decreasing = numpy.where(scheck[1:] < scheck[:-1])[0] + 1
            raise ValueError('In region ' + self.name + 'combined spacing function is '
                    + 'decreasing at indices ' + str(decreasing) + ' on contour of '
                    + 'length ' + str(len(self)) + '. It may help to increase '
                    + 'nonorthogonal_range_target or nonorthogonal_range_xpoint.')

        return new_sfunc

    def getSfuncFixedPerpSpacing(self, N, contour,
            surface_direction, lower):
        """
        Return a function s(i) giving poloidal distance as a function of index-number.
        Construct so that ds_perp/diN = d_lower at the lower end or ds_perp/diN = d_upper
        at the upper end, where s_perp is the distance normal to the vector
        'surface_direction'.
        """
        N_norm = self.poloidalSpacingParameters.N_norm
        d_lower = self.poloidalSpacingParameters.polynomial_d_lower
        d_upper = self.poloidalSpacingParameters.polynomial_d_upper

        s_of_sperp, s_perp_total = contour.interpSSperp(surface_direction)
        sperp_func = self.getPolynomialPoloidalDistanceFunc(s_perp_total, N - 1, N_norm,
                d_lower=d_lower, d_upper=d_upper)
        return lambda i: s_of_sperp(sperp_func(i)), sperp_func

    def getPolynomialPoloidalDistanceFunc(self, length, N, N_norm, *, d_lower=None,
            d_upper=None):
        """
        Return a function s(i) giving poloidal distance as a function of index-number.
        Construct s(i)=sN(iN) as a function of the normalized iN = i/N_norm so that it has the
        same form when resolution is changed. The total Ny in the grid might be a good
        choice for N_norm.
        sN(0) = 0
        sN(N/N_norm) = L
        ds/diN(0) = d_lower at iN=0
        d2s/diN2(0) = 0 if d_lower is not None
        sN(iN) = d_lower*iN for iN < 0 if d_lower is given
        ds/diN(N/N_norm) = d_upper at iN=N_norm
        d2s/diN2(N/N_norm) = 0 if d_upper is not None
        sN(iN) = L + d_upper*(iN - N/N_norm) for iN > N/N_norm if d_upper is given
        """
        if d_lower is None and d_upper is None:
            # always monotonic
            return lambda i: i*length/N
        elif d_lower is None:
            # s(iN) = a + b*iN + c*iN^2 + d*iN^3
            # s(0) = 0 = a
            # d2s/diN2(N/N_norm) = 0 = 2*c + 6*d*N/N_norm
            # c = -3*d*N/N_norm
            # ds/diN(N/N_norm) = d_upper = b + 2*c*N/N_norm + 3*d*(N/N_norm)^2
            #                            = b - 3*d*(N/N_norm)^2
            # b = d_upper + 3*d*(N/N_norm)^2
            # s(N/N_norm) = L = a + b*N/N_norm + c*(N/N_norm)^2 + d*(N/N_norm)^3
            #                 = 0 + d_upper*N/N_Norm + 3*d*(N/N_norm)^3 - 3*d*(N/N_norm)^3 + d*(N/N_norm)^3
            #               d = L*(N_norm/N)^3 - d_upper*(N_norm/N)^2
            d = length*(N_norm/N)**3 - d_upper*(N_norm/N)**2
            b = d_upper + 3*d*(N/N_norm)**2
            c = -3.*d*N/N_norm

            # check function is monotonic: gradients at beginning and end should both be
            # positive.
            # lower boundary:
            assert b >= 0., 'gradient at start should be positive'
            # upper boundary:
            assert b + 2.*c*N/N_norm + 3.*d*(N/N_norm)**2 >= 0., 'gradient at end should be positive'

            return lambda i: numpy.piecewise(i, i > N,
                    [lambda i: length + d_upper*(i - N)/N_norm,
                     lambda i: b*i/N_norm + c*(i/N_norm)**2 + d*(i/N_norm)**3])
        elif d_upper is None:
            # s(iN) = a + b*iN + c*iN^2 + d*iN^3
            # s(0) = 0 = a
            # ds/diN(0) = d_lower = b
            # d2s/diN2(0) = 0 = 2.*c
            # s(N/N_norm) = L = a + b*N/N_norm + c*(N/N_norm)^2 + d*(N/N_norm)^3
            #                 = 0 + d_lower*N/N_norm + 0 + d*(N/N_norm)^3
            # d = L*(N_norm/N)^3 - d_lower*(N_norm/N)^2
            d = length*(N_norm/N)**3 - d_lower*(N_norm/N)**2
            b = d_lower

            # check function is monotonic: gradients at beginning and end should both be
            # positive.
            # lower boundary:
            assert b >= 0., 'gradient at start should be positive'
            # upper boundary:
            assert b + 3.*d*(N/N_norm)**2 > 0., 'gradient at end should be positive'

            return lambda i: numpy.piecewise(i, i < 0.,
                    [lambda i: d_lower*i/N_norm, lambda i: b*i/N_norm + d*(i/N_norm)**3])
        else:
            # s(iN) = a + b*iN + c*iN^2 + d*iN^3 + e*iN^4 + f*iN^5
            # s(0) = 0 = a
            # ds/diN(0) = d_lower = b
            # d2s/diN2(0) = 0 = 2.*c
            # s(N/N_norm) = L = a + b*N/N_norm + c*(N/N_norm)^2 + d*(N/N_norm)^3 + e*(N/N_norm)^4 + f*(N/N_norm)^5
            #                 = 0 + d_lower*N/N_norm + 0 + d*(N/N_norm)^3 + e*(N/N_norm)^4 + f*(N/N_norm)^5
            # ds/diN(N/N_norm) = d_upper = b + 2*c*(N/N_norm) + 3*d(N/N_norm)^2 + 4*e*(N/N_norm)^3 + 5*f*(N/N_norm)^4
            #                            = d_lower + 0 + 3*d(N/N_norm)^2 + 4*e*(N/N_norm)^3 + 5*f*(N/N_norm)^4
            # d2s/diN2(N/N_norm) = 0 = 2*c + 6*d*(N/N_norm) + 12*e*(N/N_norm)^2 + 20*f*(N/N_norm)^3
            #                        = 0 + 6*d*(N/N_norm) + 12*e*(N/N_norm)^2 + 20*f*(N/N_norm)^3
            # 0 = 3*d + 6*e*N/N_norm + 10*f*(N/N_norm)^2
            # d_upper - d_lower = -2*e*(N/N_norm)^3 - 5*f*(N/N_norm)^4
            # 3*L = 3*d_lower*N/N_norm - 3*e*(N/N_norm)^4 - 7*f*(N/N_norm)^5
            # 6*L - 3*d_upper*N/N_norm + 3*d_lower*N/N_norm = 6*d_lower*N/N_norm + f*(N/N_norm)^5
            # f = 6*L*(N_norm/N)^5 - 3*(d_upper + d_lower)*(N_norm/N)^4
            f = 6.*length*(N_norm/N)**5 - 3*(d_upper + d_lower)*(N_norm/N)**4
            e = -1./2.*(d_upper - d_lower)*(N_norm/N)**3 - 5./2.*f*N/N_norm
            d = -2.*e*N/N_norm - 10./3.*f*(N/N_norm)**2
            b = d_lower

            # check function is monotonic: gradients at beginning and end should both be
            # positive. Only check the boundaries here, should really add a check that
            # gradient does not reverse in the middle somewhere...
            # lower boundary:
            assert b >= 0., 'gradient at start should be positive'
            # upper boundary:
            assert b + 3.*d*(N/N_norm)**2 + 4.*e*(N/N_norm)**3 + 5.*f*(N/N_norm)**4 >= 0., 'gradient at end should be positive'

            return lambda i: numpy.piecewise(i, [i < 0., i > N],
                    [lambda i: d_lower*i/N_norm,
                     lambda i:length + d_upper*(i - N)/N_norm,
                     lambda i: b*i/N_norm + d*(i/N_norm)**3 + e*(i/N_norm)**4 + f*(i/N_norm)**5])

    def getSqrtPoloidalDistanceFunc(self, length, N, N_norm, *, b_lower=None, a_lower=None,
            b_upper=None, a_upper=None):
        """
        Return a function s(i) giving poloidal distance as a function of index-number.
        Construct s(i)=sN(iN) as a function of the normalized iN = i/N_norm so that it has the
        same form when resolution is changed. The total Ny in the grid might be a good
        choice for N_norm.
        sN(0) = 0
        sN(N/N_norm) = L
        ds/diN(0) ~ a_lower/sqrt(iN)+b_lower at iN=0 (if a_lower not None, else no
                                                           sqrt(iN) term)
        ds/diN(N/N_norm) ~ a_upper/sqrt(N/N_norm-iN)+b_upper at iN=N_norm (if a_upper is
                                               not None, else no sqrt(N/N_norm - iN) term)

        By default a_lower=b_lower and a_upper=b_upper, unless both are
        specified explicitly
        """
        if b_lower is None and b_upper is None:
            assert a_lower is None, 'cannot set a_lower unless b_lower is set'
            assert a_upper is None, 'cannot set a_upper unless b_upper is set'
            # always monotonic
            return lambda i: i*length/N
        elif b_lower is None:
            assert a_lower is None, 'cannot set a_lower unless b_lower is set'
            if a_upper is None:
                a_upper = b_upper
            # s(iN) = -b*sqrt(N/N_norm-iN) + c + d*iN + e*(iN)^2
            # s(0) = 0 = -b*sqrt(N/N_norm) + c
            # ds/diN(N/N_norm) = b/(2*sqrt(N/N_norm-iN))+d+2*e*N/N_norm ~ a_upper/sqrt(N/N_norm-iN)+b_upper
            # b = 2*a_upper
            # d + 2*e*N/N_norm = b_upper
            # d = b_upper - 2*e*N/N_norm
            # s(N/N_norm) = L = c + d*N/N_norm + e*(N/N_norm)^2
            # L = c + b_upper*N/N_norm - 2*e*(N/N_norm)^2 + e*(N/N_norm)^2
            # e = (c + b_upper*N/N_norm - L) / (N/N_norm)^2
            b = 2.*a_upper
            c = b*numpy.sqrt(N/N_norm)
            e = (c + b_upper*N/N_norm - length) / (N/N_norm)**2
            d = b_upper - 2*e*N/N_norm

            # check function is monotonic: gradients at beginning and end should both be
            # positive.
            # lower boundary:
            assert b/(2.*numpy.sqrt(N/N_norm)) + d > 0., 'gradient at start should be positive'
            # upper boundary:
            assert b >= 0., 'sqrt part of function should be positive at end'
            assert d + 2.*e*N/N_norm >= 0., 'gradient of polynomial part should be positive at end'

            return lambda i: -b*numpy.sqrt((N-i)/N_norm) + c + d*i/N_norm + e*(i/N_norm)**2
        elif b_upper is None:
            if a_lower is None:
                a_lower = b_lower
            assert a_upper is None
            # s(iN) = a*sqrt(iN) + c + d*iN + e*iN^2
            # s(0) = 0 = c
            # ds/di(0) = a/(2*sqrt(iN))+d ~ a_lower/sqrt(iN)+b_lower
            # a = 2*a_lower
            # d = b_lower
            # s(N/N_norm) = L = a*sqrt(N/N_norm) + c + d*N/N_norm + e*(N/N_norm)^2
            a = 2.*a_lower
            d = b_lower
            e = (length - a*numpy.sqrt(N/N_norm) - d*N/N_norm)/(N/N_norm)**2

            # check function is monotonic: gradients at beginning and end should both be
            # positive.
            # lower boundary:
            assert a >= 0., 'sqrt part of function should be positive at start'
            assert d >= 0., 'gradient of polynomial part should be positive at start'
            # upper boundary:
            assert a/(2.*numpy.sqrt(N/N_norm)) + d + 2.*e*N/N_norm > 0., 'gradient at end should be positive'

            return lambda i: a*numpy.sqrt(i/N_norm) + d*i/N_norm + e*(i/N_norm)**2
        else:
            if a_lower is None:
                a_lower = b_lower
            if a_upper is None:
                a_upper = b_upper
            # s(iN) = a*sqrt(iN) - b*sqrt(N/N_norm-iN) + c + d*iN + e*iN^2 + f*iN^3
            # s(0) = 0 = -b*sqrt(N/N_norm) + c
            # c = b*sqrt(N/N_norm)
            # ds/diN(0) = a/(2*sqrt(iN))+b/(2*sqrt(N/N_norm))+d ~ a_lower/sqrt(iN)+b_lower
            # a = 2*a_lower
            # b/(2*sqrt(N/N_norm)) + d = b_lower
            # d = b_lower - b/(2*sqrt(N/N_norm)
            # ds/di(N) = b/(2*sqrt(N/N_norm-iN))+a/(2*sqrt(N))+d+2*e*N/N_norm+3*f*(N/N_norm)^2 ~ a_upper/sqrt(N/N_norm-i)+b_upper
            # b = 2*a_upper
            # a/(2*sqrt(N/N_norm) + d + 2*e*N/N_norm + 3*f*(N/N_norm)^2 = b_upper
            # e = (b_upper - a/(2*sqrt(N/N_norm)) - d)/(2*N/N_norm) - 3/2*f*N/N_norm
            # s(N/N_norm) = L = a*sqrt(N/N_norm) + c + d*N/N_norm + e*(N/N_norm)^2 + f*(N/N_norm)^3
            # L = a*sqrt(N/N_norm) + c + d*N/N_norm + (b_upper - a/(2*sqrt(N/N_norm)) - d)*N/(2*N_norm) - 3/2*f*(N/N_norm)^3 + f*(N/N_norm)^3
            # f = 2*(a*sqrt(N/N_norm) + c + d*N/(2*N_norm) + b_upper*N/(2*N_norm) - a*sqrt(N/N_norm)/4 - L)*N_norm^3/N^3
            a = 2.*a_lower
            b = 2.*a_upper
            c = b*numpy.sqrt(N/N_norm)
            d = b_lower - b/2./numpy.sqrt(N/N_norm)
            f = 2.*(a*numpy.sqrt(N/N_norm) + c + d*N/N_norm/2. + b_upper*N/N_norm/2. - a*numpy.sqrt(N/N_norm)/4. - length)*N_norm**3/N**3
            e = (b_upper - a/2./numpy.sqrt(N/N_norm) - d)*N_norm/2./N - 1.5*f*N/N_norm

            # check function is monotonic: gradients at beginning and end should both be
            # positive. Only check the boundaries here, should really add a check that
            # gradient does not reverse in the middle somewhere...
            # lower boundary:
            assert a >= 0., 'sqrt part of function should be positive at start'
            assert b/(2.*numpy.sqrt(N/N_norm)) + d >= 0., 'gradient of non-singular part should be positive at start'
            # upper boundary:
            assert b >= 0., 'sqrt part of function should be positive at end'
            assert a/(2.*numpy.sqrt(N/N_norm)) + d + 2.*e*N/N_norm + 3.*f*(N/N_norm)**2 >= 0., 'gradient of non-singular part should be positive at end'

            return lambda i: a*numpy.sqrt(i/N_norm) - b*numpy.sqrt((N-i)/N_norm) + c + d*i/N_norm + e*(i/N_norm)**2 + f*(i/N_norm)**3

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
      - self.wall: list of Point2D giving vertices of polygon representing the wall, in
        anti-clockwise order; assumed to be closed so last element and first are taken to
        be connected
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

    def wallPosition(self, s):
        """
        Get a position on the wall, where the distance along the wall is parameterized by
        0<=s<1
        """
        try:
            return Point2D(self.wallRInterp(s), self.wallZInterp(s))
        except AttributeError:
            # wall interpolation functions not created yet

            wall = deepcopy(self.wall)

            # make closed contour
            wall.append(wall[0])

            R = [p.R for p in wall]
            Z = [p.Z for p in wall]

            wallfraction = numpy.linspace(0., 1., len(wall))

            self.wallRInterp = interp1d(wallfraction, R, kind='linear', assume_sorted=True)
            self.wallZInterp = interp1d(wallfraction, Z, kind='linear', assume_sorted=True)

            return Point2D(self.wallRInterp(s), self.wallZInterp(s))

    def wallVector(self, s):
        """
        Get the vector along the wall at a point s, with the same parameterization as
        wallPosition.
        """
        try:
            return numpy.array([self.wallVectorRComponent(s), self.wallVectorZComponent(s)])
        except AttributeError:
            # wall vector interpolation functions not created yet
            Rcomponents = [self.wall[i+1].R - self.wall[i].R for i in range(len(self.wall)-1)]
            Rcomponents.append(self.wall[0].R - self.wall[-1].R)
            Rcomponents.append(self.wall[1].R - self.wall[0].R)

            Zcomponents = [self.wall[i+1].Z - self.wall[i].Z for i in range(len(self.wall)-1)]
            Zcomponents.append(self.wall[0].Z - self.wall[-1].Z)
            Zcomponents.append(self.wall[1].Z - self.wall[0].Z)

            wallfraction = numpy.linspace(0., 1., len(self.wall) + 1)

            # Vector along wall stays constant along each segment, as we assume the
            # segments are straight. Have calculated the vector at each vertex for the
            # following segment, so use 'previous' interpolation to just take the value
            # from the previous point
            self.wallVectorRComponent = interp1d(wallfraction, Rcomponents, kind='previous', assume_sorted=True)
            self.wallVectorZComponent = interp1d(wallfraction, Zcomponents, kind='previous', assume_sorted=True)

            return numpy.array([self.wallVectorRComponent(s), self.wallVectorZComponent(s)])

    def wallIntersection(self, p1, p2):
        """
        Find the intersection, if any, between the wall and the line between p1 and p2
        """
        intersect = None
        for i in range(len(self.wall) - 1):
            if intersect is None:
                intersect = find_intersection(self.wall[i], self.wall[i+1], p1, p2)
            else:
                second_intersect = find_intersection(self.wall[i], self.wall[i+1], p1, p2)
                if second_intersect is not None:
                    assert numpy.abs(intersect.R - second_intersect.R) < intersect_tolerance and numpy.abs(intersect.Z - second_intersect.Z) < intersect_tolerance, 'Multiple intersections with wall found'
        # final segment between last and first point of wall
        if intersect is None:
            intersect = find_intersection(self.wall[-1], self.wall[0], p1, p2)
        else:
            second_intersect = find_intersection(self.wall[-1], self.wall[0], p1, p2)
            if second_intersect is not None:
                assert numpy.abs(intersect.R - second_intersect.R) < 1.e-14 and numpy.abs(intersect.Z - second_intersect.Z) < 1.e-14, 'Multiple intersections with wall found'

        return intersect

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
