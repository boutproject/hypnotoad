"""
Functions for analysing an equilibrium for which an interpolating function is given for
the potential.
"""

from collections import OrderedDict
from copy import deepcopy
import warnings

import numpy
from options import Options
from scipy.optimize import minimize_scalar, brentq, root
from scipy.interpolate import interp1d

from .hypnotoad_options import HypnotoadInternalOptions, HypnotoadOptions

class SolutionError(Exception):
    """
    Solution was not found
    """
    pass

# tolerance used to try and avoid missed intersections between lines
# also if two sets of lines appear to intersect twice, only count it once if the
# distance between the intersections is less than this
intersect_tolerance = 1.e-14

def setDefault(options, name, default):
    if options[name] is None:
        options[name] = default
    return options[name]

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

def find_intersections(l1array, l2start, l2end):
    """
    Find the intersection (if there is one) between the array of lines 'l1' and the line
    'l2'.
    """
    # Copy so we don't change the inputs
    l1array = l1array.copy()
    l2start = deepcopy(l2start)
    l2end = deepcopy(l2end)

    R1array = numpy.zeros([l1array.shape[0] - 1, 2])
    R1array[:, 0] = l1array[:-1, 0]
    R1array[:, 1] = l1array[1:, 0]
    Z1array = numpy.zeros([l1array.shape[0] - 1, 2])
    Z1array[:, 0] = l1array[:-1, 1]
    Z1array[:, 1] = l1array[1:, 1]

    # for inds1, if l1 is sensible, dR1 shouldn't be too small as it's bigger than dZ1
    # l1 is Z = Z1 + dZ1/dR1 * (R - R1)
    # If the lines are parallel
    inds_a = numpy.where(
                numpy.abs(R1array[:, 0] - R1array[:, 1])
                > numpy.abs(Z1array[:, 0] - Z1array[:, 1])
            )[0]
    thisR1_a = R1array[inds_a, :]
    thisZ1_a = Z1array[inds_a, :]

    # sort inds_a points in R
    sortinds = numpy.argsort(thisR1_a, axis=1)
    na = thisR1_a.shape[0]
    thisR1_a = thisR1_a[numpy.arange(na)[:, numpy.newaxis], sortinds]
    thisZ1_a = thisZ1_a[numpy.arange(na)[:, numpy.newaxis], sortinds]

    # if l2 is sensible, dZ2 shouldn't be too small as it's bigger than dR2
    # l2 is R = R2 + dR2/dZ2 * (Z - Z2)
    inds_b = numpy.where(
                numpy.abs(R1array[:, 0] - R1array[:, 1])
                <= numpy.abs(Z1array[:, 0] - Z1array[:, 1])
            )[0]

    thisR1_b = R1array[inds_b, :]
    thisZ1_b = Z1array[inds_b, :]

    # sort inds_b points in Z
    sortinds = numpy.argsort(thisZ1_b, axis=1)
    nb = thisR1_b.shape[0]
    thisR1_b = thisR1_b[numpy.arange(nb)[:, numpy.newaxis], sortinds]
    thisZ1_b = thisZ1_b[numpy.arange(nb)[:, numpy.newaxis], sortinds]

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

        # Check intersections with 'a' lines
        #
        # If this condition is not true, lines are parallel so cannot intersect
        condition = numpy.where(
                numpy.abs((thisZ1_a[:,1]-thisZ1_a[:,0])/(thisR1_a[:,1]-thisR1_a[:,0]) - dZ2/dR2) >= 1.e-15
                )
        inds_a = inds_a[condition]
        thisR1_a = thisR1_a[condition]
        thisZ1_a = thisZ1_a[condition]

        thisdR1 = thisR1_a[:, 1] - thisR1_a[:, 0]
        thisdZ1 = thisZ1_a[:, 1] - thisZ1_a[:, 0]

        # intersection where
        # Z1 + dZ1/dR1 * (R - R1) = Z2 + dZ2/dR2 * (R - R2)
        # (dZ1/dR1 - dZ2/dR2)*R = Z2 - Z1 + dZ1/dR1*R1 - dZ2/dR2*R2
        Rcross = (Z2 - thisZ1_a[:, 0] + thisdZ1/thisdR1*thisR1_a[:, 0] - dZ2/dR2*R2) / (thisdZ1/thisdR1 - dZ2/dR2)
        intersect_inds = numpy.where(numpy.logical_and(
            Rcross >= thisR1_a[:, 0] - intersect_tolerance,
            numpy.logical_and(Rcross <= thisR1_a[:, 1] + intersect_tolerance,
            numpy.logical_and(Rcross >= R2 - intersect_tolerance,
            Rcross <= l2end.R + intersect_tolerance
        ))))
        Rintersect_a = Rcross[intersect_inds]
        Zintersect_a = thisZ1_a[:, 0][intersect_inds] + thisdZ1[intersect_inds]/thisdR1[intersect_inds] * (Rintersect_a - thisR1_a[:, 0][intersect_inds])


        # Check intersections with 'b' lines
        #
        thisdR1 = thisR1_b[:, 1] - thisR1_b[:, 0]
        thisdZ1 = thisZ1_b[:, 1] - thisZ1_b[:, 0]

        # intersection where
        # R = R1 + dR1/dZ1 * (Z2 + dZ2/dR2 * (R - R2) - Z1)
        # (1 - dR1/dZ1*dZ2/dR2) * R = R1 + dR1/dZ1 * (Z2 - dZ2/dR2*R2 - Z1)
        Rcross = (thisR1_b[:, 0] + thisdR1/thisdZ1 * (Z2 - dZ2/dR2*R2 - thisZ1_b[:, 0])) / (1. - thisdR1/thisdZ1 * dZ2/dR2)
        Zcross = Z2 + dZ2/dR2 * (Rcross - R2)
        intersect_inds = numpy.where(numpy.logical_and(
            Zcross >= thisZ1_b[:, 0] - intersect_tolerance,
            numpy.logical_and(Zcross <= thisZ1_b[:, 1] + intersect_tolerance,
            numpy.logical_and(Rcross >= R2 - intersect_tolerance,
            Rcross <= l2end.R + intersect_tolerance
        ))))
        Rintersect_b = Rcross[intersect_inds]
        Zintersect_b = Zcross[intersect_inds]
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

        # Check intersections with 'a' lines
        #
        thisdR1 = thisR1_a[:, 1] - thisR1_a[:, 0]
        thisdZ1 = thisZ1_a[:, 1] - thisZ1_a[:, 0]

        # intersection where
        # Z = Z1 + dZ1/dR1 * (R2 + dR2/dZ2 * (Z - Z2) - R1)
        # (1 - dZ1*dR2/dR1/dZ2) * Z = Z1 + dZ1/dR1 * (R2 - dR2/dZ2*Z2 - R1)
        Zcross = (thisZ1_a[:, 0] + thisdZ1/thisdR1 * (R2 - dR2/dZ2*Z2 - thisR1_a[:, 0])) / (1. - thisdZ1*dR2/(thisdR1*dZ2))
        Rcross = R2 + dR2/dZ2 * (Zcross - Z2)
        intersect_inds = numpy.where(numpy.logical_and(
            Rcross >= thisR1_a[:, 0] - intersect_tolerance,
            numpy.logical_and(Rcross <= thisR1_a[:, 1] + intersect_tolerance,
            numpy.logical_and(Zcross >= Z2 - intersect_tolerance,
            Zcross <= l2end.Z + intersect_tolerance
        ))))
        Rintersect_a = Rcross[intersect_inds]
        Zintersect_a = Zcross[intersect_inds]

        # Check intersections with 'b' lines
        #
        # If this condition is not true, lines are parallel so cannot intersect
        condition = numpy.where(
                numpy.abs(dR2/dZ2 - (thisR1_b[:, 1] - thisR1_b[:, 0])/(thisZ1_b[:, 1] - thisZ1_b[:, 0])) >= 1.e-15
                )
        inds_b = inds_b[condition]
        thisR1_b = thisR1_b[condition]
        thisZ1_b = thisZ1_b[condition]

        thisdR1 = thisR1_b[:, 1] - thisR1_b[:, 0]
        thisdZ1 = thisZ1_b[:, 1] - thisZ1_b[:, 0]

        # intersection where
        # R2 + dR2/dZ2 * (Z - Z2) = R1 + dR1/dZ1 * (Z - Z1)
        # (dR2/dZ2 - dR1*dZ1) * Z = R1 - R2 + dR2/dZ2*Z2 - dR1/dZ1*Z1
        Zcross = (thisR1_b[:, 0] - R2 + dR2/dZ2*Z2 - thisdR1/thisdZ1*thisZ1_b[:, 0]) / (dR2/dZ2 - thisdR1/thisdZ1)
        intersect_inds = numpy.where(numpy.logical_and(
            Zcross >= thisZ1_b[:, 0] - intersect_tolerance,
            numpy.logical_and(Zcross <= thisZ1_b[:, 1] + intersect_tolerance,
            numpy.logical_and(Zcross >= Z2 - intersect_tolerance,
            Zcross <= l2end.Z + intersect_tolerance
        ))))
        Zintersect_b = Zcross[intersect_inds]
        Rintersect_b = R2 + dR2/dZ2 * (Zintersect_b - Z2)

    Rintersect = numpy.concatenate([Rintersect_a, Rintersect_b])
    Zintersect = numpy.concatenate([Zintersect_a, Zintersect_b])

    if len(Rintersect) > 0 or len(Zintersect) > 0:
        return numpy.stack([Rintersect, Zintersect], axis=1)
    else:
        return None

class FineContour:
    """
    Used to give a high-resolution representation of a contour.
    Points in FineContour are uniformly spaced in poloidal distance along the contour.
    """

    options = Options(
            finecontour_Nfine = None,
            finecontour_atol = None,
            finecontour_diagnose = None,
            )

    def __init__(self, parentContour):
        self.parentContour = parentContour
        self.distance = None
        Nfine = self.options.finecontour_Nfine
        atol = self.options.finecontour_atol

        endInd = self.parentContour.endInd
        if endInd < 0:
            # endInd might be negative, which would mean relative to the end of the list,
            # but we need the actual index below
            endInd += len(self.parentContour)
        n_input = (endInd - self.parentContour.startInd + 1)

        # Extend further than will be needed in the final contour, because extrapolation
        # past the end of the fine contour is very bad.
        extend_lower_fine = 2*(self.parentContour.extend_lower * Nfine) // n_input
        extend_upper_fine = 2*(self.parentContour.extend_upper * Nfine) // n_input

        indices_fine = numpy.linspace(-extend_lower_fine,
                (Nfine - 1 + extend_upper_fine),
                Nfine + extend_lower_fine + extend_upper_fine)

        # Initial guess from interpolation of psiContour, iterate to a more accurate
        # version below.
        # Extend a copy of parentContour to make the extrapolation more stable.
        # This makes parentCopy have twice the extra points as parentContour has.
        parentCopy = parentContour.newContourFromSelf()
        parentCopy.temporaryExtend(extend_lower=parentContour.extend_lower,
                extend_upper=parentContour.extend_upper,
                ds_lower=calc_distance(parentCopy[0], parentCopy[1]),
                ds_upper=calc_distance(parentCopy[-1], parentCopy[-2]))
        interp_input, distance_estimate = parentCopy._coarseInterp()

        sfine = distance_estimate[parentCopy.endInd] / (Nfine - 1) * indices_fine

        # 2d array with size {N,2} giving the (R,Z)-positions of points on the contour
        self.positions = numpy.array(tuple(interp_input(s).as_ndarray() for s in sfine))

        self.startInd = extend_lower_fine
        self.endInd = Nfine - 1 + extend_lower_fine

        self.refine()

        self.calcDistance()

        ds = self.distance[1:] - self.distance[:-1]
        # want constant spacing, so ds has a constant value
        ds_mean = numpy.mean(ds)
        # maximum error
        ds_error = numpy.max(numpy.sqrt((ds - ds_mean)**2))

        if FineContour.options.finecontour_diagnose:
            from matplotlib import pyplot
            print('diagnosing FineContour.__init__()')
            print('extend_lower_fine', extend_lower_fine)
            print('extend_upper_fine', extend_upper_fine)
            print('ds_error', ds_error)
            count = 1

            Rpoints = self.positions[:, 0]
            Zpoints = self.positions[:, 1]
            R = numpy.linspace(Rpoints.min(), Rpoints.max(), 100)
            Z = numpy.linspace(Zpoints.min(), Zpoints.max(), 100)

            pyplot.figure()

            pyplot.subplot(131)
            pyplot.contour(R, Z, self.parentContour.psi(R[numpy.newaxis, :], Z[:, numpy.newaxis]))
            pyplot.plot(Rpoints, Zpoints, marker='x')
            pyplot.xlabel('R')
            pyplot.ylabel('Z')

            pyplot.subplot(132)
            pyplot.plot(ds)
            pyplot.ylabel('ds')

            pyplot.subplot(133)
            pyplot.plot(Rpoints, label='R')
            pyplot.plot(Zpoints, label='Z')
            pyplot.xlabel('index')
            pyplot.legend()
            pyplot.show()

        while ds_error > atol:
            sfine = self.totalDistance() / (Nfine - 1) * indices_fine

            interpFunc = self.interpFunction()

            # 2d array with size {N,2} giving the (R,Z)-positions of points on the contour
            self.positions = numpy.array(tuple(interpFunc(s).as_ndarray() for s in sfine))

            self.refine()

            self.calcDistance()

            ds = self.distance[1:] - self.distance[:-1]
            # want constant spacing, so ds has a constant value
            ds_mean = numpy.mean(ds)
            # maximum error
            ds_error = numpy.max(numpy.sqrt((ds - ds_mean)**2))

            if FineContour.options.finecontour_diagnose:
                print('iteration', count, '  ds_error', ds_error)
                count += 1

                Rpoints = self.positions[:, 0]
                Zpoints = self.positions[:, 1]
                R = numpy.linspace(Rpoints.min(), Rpoints.max(), 100)
                Z = numpy.linspace(Zpoints.min(), Zpoints.max(), 100)

                pyplot.figure()

                pyplot.subplot(131)
                pyplot.contour(R, Z, self.parentContour.psi(R[numpy.newaxis, :], Z[:, numpy.newaxis]))
                pyplot.plot(Rpoints, Zpoints, marker='x')
                pyplot.xlabel('R')
                pyplot.ylabel('Z')

                pyplot.subplot(132)
                pyplot.plot(ds)
                pyplot.ylabel('ds')

                pyplot.subplot(133)
                pyplot.plot(Rpoints, label='R')
                pyplot.plot(Zpoints, label='Z')
                pyplot.xlabel('index')
                pyplot.legend()
                pyplot.show()

    def totalDistance(self):
        return self.distance[self.endInd] - self.distance[self.startInd]

    def calcDistance(self):
        if self.distance is None:
            self.distance = numpy.zeros(self.positions.shape[0])
        deltaSquared = (self.positions[1:] - self.positions[:-1])**2
        self.distance[1:] = numpy.cumsum(numpy.sqrt(numpy.sum(deltaSquared, axis=1)))

    def interpFunction(self, *, kind='cubic'):
        distance = self.distance - self.distance[self.startInd]

        interpR = interp1d(distance, self.positions[:,0], kind=kind,
                           assume_sorted=True, fill_value='extrapolate')
        interpZ = interp1d(distance, self.positions[:,1], kind=kind,
                           assume_sorted=True, fill_value='extrapolate')
        return lambda s: Point2D(float(interpR(s)), float(interpZ(s)))

    def refine(self):
        result = numpy.zeros(self.positions.shape)

        p = self.positions[0, :]
        tangent = self.positions[1, :] - self.positions[0, :]
        result[0, :] = self.parentContour.refinePoint(Point2D(*p),
                Point2D(*tangent)).as_ndarray()
        for i in range(1, self.positions.shape[0] - 1):
            p = self.positions[i, :]
            tangent = self.positions[i+1, :] - self.positions[i-1, :]
            result[i, :] = self.parentContour.refinePoint(Point2D(*p),
                    Point2D(*tangent)).as_ndarray()
        p = self.positions[-1, :]
        tangent = self.positions[-1, :] - self.positions[-2, :]
        result[-1, :] = self.parentContour.refinePoint(Point2D(*p),
                Point2D(*tangent)).as_ndarray()

        self.positions = result

    def reverse(self):
        if self.distance is not None:
            self.distance = self.distance[-1] - self.distance[::-1]
        self.positions = self.positions[::-1, :]

        old_start = self.startInd
        n = self.positions.shape[0]
        self.startInd = n - 1 - self.endInd
        self.endInd = n - 1 - old_start

    def interpSSperp(self, vec, kind='cubic'):
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
        start_position = self.positions[self.startInd, :]

        # s_perp = (vec_perp).(r) where r is the displacement vector of each point from self[self.startInd]
        s_perp = numpy.sum((self.positions - start_position)*vec_perp[numpy.newaxis, :],
                axis=1)

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

        distance = self.distance - self.distance[self.startInd]
        s_of_sperp = interp1d(s_perp, distance, kind=kind, assume_sorted=True,
                fill_value = 'extrapolate')

        return s_of_sperp, s_perp_total

    def getDistance(self, p):
        """
        Return the distance of a point along the contour.
        Assume p is a point on the contour so has the correct psi-value.
        """
        p = p.as_ndarray()

        distance_from_points = numpy.sqrt(numpy.sum(
            (self.positions - p[numpy.newaxis, :])**2, axis=1))

        # index of closest point
        i1 = numpy.argmin(distance_from_points)
        d1 = distance_from_points[i1]

        # index of next-closest point
        if i1 + 1 >= len(distance_from_points):
            i2 = i1 - 1
        elif i1 - 1 < 0:
            i2 = 1
        elif distance_from_points[i1+1] < distance_from_points[i1-1]:
            i2 = i1 + 1
        else:
            i2 = i1 - 1
        d2 = distance_from_points[i2]

        # linearly interpolate the distance of the two closest points in the same ratio as
        # their distances from the point
        r = d2 / (d1 + d2)

        return r*self.distance[i1] + (1. - r)*self.distance[i2]

    def plot(self, *args, plotPsi=False, **kwargs):
        from matplotlib import pyplot
        Rpoints = self.positions[:, 0]
        Zpoints = self.positions[:, 1]
        if plotPsi:
            R = numpy.linspace(min(Rpoints), max(Rpoints), 100)
            Z = numpy.linspace(min(Zpoints), max(Zpoints), 100)
            pyplot.contour(R, Z, self.psi(R[numpy.newaxis, :], Z[:, numpy.newaxis]))
        pyplot.plot(Rpoints, Zpoints, *args, **kwargs)

class PsiContour:
    """
    Represents a contour as a collection of points.
    Includes methods for interpolation.
    Mostly behaves like a list
    """
    options = Options(
            refine_width = 1.e-5,
            refine_atol = 2.e-8,
            )
    def __init__(self, points, psi, psival):
        self.points = points

        self._startInd = 0
        self._endInd = len(points) - 1

        self._fine_contour = None

        self._distance = None

        # Function that evaluates the vector potential at R,Z
        self.psi = psi

        # Value of vector potential on this contour
        self.psival = psival

        # Number of boundary guard cells at either end
        # This may be set even if the contour has not been extended yet, to specify how
        # many guard cells should be added when it is - this is extra information to
        # startInd and endInd.
        self._extend_lower = 0
        self._extend_upper = 0

    @property
    def startInd(self):
        return self._startInd

    @startInd.setter
    def startInd(self, val):
        if self._startInd != val:
            # self._fine_contour needs to be recalculated if the start position changes
            self._fine_contour = None
            self._distance = None
            self._startInd = val

    @property
    def endInd(self):
        return self._endInd

    @endInd.setter
    def endInd(self, val):
        if self._endInd != val:
            # self._fine_contour needs to be recalculated if the end position changes
            self._fine_contour = None
            self._distance = None
            self._endInd = val

    @property
    def extend_lower(self):
        return self._extend_lower

    @extend_lower.setter
    def extend_lower(self, val):
        if self._extend_lower != val:
            # self._fine_contour needs to be recalculated if extend_lower changes, to add more
            # points at the lower end
            self._fine_contour = None
            self._extend_lower = val

    @property
    def extend_upper(self):
        return self._extend_upper

    @extend_upper.setter
    def extend_upper(self, val):
        if self._extend_upper != val:
            # self._fine_contour needs to be recalculated if extend_upper changes, to add more
            # points at the upper end
            self._fine_contour = None
            self._extend_upper = val

    @property
    def fine_contour(self):
        if self._fine_contour is None:
            self._fine_contour = FineContour(self)
        return self._fine_contour

    @property
    def distance(self):
        if self._distance is None:
            self._distance = [self.fine_contour.getDistance(p) for p in self]
        return self._distance

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
        self._distance = contour._distance
        self.psi = contour.psi
        self.psival = contour.psival
        self.extend_lower = contour.extend_lower
        self.extend_upper = contour.extend_upper
        self._fine_contour = contour._fine_contour

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
        if points is None:
            new_contour._fine_contour = self._fine_contour

        return new_contour

    def append(self, point):
        self._fine_contour = None
        self._distance = None
        self.points.append(point)

    def prepend(self, point):
        self._fine_contour = None
        self._distance = None
        self.points.insert(0, point)

    def insert(self, index, point):
        self._distance = None

        # Make sure index is positive, following behaviour of list.insert()
        if index < 0:
            index += len(self)
            if index < 0:
                index = 0

        self.points.insert(index, point)

        if index <= self.startInd:
            self.startInd += 1
        if index <= self.endInd:
            self.endInd += 1
        if self.endInd < 0 and index > len(self) + self.endInd:
            self.endInd -= 1

    def totalDistance(self):
        return self.distance[self.endInd] - self.distance[self.startInd]

    def reverse(self):
        self.points.reverse()
        old_start = self.startInd
        self.startInd = len(self) - 1 - self.endInd
        self.endInd = len(self) - 1 - old_start
        if self._distance is not None:
            self._distance.reverse()
            self._distance = [self.distance[0] - d for d in self.distance]
        if self._fine_contour is not None:
            self._fine_contour.reverse()

    def refine(self, *args, **kwargs):
        new = self.getRefined(*args, **kwargs)
        self.points = new.points
        self._distance = new._distance

    def refinePoint(self, p, tangent, width=None, atol=None):
        if width is None:
            width = PsiContour.options.refine_width
        if atol is None:
            atol = PsiContour.options.refine_atol

        f = lambda R,Z: self.psi(R, Z) - self.psival

        if numpy.abs(f(*p)) < atol*numpy.abs(self.psival):
            # don't need to refine
            return p

        def perpLine(w):
            # p - point through which to draw perpLine
            # tangent - vector tangent to original curve, result will be perpendicular to this
            # w - width on either side of p to draw the perpLine to
            modTangent = numpy.sqrt(tangent.R**2 + tangent.Z**2)
            perpIdentityVector = Point2D(tangent.Z/modTangent, -tangent.R/modTangent)
            return lambda s: p + 2.*(s-0.5)*w*perpIdentityVector

        converged = False
        w = width
        sp = []
        ep = []
        while not converged:
            try:
                pline = perpLine(w)
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
                pline0 = perpLine(width)
                Rbox = numpy.linspace(p.R-.1,p.R+.1,100)[numpy.newaxis,:]
                Zbox = numpy.linspace(p.Z-.1,p.Z+.1,100)[:,numpy.newaxis]
                svals = numpy.linspace(0., 1., 40)
                pyplot.figure()
                self.plot('+')
                pyplot.contour(Rbox+0.*Zbox,Zbox+0.*Rbox,self.psi(Rbox,Zbox), 200)
                pyplot.plot([pline0(s).R for s in svals], [pline0(s).Z for s in svals], 'x')
                pyplot.figure()
                pyplot.plot([f(*pline0(s)) for s in svals])
                pyplot.show()
                raise SolutionError("Could not find interval to refine point at "+str(p))

        return pline(snew)

    def getRefined(self, **kwargs):
        newpoints = []
        newpoints.append(self.refinePoint(self.points[0], self.points[1] - self.points[0],
            **kwargs))
        for i,p in enumerate(self.points[1:-1]):
            # note i+1 here is the index of point p
            newpoints.append(self.refinePoint(p, self.points[i+2] - self.points[i],
                **kwargs))
        newpoints.append(self.refinePoint(self.points[-1], self.points[-1] -
            self.points[-2], **kwargs))

        return self.newContourFromSelf(points=newpoints)

    def interpFunction(self):
        return self.fine_contour.interpFunction()

    def _coarseInterp(self, *, kind='cubic'):
        distance = [0.]
        for i in range(len(self) - 1):
            distance.append(distance[i] + calc_distance(self[i+1], self[i]))
        distance = numpy.array(numpy.float64(distance)) - distance[self.startInd]

        R = numpy.array(numpy.float64([p.R for p in self.points]))
        Z = numpy.array(numpy.float64([p.Z for p in self.points]))

        interpR = interp1d(distance, R, kind=kind,
                           assume_sorted=True, fill_value='extrapolate')
        interpZ = interp1d(distance, Z, kind=kind,
                           assume_sorted=True, fill_value='extrapolate')
        return lambda s: Point2D(interpR(s), interpZ(s)), distance

    def contourSfunc(self, kind='cubic'):
        """
        Function interpolating distance as a function of index for the current state of
        this contour. When outside [startInd, endInd], set to constant so the results
        aren't affected by extrapolation errors.
        """
        interpS = interp1d(numpy.arange(len(self), dtype=float),
                self.distance, kind=kind, assume_sorted=True, fill_value='extrapolate')
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
        return self.fine_contour.interpSSperp(vec)

    def regrid(self, *args, **kwargs):
        """
        Regrid this contour, modifying the object
        """
        self.setSelfToContour(self.getRegridded(*args, **kwargs))
        return self

    def getRegridded(self, npoints, *, width=None, atol=None, sfunc=None,
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
        if width is None:
            width = PsiContour.options.refine_width
        if atol is None:
            atol = PsiContour.options.refine_atol

        if extend_lower is not None:
            self.extend_lower = extend_lower
        if extend_upper is not None:
            self.extend_upper = extend_upper
        self.temporaryExtend(extend_lower=self.extend_lower,
                extend_upper=self.extend_upper, ds_lower=calc_distance(self[1], self[0]),
                ds_upper=calc_distance(self[-2], self[-1]))

        indices = numpy.linspace(-self.extend_lower, (npoints - 1 + self.extend_upper),
                npoints + self.extend_lower + self.extend_upper)
        if sfunc is not None:
            s = sfunc(indices)

            # offset fine_contour.interpFunction in case sfunc(0.)!=0.
            sbegin = sfunc(0.)
        else:
            s = (self.distance[self.endInd] - self.distance[self.startInd]) / (npoints - 1) * indices
            sbegin = 0.

        interp_unadjusted = self.fine_contour.interpFunction()
        interp = lambda s: interp_unadjusted(s - sbegin)

        new_contour = self.newContourFromSelf(points=[interp(x) for x in s])
        new_contour.startInd = self.extend_lower
        new_contour.endInd = len(new_contour) - 1 - self.extend_upper
        # new_contour was interpolated from a high-resolution contour, so should not need
        # a large width for refinement - use width/100. instead of 'width'
        return new_contour.getRefined(width=width/100., atol=atol)

    def temporaryExtend(self, *, extend_lower=0, extend_upper=0, ds_lower=None,
            ds_upper=None):
        """
        Add temporary guard-cell points to the beginning and/or end of a contour
        Use coarseInterp to extrapolate as using a bigger spacing gives a more stable
        extrapolation.
        """
        if extend_lower > 0:
            if ds_lower is None:
                ds = self.distance[1] - self.distance[0]
            else:
                ds = ds_lower
            for i in range(extend_lower):
                interp, distance_estimate = self._coarseInterp()
                new_point = interp(distance_estimate[0] - ds)
                self.prepend(self.refinePoint(new_point, new_point - self[0]))
                if self.startInd >= 0:
                    self.startInd += 1
                if self.endInd >= 0:
                    self.endInd += 1
        if extend_upper > 0:
            if ds_upper is None:
                ds = self.distance[-1] - self.distance[-2]
            else:
                ds = ds_upper
            for i in range(extend_upper):
                interp, distance_estimate = self._coarseInterp()
                new_point = interp(distance_estimate[-1] + ds)
                self.append(self.refinePoint(new_point, new_point - self[-1]))

    def plot(self, *args, plotPsi=False, **kwargs):
        from matplotlib import pyplot
        Rpoints = [p.R for p in self]
        Zpoints = [p.Z for p in self]
        if plotPsi:
            R = numpy.linspace(min(Rpoints), max(Rpoints), 100)
            Z = numpy.linspace(min(Zpoints), max(Zpoints), 100)
            pyplot.contour(R, Z, self.psi(R[numpy.newaxis, :], Z[:, numpy.newaxis]))
        pyplot.plot(Rpoints, Zpoints, *args, **kwargs)

class EquilibriumRegion(PsiContour):
    """
    Specialization of PsiContour for representing an equilibrium segment, which is a
    poloidal segment based around a contour (normally a segment of a separatrix). Includes
    members giving the connections to other regions and to list the X-points at the
    boundaries where the contour starts or ends.
    """

    def __init__(self, equilibrium, name, nSegments, user_options, options, *args,
            **kwargs):
        super().__init__(*args, **kwargs)
        self.equilibrium = equilibrium
        self.name = name
        self.nSegments = nSegments

        self.user_options  = user_options

        # Set up options for this object: poloidal spacing options need setting
        self.options = options.copy()

        # Set object-specific options
        assert self.options.nx is not None, 'nx must be set'
        assert self.options.ny is not None, 'ny must be set'

        # Allow options to be overridden by kwargs
        self.options = self.options.push(kwargs)

        self.setupOptions(force=False)
        self.ny_noguards = self.options.ny
        self.global_xind = 0 # 0 since EquilibriumRegion represents the contour at the separatrix

        self.xPointsAtStart = []
        self.xPointsAtEnd = []

        # Set if this segment starts on a wall, with value of vector along wall
        self.wallSurfaceAtStart = None

        # Set if this segment ends on a wall, with value of vector along wall
        self.wallSurfaceAtEnd = None

        self.connections = []
        self.psi_vals = []
        self.separatrix_radial_index = 0

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

    def setupOptions(self, *, force):
        def setoption(key, val):
            if force:
                self.options[key] = val
            else:
                setDefault(self.options, key, val)

        # Set default values depending on options.kind
        if self.options.kind.split('.')[0] == 'wall':
            setoption('sqrt_b_lower', self.user_options.target_poloidal_spacing_length)
            setoption('polynomial_d_lower', self.user_options.nonorthogonal_target_poloidal_spacing_length)
            setoption('nonorthogonal_range_lower', self.user_options.nonorthogonal_target_poloidal_spacing_range)
            setoption('nonorthogonal_range_lower_inner',
                    self.user_options.nonorthogonal_target_poloidal_spacing_range_inner)
            setoption('nonorthogonal_range_lower_outer',
                    self.user_options.nonorthogonal_target_poloidal_spacing_range_outer)
        elif self.options.kind.split('.')[0] == 'X':
            setoption('sqrt_a_lower', self.user_options.xpoint_poloidal_spacing_length)
            setoption('sqrt_b_lower', 0.)
            setoption('polynomial_d_lower', self.user_options.nonorthogonal_xpoint_poloidal_spacing_length)
            setoption('nonorthogonal_range_lower', self.user_options.nonorthogonal_xpoint_poloidal_spacing_range)
            setoption('nonorthogonal_range_lower_inner',
                    self.user_options.nonorthogonal_xpoint_poloidal_spacing_range_inner)
            setoption('nonorthogonal_range_lower_outer',
                    self.user_options.nonorthogonal_xpoint_poloidal_spacing_range_outer)
        else:
            raise ValueError('Unrecognized value before \'.\' in kind=' + str(kind))
        if self.options.kind.split('.')[1] == 'wall':
            setoption('sqrt_b_upper', self.user_options.target_poloidal_spacing_length)
            setoption('polynomial_d_upper', self.user_options.nonorthogonal_target_poloidal_spacing_length)
            setoption('nonorthogonal_range_upper', self.user_options.nonorthogonal_target_poloidal_spacing_range)
            setoption('nonorthogonal_range_upper_inner',
                    self.user_options.nonorthogonal_target_poloidal_spacing_range_inner)
            setoption('nonorthogonal_range_upper_outer',
                    self.user_options.nonorthogonal_target_poloidal_spacing_range_outer)
        elif self.options.kind.split('.')[1] == 'X':
            setoption('sqrt_a_upper', self.user_options.xpoint_poloidal_spacing_length)
            setoption('sqrt_b_upper', 0.)
            setoption('polynomial_d_upper', self.user_options.nonorthogonal_xpoint_poloidal_spacing_length)
            setoption('nonorthogonal_range_upper', self.user_options.nonorthogonal_xpoint_poloidal_spacing_range)
            setoption('nonorthogonal_range_upper_inner',
                    self.user_options.nonorthogonal_xpoint_poloidal_spacing_range_inner)
            setoption('nonorthogonal_range_upper_outer',
                    self.user_options.nonorthogonal_xpoint_poloidal_spacing_range_outer)
        else:
            raise ValueError('Unrecognized value before \'.\' in self.options.kind='
                    + str(self.options.kind))

    def copy(self):
        result = EquilibriumRegion(self.equilibrium, self.name, self.nSegments,
                self.user_options, self.options, deepcopy(self.points), self.psi,
                self.psival)
        result.xPointsAtStart = deepcopy(self.xPointsAtStart)
        result.xPointsAtEnd = deepcopy(self.xPointsAtEnd)
        result.wallSurfaceAtStart = deepcopy(self.wallSurfaceAtStart)
        result.wallSurfaceAtEnd = deepcopy(self.wallSurfaceAtEnd)
        result.connections = deepcopy(self.connections)
        result.psi_vals = deepcopy(self.psi_vals)
        result.separatrix_radial_index = self.separatrix_radial_index
        result.startInd = self.startInd
        result.endInd = self.endInd
        result.extend_lower = self.extend_lower
        result.extend_upper = self.extend_upper
        return result

    def newRegionFromPsiContour(self, contour):
        result = EquilibriumRegion(self.equilibrium, self.name, self.nSegments,
                self.user_options, self.options, contour.points, contour.psi,
                contour.psival)
        result.xPointsAtStart = deepcopy(self.xPointsAtStart)
        result.xPointsAtEnd = deepcopy(self.xPointsAtEnd)
        result.wallSurfaceAtStart = deepcopy(self.wallSurfaceAtStart)
        result.wallSurfaceAtEnd = deepcopy(self.wallSurfaceAtEnd)
        result.connections = deepcopy(self.connections)
        result.psi_vals = deepcopy(self.psi_vals)
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
            result += self.user_options.y_boundary_guards
        if self.connections[radialIndex]['upper'] is None:
            result += self.user_options.y_boundary_guards
        return result

    def nxOutsideSeparatrix(self):
        # Note: includes point at separatrix
        return 1 + sum(2*n for n in self.options.nx[self.separatrix_radial_index:])

    def nxInsideSeparatrix(self):
        # Note: also includes point at separatrix
        return 1 + sum(2*n for n in self.options.nx[:self.separatrix_radial_index])

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
            extend_lower = 2*self.user_options.y_boundary_guards
        else:
            extend_lower = 0
        if self.connections[radialIndex]['upper'] is None:
            extend_upper = 2*self.user_options.y_boundary_guards
        else:
            extend_upper = 0
        sfunc = self.getSfuncFixedSpacing(2*self.ny_noguards+1,
                self.distance[self.endInd] - self.distance[self.startInd])
        return self.newRegionFromPsiContour(super().getRegridded(2*self.ny_noguards + 1,
            extend_lower=extend_lower, extend_upper=extend_upper, sfunc=sfunc, **kwargs))

    def _checkMonotonic(self, sfunc_list, *, xind=None, total_distance=0., prefix=''):
        # Check new_sfunc is monotonically increasing
        indices = numpy.arange(-self.extend_lower,
                2*self.ny_noguards + self.extend_upper + 1, dtype=float)
        scheck = sfunc_list[0][0](indices)
        if numpy.any(scheck[1:] < scheck[:-1]):
            from matplotlib import pyplot
            print('at global xind', xind)
            pyplot.figure()
            for sfunc, label in sfunc_list:
                if sfunc is not None:
                    pyplot.plot(indices, sfunc(indices), label=label)
            pyplot.axhline(0.)
            pyplot.axhline(total_distance)
            pyplot.legend()
            pyplot.show()
            decreasing = numpy.where(scheck[1:] < scheck[:-1])[0] + 1
            raise ValueError('In region ' + self.name + 'combined spacing function is '
                    + 'decreasing at indices ' + str(decreasing) + ' on contour of '
                    + 'length ' + str(len(self)) + '. It may help to increase/decrease '
                    + prefix + 'target_poloidal_spacing_length or '
                    + prefix + 'xpoint_poloidal_spacing_length.')

    def getSfuncFixedSpacing(self, npoints, distance, *, method=None):
        if method is None:
            if self.user_options.orthogonal:
                method = self.user_options.poloidal_spacing_method
            else:
                method = 'nonorthogonal'

        if method == 'sqrt':
            if self.user_options.poloidalfunction_diagnose:
                print('in sqrt method:')
                print('N_norm =', self.options.N_norm)
                print('a_lower =', self.options.sqrt_a_lower)
                print('b_lower =', self.options.sqrt_b_lower)
                print('a_upper =', self.options.sqrt_a_upper)
                print('b_upper =', self.options.sqrt_b_upper)
            sfunc = self.getSqrtPoloidalDistanceFunc(distance, npoints-1,
                    self.options.N_norm, b_lower=self.options.sqrt_b_lower,
                    a_lower=self.options.sqrt_a_lower, b_upper=self.options.sqrt_b_upper,
                    a_upper=self.options.sqrt_a_upper)
            self._checkMonotonic([(sfunc, 'sqrt')], total_distance=distance)
        elif method == 'polynomial':
            sfunc = self.getPolynomialPoloidalDistanceFunc(distance, npoints-1,
                    self.options.N_norm, d_lower=self.options.polynomial_d_lower,
                    d_upper=self.options.polynomial_d_upper)
            self._checkMonotonic([(sfunc, 'sqrt')], total_distance=distance)
        elif method == 'nonorthogonal':
            nonorth_method = self.user_options.nonorthogonal_spacing_method
            if nonorth_method == 'poloidal_orthogonal_combined':
                return self.combineSfuncs(self, None)
            elif nonorth_method == 'perp_orthogonal_combined':
                if self.wallSurfaceAtStart is not None:
                    # surface is a wall
                    lower_surface = self.wallSurfaceAtStart
                else:
                    # Use fixed poloidal spacing when gridding the separatrix contour so that
                    # the grid spacing is the same in different regions which share a
                    # separatrix segment but have different perpendicular vectors at the
                    # X-point
                    lower_surface = None

                if self.wallSurfaceAtEnd is not None:
                    # surface is a wall
                    upper_surface = self.wallSurfaceAtEnd
                else:
                    # Use fixed poloidal spacing when gridding the separatrix contour so that
                    # the grid spacing is the same in different regions which share a
                    # separatrix segment but have different perpendicular vectors at the
                    # X-point
                    upper_surface = None

                sfunc = self.combineSfuncs(self, None, lower_surface, upper_surface)
            elif nonorth_method == 'combined':
                # Use fixed poloidal spacing when gridding the separatrix contour so that
                # the grid spacing is the same in different regions which share a
                # separatrix segment but have different perpendicular vectors at the
                # X-point
                return self.combineSfuncs(self, None)
            elif nonorth_method == 'orthogonal':
                orth_method = self.user_options.poloidal_spacing_method
                sfunc = self.getSfuncFixedSpacing(npoints, distance, method=orth_method)
            else:
                sfunc = self.getSfuncFixedSpacing(npoints, distance, method=nonorth_method)
        else:
            raise ValueError('Unrecognized option '
                             + str(self.user_options.poloidal_spacing_method)
                             + ' for poloidal spacing method')

        if self.user_options.poloidalfunction_diagnose:
            from matplotlib import pyplot
            indices = numpy.linspace(0., npoints - 1, 1000)
            pyplot.plot(indices, sfunc(indices))
            pyplot.axhline(0., color='r')
            pyplot.axhline(distance, color='r')
            pyplot.xlabel('index')
            pyplot.ylabel('s')
            pyplot.title(self.name + ' ' + method)
            pyplot.show()

        return sfunc

    def combineSfuncs(self, contour, sfunc_orthogonal, vec_lower=None, vec_upper=None):
        # this sfunc gives:
        # * - if vec_lower is None: fixed poloidal spacing at the beginning of the contour
        #   - otherwise fixed spacing perpendicular to vec_lower at the beginning of the
        #     contour
        # * - if vec_upper is None: fixed poloidal spacing at the end of the contour
        #   - otherwise fixed spacing perpendicular to vec_lower at the end of the contour
        # * Tends to orthogonal spacing far from the ends, unless sfunc_orthogonal is
        #   None, in which case it sets points so that if combineSfuncs is called again on
        #   the same contour, but with sfunc_orthogonal=contour.contourSfunc() then the
        #   same spacing is given
        if vec_lower is None:
            sfunc_fixed_lower = self.getSfuncFixedSpacing(
                    2*self.ny_noguards + 1, contour.totalDistance(), method='polynomial')
        else:
            sfunc_fixed_lower, sperp_func_lower = self.getSfuncFixedPerpSpacing(
                    2*self.ny_noguards + 1, contour, vec_lower, True)

        if vec_upper is None:
            sfunc_fixed_upper = self.getSfuncFixedSpacing(
                    2*self.ny_noguards + 1, contour.totalDistance(), method='polynomial')
        else:
            sfunc_fixed_upper, sperp_func_upper = self.getSfuncFixedPerpSpacing(
                    2*self.ny_noguards + 1, contour, vec_upper, False)

        if self.options.nonorthogonal_range_lower is not None:
            range_lower = self.options.nonorthogonal_range_lower
            range_lower_inner = self.options.nonorthogonal_range_lower_inner
            range_lower_outer = self.options.nonorthogonal_range_lower_outer
        else:
            range_lower = self.options.polynomial_d_lower
            range_lower_inner = self.options.polynomial_d_lower
            range_lower_outer = self.options.polynomial_d_lower

        if self.options.nonorthogonal_range_upper is not None:
            range_upper = self.options.nonorthogonal_range_upper
            range_upper_inner = self.options.nonorthogonal_range_upper_inner
            range_upper_outer = self.options.nonorthogonal_range_upper_outer
        else:
            range_upper = self.options.polynomial_d_upper
            range_upper_inner = self.options.polynomial_d_upper
            range_upper_outer = self.options.polynomial_d_upper

        N_norm = self.options.N_norm

        index_length = 2.*self.ny_noguards

        # Set up radial variation of weights
        if range_lower is not None:
            # this_range_lower is range_lower at separatrix, range_lower_outer at outer
            # radial boundary, range_lower_inner at inner radial boundary and has zero
            # radial derivative at the separatrix
            ix = float(contour.global_xind)
            if ix >= 0:
                xweight = (ix / (self.nxOutsideSeparatrix() - 1.))**self.user_options.nonorthogonal_radial_range_power
                this_range_lower = (1. - xweight)*range_lower + xweight*range_lower_outer
            else:
                xweight = (-ix / (self.nxInsideSeparatrix() - 1.))**self.user_options.nonorthogonal_radial_range_power
                this_range_lower = (1. - xweight)*range_lower + xweight*range_lower_inner
        if range_upper is not None:
            # this_range_upper is range_upper at separatrix, range_upper_outer at outer
            # radial boundary, range_upper_inner at inner radial boundary and has zero
            # radial derivative at the separatrix
            ix = float(contour.global_xind)
            if ix >= 0:
                xweight = (ix / (self.nxOutsideSeparatrix() - 1.))**self.user_options.nonorthogonal_radial_range_power
                this_range_upper = (1. - xweight)*range_upper + xweight*range_upper_outer
            else:
                xweight = (-ix / (self.nxInsideSeparatrix() - 1.))**self.user_options.nonorthogonal_radial_range_power
                this_range_upper = (1. - xweight)*range_upper + xweight*range_upper_inner

        if range_lower is not None and range_upper is not None:
            def new_sfunc(i):
                sfixed_lower = sfunc_fixed_lower(i)

                sfixed_upper = sfunc_fixed_upper(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_lower so it is 1. at the lower boundary and 0. at the
                # upper boundary and the gradient is zero at both boundaries
                weight_lower = numpy.piecewise(i,
                        [i < 0., i > index_length],
                        [1., 0., lambda i: numpy.exp(-(i/N_norm/this_range_lower)**2)])

                # define weight_upper so it is 1. at the upper boundary and 0. at the
                # lower boundary and the gradient is zero at both boundaries
                weight_upper = numpy.piecewise(i,
                    [i < 0., i > index_length],
                    [0., 1., lambda i: numpy.exp(-((index_length - i)/N_norm/this_range_upper)**2)])

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

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_lower so it is 1. at the lower boundary and the gradient
                # is zero at the lower boundary.
                weight_lower = numpy.piecewise(i,
                        [i < 0., i > index_length],
                        [1., 0., lambda i: numpy.exp(-(i/N_norm/this_range_lower)**2)])

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncs again for this contour
                    # with sfunc_orthogonal from self.contourSfunc() then we get the same
                    # spacing again. This is used to make the contours along the
                    # separatrix keep the same values when pushing the other contours
                    # towards orthogonal positions
                    # s = weight_lower*sfixed_lower + (1. - weight_lower)*s
                    sorth = sfixed_lower

                return (weight_lower)*sfixed_lower + (1. - weight_lower) * sorth
        elif range_upper is not None:
            def new_sfunc(i):
                sfixed_upper = sfunc_fixed_upper(i)

                if sfunc_orthogonal is None:
                    sorth = None
                else:
                    sorth = sfunc_orthogonal(i)

                # define weight_upper so it is 1. at the upper boundary and the gradient
                # is zero at the upper boundary.
                weight_upper = numpy.piecewise(i,
                    [i < 0., i > index_length],
                    [0., 1., lambda i: numpy.exp(-((index_length - i)/N_norm/this_range_upper)**2)])

                if sorth is None:
                    # Fix spacing so that if we call combineSfuncs again for this contour
                    # with sfunc_orthogonal from self.contourSfunc() then we get the same
                    # spacing again. This is used to make the contours along the
                    # separatrix keep the same values when pushing the other contours
                    # towards orthogonal positions
                    # s = weight_upper*sfixed_upper + (1. - weight_upper)*s
                    sorth = sfixed_upper

                return (weight_upper)*sfixed_upper + (1. - weight_upper) * sorth
        else:
            assert sfunc_orthogonal is not None, 'Without range_lower or range_upper, cannot use with sfunc_orthogonal=None'
            def new_sfunc(i):
                return sfunc_orthogonal(i)

        try:
            self._checkMonotonic([(new_sfunc, 'combined'), (sfunc_orthogonal, 'orthogonal'),
                (sfunc_fixed_lower, 'fixed perp lower'), (sfunc_fixed_upper, 'fixed perp upper')],
                xind=contour.global_xind, total_distance=contour.totalDistance(),
                prefix='nonorthogonal_')
        except ValueError:
            print('check lower ranges', range_lower_inner, range_lower, range_lower_outer,
                    this_range_lower)
            print('check upper ranges', range_upper_inner, range_upper, range_upper_outer,
                    this_range_upper)
            raise

        return new_sfunc

    def getSfuncFixedPerpSpacing(self, N, contour,
            surface_direction, lower):
        """
        Return a function s(i) giving poloidal distance as a function of index-number.
        Construct so that ds_perp/diN = d_lower at the lower end or ds_perp/diN = d_upper
        at the upper end, where s_perp is the distance normal to the vector
        'surface_direction'.
        """
        N_norm = self.options.N_norm

        if self.options.perp_d_lower is not None:
            d_lower = self.options.perp_d_lower
        else:
            d_lower = self.options.polynomial_d_lower

        if self.options.perp_d_upper is not None:
            d_upper = self.options.perp_d_upper
        else:
            d_upper = self.options.polynomial_d_upper

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
                a_lower = 0.
            if a_upper is None:
                a_upper = 0.
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
            if a_lower == 0.:
                # Gradient must be strictly positive as there is no positive a_lower piece
                assert b/(2.*numpy.sqrt(N/N_norm)) + d > 0., 'gradient of non-singular part should be positive at start'
            else:
                # Might be 0., so allow tolerance for small negative values due to
                # rounding errors
                assert b/(2.*numpy.sqrt(N/N_norm)) + d > -self.user_options.sfunc_checktol, 'gradient of non-singular part should be positive at start'
            # upper boundary:
            assert b >= 0., 'sqrt part of function should be positive at end'
            if a_upper == 0.:
                # Gradient must be strictly positive as there is no positive a_upper piece
                assert a/(2.*numpy.sqrt(N/N_norm)) + d + 2.*e*N/N_norm + 3.*f*(N/N_norm)**2 > 0., 'gradient of non-singular part should be positive at end'
            else:
                # Might be 0., so allow tolerance for small negative values due to
                # rounding errors
                assert a/(2.*numpy.sqrt(N/N_norm)) + d + 2.*e*N/N_norm + 3.*f*(N/N_norm)**2 > - self.user_options.sfunc_checktol, 'gradient of non-singular part should be positive at end'

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
      - self.fpolprime: psi-derivative of fpol
      - self.Rmin, self.Rmax, self.Zmin, self.Zmax: positions of the corners of a bounding
        box for the gridding
      - self.regions: OrderedDict of EquilibriumRegion objects that specify this equilibrium
      - self.wall: list of Point2D giving vertices of polygon representing the wall, in
        anti-clockwise order; assumed to be closed so last element and first are taken to
        be connected
    """
    def __init__(self, **kwargs):
        """
        Does some generic setup common to all Equilibrium derived classes.
        Note: should be called by derived class __init__() constructor after the
        user_options have been initialized.
        """

        # Set up internal options
        # '.push(kwargs)' here lets the kwargs override any values (including for
        # 'internal' options that should not need to be set by the user) set as defaults
        # from HypnotoadOptions
        self.options = HypnotoadInternalOptions.push(kwargs)

        # Set some global parameters for PsiContours and FineContours
        # Convert self.user_options to a dict so we can use it to set the values in
        # FineContour.options using 'push'
        PsiContour.options = PsiContour.options.push(dict(self.user_options))
        FineContour.options = FineContour.options.push(dict(self.user_options))

        # Set some default options
        setDefault(self.user_options, 'nonorthogonal_xpoint_poloidal_spacing_range_inner',
                self.user_options.nonorthogonal_xpoint_poloidal_spacing_range)
        setDefault(self.user_options, 'nonorthogonal_xpoint_poloidal_spacing_range_outer',
                self.user_options.nonorthogonal_xpoint_poloidal_spacing_range)
        setDefault(self.user_options, 'nonorthogonal_target_poloidal_spacing_range_inner',
                self.user_options.nonorthogonal_target_poloidal_spacing_range)
        setDefault(self.user_options, 'nonorthogonal_target_poloidal_spacing_range_outer',
                self.user_options.nonorthogonal_target_poloidal_spacing_range)

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
        assert lRegion.options.nx[lowerSegment] == uRegion.options.nx[upperSegment], 'nx should match across connection'

        lRegion.connections[lowerSegment]['upper'] = (upperRegion, upperSegment)
        uRegion.connections[upperSegment]['lower'] = (lowerRegion, lowerSegment)

    def magneticFunctionsFromGrid(self, R, Z, psiRZ):
        from .dct_interpolation import DCT_2D

        self._dct = DCT_2D(R, Z, psiRZ)

        self.psi = lambda R, Z: self._dct(R, Z)
        modGradpsiSquared = lambda R, Z: self._dct.ddR(R, Z)**2 + self._dct.ddZ(R, Z)**2
        self.f_R = lambda R, Z: self._dct.ddR(R, Z) / modGradpsiSquared(R, Z)
        self.f_Z = lambda R, Z: self._dct.ddZ(R, Z) / modGradpsiSquared(R, Z)
        self.Bp_R = lambda R, Z: self._dct.ddZ(R, Z) / R
        self.Bp_Z = lambda R, Z: -self._dct.ddR(R, Z) / R
        self.d2psidR2 = self._dct.d2dR2
        self.d2psidZ2 = self._dct.d2dZ2
        self.d2psidRdZ = self._dct.d2dRdZ

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

    def findSaddlePoint(self, p1, p2, atol=2.e-8):
        """
        Find a saddle point in the function self.psi atol is the tolerance on the position
        of the saddle point. p1, p2 are the positions of two adjacent corners of the
        square box to search for a saddle point in (the other corners p3 and p4 are taken
        to be to the right of the line p1->p2).
        """

        # Note: in this method, Point2D objects are used as displacement vectors as well
        # as as points.
        def dot(v1, v2):
            return v1.R*v2.R + v1.Z*v2.Z

        # length of sides of the square
        a = calc_distance(p1, p2)

        # unit vector along p1->p2 or p4->p3
        e1 = (p2 - p1) / a

        # unit vector along p2->p3 or p1->p4
        e2 = Point2D(e1.Z, -e1.R)

        p3 = p2 + a*e2
        p4 = p1 + a*e2

        # For the purposes of naming variables here, take p1 to be 'bottom left', p2 to
        # be 'top left', p3 to be 'top right' and p4 to be 'bottom right'
        posLeft, minLeft = self.findExtremum_1d(p1, p2)
        posTop, minTop = self.findExtremum_1d(p2, p3)
        posRight, minRight = self.findExtremum_1d(p3, p4)
        posBottom, minBottom = self.findExtremum_1d(p4, p1)

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

        extremumVert = p3
        extremumHoriz = p1

        count = 0
        while calc_distance(extremumVert, extremumHoriz) > atol:
            count = count+1

            extremumVert = vertSearch(posBottom, posTop, 0.5*atol)
            # set position along e2 direction, keep position along e1 fixed (i.e. stay on
            # left or right edge)
            deltaz = dot(extremumVert - p1, e1)
            posLeft = p1 + deltaz*e1
            posRight = p4 + deltaz*e1

            extremumHoriz = horizSearch(posLeft, posRight, 0.5*atol)
            # set position along e1 direction, keep position along e2 fixed (i.e. stay on
            # top or bottom edge)
            deltar = dot(extremumHoriz - p1, e2)
            posBottom = p1 + deltar*e2
            posTop = p2 + deltar*e2

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
        closed_wall = self.wall + [self.wall[0]]
        wallarray = numpy.array([(p.R, p.Z) for p in closed_wall])
        intersects = find_intersections(wallarray, p1, p2)
        if intersects is not None:
            intersect = Point2D(*intersects[0, :])
            assert intersects.shape[0] < 3, 'too many intersections with wall'
            if intersects.shape[0] > 1:
                second_intersect = Point2D(*intersects[1, :])
                assert numpy.abs(intersect.R - second_intersect.R) < intersect_tolerance and numpy.abs(intersect.Z - second_intersect.Z) < intersect_tolerance, 'Multiple intersections with wall found'
        else:
            intersect = None

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

    def plotPotential(self, Rmin=None, Rmax=None, Zmin=None, Zmax=None, npoints=100,
            ncontours=40, labels=True, **kwargs):
        from matplotlib import pyplot

        if Rmin is None: Rmin = self.Rmin
        if Rmax is None: Rmax = self.Rmax
        if Zmin is None: Zmin = self.Zmin
        if Zmax is None: Zmax = self.Zmax

        R = numpy.linspace(Rmin, Rmax, npoints)
        Z = numpy.linspace(Zmin, Zmax, npoints)
        ax = pyplot.axes(aspect='equal')
        contours = ax.contour(R, Z, self.psi(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T,
                ncontours, **kwargs)
        if labels:
            pyplot.clabel(contours, inline=False, fmt='%1.3g')

    def plotSeparatrix(self):
        from matplotlib import pyplot
        for region in self.regions.values():
            R = [p.R for p in region]
            Z = [p.Z for p in region]
            pyplot.scatter(R, Z, marker='x', label=region.name)

    def _getOptionsAsString(self):
        import yaml

        result = ""
        result += yaml.dump(self.equilibOptions)

        mesh_options_dict = {'Mesh':{}}
        m = mesh_options_dict['Mesh']
        for key,val in self.user_options.items():
            if val is not None:
                m[key] = str(val)

        result += yaml.dump(mesh_options_dict)

        return result

    def saveOptions(self, filename='hypnotoad_options.yaml'):
        with open(filename, 'x') as f:
            f.write(self._getOptionsAsString())

class DoubleNull(Equilibrium):
    """
    Analyse tokamak equilibrium - single-null, connected double-null or disconnected
    double-null - with regions lined up according to BOUT++ requirements.
    """

    # Add DoubleNull-specific options and default values
    user_options = HypnotoadOptions.add(
            nx_core = None,
            nx_between = None,
            nx_sol = None,
            ny_inner_lower_divertor = None,
            ny_inner_core = None,
            ny_inner_upper_divertor = None,
            ny_outer_upper_divertor = None,
            ny_outer_core = None,
            ny_outer_lower_divertor = None,
            upper_target_y_boundary_guards = None,
            psi_core = None,
            psi_sol = None,
            psi_inner_sol = None,
            )

    def __init__(self, equilibOptions, meshOptions, **kwargs):
        self.user_options = DoubleNull.user_options.push(meshOptions).push(kwargs)
        self.options = HypnotoadInternalOptions.push(kwargs)
        self.equilibOptions = equilibOptions

        raise ValueError('need to set up options here')

        setDefault(self.options, 'psi_lower_pf', self.user_options.psi_core)
        setDefault(self.options, 'psi_upper_pf', self.user_options.psi_core)

        setDefault(self.options, 'psi_inner_sol', self.options.psi_sol)
        # this option can only be set different from psi_sol in a double-null
        # configuration (i.e. if there are upper divertor legs)
        if self.options.psi_sol != self.options.psi_inner_sol:
            assert self.options.ny_inner_upper_divertor > 0 or self.options.ny_outer_upper_divertor > 0, 'inner and outer SOL should be separated by an upper divertor, i.e. topology should be double-null'

        self.psi_spacing_separatrix_multiplier = self.readOption('psi_spacing_separatrix_multiplier', None)
        if self.psi_spacing_separatrix_multiplier is not None:
            if self.options.nx_between > 0:
                raise ValueError("Cannot use psi_spacing_separatrix_multiplier when "
                                 "there are points between two separatrices - to get "
                                 "the same effect, increase the number of points in the "
                                 "between-separatrix region.")

        if ( (self.options.ny_inner_upper_divertor > 0 or self.options.ny_outer_upper_divertor > 0)
                and (self.options.ny_inner_lower_divertor > 0 or self.options.ny_inner_core > 0
                     or self.options.ny_outer_core > 0 or self.options.ny_outer_lower_divertor > 0) ):
            # there are two distinct divertor/limiter targets
            #  - the upper divertor legs are not empty, and also the other regions are not
            #    all empty
            setDefault(self.options, 'upper_target_y_boundary_guards',
                    self.options.y_boundary_guards)
        else:
            assert self.options.upper_target_y_boundary_guards is None, 'Does not make sense for upper target to have guard cells, because it does not exist'
