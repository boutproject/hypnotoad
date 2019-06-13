"""
Classes to handle Meshes and geometrical quantities for generating BOUT++ grids
"""

from copy import deepcopy
import numbers
import warnings

import numpy
from scipy.integrate import solve_ivp

from .equilibrium import calc_distance, Point2D, FineContour

class MultiLocationArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    """
    Container for arrays representing points at different cell locations
    Not all have to be filled.
    """
    _centre_array = None
    _xlow_array = None
    _ylow_array = None
    _corners_array = None

    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny

    @property
    def centre(self):
        if self._centre_array is None:
            self._centre_array = numpy.zeros([self.nx, self.ny])
        return self._centre_array

    @centre.setter
    def centre(self, value):
        if self._centre_array is None:
            self._centre_array = numpy.zeros([self.nx, self.ny])
        self._centre_array[...] = value

    @property
    def xlow(self):
        if self._xlow_array is None:
            self._xlow_array = numpy.zeros([self.nx + 1, self.ny])
        return self._xlow_array

    @xlow.setter
    def xlow(self, value):
        if self._xlow_array is None:
            self._xlow_array = numpy.zeros([self.nx + 1, self.ny])
        self._xlow_array[...] = value

    @property
    def ylow(self):
        if self._ylow_array is None:
            self._ylow_array = numpy.zeros([self.nx, self.ny + 1])
        return self._ylow_array

    @ylow.setter
    def ylow(self, value):
        if self._ylow_array is None:
            self._ylow_array = numpy.zeros([self.nx, self.ny + 1])
        self._ylow_array[...] = value

    @property
    def corners(self):
        if self._corners_array is None:
            self._corners_array = numpy.zeros([self.nx + 1, self.ny + 1])
        return self._corners_array

    @corners.setter
    def corners(self, value):
        if self._corners_array is None:
            self._corners_array = numpy.zeros([self.nx + 1, self.ny + 1])
        self._corners_array[...] = value

    # The following __array_ufunc__ implementation allows the MultiLocationArray class to
    # be handled by Numpy functions, and add, subtract, etc. like an ndarray.
    # The implementation is mostly copied from the example in
    # https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.mixins.NDArrayOperatorsMixin.html#numpy.lib.mixins.NDArrayOperatorsMixin

    # One might also consider adding the built-in list type to this
    # list, to support operations like np.add(array_like, list)
    _HANDLED_TYPES = (numpy.ndarray, numbers.Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.get('out', ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use MultiLocationArray instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle MultiLocationArray objects.
            if not isinstance(x, self._HANDLED_TYPES + (MultiLocationArray,)):
                return NotImplemented

        result = MultiLocationArray(self.nx, self.ny)

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._centre_array is not None:
            this_inputs = tuple(x._centre_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.centre if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].centre = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.centre = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._xlow_array is not None:
            this_inputs = tuple(x._xlow_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.xlow if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].xlow = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.xlow = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._ylow_array is not None:
            this_inputs = tuple(x._ylow_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.ylow if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].ylow = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.ylow = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._corners_array is not None:
            this_inputs = tuple(x._corners_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.corners if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].corners = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.corners = this_result

        return result

    def zero(self):
        # Initialise all locations, set them to zero and return the result
        self.centre = 0.
        self.xlow = 0.
        self.ylow = 0.
        self.corners = 0.
        return self

class MeshRegion:
    """
    A simple rectangular region of a Mesh, that connects to one other region (or has a
    boundary) on each edge.
    Note that these regions include cell face and boundary points, so there are
    (2nx+1)*(2ny+1) points for an nx*ny grid.
    """
    def __init__(self, meshParent, myID, equilibriumRegion, connections, radialIndex):
        self.name = equilibriumRegion.name+'('+str(radialIndex)+')'
        print('creating region', myID, '-', self.name)

        # the Mesh object that owns this MeshRegion
        self.meshParent = meshParent

        # ID that Mesh uses to keep track of its MeshRegions
        self.myID = myID

        # EquilibriumRegion representing the segment associated with this region
        self.equilibriumRegion = equilibriumRegion.copy()

        self.user_options = self.equilibriumRegion.user_options
        self.options = self.equilibriumRegion.options

        # sizes of the grid in this MeshRegion, include boundary guard cells
        self.nx = self.options.nx[radialIndex]
        self.ny = self.equilibriumRegion.ny(radialIndex)
        self.ny_noguards = self.equilibriumRegion.ny_noguards

        # psi values for radial grid
        self.psi_vals = numpy.array(self.equilibriumRegion.psi_vals[radialIndex])
        assert len(self.psi_vals) == 2*self.nx + 1, 'should be a psi value for each radial point'

        # Dictionary that specifies whether a boundary is connected to another region or
        # is an actual boundary
        self.connections = connections

        # Number of this region, counting radially outward
        self.radialIndex = radialIndex

        # Number of this region in its y-group
        self.yGroupIndex = None

        # Absolute tolerance for checking if two points are the same
        self.atol = 1.e-7

        # get points in this region
        self.contours = []
        if self.radialIndex < self.equilibriumRegion.separatrix_radial_index:
            # region is inside separatrix, so need to follow line from the last psi_val to
            # the first
            temp_psi_vals = self.psi_vals[::-1]
        else:
            temp_psi_vals = self.psi_vals

        # set sign of step in psi towards this region from primary separatrix
        if temp_psi_vals[-1] - self.equilibriumRegion.psival > 0:
            psi_sep_plus_delta = self.equilibriumRegion.psival + self.user_options.poloidal_spacing_delta_psi
        else:
            psi_sep_plus_delta = self.equilibriumRegion.psival - self.user_options.poloidal_spacing_delta_psi

        # Make vector along grad(psi) at start of equilibriumRegion
        vec_points = followPerpendicular(
                self.meshParent.equilibrium.f_R, self.meshParent.equilibrium.f_Z,
                self.equilibriumRegion[self.equilibriumRegion.startInd],
                self.equilibriumRegion.psival,
                [self.equilibriumRegion.psival, psi_sep_plus_delta],
                rtol=self.user_options.follow_perpendicular_rtol,
                atol=self.user_options.follow_perpendicular_atol)
        self.equilibriumRegion.gradPsiSurfaceAtStart = (vec_points[1].as_ndarray() - vec_points[0].as_ndarray())
        # Make vector along grad(psi) at end of equilibriumRegion
        vec_points = followPerpendicular(
                self.meshParent.equilibrium.f_R, self.meshParent.equilibrium.f_Z,
                self.equilibriumRegion[self.equilibriumRegion.endInd],
                self.equilibriumRegion.psival,
                [self.equilibriumRegion.psival, psi_sep_plus_delta],
                rtol=self.user_options.follow_perpendicular_rtol,
                atol=self.user_options.follow_perpendicular_atol)
        self.equilibriumRegion.gradPsiSurfaceAtEnd = (vec_points[1].as_ndarray() - vec_points[0].as_ndarray())

        # Calculate the perp_d_lower/perp_d_upper corresponding to d_lower/d_upper on the
        # separatrix contour
        # Use self.equilibriumRegion.fine_contour for the vector along the separatrix
        # because then the vector will not change when the grid resolution changes
        if self.equilibriumRegion.wallSurfaceAtStart is None:
            # lower end
            unit_vec_separatrix = (
                    self.equilibriumRegion.fine_contour.positions[
                        self.equilibriumRegion.fine_contour.startInd + 1, :]
                    - self.equilibriumRegion.fine_contour.positions[
                        self.equilibriumRegion.fine_contour.startInd, :])
            unit_vec_separatrix /= numpy.sqrt(numpy.sum(unit_vec_separatrix**2))
            unit_vec_surface = self.equilibriumRegion.gradPsiSurfaceAtStart
            unit_vec_surface /= numpy.sqrt(numpy.sum(unit_vec_surface**2))
            cos_angle = numpy.sum(unit_vec_separatrix*unit_vec_surface)
            # this gives abs(sin_angle), but that's OK because we only want the magnitude to
            # calculate perp_d
            sin_angle = numpy.sqrt(1. - cos_angle**2)
            self.options.set(perp_d_lower=self.options.polynomial_d_lower * sin_angle)
        if self.equilibriumRegion.wallSurfaceAtEnd is None:
            # upper end
            unit_vec_separatrix = (
                    self.equilibriumRegion.fine_contour.positions[
                        self.equilibriumRegion.fine_contour.endInd - 1, :]
                    - self.equilibriumRegion.fine_contour.positions[
                        self.equilibriumRegion.fine_contour.endInd, :])
            unit_vec_separatrix /= numpy.sqrt(numpy.sum(unit_vec_separatrix**2))
            unit_vec_surface = self.equilibriumRegion.gradPsiSurfaceAtEnd
            unit_vec_surface /= numpy.sqrt(numpy.sum(unit_vec_surface**2))
            cos_angle = numpy.sum(unit_vec_separatrix*unit_vec_surface)
            # this gives abs(sin_angle), but that's OK because we only want the magnitude to
            # calculate perp_d
            sin_angle = numpy.sqrt(1. - cos_angle**2)
            self.options.set(perp_d_upper=self.options.polynomial_d_upper * sin_angle)

        print('Following perpendicular: ' + str(1) + '/'
                + str(len(self.equilibriumRegion)), end='\r')
        perp_points = followPerpendicular(self.meshParent.equilibrium.f_R,
                self.meshParent.equilibrium.f_Z, self.equilibriumRegion[0],
                self.equilibriumRegion.psival, temp_psi_vals,
                rtol=self.user_options.follow_perpendicular_rtol,
                atol=self.user_options.follow_perpendicular_atol)
        if self.radialIndex < self.equilibriumRegion.separatrix_radial_index:
            # region is inside separatrix, so points were found from last to first
            perp_points.reverse()
        for i,point in enumerate(perp_points):
            self.contours.append(self.equilibriumRegion.newContourFromSelf(points=[point],
                psival=self.psi_vals[i]))
            self.contours[i].global_xind = self.globalXInd(i)
        for i,p in enumerate(self.equilibriumRegion[1:]):
            print('Following perpendicular: ' + str(i+2) + '/'
                    + str(len(self.equilibriumRegion)), end='\r')
            perp_points = followPerpendicular(self.meshParent.equilibrium.f_R,
                    self.meshParent.equilibrium.f_Z, p,
                    self.equilibriumRegion.psival, temp_psi_vals,
                    rtol=self.user_options.follow_perpendicular_rtol,
                    atol=self.user_options.follow_perpendicular_atol)
            if self.radialIndex < self.equilibriumRegion.separatrix_radial_index:
                perp_points.reverse()
            for j,point in enumerate(perp_points):
                self.contours[j].append(point)

        # refine the contours to make sure they are at exactly the right psi-value
        for contour in self.contours:
            contour.refine(width=self.user_options.refine_width)

        if not self.user_options.orthogonal:
            self.addPointAtWallToContours()
            self.distributePointsNonorthogonal()

    def addPointAtWallToContours(self):
        # maximum number of times to extend the contour when it has not yet hit the wall
        max_extend = 100

        # should the contour intersect a wall at the lower end?
        lower_wall = self.connections['lower'] is None

        # should the contour intersect a wall at the upper end?
        upper_wall = self.connections['upper'] is None

        # sfunc_orthogonal functions created after contour has been extended past wall (if
        # necessary) but before adding the wall point to the contour (as adding this point
        # makes the spacing of points on the contour not-smooth) and adjusted for the
        # change in distance after redefining startInd to be at the wall
        self.sfunc_orthogonal_list = []

        # find wall intersections
        def correct_sfunc_orthogonal(contour, sfunc_orthogonal_original):
            distance_at_original_start = contour.distance[contour.startInd]

            distance_at_wall = contour.distance[lower_intersect_index]

            # correct sfunc_orthogonal for the distance between the point at the lower
            # wall and the original start-point
            return lambda i: sfunc_orthogonal_original(i) + distance_at_original_start - distance_at_wall

        for i_contour, contour in enumerate(self.contours):
            print('finding wall intersections:',
                    str(i_contour+1)+'/'+str(len(self.contours)), end = '\r')

            # point where contour intersects the lower wall
            lower_intersect = None

            # index of the segment of the contour that intersects the lower wall
            lower_intersect_index = 0

            # point where contour intersects the upper wall
            upper_intersect = None

            # index of the segment of the contour that intersects the upper wall
            upper_intersect_index = -2
            if lower_wall:
                if upper_wall:
                    starti = len(contour)//2
                else:
                    starti = len(contour) - 1

                # find whether one of the segments of the contour already intersects the wall
                for i in range(starti, 0, -1):
                    lower_intersect = self.meshParent.equilibrium.wallIntersection(contour[i],
                            contour[i-1])
                    if lower_intersect is not None:
                        lower_intersect_index = i-1
                        break

                count = 0
                while lower_intersect is None:
                    # contour has not yet intersected with wall, so make it longer and try
                    # again
                    contour.temporaryExtend(extend_lower=1)
                    lower_intersect = self.meshParent.equilibrium.wallIntersection(contour[1],
                            contour[0])
                    count += 1
                    assert count < max_extend, 'extended contour too far without finding wall'

            if upper_wall:
                if lower_wall:
                    starti = len(contour//2)
                else:
                    starti = 0

                # find whether one of the segments of the contour already intersects the wall
                for i in range(starti, len(contour) - 1):
                    upper_intersect = self.meshParent.equilibrium.wallIntersection(contour[i],
                            contour[i+1])
                    if upper_intersect is not None:
                        upper_intersect_index = i
                        break

                count = 0
                while upper_intersect is None:
                    # contour has not yet intersected with wall, so make it longer and try
                    # again
                    contour.temporaryExtend(extend_upper=1)
                    upper_intersect = self.meshParent.equilibrium.wallIntersection(contour[-2],
                            contour[-1])
                    count += 1
                    assert count < max_extend, 'extended contour too far without finding wall'

            # now add points on the wall(s) to the contour
            if lower_wall:
                # need to construct a new sfunc which gives distance from the wall, not the
                # distance from the original startInd

                # this sfunc would put the points at the positions along the contour where the
                # grid would be orthogonal
                sfunc_orthogonal_original = contour.contourSfunc()

                # now make lower_intersect_index the index where the point at the wall is
                # check whether one of the points is already on the wall
                if calc_distance(contour[lower_intersect_index], lower_intersect) < self.atol:
                    pass
                elif calc_distance(contour[lower_intersect_index+1], lower_intersect) < self.atol:
                    lower_intersect_index = lower_intersect_index + 1
                else:
                    # otherwise insert a new point
                    lower_intersect_index += 1
                    contour.insert(lower_intersect_index, lower_intersect)

                # contour.contourSfunc() would put the points at the positions along the
                # contour where the grid would be orthogonal
                # need to correct sfunc_orthogonal for the distance between the point at
                # the lower wall and the original start-point
                sfunc_orthogonal = correct_sfunc_orthogonal(contour,
                        sfunc_orthogonal_original)

                # start contour from the wall
                contour.startInd = lower_intersect_index

            if upper_wall:
                if lower_wall:
                    # need to correct for point already added at lower wall
                    upper_intersect_index += 1

                # this sfunc would put the points at the positions along the contour where the
                # grid would be orthogonal
                sfunc_orthogonal = contour.contourSfunc()

                # now make upper_intersect_index the index where the point at the wall is
                # check whether one of the points is already on the wall
                if calc_distance(contour[upper_intersect_index], upper_intersect) < self.atol:
                    pass
                elif calc_distance(contour[upper_intersect_index+1], upper_intersect) < self.atol:
                    upper_intersect_index = upper_intersect_index + 1
                else:
                    # otherwise insert a new point
                    upper_intersect_index += 1
                    contour.insert(upper_intersect_index, upper_intersect)

                # end point is now at the wall
                contour.endInd = upper_intersect_index

            self.sfunc_orthogonal_list.append(sfunc_orthogonal)

            contour.refine(width=self.user_options.refine_width)

    def distributePointsNonorthogonal(self):
        # regrid the contours (which all know where the wall is)
        for i_contour, contour in enumerate(self.contours):
            print('distributing points on contour:',
                    str(i_contour+1)+'/'+str(len(self.contours)), end = '\r')

            contour_is_separatrix = (numpy.abs((contour.psival -
                self.meshParent.equilibrium.psi_sep[0]) /
                self.meshParent.equilibrium.psi_sep[0]) < 1.e-9)

            def surface_vec(lower):
                if contour_is_separatrix:
                    if lower:
                        if self.equilibriumRegion.wallSurfaceAtStart is not None:
                            return self.equilibriumRegion.wallSurfaceAtStart
                        else:
                            # Use poloidal spacing on a separatrix contour
                            return None
                    else:
                        if self.equilibriumRegion.wallSurfaceAtEnd is not None:
                            return self.equilibriumRegion.wallSurfaceAtEnd
                        else:
                            # Use poloidal spacing on a separatrix contour
                            return None

                if i_contour == 0:
                    c_in = self.contours[0]
                else:
                    # contours are being changed, but start and end points are fixed so it is OK
                    # to use contours[i_contour-1] anyway
                    c_in = self.contours[i_contour-1]
                if i_contour == len(self.contours) - 1:
                    c_out = self.contours[i_contour]
                else:
                    c_out = self.contours[i_contour+1]
                if lower:
                    p_in = c_in[c_in.startInd]
                    p_out = c_out[c_out.startInd]
                else:
                    p_in = c_in[c_in.endInd]
                    p_out = c_out[c_out.endInd]
                return [p_out.R - p_in.R, p_out.Z - p_in.Z]

            if self.user_options.nonorthogonal_spacing_method == 'orthogonal':
                warnings.warn('\'orthogonal\' option is not currently compatible with '
                        'extending grid past targets')
                sfunc = self.sfunc_orthogonal_list[i_contour]
            elif self.user_options.nonorthogonal_spacing_method == 'fixed_poloidal':
                # this sfunc gives a fixed poloidal spacing at beginning and end of contours
                sfunc = self.equilibriumRegion.getSfuncFixedSpacing(
                        2*self.ny_noguards + 1, contour.totalDistance(), method='polynomial')
            elif self.user_options.nonorthogonal_spacing_method == 'poloidal_orthogonal_combined':
                sfunc = self.equilibriumRegion.combineSfuncs(contour,
                        self.sfunc_orthogonal_list[i_contour])
            elif self.user_options.nonorthogonal_spacing_method == 'fixed_perp_lower':
                sfunc = self.equilibriumRegion.getSfuncFixedPerpSpacing(
                        2*self.ny_noguards + 1, contour, surface_vec(True), True)
            elif self.user_options.nonorthogonal_spacing_method == 'fixed_perp_upper':
                sfunc = self.equilibriumRegion.getSfuncFixedPerpSpacing(
                        2*self.ny_noguards + 1, contour, surface_vec(False), False)
            elif self.user_options.nonorthogonal_spacing_method == 'perp_orthogonal_combined':
                sfunc = self.equilibriumRegion.combineSfuncs(contour,
                        self.sfunc_orthogonal_list[i_contour],
                        surface_vec(True), surface_vec(False))
            elif self.user_options.nonorthogonal_spacing_method == 'combined':
                if self.equilibriumRegion.wallSurfaceAtStart is not None:
                    # use poloidal spacing near a wall
                    surface_vec_lower = None
                else:
                    # use perp spacing
                    surface_vec_lower = surface_vec(True)
                if self.equilibriumRegion.wallSurfaceAtEnd is not None:
                    # use poloidal spacing near a wall
                    surface_vec_upper = None
                else:
                    # use perp spacing
                    surface_vec_upper = surface_vec(False)
                sfunc = self.equilibriumRegion.combineSfuncs(contour,
                        self.sfunc_orthogonal_list[i_contour],
                        surface_vec_lower, surface_vec_upper)
            else:
                raise ValueError('Unrecognized option \'' +
                        str(self.user_options.nonorthogonal_spacing_method)
                        + '\' for nonorthogonal poloidal spacing function')

            contour.regrid(2*self.ny_noguards + 1, sfunc=sfunc,
                    width=self.user_options.refine_width,
                    extend_lower=self.equilibriumRegion.extend_lower,
                    extend_upper=self.equilibriumRegion.extend_upper)

    def globalXInd(self, i):
        """
        Get the global x-index in the set of MeshRegions connected radially to this one of
        the point with local x-index i. Define so globalXInd=0 is at the primary separatrix
        """
        if self.radialIndex >= self.equilibriumRegion.separatrix_radial_index:
            # outside separatrix
            return i + sum(2*n for n in
                    self.equilibriumRegion.options.nx[self.equilibriumRegion.separatrix_radial_index:self.radialIndex])
        else:
            # inside separatrix
            return i - sum(2*n for n in
                    self.equilibriumRegion.options.nx[self.equilibriumRegion.separatrix_radial_index:self.radialIndex:-1])

    def fillRZ(self):
        """
        Fill the Rxy, Rxy_ylow and Zxy, Zxy_ylow arrays for this region

        xlow values include the outer point, after the final cell-centre grid point
        ylow values include the upper point, above the final cell-centre grid point
        """

        self.Rxy = MultiLocationArray(self.nx, self.ny)
        self.Zxy = MultiLocationArray(self.nx, self.ny)

        self.Rxy.centre = numpy.array([[p.R for p in contour[1::2]]
            for contour in self.contours[1::2]])

        self.Rxy.ylow = numpy.array([[p.R for p in contour[0::2]]
            for contour in self.contours[1::2]])

        self.Rxy.xlow = numpy.array([[p.R for p in contour[1::2]]
            for contour in self.contours[0::2]])

        self.Zxy.centre = numpy.array( [[p.Z for p in contour[1::2]]
            for contour in self.contours[1::2]])

        self.Zxy.ylow = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in self.contours[1::2]])

        self.Zxy.xlow = numpy.array([[p.Z for p in contour[1::2]]
            for contour in self.contours[0::2]])

        self.Rxy.corners = numpy.array( [[p.R for p in contour[0::2]]
            for contour in self.contours[0::2]])
        self.Zxy.corners = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in self.contours[0::2]])

        # Fix up the corner values at the X-points. Because the PsiContour have to start
        # slightly away from the X-point in order for the integrator to go in the right
        # direction, the points that should be at the X-point will be slighly displaced,
        # and will not be consistent between regions. So replace these points with the
        # X-point position instead.
        xpoint = self.equilibriumRegion.xPointsAtStart[self.radialIndex]
        if xpoint is not None:
            self.Rxy.corners[0,0] = xpoint.R
            self.Zxy.corners[0,0] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtStart[self.radialIndex+1]
        if xpoint is not None:
            self.Rxy.corners[-1,0] = xpoint.R
            self.Zxy.corners[-1,0] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtEnd[self.radialIndex]
        if xpoint is not None:
            self.Rxy.corners[0,-1] = xpoint.R
            self.Zxy.corners[0,-1] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtEnd[self.radialIndex+1]
        if xpoint is not None:
            self.Rxy.corners[-1,-1] = xpoint.R
            self.Zxy.corners[-1,-1] = xpoint.Z

    def getRZBoundary(self):
        # Upper value of ylow array logically overlaps with the lower value in the upper
        # neighbour. They should be close, but aren't guaranteed to be identical already
        # because they were taken from separate PsiContour obects. Use the value from the
        # upper neighbour to ensure consistency.
        # Also do similarly for the corner arrays.
        # Don't need to do this for the x-boundaries, because there the PsiContour
        # objects are shared between neighbouring regions.
        #
        # This needs to be a separate method from fillRZ() so that it can be called after
        # all regions have filled their Rxy and Zxy arrays.
        if self.connections['upper'] is not None:
            up = self.getNeighbour('upper')
            self.Rxy.ylow[:,-1] = up.Rxy.ylow[:,0]
            self.Zxy.ylow[:,-1] = up.Zxy.ylow[:,0]
            self.Rxy.corners[:,-1] = up.Rxy.corners[:,0]
            self.Zxy.corners[:,-1] = up.Zxy.corners[:,0]

    def geometry(self):
        """
        Calculate geometrical quantities for this region
        """

        self.psixy = self.meshParent.equilibrium.psi(self.Rxy, self.Zxy)

        self.dx = MultiLocationArray(self.nx, self.ny)
        self.dx.centre = (self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]
        self.dx.ylow = (self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]

        if self.psi_vals[0] > self.psi_vals[-1]:
            # x-coordinate is -psixy so x always increases radially across grid
            self.bpsign = -1.
            self.xcoord = -self.psixy
        else:
            self.bpsign = 1.
            self.xcoord = self.psixy

        self.dy = MultiLocationArray(self.nx, self.ny)
        self.dy.centre = self.meshParent.dy_scalar
        self.dy.ylow = self.meshParent.dy_scalar
        self.dy.xlow = self.meshParent.dy_scalar
        self.dy.corners = self.meshParent.dy_scalar

        self.Brxy = self.meshParent.equilibrium.Bp_R(self.Rxy, self.Zxy)
        self.Bzxy = self.meshParent.equilibrium.Bp_Z(self.Rxy, self.Zxy)
        self.Bpxy = numpy.sqrt(self.Brxy**2 + self.Bzxy**2)

        # determine direction - dot Bp with Grad(y) vector
        # evaluate in 'sol' at outer radial boundary
        Bp_dot_grady = (
            self.Brxy.centre[-1, self.ny//2]
            *(self.Rxy.centre[-1, self.ny//2 + 1] - self.Rxy.centre[-1, self.ny//2 - 1])
            + self.Bzxy.centre[-1, self.ny//2]
              *(self.Zxy.centre[-1, self.ny//2 + 1] - self.Zxy.centre[-1, self.ny//2 - 1]) )
        #print(self.myID, self.psi_vals[0], self.psi_vals[1], Bp_dot_grady)
        #print(self.Brxy.centre[-1, self.ny//2], self.Bzxy.centre[-1, self.ny//2],
        #        (self.Rxy.centre[-1, self.ny//2 - 1], self.Rxy.centre[-1, self.ny//2 +
        #            1]), (self.Zxy.centre[-1, self.ny//2 - 1], self.Zxy.centre[-1, self.ny//2 + 1]))
        if Bp_dot_grady < 0.:
            print("Poloidal field is in opposite direction to Grad(theta) -> Bp negative")
            self.Bpxy = -self.Bpxy
            if self.bpsign > 0.:
                raise ValueError("Sign of Bp should be negative? (note this check will "
                        "raise an exception when bpsign was correct if you only have a "
                        "private flux region)")
        else:
            if self.bpsign < 0.:
                raise ValueError("Sign of Bp should be negative? (note this check will "
                        "raise an exception when bpsign was correct if you only have a "
                        "private flux region)")

        # Get toroidal field from poloidal current function fpol
        self.Btxy = self.meshParent.equilibrium.fpol(self.psixy) / self.Rxy

        self.Bxy = numpy.sqrt(self.Bpxy**2 + self.Btxy**2)

        self.hy = self.calcHy()

        #if not self.user_options.orthogonal:
        #    # Calculate beta (angle between x and y coordinates), used for non-orthogonal grid
        #    # Also calculate radial grid spacing
        #    self.beta, self.hrad = self.calcBeta()

        #    # eta is the polodial non-orthogonality parameter
        #    self.eta = numpy.sin(self.beta)
        #else:
        #    self.beta.centre = 0.
        #    self.eta.centre = 0.

        # variation of toroidal angle with y following a field line. Called 'pitch' in
        # Hypnotoad1 because if y was the poloidal angle then dphidy would be the pitch
        # angle.
        self.dphidy = self.hy * self.Btxy / (self.Bpxy * self.Rxy)

        self.d2phidxdy = MultiLocationArray(self.nx, self.ny)
        self.d2phidxdy.centre = self.DDX_L2C(self.dphidy.xlow)
        self.d2phidxdy.ylow = self.DDX_L2C(self.dphidy.corners, ylow=True)

    def calcMetric(self):
        """
        Calculate the metrics using geometrical information calculated in geometry().
        Needs to be a separate method as zShift can only be calculated when calcZShift()
        has been called on the MeshRegion at the beginning of the y-group. To ensure this,
        call geometry() on all regions first, then calcMetric on all regions.
        """
        if not self.user_options.shiftedmetric:
            # To implement the shiftedmetric==False case, would have to define a
            # consistent zShift=0 location for all regions, for example in the style of
            # Hypnotoad1. This could be done by a particular implementation of 'Mesh'
            # (e.g. 'BoutMesh') before calling this method. Needs to be a particular
            # implementation which knows about the topology of the grid - still not clear
            # it is possible to do consistently, e.g. in private-flux regions.
            raise ValueError("'shiftedmetric == False' not handled at present.\n"
                             "Cannot make grid for field-aligned toroidal coordinates "
                             "without making zShift consistent between all regions. "
                             "Don't know how to do this in general, and haven't "
                             "implemented the Hypnototoad1-style solution as it does not "
                             "seem consistent in the private-flux region, or the "
                             "inner-SOL of a double-null configuration.")
            # integrated shear
            self.sinty = MultiLocationArray(self.nx, self.ny)
            self.sinty.centre = self.DDX_L2C(self.zShift.xlow)
            self.sinty.ylow = self.DDX_L2C(self.zShift.corners, ylow=True)
            I = self.sinty
        else:
            # Zero integrated shear, because the coordinate system is defined locally to
            # each value of y, and defined to have no shear.
            # In this case zShift only needs to be defined consistently *along* each field
            # line - don't need to be able to take radial (x-direction) derivatives. This
            # means different (radial) regions can use different locations for where
            # zShift=0.
            I = MultiLocationArray(self.nx, self.ny).zero()

        self.g11 = (self.Rxy*self.Bpxy)**2
        self.g22 = 1./self.hy**2
        self.g33 = I*self.g11 + (self.dphidy/self.hy)**2 + 1./self.Rxy**2
        self.g12 = MultiLocationArray(self.nx, self.ny).zero()
        self.g13 = -I*self.g11
        self.g23 = -self.dphidy/self.hy**2

        self.J = self.hy / self.Bpxy

        self.g_11 = 1./self.g11 + (I*self.Rxy)**2
        self.g_22 = self.hy**2 + (self.Rxy/self.dphidy)**2
        self.g_33 = self.Rxy**2
        self.g_12 = self.Rxy**2*self.dphidy*I
        self.g_13 = self.Rxy**2*I
        self.g_23 = self.dphidy*self.Rxy**2

        # check Jacobian is OK
        one_over_sqrt_g = self.bpsign*1./numpy.sqrt(self.g11*self.g22*self.g33
                + 2.*self.g12*self.g13*self.g23 - self.g11*self.g23**2
                - self.g22*self.g13**2 - self.g33*self.g12**2)
        # ignore grid points at X-points as J should diverge there (as Bp->0)
        if self.equilibriumRegion.xPointsAtStart[self.radialIndex] is not None:
            one_over_sqrt_g.corners[0, 0] = self.J.corners[0,0]
        if self.equilibriumRegion.xPointsAtStart[self.radialIndex + 1] is not None:
            one_over_sqrt_g.corners[-1, 0] = self.J.corners[-1,0]
        if self.equilibriumRegion.xPointsAtEnd[self.radialIndex] is not None:
            one_over_sqrt_g.corners[0, -1] = self.J.corners[0, -1]
        if self.equilibriumRegion.xPointsAtEnd[self.radialIndex + 1] is not None:
            one_over_sqrt_g.corners[-1, -1] = self.J.corners[-1, -1]

        check = numpy.abs(self.J - one_over_sqrt_g) / numpy.abs(self.J) < self.user_options.geometry_rtol
        def ploterror(location):
            if location == 'centre':
                thisJ = self.J.centre
                this_one_over_sqrt_g = one_over_sqrt_g.centre
            elif location == 'ylow':
                thisJ = self.J.ylow
                this_one_over_sqrt_g = one_over_sqrt_g.ylow
            elif location == 'xlow':
                thisJ = self.J.xlow
                this_one_over_sqrt_g = one_over_sqrt_g.xlow
            elif location == 'corners':
                thisJ = self.J.corners
                this_one_over_sqrt_g = one_over_sqrt_g.corners
            else:
                raise ValueError('wrong location argument: '+str(location))
            print('rtol = ' + str(self.user_options.geometry_rtol))
            from matplotlib import pyplot
            pyplot.figure(location)
            pyplot.subplot(221)
            pyplot.pcolor(thisJ)
            pyplot.title('J')
            pyplot.colorbar()
            pyplot.subplot(222)
            pyplot.pcolor(this_one_over_sqrt_g)
            pyplot.title('1/sqrt(g)')
            pyplot.colorbar()
            pyplot.subplot(223)
            pyplot.pcolor(thisJ - this_one_over_sqrt_g)
            pyplot.title('abs difference')
            pyplot.colorbar()
            pyplot.subplot(224)
            pyplot.pcolor((thisJ - this_one_over_sqrt_g)/thisJ)
            pyplot.title('rel difference')
            pyplot.colorbar()
            pyplot.show()
        if not numpy.all(check.centre):
            ploterror('centre')
        if not numpy.all(check.ylow):
            ploterror('ylow')
        if not numpy.all(check.xlow):
            ploterror('xlow')
        if not numpy.all(check.corners):
            ploterror('corners')
        assert numpy.all(check.centre), 'Jacobian should be consistent with 1/sqrt(det(g)) calculated from the metric tensor'
        assert numpy.all(check.ylow), 'Jacobian should be consistent with 1/sqrt(det(g)) calculated from the metric tensor'
        assert numpy.all(check.xlow), 'Jacobian should be consistent with 1/sqrt(det(g)) calculated from the metric tensor'
        assert numpy.all(check.corners), 'Jacobian should be consistent with 1/sqrt(det(g)) calculated from the metric tensor'

    def calcHy(self, ylow=False):
        # hy = |Grad(theta)|
        # hy = dtheta/ds at constant psi, phi when psi and theta are orthogonal
        # approx dtheta/sqrt((R(j+1/2)-R(j-1/2))**2 + (Z(j+1/2)-Z(j-1/2)**2)
        if self.user_options.orthogonal:
            warnings.warn('need to check that this is correct for non-orthogonal grids')

        # get positions at j+/-0.5
        R = self.Rxy.ylow
        Z = self.Zxy.ylow

        hy = MultiLocationArray(self.nx, self.ny)
        hy.centre = numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2) / self.dy.centre

        # for hthe_ylow, need R, Z values from below the lower face of this region and
        # above the upper face
        R = numpy.zeros([self.nx, self.ny + 2])
        R[:,1:-1] = self.Rxy.centre
        Z = numpy.zeros([self.nx, self.ny + 2])
        Z[:,1:-1] = self.Zxy.centre
        if self.connections['lower'] is not None:
            R[:,0] = self.getNeighbour('lower').Rxy.centre[:, -1]
            Z[:,0] = self.getNeighbour('lower').Zxy.centre[:, -1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:,0] = 2.*self.Rxy.ylow[:,0] - self.Rxy.centre[:,0]
            Z[:,0] = 2.*self.Zxy.ylow[:,0] - self.Zxy.centre[:,0]
        if self.connections['upper'] is not None:
            R[:,-1] = self.getNeighbour('upper').Rxy.centre[:,0]
            Z[:,-1] = self.getNeighbour('upper').Zxy.centre[:,0]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # value will never even be passed to Mesh, since we stagger to YLOW)
            R[:,-1] = 2.*self.Rxy.ylow[:,-1] - self.Rxy.centre[:,-1]
            Z[:,-1] = 2.*self.Zxy.ylow[:,-1] - self.Zxy.centre[:,-1]

        hy.ylow =  numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2) / self.dy.ylow

        # for hthe_xlow, need R, Z values from the cell corners
        R = self.Rxy.corners
        Z = self.Zxy.corners

        hy.xlow = numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2) / self.dy.xlow

        # for hthecorners, need R, Z values from xlow
        R = numpy.zeros([self.nx+1, self.ny+2])
        Z = numpy.zeros([self.nx+1, self.ny+2])
        R[:,1:-1] = self.Rxy.xlow
        Z[:,1:-1] = self.Zxy.xlow
        if self.connections['lower'] is not None:
            R[:,0] = self.getNeighbour('lower').Rxy.xlow[:,-1]
            Z[:,0] = self.getNeighbour('lower').Zxy.xlow[:,-1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:,0] = 2.*self.Rxy.corners[:,0] - self.Rxy.xlow[:,0]
            Z[:,0] = 2.*self.Zxy.corners[:,0] - self.Zxy.xlow[:,0]
        if self.connections['upper'] is not None:
            R[:,-1] = self.getNeighbour('upper').Rxy.xlow[:,0]
            Z[:,-1] = self.getNeighbour('upper').Zxy.xlow[:,0]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # will not even be stored in Mesh, since we stagger to YLOW)
            R[:,-1] = 2.*self.Rxy.corners[:,-1] - self.Rxy.xlow[:,-1]
            Z[:,-1] = 2.*self.Zxy.corners[:,-1] - self.Zxy.xlow[:,-1]

        hy.corners = numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2) / self.dy.corners

        return hy

    def calcBeta(self, ylow=False):
        """
        beta is the angle between x and y coordinates, used for non-orthogonal grid.
        Also calculate radial grid spacing, hrad
        """
        #raise ValueError("non-orthogonal grids not calculated yet")
        warnings.warn("non-orthogonal grids not calculated yet")

        #if not ylow:
        #    # need to multiply f_R and f_Z by bpsign because we want the radially-outward
        #    # vector perpendicular to psi contours, and if bpsign is negative then psi
        #    # increases inward instead of outward so (f_R,f_Z) would be in the opposite
        #    # direction
        #    # Actually want the angle of the vector in the y-direction, i.e. (f_Z,-f_R)
        #    angle_grad_psi = numpy.arctan2(
        #            self.bpsign*self.meshParent.equilibrium.f_Z(self.Rxy, self.Zxy),
        #            -self.bpsign*self.meshParent.equilibrium.f_R(self.Rxy, self.Zxy))

        #    R = numpy.zeros([self.nx + 1, self.ny])
        #    R[:-1, :] = self.Rxy_xlow
        #    R[-1, :] = self.Rxy_extra_outer
        #    Z = numpy.zeros([self.nx + 1, self.ny])
        #    Z[:-1 :] = self.Zxy_xlow
        #    Z[-1, :] = self.Zxy_extra_outer
        #    # could calculate radial grid spacing - is it ever needed?
        #    hrad = numpy.sqrt((R[1:,:] - R[:-1,:])**2 + (Z[1:,:] - Z[:-1,:])**2)

        #    dR = R[1:,:] - R[:-1,:]
        #    dZ = Z[1:,:] - Z[:-1,:]
        #    angle_dr = numpy.arctan2(dR, dZ)
        #else:
        #    # need to multiply f_R and f_Z by bpsign because we want the radially-outward
        #    # vector perpendicular to psi contours, and if bpsign is negative then psi
        #    # increases inward instead of outward so (f_R,f_Z) would be in the opposite
        #    # direction
        #    # Actually want the angle of the vector in the y-direction, i.e. (f_Z,-f_R)
        #    angle_grad_psi = numpy.arctan2(
        #            self.bpsign*self.meshParent.equilibrium.f_Z(self.Rxy_ylow, self.Zxy_ylow),
        #            -self.bpsign*self.meshParent.equilibrium.f_R(self.Rxy_ylow, self.Zxy_ylow))

        #    # could calculate radial grid spacing - is it ever needed?
        #    ## for hrad at ylow, can use Rcorners and Zcorners
        #    hrad = numpy.sqrt((self.Rcorners[1:,:-1] - self.Rcorners[:-1,:-1])**2 +
        #                      (self.Zcorners[1:,:-1] - self.Zcorners[:-1,:-1])**2)

        #    dR = self.Rcorners[1:,:-1] - self.Rcorners[:-1,:-1]
        #    dZ = self.Zcorners[1:,:-1] - self.Zcorners[:-1,:-1]
        #    angle_dr = numpy.arctan2(dR, dZ)

        #return (angle_grad_psi - angle_dr - numpy.pi/2.), hrad

    def calcZShift(self):
        """
        Calculate zShift by integrating dphidy in y.
        """
        # Integrate using all available points - centre+ylow or xlow+corner.
        # Integrate from lower boundary on open field lines.
        # Integrate from lower side of MeshRegion with yGroupIndex=0 on closed field
        # lines.
        # Use trapezoid rule. If int_f = \int f dy
        # int_f.centre[j] = int_f.centre[j-1]
        #                   + 0.5*(f.centre[j-1] + f.ylow[j]) * (0.5*dy.centre[j-1])
        #                   + 0.5*(f.ylow[j] + f.centre[j]) * (0.5*dy.centre[j])
        #                 = i_centre[j-1] + i_ylow_upper[j]
        #                   + i_ylow_lower[j] + i_centre[j]
        # At the moment dy is a constant, but we allow for future changes with variable
        # grid-spacing in y. The cell-centre points should be half way between the
        # cell-face points, so the distance between centre[j-1] and ylow[j] is
        # 0.5*dy[j-1], and the distace between ylow[j] and centre[j] is 0.5*dy[j]
        #
        # Also
        # int_f.ylow[j] = int_f.ylow[j-1]
        #                 + 0.5*(f.ylow[j-1] + f.centre[j-1]) * (0.5*dy.centre[j-1])
        #                 + 0.5*(f.centre[j-1] + f.ylow[j]) * (0.5*dy.centre[j-1])
        #               = i_ylow_lower[j-1] + i_centre[j-1]
        #                 + i_centre[j-1] + i_ylow_upper[j]

        # Cannot just test 'connections['lower'] is not None' because periodic regions
        # always have a lower connection - requires us to give a yGroupIndex to each
        # region when creating the groups.
        if self.yGroupIndex is not 0:
            return None

        region = self
        region.zShift = MultiLocationArray(region.nx, region.ny)
        while True:
            # calculate integral for field lines with centre and ylow points
            i_centre = 0.25*numpy.cumsum(region.dphidy.centre * region.dy.centre, axis=1)
            i_ylow_lower = 0.25*numpy.cumsum(region.dphidy.ylow[:, :-1] \
                           * region.dy.centre, axis=1)
            i_ylow_upper = 0.25*numpy.cumsum(region.dphidy.ylow[:, 1:] \
                           * region.dy.centre, axis=1)

            region.zShift.centre[:,0] = region.zShift.ylow[:, 0] \
                                        + i_ylow_lower[:, 0] + i_centre[:, 0]
            region.zShift.centre[:,1:] = region.zShift.ylow[:, 0, numpy.newaxis] \
                                         + i_centre[:, :-1] + i_ylow_upper[:, :-1] \
                                         + i_ylow_lower[:, 1:] + i_centre[:, 1:]

            region.zShift.ylow[:, 1:] = region.zShift.ylow[:, 0, numpy.newaxis] \
                                        + i_ylow_lower + 2.*i_centre \
                                        + i_ylow_upper

            # repeat for field lines with xlow and corner points
            i_xlow = 0.25*numpy.cumsum(region.dphidy.xlow * region.dy.xlow, axis=1)
            i_corners_lower = 0.25*numpy.cumsum(region.dphidy.corners[:, :-1] \
                              * region.dy.xlow, axis=1)
            i_corners_upper = 0.25*numpy.cumsum(region.dphidy.corners[:, 1:] \
                              * region.dy.xlow, axis=1)

            region.zShift.xlow[:,0] = region.zShift.corners[:, 0] \
                                        + i_corners_lower[:, 0] + i_xlow[:, 0]
            region.zShift.xlow[:,1:] = region.zShift.corners[:, 0, numpy.newaxis] \
                                         + i_xlow[:, :-1] + i_corners_upper[:, :-1] \
                                         + i_corners_lower[:, 1:] + i_xlow[:, 1:]

            region.zShift.corners[:, 1:] = region.zShift.corners[:, 0, numpy.newaxis] \
                                        + i_corners_lower + 2.*i_xlow \
                                        + i_corners_upper

            next_region = region.getNeighbour('upper')
            if next_region is None:
                break
            else:
                next_region.zShift = MultiLocationArray(next_region.nx, next_region.ny)
                next_region.zShift.ylow[:, 0] = region.zShift.ylow[:, -1]
                next_region.zShift.corners[:, 0] = region.zShift.corners[:, -1]
                region = next_region

    def getNeighbour(self, face):
        if self.connections[face] is None:
            return None
        else:
            return self.meshParent.regions[self.connections[face]]

    def DDX_L2C(self, f, ylow=False):
        # x-derivative at x-centre points, calculated by 2nd order central difference from
        # x-staggered points.
        # Assume the 'xlow' quantity f has nx+1 values and includes the outer point after
        # the last cell-centre grid point.
        assert f.shape[0] == self.nx + 1, 'input field f should be at xlow or corner, so should have x-size nx+1'

        if not ylow:
            dx = self.dx.centre
        else:
            dx = self.dx.ylow
        result = (f[1:, :] - f[:-1, :]) / dx
        return result

class Mesh:
    """
    Mesh represented by a collection of connected MeshRegion objects
    """
    def __init__(self, equilibrium):
        self.user_options = equilibrium.user_options
        self.options = equilibrium.options

        # Set some global parameters for PsiContours and FineContours
        # Convert self.user_options to a dict so we can use it to set the values in
        # FineContour.options using 'push'
        PsiContour.options = PsiContour.options.push(dict(self.user_options))
        FineContour.options = FineContour.options.push(dict(self.user_options))

        self.equilibrium = equilibrium

        # Generate MeshRegion object for each section of the mesh
        self.regions = {}

        # Make consecutive numbering scheme for regions
        regionlist = []
        self.region_lookup = {}
        for reg_name,eq_reg in equilibrium.regions.items():
            for i in range(eq_reg.nSegments):
                region_number = len(regionlist)
                regionlist.append((reg_name, i))
                self.region_lookup[(reg_name, i)] = region_number

        # Get connections between regions
        self.connections = {}
        for region_id,(eq_reg,i) in enumerate(regionlist):
            self.connections[region_id] = {}
            region = equilibrium.regions[eq_reg]
            c = region.connections[i]
            for key, val in c.items():
                if val is not None:
                    self.connections[region_id][key] = self.region_lookup[val]
                else:
                    self.connections[region_id][key] = None

        self.makeRegions()

    def makeRegions(self):
        for eq_region in self.equilibrium.regions.values():
            for i in range(eq_region.nSegments):
                region_id = self.region_lookup[(eq_region.name,i)]
                eq_region_with_boundaries = eq_region.getRegridded(radialIndex=i,
                        width=self.user_options.refine_width)
                self.regions[region_id] = MeshRegion(self, region_id,
                        eq_region_with_boundaries, self.connections[region_id], i)

        # create groups that connect in x
        self.x_groups = []
        region_set = set(self.regions.values())
        while region_set:
            for region in region_set:
                if region.connections['inner'] is None:
                    break
            group = []
            while True:
                group.append(region)
                region_set.remove(region)
                region = region.getNeighbour('outer')
                if region is None or group.count(region) > 0:
                    # reached boundary or have all regions in a periodic group
                    break
            self.x_groups.append(group)

        # create groups that connect in y
        self.y_groups = []
        region_set = set(self.regions.values())
        while region_set:
            for region in region_set:
                if region.connections['lower'] is None:
                    break
                # note, if no region with connections['lower']=None is found, then some
                # arbitrary region will be 'region' after this loop. This is OK, as this
                # region must be part of a periodic group, which we will handle.
            group = []
            while True:
                assert region.yGroupIndex == None, 'region should not have been added to any yGroup before'
                region.yGroupIndex = len(group)
                group.append(region)
                region_set.remove(region)
                region = region.getNeighbour('upper')
                if region is None or group.count(region) > 0:
                    # reached boundary or have all regions in a periodic group
                    break
            self.y_groups.append(group)

    def redistributePoints(self, **kwargs):
        warnings.warn('It is not recommended to use Mesh.redistributePoints() for '
                '\'production\' output. Suggest saving the final settings to a .yaml '
                'file and creating the \'production\' grid non-interactively to ensure '
                'reproducibility.')

        self.user_options.set(**kwargs)

        assert not self.user_options.orthogonal, 'redistributePoints would do nothing for an orthogonal grid.'
        for region in self.regions.values():
            print('redistributing', region.name)
            region.equilibriumRegion.setupOptions(force=True)
            region.distributePointsNonorthogonal()

    def geometry(self):
        """
        Calculate geometrical quantities for BOUT++
        """
        for region in self.regions.values():
            region.fillRZ()
        for region in self.regions.values():
            region.getRZBoundary()
        for region in self.regions.values():
            region.geometry()
        for region in self.regions.values():
            region.calcZShift()
        for region in self.regions.values():
            region.calcMetric()

    def plotPoints(self, xlow=False, ylow=False, corners=False):
        from matplotlib import pyplot
        from cycler import cycle

        colors = cycle(pyplot.rcParams['axes.prop_cycle'].by_key()['color'])
        for region in self.regions.values():
            c = next(colors)
            pyplot.scatter(region.Rxy.centre, region.Zxy.centre, marker='x', c=c,
                    label=region.myID)
            if xlow:
                pyplot.scatter(region.Rxy.xlow, region.Zxy.xlow, marker='1', c=c)
            if ylow:
                pyplot.scatter(region.Rxy.ylow, region.Zxy.ylow, marker='2', c=c)
            if corners:
                pyplot.scatter(region.Rxy.corners, region.Zxy.corners, marker='+', c=c)
        pyplot.legend()

    def plotPotential(self, *args, **kwargs):
        """
        Plot the flux function psi. Passes through to self.equilibrium.plotPotential.
        """
        return self.equilibrium.plotPotential(*args, **kwargs)

def followPerpendicular(f_R, f_Z, p0, A0, Avals, rtol=2.e-8, atol=1.e-8):
    """
    Follow a line perpendicular to Bp from point p0 until magnetic potential A_target is
    reached.
    """
    f = lambda A,x: (f_R(x[0], x[1]), f_Z(x[0], x[1]))
    Arange = (A0, Avals[-1])
    solution = solve_ivp(f, Arange, tuple(p0), t_eval=Avals, rtol=rtol, atol=atol,
            vectorized=True)

    return [Point2D(*p) for p in solution.y.T]

class BoutMesh(Mesh):
    """
    Mesh quantities to be written to a grid file for BOUT++

    Requires that the MeshRegion members fit together into a global logically-rectangular
    Mesh, with the topology assumed by BOUT++ (allowing complexity up to
    disconnected-double-null).

    For compatibility with BOUT++, the regions in the OrderedDict equilibrium.regions must
    be in the order: inner_lower_divertor, inner_core, inner_upper_divertor,
    outer_upper_divertor, outer_core, outer_lower_divertor. This ensures the correct
    positioning in the global logically rectangular grid. Regions are allowed to not be
    present (if they would have size 0).
    """
    def __init__(self, equilibrium, *args, **kwargs):

        super().__init__(equilibrium, *args, **kwargs)

        # nx, ny both include boundary guard cells
        eq_region0 = next(iter(self.equilibrium.regions.values()))
        self.nx = sum(eq_region0.options.nx)

        self.ny = sum(r.ny(0) for r in self.equilibrium.regions.values())

        self.ny_noguards = sum(r.ny_noguards for r in self.equilibrium.regions.values())

        # Keep ranges of global indices for each region, separately from the MeshRegions,
        # because we don't want MeshRegion objects to depend on global indices
        assert all([r.options.nx == eq_region0.options.nx for r in self.equilibrium.regions.values()]), 'all regions should have same set of x-grid sizes to be compatible with a global, logically-rectangular grid'
        x_sizes = [0] + list(eq_region0.options.nx)
        x_startinds = numpy.cumsum(x_sizes)
        x_regions = tuple(slice(x_startinds[i], x_startinds[i+1], None)
                     for i in range(len(x_startinds)-1))
        y_total = 0
        y_regions = {}
        for regname, region in self.equilibrium.regions.items():
            # all segments must have the same ny, i.e. same number of y-boundary guard
            # cells
            this_ny = region.ny(0)
            assert all(region.ny(i) == this_ny for i in range(region.nSegments)), 'all radial segments in an equilibrium-region must have the same ny (i.e. same number of boundary guard cells) to be compatible with a global, logically-rectangular grid'

            y_total_new = y_total + this_ny
            reg_slice = slice(y_total, y_total_new, None)
            y_total = y_total_new
            y_regions[regname] = reg_slice

        self.region_indices = {}
        for reg_name in self.equilibrium.regions:
            for i in range(len(x_regions)):
                self.region_indices[self.region_lookup[(reg_name, i)]] = numpy.index_exp[
                        x_regions[i], y_regions[reg_name]]

        # constant spacing in y for now
        self.dy_scalar = 2.*numpy.pi / self.ny_noguards

    def geometry(self):
        # Call geometry() method of base class
        super().geometry()

        def addFromRegion(f, f_region, regionID):
            if f_region._centre_array is not None:
                f.centre[self.region_indices[regionID]] = f_region.centre
            if f_region._xlow_array is not None:
                f.xlow[self.region_indices[regionID]] = f_region.xlow[:-1,:]
            if f_region._ylow_array is not None:
                f.ylow[self.region_indices[regionID]] = f_region.ylow[:,:-1]
            if f_region._corners_array is not None:
                f.corners[self.region_indices[regionID]] = f_region.corners[:-1,:-1]

        self.Rxy = MultiLocationArray(self.nx, self.ny)
        self.Zxy = MultiLocationArray(self.nx, self.ny)
        self.psixy = MultiLocationArray(self.nx, self.ny)
        self.dx = MultiLocationArray(self.nx, self.ny)
        self.dy = MultiLocationArray(self.nx, self.ny)
        self.Brxy = MultiLocationArray(self.nx, self.ny)
        self.Bzxy = MultiLocationArray(self.nx, self.ny)
        self.Bpxy = MultiLocationArray(self.nx, self.ny)
        self.Btxy = MultiLocationArray(self.nx, self.ny)
        self.Bxy = MultiLocationArray(self.nx, self.ny)
        self.hy = MultiLocationArray(self.nx, self.ny)
        #if not self.user_options.orthogonal:
        #    self.beta = MultiLocationArray(self.nx, self.ny)
        #    self.eta = MultiLocationArray(self.nx, self.ny)
        self.dphidy = MultiLocationArray(self.nx, self.ny)
        self.d2phidxdy = MultiLocationArray(self.nx, self.ny)
        self.zShift = MultiLocationArray(self.nx, self.ny)
        if not self.user_options.shiftedmetric:
            self.sinty = MultiLocationArray(self.nx, self.ny)
        self.g11 = MultiLocationArray(self.nx, self.ny)
        self.g22 = MultiLocationArray(self.nx, self.ny)
        self.g33 = MultiLocationArray(self.nx, self.ny)
        self.g12 = MultiLocationArray(self.nx, self.ny)
        self.g13 = MultiLocationArray(self.nx, self.ny)
        self.g23 = MultiLocationArray(self.nx, self.ny)
        self.J = MultiLocationArray(self.nx, self.ny)
        self.g_11 = MultiLocationArray(self.nx, self.ny)
        self.g_22 = MultiLocationArray(self.nx, self.ny)
        self.g_33 = MultiLocationArray(self.nx, self.ny)
        self.g_12 = MultiLocationArray(self.nx, self.ny)
        self.g_13 = MultiLocationArray(self.nx, self.ny)
        self.g_23 = MultiLocationArray(self.nx, self.ny)

        for region in self.regions.values():
            addFromRegion(self.Rxy, region.Rxy, region.myID)
            addFromRegion(self.Zxy, region.Zxy, region.myID)
            addFromRegion(self.psixy, region.psixy, region.myID)
            addFromRegion(self.dx, region.dx, region.myID)
            addFromRegion(self.dy, region.dy, region.myID)
            addFromRegion(self.Brxy, region.Brxy, region.myID)
            addFromRegion(self.Bzxy, region.Bzxy, region.myID)
            addFromRegion(self.Bpxy, region.Bpxy, region.myID)
            addFromRegion(self.Btxy, region.Btxy, region.myID)
            addFromRegion(self.Bxy, region.Bxy, region.myID)
            addFromRegion(self.hy, region.hy, region.myID)
            #if not self.user_options.orthogonal:
            #    addFromRegion(self.beta, region.beta, region.myID)
            #    addFromRegion(self.eta, region.eta, region.myID)
            addFromRegion(self.dphidy, region.dphidy, region.myID)
            addFromRegion(self.d2phidxdy, region.d2phidxdy, region.myID)
            addFromRegion(self.zShift, region.zShift, region.myID)
            if not self.user_options.shiftedmetric:
                addFromRegion(self.sinty, region.sinty, region.myID)
            addFromRegion(self.g11, region.g11, region.myID)
            addFromRegion(self.g22, region.g22, region.myID)
            addFromRegion(self.g33, region.g33, region.myID)
            addFromRegion(self.g12, region.g12, region.myID)
            addFromRegion(self.g13, region.g13, region.myID)
            addFromRegion(self.g23, region.g23, region.myID)
            addFromRegion(self.J, region.J, region.myID)
            addFromRegion(self.g_11, region.g_11, region.myID)
            addFromRegion(self.g_22, region.g_22, region.myID)
            addFromRegion(self.g_33, region.g_33, region.myID)
            addFromRegion(self.g_12, region.g_12, region.myID)
            addFromRegion(self.g_13, region.g_13, region.myID)
            addFromRegion(self.g_23, region.g_23, region.myID)

    def writeArray(self, name, array, f):
        f.write(name, array.centre)
        f.write(name+'_ylow', array.ylow[:, :-1])

    def writeGridfile(self, filename):
        from boututils.datafile import DataFile

        with DataFile(filename, create=True, format='NETCDF4') as f:
            f.write('nx', self.nx)
            # ny for BOUT++ excludes boundary guard cells
            f.write('ny', self.ny_noguards)
            f.write('y_boundary_guards', self.user_options.y_boundary_guards)
            self.writeArray('Rxy', self.Rxy, f)
            self.writeArray('Zxy', self.Zxy, f)
            self.writeArray('psixy', self.psixy, f)
            self.writeArray('dx', self.dx, f)
            self.writeArray('dy', self.dy, f)
            self.writeArray('Bpxy', self.Bpxy, f)
            self.writeArray('Btxy', self.Btxy, f)
            self.writeArray('Bxy', self.Bxy, f)
            self.writeArray('hthe', self.hy, f)
            #if not self.user_options.orthogonal:
            #    self.writeArray('beta', self.beta, f)
            #    self.writeArray('eta', self.eta, f)
            self.writeArray('dphidy', self.dphidy, f)
            self.writeArray('zShift', self.dphidy, f)

            # Haven't checked this is exactly the quantity needed by BOUT++...
            # ShiftTorsion is only used in Curl operator - Curl is rarely used.
            self.writeArray('ShiftTorsion', self.d2phidxdy, f)

            # I think IntShiftTorsion should be the same as sinty in Hypnotoad1.
            # IntShiftTorsion should never be used. It is only for some 'BOUT-06 style
            # differencing'. IntShiftTorsion is not written by Hypnotoad1, so don't write
            # here. /JTO 19/5/2019
            if not self.user_options.shiftedmetric:
                self.writeArray('sinty', self.sinty, f)

            self.writeArray('g11', self.g11, f)
            self.writeArray('g22', self.g22, f)
            self.writeArray('g33', self.g33, f)
            self.writeArray('g12', self.g12, f)
            self.writeArray('g13', self.g13, f)
            self.writeArray('g23', self.g23, f)
            self.writeArray('J', self.J, f)
            self.writeArray('g_11', self.g_11, f)
            self.writeArray('g_22', self.g_22, f)
            self.writeArray('g_33', self.g_33, f)
            self.writeArray('g_12', self.g_12, f)
            self.writeArray('g_13', self.g_13, f)
            self.writeArray('g_23', self.g_23, f)

            f.write('hypnotoad_inputs', self.equilibrium._getOptionsAsString())

    def plot2D(self, f, title=None):
        from matplotlib import pyplot

        try:
            vmin = f.min()
            vmax = f.max()
            if vmin == vmax:
                vmin -= 0.1
                vmax += 0.1

            for region, indices in zip(self.regions.values(), self.region_indices.values()):
                pyplot.pcolor(region.Rxy.corners, region.Zxy.corners, f[indices],
                              vmin=vmin, vmax=vmax)

            pyplot.colorbar()
        except NameError:
            raise NameError('Some variable has not been defined yet: have you called Mesh.geometry()?')

    def saveOptions(self, filename='hypnotoad_options.yaml'):
        self.equilibrium.saveOptions(filename)
