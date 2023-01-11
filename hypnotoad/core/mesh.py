# Copyright 2019 J.T. Omotani
#
# Contact John Omotani john.omotani@ukaea.uk
#
# This file is part of Hypnotoad 2.
#
# Hypnotoad 2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Hypnotoad 2 is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# Hypnotoad 2.  If not, see <http://www.gnu.org/licenses/>.

"""
Classes to handle Meshes and geometrical quantities for generating BOUT++ grids
"""

from copy import deepcopy
import re
import sys
import warnings
import yaml

import numpy
from optionsfactory import OptionsFactory, WithMeta
from optionsfactory.checks import (
    NoneType,
    is_non_negative,
    is_positive,
)
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import interp1d

from boututils.boutarray import BoutArray
from boututils.run_wrapper import shell_safe

from .equilibrium import (
    calc_distance,
    Equilibrium,
    EquilibriumRegion,
    Point2D,
    PsiContour,
)
from .multilocationarray import MultiLocationArray
from ..utils.utils import module_versions_formatted
from ..utils.parallel_map import ParallelMap
from ..__version__ import get_versions


class MeshRegion:
    """
    Collection of :class:`PsiContour <hypnotoad.core.equilibrium.PsiContour>` objects
    representing a logically rectangular sub-region of the grid.

    Each side of the ``MeshRegion`` is either a grid boundary or connects to one other
    ``MeshRegion``. The ``MeshRegion`` includes points on the boundary, whose positions
    are shared with the boundary points of the neighbouring ``MeshRegion``, if there is
    one.

    Includes methods for calculating R-Z positions, magnetic field, and geometric
    quantities (for the standard BOUT++ locally field aligned coordinate system) on the
    points in the ``MeshRegion``.

    Developers, see :ref:`developer/meshregion:MeshRegion notes`.
    """

    user_options_factory = OptionsFactory(
        EquilibriumRegion.user_options_factory,
        shiftedmetric=WithMeta(
            True,
            doc="Is grid generated for paralleltransform=ShiftedMetric?",
            value_type=bool,
        ),
        curvature_type=WithMeta(
            "curl(b/B)",
            doc="Expression used to calculate curvature operator 'bxcv'",
            value_type=str,
            allowed=["curl(b/B)", "curl(b/B) with x-y derivatives", "bxkappa"],
        ),
        curvature_smoothing=WithMeta(
            None,
            doc=(
                "Smoothing for components of bxcv output: None - no smoothing; "
                "'smoothnl' - non-linear smoothing using the algorithm from IDL "
                "hypnotoad."
            ),
            value_type=[str, NoneType],
            allowed=[None, "smoothnl"],
        ),
        follow_perpendicular_rtol=WithMeta(
            2.0e-8,
            doc="Relative tolerance for following Grad(psi)",
            value_type=[float, int],
            check_all=is_non_negative,
        ),
        follow_perpendicular_atol=WithMeta(
            1.0e-8,
            doc="Absolute tolerance for following Grad(psi)",
            value_type=[float, int],
            check_all=is_non_negative,
        ),
        geometry_rtol=WithMeta(
            1.0e-10,
            doc=(
                "Tolerance for checking identities on calculated geometrical quantities "
                "(for example the Jacobian)"
            ),
            value_type=[float, int],
            check_all=is_positive,
        ),
        cap_Bp_ylow_xpoint=WithMeta(
            False,
            doc=(
                "Fudge to get rid of minimum in Bpxy.ylow field, because it can cause "
                "large spikes in some metric coefficients, which may cause numerical "
                "problems in simulations"
            ),
            value_type=bool,
        ),
        wall_point_exclude_radius=WithMeta(
            1.0e-3,
            doc=(
                "When adding the wall point to a contour, this is the radius within "
                "which existing contour points are removed as too close to the wall "
                "point (points being too close might lead to bad interpolation). "
                "This parameter should not usually need to be changed."
            ),
            value_type=float,
            check_all=is_positive,
        ),
    )

    def __init__(
        self,
        meshParent,
        myID,
        equilibriumRegion,
        connections,
        radialIndex,
        settings,
        parallel_map,
    ):

        self.user_options = self.user_options_factory.create(settings)

        self.name = equilibriumRegion.name + "(" + str(radialIndex) + ")"
        print("creating region", myID, "-", self.name, flush=True)

        # the Mesh object that owns this MeshRegion
        self.meshParent = meshParent

        # ID that Mesh uses to keep track of its MeshRegions
        self.myID = myID

        # EquilibriumRegion representing the segment associated with this region
        self.equilibriumRegion = equilibriumRegion.copy()

        # check settings were consistent
        for key, value in self.equilibriumRegion.user_options.items():
            if value != self.user_options[key]:
                raise ValueError(
                    f"{key} not consistent in MeshRegion user_options and "
                    f"EquilibriumRegion user_options"
                )

        # sizes of the grid in this MeshRegion, include boundary guard cells
        self.nx = self.equilibriumRegion.nx[radialIndex]
        self.ny = self.equilibriumRegion.ny(radialIndex)
        self.ny_noguards = self.equilibriumRegion.ny_noguards

        # psi values for radial grid
        self.psi_vals = numpy.array(self.equilibriumRegion.psi_vals[radialIndex])
        if len(self.psi_vals) != 2 * self.nx + 1:
            raise ValueError("should be a psi value for each radial point")

        # Dictionary that specifies whether a boundary is connected to another region or
        # is an actual boundary
        self.connections = connections

        # Number of this region, counting radially outward
        self.radialIndex = radialIndex

        # Number of this region in its y-group
        self.yGroupIndex = None

        # Absolute tolerance for checking if two points are the same
        self.atol = 1.0e-7

        # Create ParallelMap object for parallelising loops
        self.parallel_map = parallel_map

        # get points in this region
        self.contours = []
        if self.radialIndex < self.equilibriumRegion.separatrix_radial_index:
            # region is inside separatrix, so need to follow line from the last psi_val
            # to the first
            temp_psi_vals = self.psi_vals[::-1]
        else:
            temp_psi_vals = self.psi_vals

        # Make vector along grad(psi) at start of equilibriumRegion Here we assume that
        # the equilibriumRegion at a separatrix at the beginning and end, but not
        # necessarily in between.  This is to handle disconnected double null
        start_point = self.equilibriumRegion[self.equilibriumRegion.startInd]
        start_psi = self.equilibriumRegion.psi(*start_point)

        # set sign of step in psi towards this region from primary separatrix at start of
        # region
        if temp_psi_vals[-1] - start_psi > 0:
            start_psi_sep_plus_delta = (
                start_psi
                + self.equilibriumRegion.equilibrium.poloidal_spacing_delta_psi
            )
        else:
            start_psi_sep_plus_delta = (
                start_psi
                - self.equilibriumRegion.equilibrium.poloidal_spacing_delta_psi
            )

        vec_points = followPerpendicular(
            0,
            start_point,
            start_psi,
            f_R=self.meshParent.equilibrium.f_R,
            f_Z=self.meshParent.equilibrium.f_Z,
            psivals=[start_psi, start_psi_sep_plus_delta],
            rtol=self.user_options.follow_perpendicular_rtol,
            atol=self.user_options.follow_perpendicular_atol,
        )
        self.equilibriumRegion.gradPsiSurfaceAtStart = (
            vec_points[1].as_ndarray() - vec_points[0].as_ndarray()
        )

        # Make vector along grad(psi) at end of equilibriumRegion
        end_point = self.equilibriumRegion[self.equilibriumRegion.endInd]
        end_psi = self.equilibriumRegion.psi(*end_point)

        # set sign of step in psi towards this region from primary separatrix at end of
        # region
        if temp_psi_vals[-1] - end_psi > 0:
            end_psi_sep_plus_delta = (
                end_psi + self.equilibriumRegion.equilibrium.poloidal_spacing_delta_psi
            )
        else:
            end_psi_sep_plus_delta = (
                end_psi - self.equilibriumRegion.equilibrium.poloidal_spacing_delta_psi
            )

        vec_points = followPerpendicular(
            -1,
            end_point,
            end_psi,
            f_R=self.meshParent.equilibrium.f_R,
            f_Z=self.meshParent.equilibrium.f_Z,
            psivals=[end_psi, end_psi_sep_plus_delta],
            rtol=self.user_options.follow_perpendicular_rtol,
            atol=self.user_options.follow_perpendicular_atol,
        )
        self.equilibriumRegion.gradPsiSurfaceAtEnd = (
            vec_points[1].as_ndarray() - vec_points[0].as_ndarray()
        )

        # Calculate the angles for perp_d_lower/perp_d_upper corresponding to
        # d_lower/d_upper on the separatrix contour
        # Use self.equilibriumRegion.fine_contour for the vector along the separatrix
        # because then the vector will not change when the grid resolution changes
        fine_contour = self.equilibriumRegion.get_fine_contour(
            psi=self.equilibriumRegion.psi
        )
        if self.equilibriumRegion.wallSurfaceAtStart is None:
            # lower end
            unit_vec_separatrix = (
                fine_contour.positions[fine_contour.startInd + 1, :]
                - fine_contour.positions[fine_contour.startInd, :]
            )
            unit_vec_separatrix /= numpy.sqrt(numpy.sum(unit_vec_separatrix**2))
            unit_vec_surface = self.equilibriumRegion.gradPsiSurfaceAtStart
            unit_vec_surface /= numpy.sqrt(numpy.sum(unit_vec_surface**2))
            cos_angle = numpy.sum(unit_vec_separatrix * unit_vec_surface)
            # this gives abs(sin_angle), but that's OK because we only want the magnitude
            # to calculate perp_d
            self.equilibriumRegion.sin_angle_at_start = numpy.sqrt(1.0 - cos_angle**2)
        if self.equilibriumRegion.wallSurfaceAtEnd is None:
            # upper end
            unit_vec_separatrix = (
                fine_contour.positions[fine_contour.endInd - 1, :]
                - fine_contour.positions[fine_contour.endInd, :]
            )
            unit_vec_separatrix /= numpy.sqrt(numpy.sum(unit_vec_separatrix**2))
            unit_vec_surface = self.equilibriumRegion.gradPsiSurfaceAtEnd
            unit_vec_surface /= numpy.sqrt(numpy.sum(unit_vec_surface**2))
            cos_angle = numpy.sum(unit_vec_separatrix * unit_vec_surface)
            # this gives abs(sin_angle), but that's OK because we only want the magnitude
            # to calculate perp_d
            self.equilibriumRegion.sin_angle_at_end = numpy.sqrt(1.0 - cos_angle**2)

        perp_points_list = self.parallel_map(
            followPerpendicular,
            zip(
                range(len(self.equilibriumRegion)),
                self.equilibriumRegion,
                [self.equilibriumRegion.psi(*p) for p in self.equilibriumRegion],
            ),
            psivals=temp_psi_vals,
            rtol=self.user_options.follow_perpendicular_rtol,
            atol=self.user_options.follow_perpendicular_atol,
        )
        if self.radialIndex < self.equilibriumRegion.separatrix_radial_index:
            for perp_points in perp_points_list:
                perp_points.reverse()

        for i, point in enumerate(perp_points_list[0]):
            self.contours.append(
                self.equilibriumRegion.newContourFromSelf(
                    points=[point], psival=self.psi_vals[i]
                )
            )
            self.contours[i].global_xind = self.globalXInd(i)
        for perp_points in perp_points_list[1:]:
            for i, point in enumerate(perp_points):
                self.contours[i].append(point)

        # refine the contours to make sure they are at exactly the right psi-value
        self.contours = self.parallel_map(
            PsiContour.refine,
            ((c,) for c in self.contours),
            width=self.user_options.refine_width,
        )

        if not self.user_options.orthogonal:
            self.addPointAtWallToContours()
            self.distributePointsNonorthogonal()

    def addPointAtWallToContours(self):
        # maximum number of times to extend the contour when it has not yet hit the wall
        max_extend = 100

        # should the contour intersect a wall at the lower end?
        lower_wall = self.connections["lower"] is None

        # should the contour intersect a wall at the upper end?
        upper_wall = self.connections["upper"] is None

        map_result = self.parallel_map(
            _find_intersection,
            enumerate(self.contours),
            lower_wall=lower_wall,
            upper_wall=upper_wall,
            max_extend=max_extend,
            atol=self.atol,
            refine_width=self.user_options.refine_width,
        )

        self.contours = [x[0] for x in map_result]
        intersect_info_list = [x[1:] for x in map_result]

        # sfunc_orthogonal functions created after contour has been extended past wall
        # (if necessary) but before adding the wall point to the contour (as adding this
        # point makes the spacing of points on the contour not-smooth) and adjusted for
        # the change in distance after redefining startInd to be at the wall
        self.sfunc_orthogonal_list = [
            contour.contourSfunc(psi=self.equilibriumRegion.psi)
            for contour in self.contours
        ]

        for i, contour, intersect_info in zip(
            range(len(self.contours)), self.contours, intersect_info_list
        ):
            (
                lower_intersect_index,
                lower_intersect,
                upper_intersect_index,
                upper_intersect,
            ) = intersect_info
            # now add points on the wall(s) to the contour
            if lower_wall:
                # need to construct a new sfunc which gives distance from the wall, not
                # the distance from the original startInd

                # now make lower_intersect_index the index where the point at the wall
                # is.
                # check whether one of the points is close to the wall point and remove
                # if it is - use a pretty loose check because points apart from
                # the wall point are about to be redistributed anyway, so does
                # not matter if we remove one.
                if (
                    calc_distance(contour[lower_intersect_index], lower_intersect)
                    < self.user_options.wall_point_exclude_radius
                ):
                    contour.replace(lower_intersect_index, lower_intersect)
                elif (
                    calc_distance(contour[lower_intersect_index + 1], lower_intersect)
                    < self.user_options.wall_point_exclude_radius
                ):
                    lower_intersect_index = lower_intersect_index + 1
                    contour.replace(lower_intersect_index, lower_intersect)
                else:
                    # otherwise insert a new point
                    lower_intersect_index += 1
                    contour.insert(lower_intersect_index, lower_intersect)
                    if upper_wall and upper_intersect_index >= 0:
                        # need to correct for point already added at lower wall
                        upper_intersect_index += 1

                def correct_sfunc_orthogonal_and_set_startInd(
                    contour, sfunc_orthogonal_original
                ):
                    # Shift sfunc_orthogonal so that the distance at the endInd
                    # value is exactly that calculated using the *new*
                    # FineContour
                    original_total_distance = sfunc_orthogonal_original(
                        float(contour.endInd)
                    )

                    # start contour from the wall
                    contour.startInd = lower_intersect_index
                    # Force reset of contour._fine_contour, etc. even if
                    # lower_intersect_index happens to be the same as previous startInd
                    contour._reset_cached()

                    # Calculate total distance using a new FineContour, created at this
                    # point
                    new_total_distance = contour.totalDistance(
                        psi=self.equilibriumRegion.psi
                    )

                    return (
                        lambda i: sfunc_orthogonal_original(i)
                        - original_total_distance
                        + new_total_distance
                    )

                # original sfuncs_orthogonal would put the points at the positions
                # along the contour where the grid would be orthogonal need to
                # correct sfunc_orthogonal for the distance between the point at
                # the lower wall and the original start-point
                self.sfunc_orthogonal_list[
                    i
                ] = correct_sfunc_orthogonal_and_set_startInd(
                    contour,
                    self.sfunc_orthogonal_list[i],
                )

            if upper_wall:
                # now make upper_intersect_index the index where the point at the wall is
                # check whether one of the points is close to the wall point and remove
                # if it is - use a pretty loose check because points apart from the
                # wall point are about to be redistributed anyway, so does not matter if
                # we remove one.
                if (
                    calc_distance(contour[upper_intersect_index], upper_intersect)
                    < self.user_options.wall_point_exclude_radius
                ):
                    contour.replace(upper_intersect_index, upper_intersect)
                elif (
                    calc_distance(contour[upper_intersect_index + 1], upper_intersect)
                    < self.user_options.wall_point_exclude_radius
                ):
                    upper_intersect_index = upper_intersect_index + 1
                    contour.replace(upper_intersect_index, upper_intersect)
                else:
                    # otherwise insert a new point
                    contour.insert(upper_intersect_index + 1, upper_intersect)
                    if upper_intersect_index >= 0:
                        # If upper_intersect_index is negative, then it still points to
                        # the correct element after insertion
                        upper_intersect_index += 1

                # end point is now at the wall
                contour.endInd = upper_intersect_index
                # Force reset of contour._fine_contour, etc. even if
                # upper_intersect_index happens to be the same as previous endInd
                contour._reset_cached()
                # Note: no need to adjust sfunc_orthogonal[i] because the startInd (at
                # an X-point) is still at distance=0.0

        self.contours = self.parallel_map(_refine_extend, enumerate(self.contours))

    def distributePointsNonorthogonal(self, nonorthogonal_settings=None):
        if nonorthogonal_settings is not None:
            self.equilibriumRegion.resetNonorthogonalOptions(nonorthogonal_settings)

        def surface_vec(i_contour, contour, lower):
            psi_sep = self.meshParent.equilibrium.psi_sep[0]
            contour_is_separatrix = (
                numpy.abs((contour.psival - psi_sep) / psi_sep) < 1.0e-9
            )

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
                c_in = self.contours[i_contour - 1]
            if i_contour == len(self.contours) - 1:
                c_out = self.contours[i_contour]
            else:
                c_out = self.contours[i_contour + 1]
            if lower:
                p_in = c_in[c_in.startInd]
                p_out = c_out[c_out.startInd]
            else:
                p_in = c_in[c_in.endInd]
                p_out = c_out[c_out.endInd]
            return [p_out.R - p_in.R, p_out.Z - p_in.Z]

        # surface_vecs_lower = [
        #    surface_vec(i, c, True) for i, c in enumerate(self.contours)
        # ]
        # surface_vecs_upper = [
        #    surface_vec(i, c, False) for i, c in enumerate(self.contours)
        # ]

        def get_sfunc(i_contour, contour, sfunc_orthogonal):
            surface_vec_lower = surface_vec(i_contour, contour, True)
            surface_vec_upper = surface_vec(i_contour, contour, False)
            if (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "orthogonal"
            ):
                warnings.warn(
                    "'orthogonal' option is not currently compatible with "
                    "extending grid past targets"
                )
                return sfunc_orthogonal
            elif (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "fixed_poloidal"
            ):
                # this sfunc gives a fixed poloidal spacing at beginning and end of
                # contours
                return self.equilibriumRegion.getSfuncFixedSpacing(
                    2 * self.ny_noguards + 1,
                    contour.totalDistance(),
                    method="monotonic",
                )
            elif (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "poloidal_orthogonal_combined"
            ):
                return self.equilibriumRegion.combineSfuncs(contour, sfunc_orthogonal)
            elif (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "fixed_perp_lower"
            ):
                return self.equilibriumRegion.getSfuncFixedPerpSpacing(
                    2 * self.ny_noguards + 1, contour, surface_vec_lower, True
                )
            elif (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "fixed_perp_upper"
            ):
                return self.equilibriumRegion.getSfuncFixedPerpSpacing(
                    2 * self.ny_noguards + 1, contour, surface_vec_upper, False
                )
            elif (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "perp_orthogonal_combined"
            ):
                return self.equilibriumRegion.combineSfuncs(
                    contour,
                    sfunc_orthogonal,
                    surface_vec_lower,
                    surface_vec_upper,
                )
            elif (
                self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method
                == "combined"
            ):
                if self.equilibriumRegion.wallSurfaceAtStart is not None:
                    # use poloidal spacing near a wall
                    surface_vec_lower = None
                else:
                    # use perp spacing
                    surface_vec_lower = surface_vec_lower
                if self.equilibriumRegion.wallSurfaceAtEnd is not None:
                    # use poloidal spacing near a wall
                    surface_vec_upper = None
                else:
                    # use perp spacing
                    surface_vec_upper = surface_vec_upper
                return self.equilibriumRegion.combineSfuncs(
                    contour,
                    sfunc_orthogonal,
                    surface_vec_lower,
                    surface_vec_upper,
                )
            else:
                raise ValueError(
                    "Unrecognized option '"
                    + str(
                        self.equilibriumRegion.nonorthogonal_options.nonorthogonal_spacing_method  # noqa: E501
                    )
                    + "' for nonorthogonal poloidal spacing function"
                )

        # sfunc_list = [
        #    get_sfunc(i, c, sfunc_orth)
        #    for i, c, sfunc_orth in zip(
        #        range(len(self.contours)), self.contours, self.sfunc_orthogonal_list
        #    )
        # ]

        # regrid the contours (which all know where the wall is)
        for i, c, sfunc_orth in zip(
            range(len(self.contours)), self.contours, self.sfunc_orthogonal_list
        ):
            print("Regridding:", i, flush=True, end="\r")
            c.regrid(
                2 * self.ny_noguards + 1,
                psi=self.equilibriumRegion.psi,
                sfunc=get_sfunc(i, c, sfunc_orth),
                extend_lower=self.equilibriumRegion.extend_lower,
                extend_upper=self.equilibriumRegion.extend_upper,
                refine=False,
            )

        # refine the regridded contours in parallel
        self.contours = self.parallel_map(
            PsiContour.refine,
            ((c,) for c in self.contours),
        )

    def globalXInd(self, i):
        """
        Get the global x-index in the set of MeshRegions connected radially to this one
        of the point with local x-index i. Define so globalXInd=0 is at the primary
        separatrix
        """
        if self.radialIndex >= self.equilibriumRegion.separatrix_radial_index:
            # outside separatrix
            return i + sum(
                2 * n
                for n in self.equilibriumRegion.nx[
                    self.equilibriumRegion.separatrix_radial_index : self.radialIndex
                ]
            )
        else:
            # inside separatrix
            return i - sum(
                2 * n
                for n in self.equilibriumRegion.nx[
                    self.equilibriumRegion.separatrix_radial_index : self.radialIndex : -1  # noqa: E501
                ]
            )

    def fillRZ(self):
        """
        Fill the Rxy, Rxy_ylow and Zxy, Zxy_ylow arrays for this region

        xlow values include the outer point, after the final cell-centre grid point
        ylow values include the upper point, above the final cell-centre grid point
        """

        self.Rxy = MultiLocationArray(self.nx, self.ny)
        self.Zxy = MultiLocationArray(self.nx, self.ny)

        self.Rxy.centre = numpy.array(
            [[p.R for p in contour[1::2]] for contour in self.contours[1::2]]
        )

        self.Rxy.ylow = numpy.array(
            [[p.R for p in contour[0::2]] for contour in self.contours[1::2]]
        )

        self.Rxy.xlow = numpy.array(
            [[p.R for p in contour[1::2]] for contour in self.contours[0::2]]
        )

        self.Zxy.centre = numpy.array(
            [[p.Z for p in contour[1::2]] for contour in self.contours[1::2]]
        )

        self.Zxy.ylow = numpy.array(
            [[p.Z for p in contour[0::2]] for contour in self.contours[1::2]]
        )

        self.Zxy.xlow = numpy.array(
            [[p.Z for p in contour[1::2]] for contour in self.contours[0::2]]
        )

        self.Rxy.corners = numpy.array(
            [[p.R for p in contour[0::2]] for contour in self.contours[0::2]]
        )
        self.Zxy.corners = numpy.array(
            [[p.Z for p in contour[0::2]] for contour in self.contours[0::2]]
        )

        # Fix up the corner values at the X-points. Because the PsiContour have to start
        # slightly away from the X-point in order for the integrator to go in the right
        # direction, the points that should be at the X-point will be slighly displaced,
        # and will not be consistent between regions. So replace these points with the
        # X-point position instead.
        xpoint = self.equilibriumRegion.xPointsAtStart[self.radialIndex]
        if xpoint is not None:
            self.Rxy.corners[0, 0] = xpoint.R
            self.Zxy.corners[0, 0] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtStart[self.radialIndex + 1]
        if xpoint is not None:
            self.Rxy.corners[-1, 0] = xpoint.R
            self.Zxy.corners[-1, 0] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtEnd[self.radialIndex]
        if xpoint is not None:
            self.Rxy.corners[0, -1] = xpoint.R
            self.Zxy.corners[0, -1] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtEnd[self.radialIndex + 1]
        if xpoint is not None:
            self.Rxy.corners[-1, -1] = xpoint.R
            self.Zxy.corners[-1, -1] = xpoint.Z

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
        if self.connections["upper"] is not None:
            up = self.getNeighbour("upper")
            self.Rxy.ylow[:, -1] = up.Rxy.ylow[:, 0]
            self.Zxy.ylow[:, -1] = up.Zxy.ylow[:, 0]
            self.Rxy.corners[:, -1] = up.Rxy.corners[:, 0]
            self.Zxy.corners[:, -1] = up.Zxy.corners[:, 0]

    def calcDistances(self):
        """
        Calculate the distances for all PsiContours in the region
        """
        if self.user_options.orthogonal:
            # Distances already calculated in non-orthogonal case
            self.contours = self.parallel_map(
                _calc_contour_distance, enumerate(self.contours)
            )

    def geometry1(self):
        """
        Calculate geometrical quantities for this region
        """
        self.psixy = self.meshParent.equilibrium.psi(self.Rxy, self.Zxy)

        self.dx = MultiLocationArray(self.nx, self.ny)
        self.dx.centre = (self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]
        self.dx.ylow = (self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]

        if self.psi_vals[0] > self.psi_vals[-1]:
            # x-coordinate is -psixy so x always increases radially across grid
            self.bpsign = -1.0
            self.xcoord = -self.psixy
        else:
            self.bpsign = 1.0
            self.xcoord = self.psixy

        self.dy = MultiLocationArray(self.nx, self.ny)
        self.dy.centre = self.meshParent.dy_scalar
        self.dy.ylow = self.meshParent.dy_scalar
        self.dy.xlow = self.meshParent.dy_scalar
        self.dy.corners = self.meshParent.dy_scalar

        self.Brxy = self.meshParent.equilibrium.Bp_R(self.Rxy, self.Zxy)
        self.Bzxy = self.meshParent.equilibrium.Bp_Z(self.Rxy, self.Zxy)
        self.Bpxy = numpy.sqrt(self.Brxy**2 + self.Bzxy**2)

        self.calcPoloidalDistance()

        if hasattr(
            self.meshParent.equilibrium.regions[self.equilibriumRegion.name], "pressure"
        ):
            self.pressure = self.meshParent.equilibrium.regions[
                self.equilibriumRegion.name
            ].pressure(self.psixy)

        # determine direction - dot Bp with Grad(y) vector
        # evaluate in 'sol' at outer radial boundary
        Bp_dot_grady = self.Brxy.centre[-1, self.ny // 2] * (
            self.Rxy.centre[-1, self.ny // 2 + 1]
            - self.Rxy.centre[-1, self.ny // 2 - 1]
        ) + self.Bzxy.centre[-1, self.ny // 2] * (
            self.Zxy.centre[-1, self.ny // 2 + 1]
            - self.Zxy.centre[-1, self.ny // 2 - 1]
        )
        # print(self.myID, self.psi_vals[0], self.psi_vals[1], Bp_dot_grady)
        # print(self.Brxy.centre[-1, self.ny//2], self.Bzxy.centre[-1, self.ny//2],
        #        (self.Rxy.centre[-1, self.ny//2 - 1], self.Rxy.centre[-1, self.ny//2 +
        #            1]), (self.Zxy.centre[-1, self.ny//2 - 1],
        #                  self.Zxy.centre[-1, self.ny//2 + 1]))
        if Bp_dot_grady < 0.0:
            print(
                "Poloidal field is in opposite direction to Grad(theta) -> Bp negative"
            )
            self.Bpxy = -self.Bpxy
            if self.bpsign > 0.0:
                raise ValueError(
                    "Sign of Bp should be negative? (note this check will raise an "
                    "exception when bpsign was correct if you only have a private flux "
                    "region)"
                )
        else:
            if self.bpsign < 0.0:
                raise ValueError(
                    "Sign of Bp should be negative? (note this check will raise an "
                    "exception when bpsign was correct if you only have a private flux "
                    "region)"
                )

        # Get toroidal field from poloidal current function fpol
        self.Btxy = self.meshParent.equilibrium.fpol(self.psixy) / self.Rxy

        self.Bxy = numpy.sqrt(self.Bpxy**2 + self.Btxy**2)

    def geometry2(self):
        """
        Continuation of geometry1(), but needs neighbours to have calculated Bp so called
        after all regions call geometry1().
        """

        if self.user_options.cap_Bp_ylow_xpoint:
            # (By default do _not_ do this)
            # Get rid of minimum in Bpxy.ylow field, because it can cause large spikes in
            # some metric coefficients, which may cause numerical problems in simulations
            self.capBpYlowXpoint()

        self.hy = self.calcHy()

        if not self.user_options.orthogonal:
            # Calculate beta (angle between e_x and Grad(x), also the angle between e_y
            # and Grad(y)), used for non-orthogonal grid
            self.calcBeta()
        # else:
        #    self.beta.centre = 0.
        #    self.eta.centre = 0.

        # variation of toroidal angle with y following a field line. Called 'pitch' in
        # Hypnotoad1 because if y was the poloidal angle then dphidy would be the pitch
        # angle.
        self.dphidy = self.hy * self.Btxy / (self.Bpxy * self.Rxy)

    def capBpYlowXpoint(self):
        if self.equilibriumRegion.xPointsAtStart[self.radialIndex] is not None:
            # Choose a minumum Bp as the average of the two values of Bpxy.centre
            # nearest to the X-point
            Bp_min = min(
                self.Bpxy.centre[0, 0], self.getNeighbour("lower").Bpxy.centre[0, -1]
            )
            for i in range(self.nx):
                if self.Bpxy.ylow[i, 0] < Bp_min:
                    self.Bpxy.ylow[i, 0] = Bp_min
                else:
                    break
        if self.equilibriumRegion.xPointsAtStart[self.radialIndex + 1] is not None:
            # Choose a minumum Bp as the average of the two values of Bpxy.centre
            # nearest to the X-point
            Bp_min = min(
                self.Bpxy.centre[-1, 0], self.getNeighbour("lower").Bpxy.centre[-1, -1]
            )
            for i in range(self.nx):
                if self.Bpxy.ylow[-i - 1, 0] < Bp_min:
                    self.Bpxy.ylow[-i - 1, 0] = Bp_min
                else:
                    break
        if self.equilibriumRegion.xPointsAtEnd[self.radialIndex] is not None:
            # Choose a minumum Bp as the average of the two values of Bpxy.centre
            # nearest to the X-point
            Bp_min = min(
                self.Bpxy.centre[0, -1], self.getNeighbour("upper").Bpxy.centre[0, 0]
            )
            for i in range(self.nx):
                if self.Bpxy.ylow[i, -1] < Bp_min:
                    self.Bpxy.ylow[i, -1] = Bp_min
                else:
                    break
        if self.equilibriumRegion.xPointsAtEnd[self.radialIndex + 1] is not None:
            # Choose a minumum Bp as the average of the two values of Bpxy.centre
            # nearest to the X-point
            Bp_min = min(
                self.Bpxy.centre[-1, -1], self.getNeighbour("upper").Bpxy.centre[-1, 0]
            )
            for i in range(self.nx):
                if self.Bpxy.ylow[-i - 1, -1] < Bp_min:
                    self.Bpxy.ylow[-i - 1, -1] = Bp_min
                else:
                    break

    def calcMetric(self):
        """
        Calculate the metrics using geometrical information calculated in geometry1() and
        geometry2().
        Needs to be a separate method as zShift can only be calculated when calcZShift()
        has been called on the MeshRegion at the beginning of the y-group. To ensure
        this, call geometry1() and geometry2() on all regions first, then calcMetric on
        all regions.
        """
        if not self.user_options.shiftedmetric:
            # To implement the shiftedmetric==False case, would have to define a
            # consistent zShift=0 location for all regions, for example in the style of
            # Hypnotoad1. This could be done by a particular implementation of 'Mesh'
            # (e.g. 'BoutMesh') before calling this method. Needs to be a particular
            # implementation which knows about the topology of the grid - still not clear
            # it is possible to do consistently, e.g. in private-flux regions.
            raise ValueError(
                "'shiftedmetric == False' not handled at present.\n"
                "Cannot make grid for field-aligned toroidal coordinates without making "
                "zShift consistent between all regions. Don't know how to do this in "
                "general, and haven't implemented the Hypnototoad1-style solution as it "
                "does not seem consistent in the private-flux region, or the inner-SOL "
                "of a double-null configuration."
            )
            # integrated shear
            self.sinty = self.DDX("zShift")
            self.I = self.sinty
        else:
            # Zero integrated shear, because the coordinate system is defined locally to
            # each value of y, and defined to have no shear.
            # In this case zShift only needs to be defined consistently *along* each
            # field line - don't need to be able to take radial (x-direction)
            # derivatives. This means different (radial) regions can use different
            # locations for where zShift=0.
            self.I = MultiLocationArray(self.nx, self.ny).zero()

        # Here ShiftTorsion = d2phidxdy
        # Haven't checked this is exactly the quantity needed by BOUT++...
        # ShiftTorsion is only used in Curl operator - Curl is rarely used.
        self.ShiftTorsion = self.DDX("#dphidy")

        if self.user_options.orthogonal:
            self.g11 = (self.Rxy * self.Bpxy) ** 2
            self.g22 = 1.0 / self.hy**2
            self.g33 = (
                self.I * self.g11 + (self.dphidy / self.hy) ** 2 + 1.0 / self.Rxy**2
            )
            self.g12 = MultiLocationArray(self.nx, self.ny).zero()
            self.g13 = -self.I * self.g11
            self.g23 = -self.dphidy / self.hy**2

            self.J = self.hy / self.Bpxy

            self.g_11 = 1.0 / self.g11 + (self.I * self.Rxy) ** 2
            self.g_22 = self.hy**2 + (self.Rxy * self.dphidy) ** 2
            self.g_33 = self.Rxy**2
            self.g_12 = self.Rxy**2 * self.dphidy * self.I
            self.g_13 = self.Rxy**2 * self.I
            self.g_23 = self.dphidy * self.Rxy**2
        else:
            self.g11 = (self.Rxy * self.Bpxy) ** 2
            self.g22 = 1.0 / (self.hy * self.cosBeta) ** 2
            self.g33 = (
                1.0 / self.Rxy**2
                + (self.Rxy * self.Bpxy * self.I) ** 2
                + (self.dphidy / (self.hy * self.cosBeta)) ** 2
                + 2.0
                * self.Rxy
                * self.Bpxy
                * self.I
                * self.dphidy
                * self.tanBeta
                / self.hy
            )
            self.g12 = self.Rxy * numpy.abs(self.Bpxy) * self.tanBeta / self.hy
            self.g13 = (
                -self.Rxy * self.Bpxy * self.dphidy * self.tanBeta / self.hy
                - self.I * (self.Rxy * self.Bpxy) ** 2
            )
            self.g23 = (
                -self.bpsign * self.dphidy / (self.hy * self.cosBeta) ** 2
                - self.Rxy * numpy.abs(self.Bpxy) * self.I * self.tanBeta / self.hy
            )

            self.J = self.hy / self.Bpxy

            self.g_11 = (
                1.0 / (self.Rxy * self.Bpxy * self.cosBeta) ** 2
                + (self.I * self.Rxy) ** 2
            )
            self.g_22 = self.hy**2 + (self.dphidy * self.Rxy) ** 2
            self.g_33 = self.Rxy**2
            self.g_12 = (
                self.bpsign * self.I * self.dphidy * self.Rxy**2
                - self.hy * self.tanBeta / (self.Rxy * numpy.abs(self.Bpxy))
            )
            self.g_13 = self.I * self.Rxy**2
            self.g_23 = self.bpsign * self.dphidy * self.Rxy**2

        # check Jacobian is OK
        Jcheck = (
            self.bpsign
            * 1.0
            / numpy.sqrt(
                self.g11 * self.g22 * self.g33
                + 2.0 * self.g12 * self.g13 * self.g23
                - self.g11 * self.g23**2
                - self.g22 * self.g13**2
                - self.g33 * self.g12**2
            )
        )
        # ignore grid points at X-points as J should diverge there (as Bp->0)
        if Jcheck._corners_array is not None:
            # If Jcheck was not calculated at the corners location, no check is needed.
            # Skip these fixes because they would initialise Jcheck.corners, which we do
            # not want to do.
            if self.equilibriumRegion.xPointsAtStart[self.radialIndex] is not None:
                Jcheck.corners[0, 0] = self.J.corners[0, 0]
            if self.equilibriumRegion.xPointsAtStart[self.radialIndex + 1] is not None:
                Jcheck.corners[-1, 0] = self.J.corners[-1, 0]
            if self.equilibriumRegion.xPointsAtEnd[self.radialIndex] is not None:
                Jcheck.corners[0, -1] = self.J.corners[0, -1]
            if self.equilibriumRegion.xPointsAtEnd[self.radialIndex + 1] is not None:
                Jcheck.corners[-1, -1] = self.J.corners[-1, -1]

        check = (
            numpy.abs(self.J - Jcheck) / numpy.abs(self.J)
            < self.user_options.geometry_rtol
        )

        def ploterror(location):
            if location == "centre":
                thisJ = self.J.centre
                this_one_over_sqrt_g = Jcheck.centre
            elif location == "ylow":
                thisJ = self.J.ylow
                this_one_over_sqrt_g = Jcheck.ylow
            elif location == "xlow":
                thisJ = self.J.xlow
                this_one_over_sqrt_g = Jcheck.xlow
            elif location == "corners":
                thisJ = self.J.corners
                this_one_over_sqrt_g = Jcheck.corners
            else:
                raise ValueError("wrong location argument: " + str(location))
            print(self.name, "rtol = " + str(self.user_options.geometry_rtol))
            from matplotlib import pyplot

            pyplot.figure(location)
            pyplot.subplot(221)
            pyplot.pcolor(thisJ)
            pyplot.title("J")
            pyplot.colorbar()
            pyplot.subplot(222)
            pyplot.pcolor(this_one_over_sqrt_g)
            pyplot.title("1/sqrt(g)")
            pyplot.colorbar()
            pyplot.subplot(223)
            pyplot.pcolor(thisJ - this_one_over_sqrt_g)
            pyplot.title("abs difference")
            pyplot.colorbar()
            pyplot.subplot(224)
            pyplot.pcolor((thisJ - this_one_over_sqrt_g) / thisJ)
            pyplot.title("rel difference")
            pyplot.colorbar()
            pyplot.show()

        if not numpy.all(check.centre):
            ploterror("centre")
            raise ValueError(
                f"Geometry: Jacobian at centre should be consistent with "
                f"1/sqrt(det(g)) calculated from the metric tensor. If the plot "
                f"looks OK, you may want to increase the value of "
                f"geometry_rtol={self.user_options.geometry_rtol}"
            )

        if not numpy.all(check.ylow):
            ploterror("ylow")
            raise ValueError(
                f"Geometry: Jacobian at ylow should be consistent with "
                f"1/sqrt(det(g)) calculated from the metric tensor. If the plot "
                f"looks OK, you may want to increase the value of "
                f"geometry_rtol={self.user_options.geometry_rtol}"
            )

        if check._xlow_array is not None:
            if not numpy.all(check.xlow):
                ploterror("xlow")
                raise ValueError(
                    f"Geometry: Jacobian at xlow should be consistent with "
                    f"1/sqrt(det(g)) calculated from the metric tensor. If the "
                    f"plot looks OK, you may want to increase the value of "
                    f"geometry_rtol={self.user_options.geometry_rtol}"
                )
        if check._corners_array is not None:
            if not numpy.all(check.corners):
                ploterror("corners")
                raise ValueError(
                    f"Geometry: Jacobian at corners should be consistent with "
                    f"1/sqrt(det(g)) calculated from the metric tensor. If the "
                    f"plot looks OK, you may want to increase the value of "
                    f"geometry_rtol={self.user_options.geometry_rtol}"
                )

        # curvature terms
        self.calc_curvature()

    def calc_curvature(self):
        """
        Calculate curvature components. Note that curl_bOverB_x, curl_bOverB_y, and
        curl_bOverB_z are contravariant components - Curl(b/B)^x, Curl(b/B)^y, and
        Curl(b/B)^z - despite the slightly misleading variable names.
        """
        if self.user_options.curvature_type == "curl(b/B) with x-y derivatives":
            if not self.user_options.orthogonal:
                raise ValueError(
                    'curvature_type = "curl(b/B) with x-y derivatives" does not '
                    "support non-orthogonal grids."
                )

            # calculate curl on x-y grid
            self.curl_bOverB_x = (
                -2.0
                * self.bpsign
                * self.Bpxy
                * self.Btxy
                * self.Rxy
                / (self.hy * self.Bxy**3)
                * self.DDY("#Bxy")
            )
            self.curl_bOverB_y = (
                -self.bpsign * self.Bpxy / self.hy * self.DDX("#Btxy*#Rxy/#Bxy**2")
            )
            self.curl_bOverB_z = (
                self.Bpxy**3 / (self.hy * self.Bxy**2) * self.DDX("#hy/#Bpxy")
                - self.Btxy * self.Rxy / self.Bxy**2 * self.DDX("#Btxy/#Rxy")
                - self.I * self.curl_bOverB_x
            )
            self.bxcvx = self.Bxy / 2.0 * self.curl_bOverB_x
            self.bxcvy = self.Bxy / 2.0 * self.curl_bOverB_y
            self.bxcvz = self.Bxy / 2.0 * self.curl_bOverB_z
        elif self.user_options.curvature_type == "curl(b/B)":
            # Calculate Curl(b/B) in R-Z, then project onto x-y-z components
            # This calculates contravariant components of a curvature vector

            equilib = self.meshParent.equilibrium
            psi = equilib.psi

            def fpol(R, Z):
                return equilib.fpol(psi(R, Z))

            def fpolprime(R, Z):
                return equilib.fpolprime(psi(R, Z))

            # BR = dpsi/dZ / R
            BR = equilib.Bp_R
            # BZ = -dpsi/dR / R
            BZ = equilib.Bp_Z
            Bzeta = equilib.Bzeta
            B2 = equilib.B2
            dBzetadR = equilib.dBzetadR
            dBzetadZ = equilib.dBzetadZ
            dBRdZ = equilib.dBRdZ
            dBZdR = equilib.dBZdR
            dB2dR = equilib.dB2dR
            dB2dZ = equilib.dB2dZ

            # In cylindrical coords
            # curl(A) = (1/R*d(AZ)/dzeta - d(Azeta)/dZ)  * Rhat
            #           + 1/R*(d(R Azeta)/dR - d(AR)/dzeta) * Zhat
            #           + (d(AR)/dZ - d(AZ)/dR) * zetahat
            # Where AR, AZ and Azeta are the components on a basis of unit vectors,
            # i.e. AR = A.Rhat; AZ = A.Zhat; Azeta = A.zetahat
            # https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates,
            #
            # curl(b/B) = curl((BR/B2), (BZ/B2), (Bzeta/B2))
            # curl(b/B)_Rhat = 1/R d(BZ/B2)/dzeta - d(Bzeta/B2)/dZ
            #                = 1/(R*B2)*d(BZ)/dzeta - BZ/(R*B4)*d(B2)/dzeta
            #                  - 1/B2*d(Bzeta)/dZ + Bzeta/B4*d(B2)/dZ
            #                = -1/B2*d(Bzeta)/dZ + Bzeta/B4*d(B2)/dZ
            # curl(b/B)_Zhat = 1/R * (d(R Bzeta/B2)/dR - d(BR/B2)/dzeta)
            #                = Bzeta/(R*B2) + 1/B2*d(Bzeta)/dR - Bzeta/B4*d(B2)/dR
            #                  - 1/(R*B2)*d(BR)/dzeta + BR/(R*B4)*d(B2)/dzeta
            #                = Bzeta/(R*B2) + 1/B2*d(Bzeta)/dR - Bzeta/B4*d(B2)/dR
            # curl(b/B)_zetahat = d(BR/B2)/dZ - d(BZ/B2)/dR
            #                   = 1/B2*d(BR)/dZ - BR/B4*d(B2)/dZ
            #                     - 1/B2*d(BZ)/dR + BZ/B4*d(B2)/dR
            # remembering d/dzeta=0 for axisymmetric equilibrium
            def curl_bOverB_Rhat(R, Z):
                return -dBzetadZ(R, Z) / B2(R, Z) + Bzeta(R, Z) / B2(R, Z) ** 2 * dB2dZ(
                    R, Z
                )

            def curl_bOverB_Zhat(R, Z):
                return (
                    Bzeta(R, Z) / (R * B2(R, Z))
                    + dBzetadR(R, Z) / B2(R, Z)
                    - Bzeta(R, Z) / B2(R, Z) ** 2 * dB2dR(R, Z)
                )

            def curl_bOverB_zetahat(R, Z):
                return (
                    dBRdZ(R, Z) / B2(R, Z)
                    - BR(R, Z) / B2(R, Z) ** 2 * dB2dZ(R, Z)
                    - dBZdR(R, Z) / B2(R, Z)
                    + BZ(R, Z) / B2(R, Z) ** 2 * dB2dR(R, Z)
                )

            # We want to output contravariant components of Curl(b/B) in the
            # locally field-aligned coordinate system.
            # The contravariant components of an arbitrary vector A are
            # A^x = A.Grad(x)
            # A^y = A.Grad(y)
            # A^z = A.Grad(z)

            # Grad in cylindrical coordinates is
            # Grad(f) = df/dR Rhat + 1/R df/dzeta zetahat + df/dZ Zhat
            # https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates,

            # x = psi - psi_min
            # dpsi/dR = -R*BZ
            # dpsi/dZ = R*BR
            # => Grad(x) = (dpsi/dR, 0, dpsi/dZ).(Rhat, zetahat, Zhat)
            # => Grad(x) = (-R BZ, 0, R BR).(Rhat, zetahat, Zhat)
            self.curl_bOverB_x = curl_bOverB_Rhat(self.Rxy, self.Zxy) * (
                -self.Rxy * BZ(self.Rxy, self.Zxy)
            ) + curl_bOverB_Zhat(self.Rxy, self.Zxy) * (
                self.Rxy * BR(self.Rxy, self.Zxy)
            )

            if self.user_options.orthogonal:
                # Grad(y) = (BR, 0, BZ)/(hy Bp)
                self.curl_bOverB_y = (
                    curl_bOverB_Rhat(self.Rxy, self.Zxy) * BR(self.Rxy, self.Zxy)
                    + curl_bOverB_Zhat(self.Rxy, self.Zxy) * BZ(self.Rxy, self.Zxy)
                ) / (self.Bpxy * self.hy)
            else:
                # Grad(y) = (d_Z, 0, -d_R)/(hy*cosBeta)
                #         = (BR*cosBeta-BZ*sinBeta, 0, BZ*cosBeta+BR*sinBeta)
                #           /(Bp*hy*cosBeta)
                #         = (BR-BZ*tanBeta, 0, BZ+BR*tanBeta)/(Bp*hy)
                self.curl_bOverB_y = (
                    curl_bOverB_Rhat(self.Rxy, self.Zxy)
                    * (BR(self.Rxy, self.Zxy) - BZ(self.Rxy, self.Zxy) * self.tanBeta)
                    + curl_bOverB_Zhat(self.Rxy, self.Zxy)
                    * (BZ(self.Rxy, self.Zxy) + BR(self.Rxy, self.Zxy) * self.tanBeta)
                ) / (self.Bpxy * self.hy)

            # Grad(z) = Grad(zeta) - Bt*hy/(Bp*R)*Grad(y) - I*Grad(x)
            # Grad(z) = (0, 1/R, 0) - Bt*hy/(Bp*R)*Grad(y) - I*Grad(x)
            self.curl_bOverB_z = (
                curl_bOverB_zetahat(self.Rxy, self.Zxy) / self.Rxy
                - self.Btxy * self.hy / (self.Bpxy * self.Rxy) * self.curl_bOverB_y
                - self.I * self.curl_bOverB_x
            )

            # bxcv is calculated this way for backward compatibility with Hypnotoad.
            # bxcv stands for 'b x kappa' where kappa is the field-line curvature, which
            # is not exactly equivalent to the result here, but this is how Hypnotoad
            # passed 'curvature' calculated as curl(b/B)
            self.bxcvx = self.Bxy / 2.0 * self.curl_bOverB_x
            self.bxcvy = self.Bxy / 2.0 * self.curl_bOverB_y
            self.bxcvz = self.Bxy / 2.0 * self.curl_bOverB_z
        elif self.user_options.curvature_type == "bxkappa":
            raise ValueError("bxkappa form of curvature not implemented yet")
            self.bxcvx = float("nan")
            self.bxcvy = float("nan")
            self.bxcvz = float("nan")
        else:
            raise ValueError(
                "Unrecognized option '"
                + str(self.user_options.curvature_type)
                + "' for curvature type"
            )

    def calcHy(self):
        # hy = 1/|Grad(theta)|
        # hy = 1/(dtheta/ds) at constant psi, phi when psi and theta are orthogonal
        # dtheta/ds \approx dtheta/sqrt((R(j+1/2)-R(j-1/2))**2 + (Z(j+1/2)-Z(j-1/2)**2)
        if not self.user_options.orthogonal:
            warnings.warn("need to check that this is correct for non-orthogonal grids")

        hy = MultiLocationArray(self.nx, self.ny)
        # contours have accurately calculated distances
        # calculate distances between j+/-0.5
        for i in range(self.nx):
            print(f"{self.name} calcHy {i} / {2 * self.nx + 1}", end="\r", flush=True)
            d = numpy.array(
                self.contours[2 * i + 1].get_distance(psi=self.equilibriumRegion.psi)
            )
            hy.centre[i, :] = d[2::2] - d[:-2:2]
            hy.ylow[i, 1:-1] = d[3:-1:2] - d[1:-3:2]
            if self.connections["lower"] is not None:
                cbelow = self.getNeighbour("lower").contours[2 * i + 1]
                dbelow = cbelow.get_distance(psi=self.equilibriumRegion.psi)
                hy.ylow[i, 0] = d[1] - d[0] + dbelow[-1] - dbelow[-2]
            else:
                # no region below, so estimate distance to point before '0' as the same
                # as from '0' to '1'
                hy.ylow[i, 0] = 2.0 * (d[1] - d[0])
            if self.connections["upper"] is not None:
                cabove = self.getNeighbour("upper").contours[2 * i + 1]
                dabove = cabove.get_distance(psi=self.equilibriumRegion.psi)
                hy.ylow[i, -1] = d[-1] - d[-2] + dabove[1] - dabove[0]
            else:
                # no region below, so estimate distance to point before '0' as the same
                # as from '0' to '1'
                hy.ylow[i, -1] = 2.0 * (d[-1] - d[-2])

        for i in range(self.nx + 1):
            print(
                f"{self.name} calcHy {i + self.nx} / {2 * self.nx + 1}",
                end="\r",
                flush=True,
            )
            d = numpy.array(
                self.contours[2 * i].get_distance(psi=self.equilibriumRegion.psi)
            )
            hy.xlow[i, :] = d[2::2] - d[:-2:2]
            hy.corners[i, 1:-1] = d[3:-1:2] - d[1:-3:2]
            if self.connections["lower"] is not None:
                cbelow = self.getNeighbour("lower").contours[2 * i]
                dbelow = cbelow.get_distance(psi=self.equilibriumRegion.psi)
                hy.corners[i, 0] = d[1] - d[0] + dbelow[-1] - dbelow[-2]
            else:
                # no region below, so estimate distance to point before '0' as the same
                # as from '0' to '1'
                hy.corners[i, 0] = 2.0 * (d[1] - d[0])
            if self.connections["upper"] is not None:
                cabove = self.getNeighbour("upper").contours[2 * i]
                dabove = cabove.get_distance(psi=self.equilibriumRegion.psi)
                hy.corners[i, -1] = d[-1] - d[-2] + dabove[1] - dabove[0]
            else:
                # no region below, so estimate distance to point before '0' as the same
                # as from '0' to '1'
                hy.corners[i, -1] = 2.0 * (d[-1] - d[-2])

        hy /= self.dy

        def negative_indices(a):
            xinds, yinds = numpy.where(a < 0.0)
            return list(zip(list(xinds), list(yinds)))

        if not numpy.all(hy.centre > 0.0):
            raise ValueError(
                f"hy.centre should always be positive. Negative values found in "
                f"region '{self.name}' at (x,y) indices {negative_indices(hy.centre)}"
            )
        if not numpy.all(hy.xlow > 0.0):
            raise ValueError(
                f"hy.xlow should always be positive. Negative values found in "
                f"region '{self.name}' at (x,y) indices {negative_indices(hy.xlow)}"
            )
        if not numpy.all(hy.ylow > 0.0):
            raise ValueError(
                f"hy.ylow should always be positive. Negative values found in "
                f"region '{self.name}' at (x,y) indices {negative_indices(hy.ylow)}"
            )
        if not numpy.all(hy.corners > 0.0):
            raise ValueError(
                f"hy.corners should always be positive. Negative values found in "
                f"region '{self.name}' at (x,y) indices {negative_indices(hy.corners)}"
            )
        if not numpy.all(hy.xlow > 0.0):
            raise ValueError("hy.xlow should always be positive")
        if not numpy.all(hy.ylow > 0.0):
            raise ValueError("hy.ylow should always be positive")
        if not numpy.all(hy.corners > 0.0):
            raise ValueError("hy.corners should always be positive")

        return hy

    def calcBeta(self):
        """
        Calculate beta (angle between e_x and Grad(x), also the angle between e_y and
        Grad(y)), used for non-orthogonal grid
        """

        self.cosBeta = MultiLocationArray(self.nx, self.ny)
        self.sinBeta = MultiLocationArray(self.nx, self.ny)
        self.tanBeta = MultiLocationArray(self.nx, self.ny)

        # for centre points
        ###################

        # vector from i-1/2 to i+1/2
        delta_x = [
            self.Rxy.xlow[1:, :] - self.Rxy.xlow[:-1, :],
            self.Zxy.xlow[1:, :] - self.Zxy.xlow[:-1, :],
        ]
        # normalise to 1
        mod_delta_x = numpy.sqrt(delta_x[0] ** 2 + delta_x[1] ** 2)
        delta_x[0] /= mod_delta_x
        delta_x[1] /= mod_delta_x

        # vector in the Grad(psi) direction
        delta_psi = [
            self.meshParent.equilibrium.f_R(self.Rxy.centre, self.Zxy.centre),
            self.meshParent.equilibrium.f_Z(self.Rxy.centre, self.Zxy.centre),
        ]
        # normalise to 1
        mod_delta_psi = numpy.sqrt(delta_psi[0] ** 2 + delta_psi[1] ** 2)
        delta_psi[0] /= mod_delta_psi
        delta_psi[1] /= mod_delta_psi

        # cosBeta = delta_x.delta_psi
        self.cosBeta.centre = delta_x[0] * delta_psi[0] + delta_x[1] * delta_psi[1]

        # Rotate delta_psi 90 degrees clockwise gives unit vector in e_y direction
        delta_y = [delta_psi[1], -delta_psi[0]]

        # sin(beta) = cos(pi/2 - beta) = e_x_hat.e_y_hat = delta_x.delta_y
        self.sinBeta.centre = delta_x[0] * delta_y[0] + delta_x[1] * delta_y[1]

        # for ylow points
        #################

        # vector from i-1/2 to i+1/2
        delta_x = [
            self.Rxy.corners[1:, :] - self.Rxy.corners[:-1, :],
            self.Zxy.corners[1:, :] - self.Zxy.corners[:-1, :],
        ]
        # normalise to 1
        mod_delta_x = numpy.sqrt(delta_x[0] ** 2 + delta_x[1] ** 2)
        delta_x[0] /= mod_delta_x
        delta_x[1] /= mod_delta_x

        # unit vector in the Grad(psi) direction
        delta_psi = [
            self.meshParent.equilibrium.f_R(self.Rxy.ylow, self.Zxy.ylow),
            self.meshParent.equilibrium.f_Z(self.Rxy.ylow, self.Zxy.ylow),
        ]
        # normalise to 1
        mod_delta_psi = numpy.sqrt(delta_psi[0] ** 2 + delta_psi[1] ** 2)
        delta_psi[0] /= mod_delta_psi
        delta_psi[1] /= mod_delta_psi

        # cosBeta = delta_x.delta_psi
        self.cosBeta.ylow = delta_x[0] * delta_psi[0] + delta_x[1] * delta_psi[1]

        # Rotate delta_psi 90 degrees clockwise gives unit vector in e_y direction
        delta_y = [delta_psi[1], -delta_psi[0]]

        # sin(beta) = cos(pi/2 - beta) = e_x.e_y = delta_x.delta_y
        self.sinBeta.ylow = delta_x[0] * delta_y[0] + delta_x[1] * delta_y[1]

        self.tanBeta = self.sinBeta / self.cosBeta

    def calcZShift(self):
        """
        Calculate zShift by integrating along FineContours.
        """
        # Step in toroidal direction for a step of length ds in the poloidal
        # plane along a FineContour is ds*Bt/Bp.
        # The change in toroidal angle is therefore (ds Bt)/(R Bp).
        # The toroidal angle along the FineContour can be calculated by
        # integrating using scipy.integrate.cumulative_trapezoid

        # Cannot just test 'connections['lower'] is not None' because periodic regions
        # always have a lower connection - requires us to give a yGroupIndex to each
        # region when creating the groups.
        if self.yGroupIndex != 0:
            return None

        region = self
        region.zShift = MultiLocationArray(region.nx, region.ny)
        region.zShift.centre[:, :] = 0.0
        region.zShift.ylow[:, :] = 0.0
        region.zShift.xlow[:, :] = 0.0
        region.zShift.corners[:, :] = 0.0
        while True:
            print("calcZShift", region.name, end="\r", flush=True)
            for i, contour in enumerate(region.contours):
                fine_contour = contour.get_fine_contour(psi=self.equilibriumRegion.psi)
                fine_distance = fine_contour.distance

                def integrand_func(R, Z):
                    Bt = (
                        self.meshParent.equilibrium.fpol(
                            self.meshParent.equilibrium.psi(R, Z)
                        )
                        / R
                    )
                    Bp = numpy.sqrt(
                        self.meshParent.equilibrium.Bp_R(R, Z) ** 2
                        + self.meshParent.equilibrium.Bp_Z(R, Z) ** 2
                    )
                    return Bt / (R * Bp)

                integrand = integrand_func(
                    fine_contour.positions[:, 0],
                    fine_contour.positions[:, 1],
                )

                zShift_fine = cumulative_trapezoid(
                    integrand, x=fine_distance, initial=0.0
                )

                # Make sure zShift_fine starts at the 'startInd' of the
                # contour/fine_contour
                zShift_fine[:] -= zShift_fine[fine_contour.startInd]

                zShift_interpolator = interp1d(
                    fine_distance, zShift_fine, kind="linear", assume_sorted=True
                )
                zShift_contour = zShift_interpolator(
                    contour.get_distance(psi=self.equilibriumRegion.psi)
                )

                if i % 2 == 0:
                    # xlow and corners
                    xind = i // 2
                    region.zShift.corners[xind, :] += zShift_contour[::2]
                    region.zShift.xlow[xind, :] += zShift_contour[1::2]
                else:
                    # centre and ylow
                    xind = i // 2
                    region.zShift.ylow[xind, :] += zShift_contour[::2]
                    region.zShift.centre[xind, :] += zShift_contour[1::2]

            next_region = region.getNeighbour("upper")
            if (next_region is None) or (next_region is self):
                # Note: If periodic, next_region is self (back to start)
                break
            else:
                # Initialise with values at the lower y-boundary of next_region
                next_region.zShift = MultiLocationArray(next_region.nx, next_region.ny)
                next_region.zShift.centre[:, :] = region.zShift.ylow[
                    :, -1, numpy.newaxis
                ]
                next_region.zShift.ylow[:, :] = region.zShift.ylow[:, -1, numpy.newaxis]
                next_region.zShift.xlow[:, :] = region.zShift.corners[
                    :, -1, numpy.newaxis
                ]
                next_region.zShift.corners[:, :] = region.zShift.corners[
                    :, -1, numpy.newaxis
                ]
                region = next_region

        # Calculate ShiftAngle for closed field line regions
        self.ShiftAngle = MultiLocationArray(self.nx, 1)
        if self.connections["lower"] is not None:
            # This is a periodic region (we already checked that the self.yGroupIndex is
            # 0).
            # 'region' is the last region in the y-group
            self.ShiftAngle.centre = (
                region.zShift.ylow[:, -1] - self.zShift.ylow[:, 0]
            ).reshape((-1, 1))
            self.ShiftAngle.xlow = (
                region.zShift.corners[:, -1] - self.zShift.corners[:, 0]
            ).reshape((-1, 1))

    def calcPoloidalDistance(self):
        """
        Calculate poloidal distance by following contours between regions.
        """
        # Cannot just test 'connections['lower'] is not None' because periodic regions
        # always have a lower connection - requires us to give a yGroupIndex to each
        # region when creating the groups.
        if self.yGroupIndex != 0:
            return None

        region = self
        region.poloidal_distance = MultiLocationArray(region.nx, region.ny)
        region.poloidal_distance.centre = 0.0
        region.poloidal_distance.ylow = 0.0
        region.poloidal_distance.xlow = 0.0
        region.poloidal_distance.corners = 0.0

        # Initialise so that distance counts from the lower wall (for SOL/PFR) or wall
        # (for core)
        for i in range(self.nx):
            c = region.contours[2 * i + 1]
            # Cell-centre points
            region.poloidal_distance.centre[i, :] -= c.get_distance(
                psi=self.meshParent.equilibrium.psi
            )[c.startInd]
            # ylow points
            region.poloidal_distance.ylow[i, :] -= c.get_distance(
                psi=self.meshParent.equilibrium.psi
            )[c.startInd]
        for i in range(self.nx + 1):
            c = region.contours[2 * i]
            # Cell-centre points
            region.poloidal_distance.xlow[i, :] -= c.get_distance(
                psi=self.meshParent.equilibrium.psi
            )[c.startInd]
            # ylow points
            region.poloidal_distance.corners[i, :] -= c.get_distance(
                psi=self.meshParent.equilibrium.psi
            )[c.startInd]

        # Get distances from contours
        while True:
            for i in range(self.nx):
                c = region.contours[2 * i + 1]
                # Cell-centre points
                region.poloidal_distance.centre[i, :] += c.get_distance(
                    psi=self.meshParent.equilibrium.psi
                )[1::2]
                # ylow points
                region.poloidal_distance.ylow[i, :] += c.get_distance(
                    psi=self.meshParent.equilibrium.psi
                )[::2]
            for i in range(self.nx + 1):
                c = region.contours[2 * i]
                # Cell-centre points
                region.poloidal_distance.xlow[i, :] += c.get_distance(
                    psi=self.meshParent.equilibrium.psi
                )[1::2]
                # ylow points
                region.poloidal_distance.corners[i, :] += c.get_distance(
                    psi=self.meshParent.equilibrium.psi
                )[::2]

            next_region = region.getNeighbour("upper")
            if (next_region is None) or (next_region is self):
                # Note: If periodic, next_region is self (back to start)
                break
            else:
                # Initialise with values at the lower y-boundary of next_region
                next_region.poloidal_distance = MultiLocationArray(
                    next_region.nx, next_region.ny
                )
                next_region.poloidal_distance.centre[
                    :, :
                ] = region.poloidal_distance.ylow[:, -1, numpy.newaxis]
                next_region.poloidal_distance.ylow[
                    :, :
                ] = region.poloidal_distance.ylow[:, -1, numpy.newaxis]
                next_region.poloidal_distance.xlow[
                    :, :
                ] = region.poloidal_distance.corners[:, -1, numpy.newaxis]
                next_region.poloidal_distance.corners[
                    :, :
                ] = region.poloidal_distance.corners[:, -1, numpy.newaxis]
                region = next_region

        # Save total poloidal distance in core
        self.total_poloidal_distance = MultiLocationArray(region.nx, 1)
        if self.connections["lower"] is not None:
            # This is a periodic region (we already checked that the self.yGroupIndex is
            # 0).
            # 'region' is the last region in the y-group
            self.total_poloidal_distance.centre[:, 0] = region.poloidal_distance.ylow[
                :, -1
            ]
            self.total_poloidal_distance.xlow[:, 0] = region.poloidal_distance.corners[
                :, -1
            ]

    def getNeighbour(self, face):
        if self.connections[face] is None:
            return None
        else:
            return self.meshParent.regions[self.connections[face]]

    def _eval_from_region(self, expr, region=None, component=None):
        # Utility routine to evaluate an expression using different MeshRegions
        # Names of fields belonging to the MeshRegion are indicated by a '#' in expr,
        # e.g.  if 'foo' and 'bar' are two member variables, we could have
        # expr='#foo + #bar'

        if region is None:
            region_string = "self"
        else:
            region_string = "self.getNeighbour('" + region + "')"

        if component is None:
            component = ""
        else:
            component = "." + component

        # replace the name of the field with an expression to get that field from the
        # MeshRegion 'region'
        expr = re.sub("#(\\w+)", region_string + ".__dict__['\\1']" + component, expr)

        return eval(expr)

    def DDX(self, expr):
        # x-derivative of a MultiLocationArray, calculated with 2nd order central
        # differences

        f = self._eval_from_region(expr)

        result = MultiLocationArray(self.nx, self.ny)

        if f.xlow is not None:
            result.centre[...] = (f.xlow[1:, :] - f.xlow[:-1, :]) / self.dx.centre
        else:
            warnings.warn(
                "No xlow field available to calculate DDX(" + expr + ").centre"
            )
        if f.corners is not None:
            result.ylow[...] = (f.corners[1:, :] - f.corners[:-1, :]) / self.dx.ylow
        else:
            warnings.warn(
                "No corners field available to calculate DDX(" + expr + ").ylow"
            )

        if f.centre is not None:
            result.xlow[1:-1, :] = (f.centre[1:, :] - f.centre[:-1, :]) / self.dx.xlow[
                1:-1, :
            ]
            if self.connections["inner"] is not None:
                f_inner = self._eval_from_region(expr, "inner", "centre[-1, :]")
                result.xlow[0, :] = (f.centre[0, :] - f_inner) / self.dx.xlow[0, :]
            else:
                result.xlow[0, :] = (f.centre[0, :] - f.xlow[0, :]) / (
                    self.dx.xlow[0, :] / 2.0
                )
            if self.connections["outer"] is not None:
                f_outer = self._eval_from_region(expr, "outer", "centre[0, :]")
                result.xlow[-1, :] = (f_outer - f.centre[-1, :]) / self.dx.xlow[-1, :]
            else:
                result.xlow[-1, :] = (f.xlow[-1, :] - f.centre[-1, :]) / (
                    self.dx.xlow[-1, :] / 2.0
                )
        else:
            warnings.warn(
                "No centre field available to calculate DDX(" + expr + ").xlow"
            )

        if f.ylow is not None:
            result.corners[1:-1, :] = (
                f.ylow[1:, :] - f.ylow[:-1, :]
            ) / self.dx.corners[1:-1, :]
            if self.connections["inner"] is not None:
                f_inner = self._eval_from_region(expr, "inner", "ylow[-1, :]")
                result.corners[0, :] = (f.ylow[0, :] - f_inner) / self.dx.corners[0, :]
            else:
                result.corners[0, :] = (f.ylow[0, :] - f.corners[0, :]) / (
                    self.dx.corners[0, :] / 2.0
                )
            if self.connections["outer"] is not None:
                f_outer = self._eval_from_region(expr, "outer", "ylow[0, :]")
                result.corners[-1, :] = (f_outer - f.ylow[-1, :]) / self.dx.corners[
                    -1, :
                ]
            else:
                result.corners[-1, :] = (f.corners[-1, :] - f.ylow[-1, :]) / (
                    self.dx.corners[-1, :] / 2.0
                )
        else:
            warnings.warn(
                "No ylow field available to calculate DDX(" + expr + ").corners"
            )

        return result

    def DDY(self, expr):
        # y-derivative of a MultiLocationArray, calculated with 2nd order central
        # differences
        f = self._eval_from_region(expr)

        result = MultiLocationArray(self.nx, self.ny)

        if f.ylow is not None:
            result.centre[...] = (f.ylow[:, 1:] - f.ylow[:, :-1]) / self.dy.centre
        else:
            warnings.warn(
                "No ylow field available to calculate DDY(" + expr + ").centre"
            )
        if f.corners is not None:
            result.xlow[...] = (f.corners[:, 1:] - f.corners[:, :-1]) / self.dy.xlow
        else:
            warnings.warn(
                "No corners field available to calculate DDY(" + expr + ").xlow"
            )

        if f.centre is not None:
            result.ylow[:, 1:-1] = (f.centre[:, 1:] - f.centre[:, :-1]) / self.dy.ylow[
                :, 1:-1
            ]
            if self.connections["lower"] is not None:
                f_lower = self._eval_from_region(expr, "lower", "centre[:, -1]")
                result.ylow[:, 0] = (f.centre[:, 0] - f_lower) / self.dy.ylow[:, 0]
            else:
                result.ylow[:, 0] = (f.centre[:, 0] - f.ylow[:, 0]) / (
                    self.dy.ylow[:, 0] / 2.0
                )
            if self.connections["upper"] is not None:
                f_upper = self._eval_from_region(expr, "upper", "centre[:, 0]")
                result.ylow[:, -1] = (f_upper - f.centre[:, -1]) / self.dy.ylow[:, -1]
            else:
                result.ylow[:, -1] = (f.ylow[:, -1] - f.centre[:, -1]) / (
                    self.dy.ylow[:, -1] / 2.0
                )
        else:
            warnings.warn(
                "No centre field available to calculate DDY(" + expr + ").ylow"
            )

        if f.xlow is not None:
            result.corners[:, 1:-1] = (
                f.xlow[:, 1:] - f.xlow[:, :-1]
            ) / self.dy.corners[:, 1:-1]
            if self.connections["lower"] is not None:
                f_lower = self._eval_from_region(expr, "lower", "xlow[:, -1]")
                result.corners[:, 0] = (f.xlow[:, 0] - f_lower) / self.dy.corners[:, 0]
            else:
                result.corners[:, 0] = (f.xlow[:, 0] - f.corners[:, 0]) / (
                    self.dy.corners[:, 0] / 2.0
                )
            if self.connections["upper"] is not None:
                f_upper = self._eval_from_region(expr, "upper", "xlow[:, 0]")
                result.corners[:, -1] = (f_upper - f.xlow[:, -1]) / self.dy.corners[
                    :, -1
                ]
            else:
                result.corners[:, -1] = (f.corners[:, -1] - f.xlow[:, -1]) / (
                    self.dy.corners[:, -1] / 2.0
                )
        else:
            warnings.warn(
                "No xlow field available to calculate DDY(" + expr + ").corners"
            )

        return result

    def smoothnl_inner1(self, varname):
        f = getattr(self, varname)
        if self.connections["inner"] is not None:
            f_inner = getattr(self.getNeighbour("inner"), varname)
        else:
            f_inner = None
        if self.connections["outer"] is not None:
            f_outer = getattr(self.getNeighbour("outer"), varname)
        else:
            f_outer = None
        if self.connections["lower"] is not None:
            f_lower = getattr(self.getNeighbour("lower"), varname)
        else:
            f_lower = None
        if self.connections["upper"] is not None:
            f_upper = getattr(self.getNeighbour("upper"), varname)
        else:
            f_upper = None

        mxn = MultiLocationArray(self.nx, self.ny)
        myn = MultiLocationArray(self.nx, self.ny)

        if f.centre is not None:
            dxm = numpy.zeros(f.centre.shape)
            dxp = numpy.zeros(f.centre.shape)
            dym = numpy.zeros(f.centre.shape)
            dyp = numpy.zeros(f.centre.shape)

            dxm[1:-1, :] = f.centre[1:-1, :] - f.centre[:-2, :]
            dxp[1:-1, :] = f.centre[2:, :] - f.centre[1:-1, :]

            if f_inner is not None:
                dxm[0, :] = f.centre[0, :] - f_inner.centre[-1, :]
            if f_outer is not None:
                dxm[-1, :] = f_outer.centre[0, :] - f.centre[-1, :]

            dym[:, 1:-1] = f.centre[:, 1:-1] - f.centre[:, :-2]
            dyp[:, 1:-1] = f.centre[:, 2:] - f.centre[:, 1:-1]

            if f_lower is not None:
                dym[:, 0] = f.centre[:, 0] - f_lower.centre[:, -1]
            if f_upper is not None:
                dym[:, -1] = f_upper.centre[:, 0] - f.centre[:, -1]

            mxn.centre = 0.5 * (abs(dxm) + abs(dxp))
            myn.centre = 0.5 * (abs(dym) + abs(dyp))

        # Note indexing of staggered fields from neighbouring regions - staggered points
        # overlap at the boundaries.
        # Don't set outer or upper boundary points so we can sum the results without
        # duplicates

        if f.xlow is not None:
            dxm = numpy.zeros(f.xlow.shape)
            dxp = numpy.zeros(f.xlow.shape)
            dym = numpy.zeros(f.xlow.shape)
            dyp = numpy.zeros(f.xlow.shape)

            dxm[1:-2, :] = f.xlow[1:-2, :] - f.xlow[:-3, :]
            dxp[1:-2, :] = f.xlow[2:-1, :] - f.xlow[1:-2, :]

            if f_inner is not None:
                dxm[0, :] = f.xlow[0, :] - f_inner.xlow[-2, :]
            if f_outer is not None:
                dxm[-2, :] = f_outer.xlow[0, :] - f.xlow[-2, :]

            dym[:, 1:-1] = f.xlow[:, 1:-1] - f.xlow[:, :-2]
            dyp[:, 1:-1] = f.xlow[:, 2:] - f.xlow[:, 1:-1]

            if f_lower is not None:
                dym[:, 0] = f.xlow[:, 0] - f_lower.xlow[:, -1]
            if f_upper is not None:
                dym[:, -1] = f_upper.xlow[:, 0] - f.xlow[:, -1]

            mxn.xlow = 0.5 * (abs(dxm) + abs(dxp))
            myn.xlow = 0.5 * (abs(dym) + abs(dyp))

        if f.ylow is not None:
            dxm = numpy.zeros(f.ylow.shape)
            dxp = numpy.zeros(f.ylow.shape)
            dym = numpy.zeros(f.ylow.shape)
            dyp = numpy.zeros(f.ylow.shape)

            dxm[1:-1, :] = f.ylow[1:-1, :] - f.ylow[:-2, :]
            dxp[1:-1, :] = f.ylow[2:, :] - f.ylow[1:-1, :]

            if f_inner is not None:
                dxm[0, :] = f.ylow[0, :] - f_inner.ylow[-1, :]
            if f_outer is not None:
                dxm[-1, :] = f_outer.ylow[0, :] - f.ylow[-1, :]

            dym[:, 1:-2] = f.ylow[:, 1:-2] - f.ylow[:, :-3]
            dyp[:, 1:-2] = f.ylow[:, 2:-1] - f.ylow[:, 1:-2]

            if f_lower is not None:
                dym[:, 0] = f.ylow[:, 0] - f_lower.ylow[:, -2]
            if f_upper is not None:
                dym[:, -2] = f_upper.ylow[:, 0] - f.ylow[:, -2]

            mxn.ylow = 0.5 * (abs(dxm) + abs(dxp))
            myn.ylow = 0.5 * (abs(dym) + abs(dyp))

        if f.corners is not None:
            dxm = numpy.zeros(f.corners.shape)
            dxp = numpy.zeros(f.corners.shape)
            dym = numpy.zeros(f.corners.shape)
            dyp = numpy.zeros(f.corners.shape)

            dxm[1:-2, :-1] = f.corners[1:-2, :-1] - f.corners[:-3, :-1]
            dxp[1:-2, :-1] = f.corners[2:-1, :-1] - f.corners[1:-2, :-1]

            if f_inner is not None:
                dxm[0, :-1] = f.corners[0, :-1] - f_inner.corners[-2, :-1]
            if f_outer is not None:
                dxm[-2, :-1] = f_outer.corners[0, :-1] - f.corners[-2, :-1]

            dym[:-1, 1:-2] = f.corners[:-1, 1:-2] - f.corners[:-1, :-3]
            dyp[:-1, 1:-2] = f.corners[:-1, 2:-1] - f.corners[:-1, 1:-2]

            if f_lower is not None:
                dym[:-1, 0] = f.corners[:-1, 0] - f_lower.corners[:-1, -2]
            if f_upper is not None:
                dym[:-1, -2] = f_upper.corners[:-1, 0] - f.corners[:-1, -2]

            mxn.corners = 0.5 * (abs(dxm) + abs(dxp))
            myn.corners = 0.5 * (abs(dym) + abs(dyp))

        return mxn, myn

    def smoothnl_inner2(self, varname, markx, marky):
        tmp = getattr(self, varname).copy()
        if self.connections["inner"] is not None:
            tmp_inner = getattr(self.getNeighbour("inner"), varname).copy()
        else:
            tmp_inner = None
        if self.connections["outer"] is not None:
            tmp_outer = getattr(self.getNeighbour("outer"), varname).copy()
        else:
            tmp_outer = None
        if self.connections["lower"] is not None:
            tmp_lower = getattr(self.getNeighbour("lower"), varname).copy()
        else:
            tmp_lower = None
        if self.connections["upper"] is not None:
            tmp_upper = getattr(self.getNeighbour("upper"), varname).copy()
        else:
            tmp_upper = None

        # Smooth the smoothing mask
        def smooth_mask(mark):
            result = MultiLocationArray(self.nx, self.ny)
            if mark.centre is not None:
                result.centre = 0.1 * (
                    mark.centre[1:-1, 1:-1]
                    + mark.centre[0:-2, 1:-1]
                    + mark.centre[2:, 1:-1]
                    + mark.centre[1:-1, :-2]
                    + mark.centre[1:-1, 2:]
                )
            if mark.xlow is not None:
                result.xlow = 0.1 * (
                    mark.xlow[1:-1, 1:-1]
                    + mark.xlow[0:-2, 1:-1]
                    + mark.xlow[2:, 1:-1]
                    + mark.xlow[1:-1, :-2]
                    + mark.xlow[1:-1, 2:]
                )
            if mark.ylow is not None:
                result.ylow = 0.1 * (
                    mark.ylow[1:-1, 1:-1]
                    + mark.ylow[0:-2, 1:-1]
                    + mark.ylow[2:, 1:-1]
                    + mark.ylow[1:-1, :-2]
                    + mark.ylow[1:-1, 2:]
                )
            if mark.corners is not None:
                result.corners = 0.1 * (
                    mark.corners[1:-1, 1:-1]
                    + mark.corners[0:-2, 1:-1]
                    + mark.corners[2:, 1:-1]
                    + mark.corners[1:-1, :-2]
                    + mark.corners[1:-1, 2:]
                )

            return result

        mx = smooth_mask(markx)
        my = smooth_mask(marky)

        if tmp.centre is not None:
            tmp.centre[1:-1, 1:-1] = (
                (1.0 - mx.centre[1:-1, 1:-1] - my.centre[1:-1, 1:-1])
                * tmp.centre[1:-1, 1:-1]
                + mx.centre[1:-1, 1:-1]
                * 0.5
                * (tmp.centre[:-2, 1:-1] + tmp.centre[2:, 1:-1])
                + my.centre[1:-1, 1:-1]
                * 0.5
                * (tmp.centre[1:-1, :-2] + tmp.centre[1:-1, 2:])
            )
            if tmp_inner is not None:
                tmp.centre[0, 1:-1] = (
                    (1.0 - mx.centre[0, 1:-1] - my.centre[0, 1:-1])
                    * tmp.centre[0, 1:-1]
                    + mx.centre[0, 1:-1]
                    * 0.5
                    * (tmp_inner.centre[-1, 1:-1] + tmp.centre[1, 1:-1])
                    + my.centre[0, 1:-1]
                    * 0.5
                    * (tmp.centre[0, :-2] + tmp.centre[0, 2:])
                )
            else:
                tmp.centre[1, 1:-1] = tmp.centre[1, 1:-1]
            if tmp_outer is not None:
                tmp.centre[-1, 1:-1] = (
                    (1.0 - mx.centre[-1, 1:-1] - my.centre[-1, 1:-1])
                    * tmp.centre[-1, 1:-1]
                    + mx.centre[-1, 1:-1]
                    * 0.5
                    * (tmp.centre[-2, 1:-1] + tmp_outer.centre[0, 1:-1])
                    + my.centre[-1, 1:-1]
                    * 0.5
                    * (tmp.centre[-1, :-2] + tmp.centre[-1, 2:])
                )
            else:
                tmp.centre[-1, 1:-1] = tmp.centre[-2, 1:-1]
            if tmp_lower is not None:
                tmp.centre[1:-1, 0] = (
                    (1.0 - mx.centre[1:-1, 0] - my.centre[1:-1, 0])
                    * tmp.centre[1:-1, 0]
                    + mx.centre[1:-1, 0]
                    * 0.5
                    * (tmp.centre[:-2, 0] + tmp.centre[2:, 0])
                    + my.centre[1:-1, 0]
                    * 0.5
                    * (tmp_lower.centre[1:-1, -1] + tmp.centre[1:-1, 1])
                )
            if tmp_upper is not None:
                tmp.centre[1:-1, -1] = (
                    (1.0 - mx.centre[1:-1, -1] - my.centre[1:-1, -1])
                    * tmp.centre[1:-1, -1]
                    + mx.centre[1:-1, -1]
                    * 0.5
                    * (tmp.centre[:-2, -1] + tmp.centre[2:, -1])
                    + my.centre[1:-1, -1]
                    * 0.5
                    * (tmp.centre[1:-1, -2] + tmp_upper.centre[1:-1, 0])
                )
        # Note indexing of staggered fields from neighbouring regions - staggered points
        # overlap at the boundaries.

        if tmp.xlow is not None:
            tmp.xlow[1:-1, 1:-1] = (
                (1.0 - mx.xlow[1:-1, 1:-1] - my.xlow[1:-1, 1:-1]) * tmp.xlow[1:-1, 1:-1]
                + mx.xlow[1:-1, 1:-1] * 0.5 * (tmp.xlow[:-2, 1:-1] + tmp.xlow[2:, 1:-1])
                + my.xlow[1:-1, 1:-1] * 0.5 * (tmp.xlow[1:-1, :-2] + tmp.xlow[1:-1, 2:])
            )
            if tmp_inner is not None:
                tmp.xlow[0, 1:-1] = (
                    (1.0 - mx.xlow[0, 1:-1] - my.xlow[0, 1:-1]) * tmp.xlow[0, 1:-1]
                    + mx.xlow[0, 1:-1]
                    * 0.5
                    * (tmp_inner.xlow[-2, 1:-1] + tmp.xlow[1, 1:-1])
                    + my.xlow[0, 1:-1] * 0.5 * (tmp.xlow[0, :-2] + tmp.xlow[0, 2:])
                )
            else:
                tmp.xlow[0, 1:-1] = tmp.xlow[1, 1:-1]
            if tmp_outer is not None:
                tmp.xlow[-1, 1:-1] = (
                    (1.0 - mx.xlow[-1, 1:-1] - my.xlow[-1, 1:-1]) * tmp.xlow[-1, 1:-1]
                    + mx.xlow[-1, 1:-1]
                    * 0.5
                    * (tmp.xlow[-2, 1:-1] + tmp_outer.xlow[1, 1:-1])
                    + my.xlow[-1, 1:-1] * 0.5 * (tmp.xlow[-1, :-2] + tmp.xlow[-1, 2:])
                )
            else:
                tmp.xlow[-1, 1:-1] = tmp.xlow[-2, 1:-1]
            if tmp_lower is not None:
                tmp.xlow[1:-1, 0] = (
                    (1.0 - mx.xlow[1:-1, 0] - my.xlow[1:-1, 0]) * tmp.xlow[1:-1, 0]
                    + mx.xlow[1:-1, 0] * 0.5 * (tmp.xlow[:-2, 0] + tmp.xlow[2:, 0])
                    + my.xlow[1:-1, 0]
                    * 0.5
                    * (tmp_lower.xlow[1:-1, -1] + tmp.xlow[1:-1, 1])
                )
            if tmp_upper is not None:
                tmp.xlow[1:-1, -1] = (
                    (1.0 - mx.xlow[1:-1, -1] - my.xlow[1:-1, -1]) * tmp.xlow[1:-1, -1]
                    + mx.xlow[1:-1, -1] * 0.5 * (tmp.xlow[:-2, -1] + tmp.xlow[2:, -1])
                    + my.xlow[1:-1, -1]
                    * 0.5
                    * (tmp.xlow[1:-1, -2] + tmp_upper.xlow[1:-1, 0])
                )
        if tmp.ylow is not None:
            tmp.ylow[1:-1, 1:-1] = (
                (1.0 - mx.ylow[1:-1, 1:-1] - my.ylow[1:-1, 1:-1]) * tmp.ylow[1:-1, 1:-1]
                + mx.ylow[1:-1, 1:-1] * 0.5 * (tmp.ylow[:-2, 1:-1] + tmp.ylow[2:, 1:-1])
                + my.ylow[1:-1, 1:-1] * 0.5 * (tmp.ylow[1:-1, :-2] + tmp.ylow[1:-1, 2:])
            )
            if tmp_inner is not None:
                tmp.ylow[0, 1:-1] = (
                    (1.0 - mx.ylow[0, 1:-1] - my.ylow[0, 1:-1]) * tmp.ylow[0, 1:-1]
                    + mx.ylow[0, 1:-1]
                    * 0.5
                    * (tmp_inner.ylow[-1, 1:-1] + tmp.ylow[1, 1:-1])
                    + my.ylow[0, 1:-1] * 0.5 * (tmp.ylow[0, :-2] + tmp.ylow[0, 2:])
                )
            else:
                tmp.ylow[0, 1:-1] = tmp.ylow[1, 1:-1]
            if tmp_outer is not None:
                tmp.ylow[-1, 1:-1] = (
                    (1.0 - mx.ylow[-1, 1:-1] - my.ylow[-1, 1:-1]) * tmp.ylow[-1, 1:-1]
                    + mx.ylow[-1, 1:-1]
                    * 0.5
                    * (tmp.ylow[-2, 1:-1] + tmp_outer.ylow[0, 1:-1])
                    + my.ylow[-1, 1:-1] * 0.5 * (tmp.ylow[-1, :-2] + tmp.ylow[-1, 2:])
                )
            else:
                tmp.ylow[-1, 1:-1] = tmp.ylow[-2, 1:-1]
            if tmp_lower is not None:
                tmp.ylow[1:-1, 0] = (
                    (1.0 - mx.ylow[1:-1, 0] - my.ylow[1:-1, 0]) * tmp.ylow[1:-1, 0]
                    + mx.ylow[1:-1, 0] * 0.5 * (tmp.ylow[:-2, 0] + tmp.ylow[2:, 0])
                    + my.ylow[1:-1, 0]
                    * 0.5
                    * (tmp_lower.ylow[1:-1, -2] + tmp.ylow[1:-1, 1])
                )
            if tmp_upper is not None:
                tmp.ylow[1:-1, -1] = (
                    (1.0 - mx.ylow[1:-1, -1] - my.ylow[1:-1, -1]) * tmp.ylow[1:-1, -1]
                    + mx.ylow[1:-1, -1] * 0.5 * (tmp.ylow[:-2, -1] + tmp.ylow[2:, -1])
                    + my.ylow[1:-1, -1]
                    * 0.5
                    * (tmp.ylow[1:-1, -2] + tmp_upper.ylow[1:-1, 1])
                )
        if tmp.corners is not None:
            tmp.corners[1:-1, 1:-1] = (
                (1.0 - mx.corners[1:-1, 1:-1] - my.corners[1:-1, 1:-1])
                * tmp.corners[1:-1, 1:-1]
                + mx.corners[1:-1, 1:-1]
                * 0.5
                * (tmp.corners[:-2, 1:-1] + tmp.corners[2:, 1:-1])
                + my.corners[1:-1, 1:-1]
                * 0.5
                * (tmp.corners[1:-1, :-2] + tmp.corners[1:-1, 2:])
            )
            if tmp_inner is not None:
                tmp.corners[0, 1:-1] = (
                    (1.0 - mx.corners[0, 1:-1] - my.corners[0, 1:-1])
                    * tmp.corners[0, 1:-1]
                    + mx.corners[0, 1:-1]
                    * 0.5
                    * (tmp_inner.corners[-2, 1:-1] + tmp.corners[1, 1:-1])
                    + my.corners[0, 1:-1]
                    * 0.5
                    * (tmp.corners[0, :-2] + tmp.corners[0, 2:])
                )
            else:
                tmp.corners[0, 1:-1] = tmp.corners[1, 1:-1]
            if tmp_outer is not None:
                tmp.corners[-1, 1:-1] = (
                    (1.0 - mx.corners[-1, 1:-1] - my.corners[-1, 1:-1])
                    * tmp.corners[-1, 1:-1]
                    + mx.corners[-1, 1:-1]
                    * 0.5
                    * (tmp.corners[-2, 1:-1] + tmp_outer.corners[1, 1:-1])
                    + my.corners[-1, 1:-1]
                    * 0.5
                    * (tmp.corners[-1, :-2] + tmp.corners[-1, 2:])
                )
            else:
                tmp.corners[-1, 1:-1] = tmp.corners[-2, 1:-1]
            if tmp_lower is not None:
                tmp.corners[1:-1, 0] = (
                    (1.0 - mx.corners[1:-1, 0] - my.corners[1:-1, 0])
                    * tmp.corners[1:-1, 0]
                    + mx.corners[1:-1, 0]
                    * 0.5
                    * (tmp.corners[:-2, 0] + tmp.corners[2:, 0])
                    + my.corners[1:-1, 0]
                    * 0.5
                    * (tmp_lower.corners[1:-1, -2] + tmp.corners[1:-1, 1])
                )
            if tmp_upper is not None:
                tmp.corners[1:-1, -1] = (
                    (1.0 - mx.corners[1:-1, -1] - my.corners[1:-1, -1])
                    * tmp.corners[1:-1, -1]
                    + mx.corners[1:-1, -1]
                    * 0.5
                    * (tmp.corners[:-2, -1] + tmp.corners[2:, -1])
                    + my.corners[1:-1, -1]
                    * 0.5
                    * (tmp.corners[1:-1, -2] + tmp_upper.corners[1:-1, 1])
                )

        diff = abs(tmp - getattr(self, varname))
        changes = []
        if diff.centre is not None:
            changes.append(numpy.max(diff.centre))
        if diff.xlow is not None:
            changes.append(numpy.max(diff.xlow))
        if diff.ylow is not None:
            changes.append(numpy.max(diff.ylow))
        if diff.corners is not None:
            changes.append(numpy.max(diff.corners))
        change = max(changes)

        return tmp, change
        setattr(self, varname, tmp)


def _calc_contour_distance(i, c, *, psi, **kwargs):
    print(
        f"Calculating contour distances: {i + 1}",
        end="\r",
        flush=True,
    )
    c.get_distance(psi=psi)
    return c


def _find_intersection(
    i_contour,
    contour,
    *,
    psi,
    equilibrium,
    lower_wall,
    upper_wall,
    max_extend,
    atol,
    refine_width,
    **kwargs,
):
    print(
        f"finding wall intersections: {i_contour + 1}",
        end="\r",
        flush=True,
    )

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
            starti = len(contour) // 2
        else:
            starti = len(contour) - 1

        # find whether one of the segments of the contour already intersects the
        # wall
        for i in range(starti, 0, -1):
            coarse_lower_intersect = equilibrium.wallIntersection(
                contour[i], contour[i - 1]
            )
            if coarse_lower_intersect is not None:
                lower_intersect_index = i - 1
                break

        count = 0
        distance = contour.get_distance(psi=psi)
        ds_extend = distance[1] - distance[0]
        while coarse_lower_intersect is None:
            # contour has not yet intersected with wall, so make it longer and
            # try again
            contour.temporaryExtend(psi=psi, extend_lower=1, ds_lower=ds_extend)
            coarse_lower_intersect = equilibrium.wallIntersection(
                contour[1], contour[0]
            )
            count += 1
            if count >= max_extend:
                print("Contour intersection failed")

                import matplotlib.pyplot as plt

                contour.plot(color="b", psi=equilibrium.psi)

                plt.plot(
                    [contour[0].R, contour[1].R],
                    [contour[0].Z, contour[1].Z],
                    color="r",
                    linewidth=3,
                )

                equilibrium.plotWall()

                plt.show()
                raise RuntimeError("extended contour too far without finding wall")
        # refine lower_intersect by finding the intersection from FineContour
        # points
        #
        # first find nearest FineContour points
        fine_contour = contour.get_fine_contour(psi=psi)
        d = fine_contour.getDistance(coarse_lower_intersect)
        i_fine = numpy.searchsorted(fine_contour.distance, d)
        # Intersection should be between i_fine-1 and i_fine, but check
        # intervals on either side if necessary
        lower_intersect = equilibrium.wallIntersection(
            Point2D(*fine_contour.positions[i_fine - 1]),
            Point2D(*fine_contour.positions[i_fine]),
        )
        if lower_intersect is None:
            lower_intersect = equilibrium.wallIntersection(
                Point2D(*fine_contour.positions[i_fine - 2]),
                Point2D(*fine_contour.positions[i_fine - 1]),
            )
            if lower_intersect is not None:
                i_fine = i_fine - 1
        if lower_intersect is None:
            lower_intersect = equilibrium.wallIntersection(
                Point2D(*fine_contour.positions[i_fine]),
                Point2D(*fine_contour.positions[i_fine + 1]),
            )
            if lower_intersect is not None:
                i_fine = i_fine + 1
        if lower_intersect is None:
            raise ValueError(
                "Did not find lower_intersect from fine_contour, even though "
                "intersection was found on contour."
            )
        # Further refine, to ensure wall point is at correct psi
        lower_intersect = contour.refinePoint(
            lower_intersect,
            Point2D(
                *(fine_contour.positions[i_fine] - fine_contour.positions[i_fine - 1])
            ),
            psi=psi,
        )

    if upper_wall:
        if lower_wall:
            starti = len(contour // 2)
        else:
            starti = 0

        # find whether one of the segments of the contour already intersects the
        # wall
        for i in range(starti, len(contour) - 1):
            coarse_upper_intersect = equilibrium.wallIntersection(
                contour[i], contour[i + 1]
            )
            if coarse_upper_intersect is not None:
                upper_intersect_index = i
                break

        count = 0
        distance = contour.get_distance(psi=psi)
        ds_extend = distance[-1] - distance[-2]
        while coarse_upper_intersect is None:
            # contour has not yet intersected with wall, so make it longer and
            # try again
            contour.temporaryExtend(psi=psi, extend_upper=1, ds_upper=ds_extend)
            coarse_upper_intersect = equilibrium.wallIntersection(
                contour[-2], contour[-1]
            )
            count += 1
            if count >= max_extend:
                print("Contour intersection failed")

                import matplotlib.pyplot as plt

                plt.plot([p.R for p in contour], [p.Z for p in contour], color="b")

                plt.plot(
                    [contour[-2].R, contour[-1].R],
                    [contour[-2].Z, contour[-1].Z],
                    color="r",
                    linewidth=3,
                )

                equilibrium.plotWall(color="k")

                plt.show()
                raise RuntimeError("extended contour too far without finding wall")
        # refine upper_intersect by finding the intersection from FineContour
        # points
        #
        # first find nearest FineContour points
        fine_contour = contour.get_fine_contour(psi=psi)
        d = fine_contour.getDistance(coarse_upper_intersect)
        i_fine = numpy.searchsorted(fine_contour.distance, d)
        # Intersection should be between i_fine-1 and i_fine, but check
        # intervals on either side if necessary
        upper_intersect = equilibrium.wallIntersection(
            Point2D(*fine_contour.positions[i_fine - 1]),
            Point2D(*fine_contour.positions[i_fine]),
        )
        if upper_intersect is None:
            upper_intersect = equilibrium.wallIntersection(
                Point2D(*fine_contour.positions[i_fine - 2]),
                Point2D(*fine_contour.positions[i_fine - 1]),
            )
            if upper_intersect is not None:
                i_fine = i_fine - 1
        if upper_intersect is None:
            upper_intersect = equilibrium.wallIntersection(
                Point2D(*fine_contour.positions[i_fine]),
                Point2D(*fine_contour.positions[i_fine + 1]),
            )
            if upper_intersect is not None:
                i_fine = i_fine + 1
        if upper_intersect is None:
            raise ValueError(
                "Did not find upper_intersect from fine_contour, even though "
                "intersection was found on contour."
            )
        # Further refine, to ensure wall point is at correct psi
        upper_intersect = contour.refinePoint(
            upper_intersect,
            Point2D(
                *(fine_contour.positions[i_fine] - fine_contour.positions[i_fine - 1])
            ),
            psi=psi,
        )

    # Create FineContour for result, so that this is done in parallel
    contour.get_fine_contour(psi=psi)

    return (
        contour,
        lower_intersect_index,
        lower_intersect,
        upper_intersect_index,
        upper_intersect,
    )


def _refine_extend(i, contour, *, psi, **kwargs):
    print(
        "refine and extend",
        i,
        end="\r",
        flush=True,
    )
    contour.refine(psi=psi)
    contour.checkFineContourExtend(psi=psi)
    return contour


def regrid_contours(
    i_contour,
    contour,
    sfunc,
    *,
    psi,
    ny_noguards,
    extend_lower,
    extend_upper,
    refine_width,
    **kwargs,
):
    print(
        f"distributing points on contour: {i_contour + 1}",
        end="\r",
        flush=True,
    )

    contour.regrid(
        2 * ny_noguards + 1,
        psi=psi,
        sfunc=sfunc,
        width=refine_width,
        extend_lower=extend_lower,
        extend_upper=extend_upper,
    )

    return contour


class Mesh:
    """
    Collection of MeshRegion objects representing the entire grid.
    """

    user_options_factory = OptionsFactory(
        # Include settings for member Equilibrium object
        Equilibrium.user_options_factory,
        # Include settings for member MeshRegion objects
        MeshRegion.user_options_factory,
        number_of_processors=WithMeta(
            1,
            doc=(
                "Number of processors to use for parallelisable loops. Note overhead "
                "of the parallelisation is fairly high (probably due to cost of "
                "pickling psi interpolation function, etc.), so increasing the number "
                "of processors above 1 is only useful for relatively large grids, or "
                "large values of finecontour_Nfine. This option may even slow things "
                "down on desktop/laptop computers which can boost clock-speed for "
                "single-core operation."
            ),
            value_type=int,
            check_all=is_positive,
        ),
    )

    def __init__(self, equilibrium, settings):
        """
        Parameters
        ----------
        equilibrium : Equilibrium
            Used to generate the grid
        settings : dict
            Non-default values to use to generate the grid. Must be consistent with the
            ones that were used to create the equilibrium
        """
        self.equilibrium = equilibrium

        self.user_options = self.user_options_factory.create(settings)
        # Check settings didn't change since equilibrium was created
        for key in self.equilibrium.user_options:
            if (key in self.user_options) and (
                self.equilibrium.user_options[key] != self.user_options[key]
            ):
                raise ValueError(
                    f"Setting {key} has been changed since equilibrium was created."
                    f"Re-create or re-load the equilibrium with the current settings."
                )

        # Print the table of options
        print(self.user_options.as_table(), flush=True)

        versions = get_versions()
        self.version = versions["version"]
        self.git_hash = versions["full-revisionid"]
        self.git_diff = None

        if versions["dirty"]:
            # There are changes from the last commit, get git diff

            from pathlib import Path
            from hypnotoad.__init__ import __file__ as hypnotoad_init_file

            hypnotoad_path = Path(hypnotoad_init_file).parent

            retval, self.git_diff = shell_safe(
                "cd " + str(hypnotoad_path) + "&& git diff", pipe=True
            )
            self.git_diff = self.git_diff.strip()

        # Generate MeshRegion object for each section of the mesh
        self.regions = {}

        # Make consecutive numbering scheme for regions
        regionlist = []
        self.region_lookup = {}
        for reg_name, eq_reg in equilibrium.regions.items():
            for i in range(eq_reg.nSegments):
                region_number = len(regionlist)
                regionlist.append((reg_name, i))
                self.region_lookup[(reg_name, i)] = region_number

        # Get connections between regions
        self.connections = {}
        for region_id, (eq_reg, i) in enumerate(regionlist):
            self.connections[region_id] = {}
            region = equilibrium.regions[eq_reg]
            c = region.connections[i]
            for key, val in c.items():
                if val is not None:
                    self.connections[region_id][key] = self.region_lookup[val]
                else:
                    self.connections[region_id][key] = None

        parallel_map = ParallelMap(
            self.user_options.number_of_processors,
            equilibrium=self.equilibrium,
        )

        self.makeRegions(parallel_map)

    def makeRegions(self, parallel_map):
        for eq_region in self.equilibrium.regions.values():
            for i in range(eq_region.nSegments):
                region_id = self.region_lookup[(eq_region.name, i)]
                eq_region_with_boundaries = eq_region.getRegridded(
                    radialIndex=i,
                    psi=self.equilibrium.psi,
                    width=self.user_options.refine_width,
                )
                self.regions[region_id] = MeshRegion(
                    self,
                    region_id,
                    eq_region_with_boundaries,
                    self.connections[region_id],
                    i,
                    self.user_options,
                    parallel_map,
                )

        # create groups that connect in x
        self.x_groups = []
        region_set = set(self.regions.values())
        while region_set:
            for region in region_set:
                if region.connections["inner"] is None:
                    break
            group = []
            while True:
                group.append(region)
                region_set.remove(region)
                region = region.getNeighbour("outer")
                if region is None or group.count(region) > 0:
                    # reached boundary or have all regions in a periodic group
                    break
            self.x_groups.append(group)

        # create groups that connect in y
        self.y_groups = []
        region_list = list(self.regions.values())

        # Once 'region_list' is empty, all regions have been added to some group
        while region_list:
            for i, first_region in enumerate(region_list):
                if first_region.connections["lower"] is None:
                    # Found a region with a lower boundary - start stepping through
                    # y-connections from here
                    break
                # note, if no region with connections['lower']=None is found, then some
                # arbitrary region will be 'first_region' after this loop. This is OK,
                # as this region must be part of a periodic group, which we will handle.

            # Find all the regions connected in the y-direction to 'first_region' and
            # add them to 'group'. Remove them from 'region_list' since each region can
            # only be in one group.
            group = []
            next_region = first_region
            while True:
                if next_region.yGroupIndex is not None:
                    raise ValueError(
                        "region should not have been added to any yGroup before"
                    )
                next_region.yGroupIndex = len(group)
                group.append(next_region)
                region_list.pop(i)

                next_region = next_region.getNeighbour("upper")
                if next_region is None or group.count(next_region) > 0:
                    # reached boundary or have all regions in a periodic group
                    break
                # index of 'next_region' in 'region_list', so we can remove
                # 'next_region' after adding to 'group' in the next step of the loop
                i = region_list.index(next_region)

            self.y_groups.append(group)

    def redistributePoints(self, nonorthogonal_settings):
        warnings.warn(
            "It is not recommended to use Mesh.redistributePoints() for 'production' "
            "output. Suggest saving the final settings to a .yaml file and creating the "
            "'production' grid non-interactively to ensure reproducibility."
        )

        self.equilibrium.resetNonorthogonalOptions(nonorthogonal_settings)

        if self.user_options.orthogonal:
            raise ValueError(
                "redistributePoints would do nothing for an orthogonal grid."
            )
        for region in self.regions.values():
            print("redistributing", region.name, flush=True)
            region.distributePointsNonorthogonal(nonorthogonal_settings)

    def calculateRZ(self):
        """
        Create arrays with R and Z values of all points in the grid
        """
        print("Get RZ values", flush=True)
        for region in self.regions.values():
            region.fillRZ()
        for region in self.regions.values():
            region.getRZBoundary()

    def geometry(self):
        """
        Calculate geometrical quantities for BOUT++
        """
        for region in self.regions.values():
            if not hasattr(region, "Rxy") or not hasattr(region, "Zxy"):
                # R and Z arrays need calculating
                self.calculateRZ()
                break
        print("Calculate geometry", flush=True)
        for region in self.regions.values():
            print("Distances", region.name, flush=True)
            region.calcDistances()
        for region in self.regions.values():
            print("1", region.name, flush=True)
            region.geometry1()
        for region in self.regions.values():
            print("2", region.name, end="\r", flush=True)
            region.geometry2()
        print("Calculate zShift", flush=True)
        for region in self.regions.values():
            print(region.name, end="\r", flush=True)
            region.calcZShift()
        print("Calculate Metric", flush=True)
        for region in self.regions.values():
            print(region.name, end="\r", flush=True)
            region.calcMetric()

        if self.user_options.curvature_smoothing == "smoothnl":
            # Nonlinear smoothing. Tries to smooth only regions with large changes in
            # gradient.
            # Smooth {bxcvx,bxcvy,bxcvz} and {curl_bOverB_x,curl_bOverB_y,curl_bOverB_z}
            # separately (not consistently with each other).
            if not self.user_options.shiftedmetric:
                # If shiftedmetric==False, would need to follow IDL hypnotoad and:
                #  - calculate bz = bxcvz + I*bxcvx
                #  - smooth bxcvx, bxcvy, and bz
                #  - set bxcvz = bz - I * bxcvx
                # and similarly for curl_bOverB_z
                raise ValueError(
                    "shiftedmetric==False not handled in "
                    "curvature_smoothing=='smoothnl'. Non-zero I requires bxcvx and "
                    "bxcvz to be smoothed consistently"
                )
            self.smoothnl("bxcvx")
            self.smoothnl("bxcvy")
            self.smoothnl("bxcvz")
            self.smoothnl("curl_bOverB_x")
            self.smoothnl("curl_bOverB_y")
            self.smoothnl("curl_bOverB_z")

    def smoothnl(self, varname):
        """
        Smoothing algorithm copied from IDL hypnotoad
        https://github.com/boutproject/BOUT-dev/blob/v4.3.2/tools/tokamak_grids/gridgen/smooth_nl.pro  # noqa: E501
        """
        npoints_centre = sum(region.nx * region.ny for region in self.regions.values())
        npoints_xlow = sum(
            (region.nx + (1 if region.connections["outer"] is None else 0)) * region.ny
            for region in self.regions.values()
        )
        npoints_ylow = sum(
            region.nx * (region.ny + (1 if region.connections["upper"] is None else 0))
            for region in self.regions.values()
        )
        npoints_corners = sum(
            (region.nx + (1 if region.connections["outer"] is None else 0))
            * (region.ny + (1 if region.connections["upper"] is None else 0))
            for region in self.regions.values()
        )

        region0 = list(self.regions.values())[0]
        mxn = {}
        myn = {}
        # markx and marky include guard cells, to avoid needing to add them as members
        # to MeshRegion objects (if they did not have guard cells, we would need to be
        # able to get them from self.getNeighbour("inner"), etc. in
        # MeshRegion.smoothnl_inner2()).
        markx = {
            region_name: MultiLocationArray(region.nx + 2, region.ny + 2)
            for region_name, region in self.regions.items()
        }
        marky = {
            region_name: MultiLocationArray(region.nx + 2, region.ny + 2)
            for region_name, region in self.regions.items()
        }
        for i in range(50):
            for region_name, region in self.regions.items():
                mxn[region_name], myn[region_name] = region.smoothnl_inner1(varname)

            if getattr(region0, varname).centre is not None:
                mean_mxn = sum(x.centre.sum() for x in mxn.values()) / npoints_centre
                mean_myn = sum(x.centre.sum() for x in myn.values()) / npoints_centre
                for region_name, region in self.regions.items():
                    this_markx = 0.5 * mxn[region_name].centre / mean_mxn
                    this_markx = numpy.where(this_markx < 1.0, this_markx, 1.0)
                    markx[region_name].centre[1:-1, 1:-1] = this_markx

                    this_marky = 0.5 * myn[region_name].centre / mean_myn
                    this_marky = numpy.where(this_marky < 1.0, this_marky, 1.0)
                    marky[region_name].centre[1:-1, 1:-1] = this_marky

                    if region.connections["inner"] is not None:
                        markx[region.connections["inner"]].centre[
                            -1, 1:-1
                        ] = this_markx[0, :]
                        marky[region.connections["inner"]].centre[
                            -1, 1:-1
                        ] = this_marky[0, :]
                    if region.connections["outer"] is not None:
                        markx[region.connections["outer"]].centre[0, 1:-1] = this_markx[
                            -1, :
                        ]
                        marky[region.connections["outer"]].centre[0, 1:-1] = this_marky[
                            -1, :
                        ]
                    if region.connections["lower"] is not None:
                        markx[region.connections["lower"]].centre[
                            1:-1, -1
                        ] = this_markx[:, 0]
                        marky[region.connections["lower"]].centre[
                            1:-1, -1
                        ] = this_marky[:, 0]
                    if region.connections["upper"] is not None:
                        markx[region.connections["upper"]].centre[1:-1, 0] = this_markx[
                            :, -1
                        ]
                        marky[region.connections["upper"]].centre[1:-1, 0] = this_marky[
                            :, -1
                        ]
                    if numpy.any(markx[region_name].centre[1:-1, 1:-1] > 1.0):
                        raise ValueError(f"{markx[region_name].centre[1:-1, 1:-1]}")
                    if numpy.any(markx[region_name].centre[1:-1, -1] > 1.0):
                        raise ValueError(f"{markx[region_name].centre[1:-1, -1]}")

            if getattr(region0, varname).xlow is not None:
                mean_mxn = sum(x.xlow.sum() for x in mxn.values()) / npoints_xlow
                mean_myn = sum(x.xlow.sum() for x in myn.values()) / npoints_xlow
                for region_name, region in self.regions.items():
                    this_markx = 0.5 * mxn[region_name].xlow / mean_mxn
                    this_markx = numpy.where(this_markx < 1.0, this_markx, 1.0)
                    markx[region_name].xlow[1:-1, 1:-1] = this_markx

                    this_marky = 0.5 * myn[region_name].xlow / mean_myn
                    this_marky = numpy.where(this_marky < 1.0, this_marky, 1.0)
                    marky[region_name].xlow[1:-1, 1:-1] = this_marky

                    if region.connections["inner"] is not None:
                        markx[region.connections["inner"]].xlow[-1, 1:-1] = this_markx[
                            0, :
                        ]
                        marky[region.connections["inner"]].xlow[-1, 1:-1] = this_marky[
                            0, :
                        ]
                    if region.connections["outer"] is not None:
                        markx[region.connections["outer"]].xlow[0, 1:-1] = this_markx[
                            -1, :
                        ]
                        marky[region.connections["outer"]].xlow[0, 1:-1] = this_marky[
                            -1, :
                        ]
                    if region.connections["lower"] is not None:
                        markx[region.connections["lower"]].xlow[1:-1, -1] = this_markx[
                            :, 0
                        ]
                        marky[region.connections["lower"]].xlow[1:-1, -1] = this_marky[
                            :, 0
                        ]
                    if region.connections["upper"] is not None:
                        markx[region.connections["upper"]].xlow[1:-1, 0] = this_markx[
                            :, -1
                        ]
                        marky[region.connections["upper"]].xlow[1:-1, 0] = this_marky[
                            :, -1
                        ]

            if getattr(region0, varname).ylow is not None:
                mean_mxn = sum(x.ylow.sum() for x in mxn.values()) / npoints_ylow
                mean_myn = sum(x.ylow.sum() for x in myn.values()) / npoints_ylow
                for region_name, region in self.regions.items():
                    this_markx = 0.5 * mxn[region_name].ylow / mean_mxn
                    this_markx = numpy.where(this_markx < 1.0, this_markx, 1.0)
                    markx[region_name].ylow[1:-1, 1:-1] = this_markx

                    this_marky = 0.5 * myn[region_name].ylow / mean_myn
                    this_marky = numpy.where(this_marky < 1.0, this_marky, 1.0)
                    marky[region_name].ylow[1:-1, 1:-1] = this_marky

                    if region.connections["inner"] is not None:
                        markx[region.connections["inner"]].ylow[-1, 1:-1] = this_markx[
                            0, :
                        ]
                        marky[region.connections["inner"]].ylow[-1, 1:-1] = this_marky[
                            0, :
                        ]
                    if region.connections["outer"] is not None:
                        markx[region.connections["outer"]].ylow[0, 1:-1] = this_markx[
                            -1, :
                        ]
                        marky[region.connections["outer"]].ylow[0, 1:-1] = this_marky[
                            -1, :
                        ]
                    if region.connections["lower"] is not None:
                        markx[region.connections["lower"]].ylow[1:-1, -1] = this_markx[
                            :, 0
                        ]
                        marky[region.connections["lower"]].ylow[1:-1, -1] = this_marky[
                            :, 0
                        ]
                    if region.connections["upper"] is not None:
                        markx[region.connections["upper"]].ylow[1:-1, 0] = this_markx[
                            :, -1
                        ]
                        marky[region.connections["upper"]].ylow[1:-1, 0] = this_marky[
                            :, -1
                        ]

            if getattr(region0, varname).corners is not None:
                mean_mxn = sum(x.corners.sum() for x in mxn.values()) / npoints_corners
                mean_myn = sum(x.corners.sum() for x in myn.values()) / npoints_corners
                for region_name, region in self.regions.items():
                    this_markx = 0.5 * mxn[region_name].corners / mean_mxn
                    this_markx = numpy.where(this_markx < 1.0, this_markx, 1.0)
                    markx[region_name].corners[1:-1, 1:-1] = this_markx

                    this_marky = 0.5 * myn[region_name].corners / mean_myn
                    this_marky = numpy.where(this_marky < 1.0, this_marky, 1.0)
                    marky[region_name].corners[1:-1, 1:-1] = this_marky

                    if region.connections["inner"] is not None:
                        markx[region.connections["inner"]].corners[
                            -1, 1:-1
                        ] = this_markx[0, :]
                        marky[region.connections["inner"]].corners[
                            -1, 1:-1
                        ] = this_marky[0, :]
                    if region.connections["outer"] is not None:
                        markx[region.connections["outer"]].corners[
                            0, 1:-1
                        ] = this_markx[-1, :]
                        marky[region.connections["outer"]].corners[
                            0, 1:-1
                        ] = this_marky[-1, :]
                    if region.connections["lower"] is not None:
                        markx[region.connections["lower"]].corners[
                            1:-1, -1
                        ] = this_markx[:, 0]
                        marky[region.connections["lower"]].corners[
                            1:-1, -1
                        ] = this_marky[:, 0]
                    if region.connections["upper"] is not None:
                        markx[region.connections["upper"]].corners[
                            1:-1, 0
                        ] = this_markx[:, -1]
                        marky[region.connections["upper"]].corners[
                            1:-1, 0
                        ] = this_marky[:, -1]

            changes = []
            tmp = {}
            for region_name, region in self.regions.items():
                this_tmp, change = region.smoothnl_inner2(
                    varname, markx[region_name], marky[region_name]
                )
                changes.append(change)
                tmp[region_name] = this_tmp

            change = max(changes)

            # Need to update the variables after calculating all tmp values because the
            # variable values from neighbouring regions are used in
            # region.smoothnl_inner2()
            for region_name, region in self.regions.items():
                setattr(region, varname, tmp[region_name])

            print(f"Smoothing {varname} {i}: change={change}", flush=True)

            if change < 1.0e-3:
                break

    def plotGridLines(self, **kwargs):
        from matplotlib import pyplot
        from cycler import cycle

        colors = cycle(pyplot.rcParams["axes.prop_cycle"].by_key()["color"])

        for region in self.regions.values():
            c = next(colors)
            label = region.myID
            for i in range(region.nx):
                pyplot.plot(
                    region.Rxy.centre[i, :],
                    region.Zxy.centre[i, :],
                    c=c,
                    label=label,
                    **kwargs,
                )
                label = None
            label = region.myID
            for j in range(region.ny):
                pyplot.plot(
                    region.Rxy.centre[:, j],
                    region.Zxy.centre[:, j],
                    c=c,
                    label=None,
                    **kwargs,
                )
                label = None
        l = pyplot.legend()
        l.set_draggable(True)

    def plotPoints(
        self,
        xlow=False,
        ylow=False,
        corners=False,
        markers=None,
        ax=None,
        plot_types="scatter",
        legend=True,
        **kwargs,
    ):
        from matplotlib import pyplot
        from cycler import cycle

        if isinstance(plot_types, str):
            plot_types = [plot_types]

        colors = cycle(pyplot.rcParams["axes.prop_cycle"].by_key()["color"])

        if markers is None:
            markers = ["x"]
            if xlow:
                markers.append("1")
            if ylow:
                markers.append("2")
            if corners:
                markers.append("+")
        try:
            markers[0]
        except TypeError:
            markers = list(markers)

        if ax is None:
            fig, ax = pyplot.subplots(1)
        else:
            fig = ax.figure

        for region in self.regions.values():
            c = next(colors)

            if "scatter" in plot_types:
                m = iter(markers)
                ax.scatter(
                    region.Rxy.centre,
                    region.Zxy.centre,
                    marker=next(m),
                    c=c,
                    label=region.myID,
                    **kwargs,
                )
                if xlow:
                    ax.scatter(
                        region.Rxy.xlow, region.Zxy.xlow, marker=next(m), c=c, **kwargs
                    )
                if ylow:
                    ax.scatter(
                        region.Rxy.ylow, region.Zxy.ylow, marker=next(m), c=c, **kwargs
                    )
                if corners:
                    ax.scatter(
                        region.Rxy.corners,
                        region.Zxy.corners,
                        marker=next(m),
                        c=c,
                        **kwargs,
                    )
            if "radial" in plot_types:
                R = numpy.empty([2 * region.nx + 1, region.ny])
                R[1::2, :] = region.Rxy.centre
                R[::2, :] = region.Rxy.xlow
                Z = numpy.empty([2 * region.nx + 1, region.ny])
                Z[1::2, :] = region.Zxy.centre
                Z[::2, :] = region.Zxy.xlow
                lines = ax.plot(R, Z, linestyle="-", c=c)
                lines[0].set_label(region.myID)
                if ylow:
                    R = numpy.empty([2 * region.nx + 1, region.ny + 1])
                    R[1::2, :] = region.Rxy.ylow
                    R[::2, :] = region.Rxy.corners
                    Z = numpy.empty([2 * region.nx + 1, region.ny + 1])
                    Z[1::2, :] = region.Zxy.ylow
                    Z[::2, :] = region.Zxy.corners
                    ax.plot(R, Z, linestyle="--", c=c)
            if "poloidal" in plot_types:
                R = numpy.empty([region.nx, 2 * region.ny + 1])
                R[:, 1::2] = region.Rxy.centre
                R[:, ::2] = region.Rxy.ylow
                Z = numpy.empty([region.nx, 2 * region.ny + 1])
                Z[:, 1::2] = region.Zxy.centre
                Z[:, ::2] = region.Zxy.ylow
                lines = ax.plot(R.T, Z.T, linestyle="-", c=c, label=region.myID)
                lines[0].set_label(region.myID)
                if ylow:
                    R = numpy.empty([region.nx + 1, 2 * region.ny + 1])
                    R[:, 1::2] = region.Rxy.xlow
                    R[:, ::2] = region.Rxy.corners
                    Z = numpy.empty([region.nx + 1, 2 * region.ny + 1])
                    Z[:, 1::2] = region.Zxy.xlow
                    Z[:, ::2] = region.Zxy.corners
                    ax.plot(R.T, Z.T, linestyle="--", c=c)

        if legend:
            l = ax.legend()
            l.set_draggable(True)

        return fig, ax

    def plotPotential(self, *args, **kwargs):
        """
        Plot the flux function psi. Passes through to self.equilibrium.plotPotential.
        """
        return self.equilibrium.plotPotential(*args, **kwargs)


def followPerpendicular(
    i, p0, psi0, *, f_R, f_Z, psivals, rtol=2.0e-8, atol=1.0e-8, **kwargs
):
    """
    Follow a line perpendicular to Bp from point p0 until psi_target is reached.
    """
    if i is not None:
        print(f"Following perpendicular: {i + 1}", end="\r", flush=True)

    # psi0 might be in somewhere in the range of psivals, rather than at one end
    if min(psivals) < psi0 < max(psivals):
        # Integrate in each direction, then put together
        # Partition into left and right halves
        if psivals[0] < psi0:
            left = [psi for psi in psivals if psi < psi0]
            right = [psi for psi in psivals if psi >= psi0]
        else:
            left = [psi for psi in psivals if psi >= psi0]
            right = [psi for psi in psivals if psi < psi0]

        return followPerpendicular(
            None, p0, psi0, f_R=f_R, f_Z=f_Z, psivals=left[::-1], rtol=rtol, atol=atol
        )[::-1] + followPerpendicular(
            None, p0, psi0, f_R=f_R, f_Z=f_Z, psivals=right, rtol=rtol, atol=atol
        )

    if abs(psivals[-1] - psi0) < abs(psivals[0] - psi0):
        # Closer at the end than the start -> Reverse
        return followPerpendicular(
            None,
            p0,
            psi0,
            f_R=f_R,
            f_Z=f_Z,
            psivals=psivals[::-1],
            rtol=rtol,
            atol=atol,
        )[::-1]
    psivals = psivals.copy()

    def f(psi, x):
        return (f_R(x[0], x[1]), f_Z(x[0], x[1]))

    psirange = (psi0, psivals[-1])
    # make sure rounding errors do not cause exception:
    if psirange[1] - psirange[0] > 0:
        # psi increasing in this interval
        if psivals[0] < psirange[0] and psirange[0] - psivals[0] < 1.0e-15 * numpy.abs(
            psirange[0]
        ):
            # rounding error present, reset psivals[0]
            psivals[0] = psirange[0]
    else:
        # psi decreasing in this interval
        if psivals[0] > psirange[0] and psivals[0] - psirange[0] < 1.0e-15 * numpy.abs(
            psirange[0]
        ):
            # rounding error present, reset psivals[0]
            psivals[0] = psirange[0]
    try:
        solution = solve_ivp(
            f,
            psirange,
            tuple(p0),
            t_eval=psivals,
            rtol=rtol,
            atol=atol,
            vectorized=True,
        )
    except ValueError:
        print(psirange, psivals)
        raise

    return [Point2D(*p) for p in solution.y.T]


class BoutMesh(Mesh):
    """
    Implementation of :class:`Mesh <hypnotoad.core.mesh.Mesh>` for BOUT++ grids.

    Handles writing of the grid file in the format expected by BOUT++, including
    creation of global arrays collected from the sub-regions contained in the
    :class:`MeshRegion <hypnotoad.core.mesh.MeshRegion>` objects.

    ``BoutMesh`` requires that the MeshRegion members fit together into a global
    logically-rectangular Mesh, with one of the topologies supported by BOUT++ (slab,
    limiter, single null, connected double null, or disconnected double null).
    """

    user_options_factory = Mesh.user_options_factory.add(
        # BoutMesh-specific options
        ###########################
    )

    def __init__(self, equilibrium, settings):

        super().__init__(equilibrium, settings)

        # nx, ny both include boundary guard cells
        eq_region0 = next(iter(self.equilibrium.regions.values()))
        self.nx = sum(eq_region0.nx)

        self.ny = sum(r.ny(0) for r in self.equilibrium.regions.values())

        self.ny_noguards = sum(r.ny_noguards for r in self.equilibrium.regions.values())
        self.ny_core = sum(
            r.ny_noguards for r in self.equilibrium.regions.values() if r.kind == "X.X"
        )

        self.fields_to_output = []
        self.arrayXDirection_to_output = []

        # Keep ranges of global indices for each region, separately from the MeshRegions,
        # because we don't want MeshRegion objects to depend on global indices
        if not all([r.nx == eq_region0.nx for r in self.equilibrium.regions.values()]):
            raise ValueError(
                "all regions should have same set of x-grid sizes to be compatible "
                "with a global, logically-rectangular grid"
            )
        x_sizes = [0] + list(eq_region0.nx)

        # Note: x_startinds includes the end: self.x_startinds[-1] = nx
        self.x_startinds = numpy.cumsum(x_sizes)
        x_regions = tuple(
            slice(self.x_startinds[i], self.x_startinds[i + 1], None)
            for i in range(len(self.x_startinds) - 1)
        )
        y_total = 0
        y_regions = {}
        self.y_regions_noguards = []
        for regname, region in self.equilibrium.regions.items():
            # all segments must have the same ny, i.e. same number of y-boundary guard
            # cells
            this_ny = region.ny(0)
            if not all(region.ny(i) == this_ny for i in range(region.nSegments)):
                raise ValueError(
                    "all radial segments in an equilibrium-region must have the same "
                    "ny (i.e.  same number of boundary guard cells) to be compatible "
                    "with a global, logically-rectangular grid"
                )

            y_total_new = y_total + this_ny
            self.y_regions_noguards.append(region.ny_noguards)
            reg_slice = slice(y_total, y_total_new, None)
            y_total = y_total_new
            y_regions[regname] = reg_slice

        self.region_indices = {}
        for reg_name in self.equilibrium.regions:
            for i in range(len(x_regions)):
                self.region_indices[
                    self.region_lookup[(reg_name, i)]
                ] = numpy.index_exp[x_regions[i], y_regions[reg_name]]

        # constant spacing in y for now
        if self.ny_core > 0:
            # If there is a core region, set dy consistent with 0<=y<2pi in the core
            self.dy_scalar = 2.0 * numpy.pi / self.ny_core
        else:
            # No core region, set dy consistent with 0<=y<2pi in whole domain
            self.dy_scalar = 2.0 * numpy.pi / self.ny_noguards

    def geometry(self):
        # Call geometry() method of base class
        super().geometry()

        def addFromRegions(name):
            # Collect a 2d field from the regions
            self.fields_to_output.append(name)
            f = MultiLocationArray(self.nx, self.ny)
            self.__dict__[name] = f
            f.attributes = next(iter(self.regions.values())).__dict__[name].attributes
            for region in self.regions.values():
                f_region = region.__dict__[name]

                if f.attributes != f_region.attributes:
                    raise ValueError(
                        "attributes of a field must be set consistently in every region"
                    )
                if f_region._centre_array is not None:
                    f.centre[self.region_indices[region.myID]] = f_region.centre
                if f_region._xlow_array is not None:
                    f.xlow[self.region_indices[region.myID]] = f_region.xlow[:-1, :]
                if f_region._ylow_array is not None:
                    f.ylow[self.region_indices[region.myID]] = f_region.ylow[:, :-1]
                if f_region._corners_array is not None:
                    f.corners[self.region_indices[region.myID]] = f_region.corners[
                        :-1, :-1
                    ]

            # Set 'bout_type' so it gets saved in the grid file
            f.attributes["bout_type"] = "Field2D"

        def addFromRegionsXArray(name):
            # Collects 1d arrays, defined on a grid in the x-direction (no y-variation)
            # Data taken from the first region in each y-group
            self.arrayXDirection_to_output.append(name)
            f = MultiLocationArray(self.nx, 1)
            f.centre[...] = float("nan")
            f.xlow[...] = float("nan")
            self.__dict__[name] = f
            f.attributes = self.y_groups[0][0].__dict__[name].attributes
            for y_group in self.y_groups:
                # Get values from first region in each y_group
                region = y_group[0]

                f_region = region.__dict__[name]

                if f.attributes != f_region.attributes:
                    raise ValueError(
                        "attributes of a field must be set consistently in every region"
                    )
                if f_region._centre_array is not None:
                    f.centre[self.region_indices[region.myID][0], :] = f_region.centre
                if f_region._xlow_array is not None:
                    f.xlow[self.region_indices[region.myID]] = f_region.xlow[:-1, :]
                if f_region._ylow_array is not None:
                    raise ValueError("Cannot have an x-direction array at ylow")
                if f_region._corners_array is not None:
                    raise ValueError("Cannot have an x-direction array at corners")

            # Set 'bout_type' so it gets saved in the grid file
            f.attributes["bout_type"] = "ArrayX"

        addFromRegions("Rxy")
        addFromRegions("Zxy")
        addFromRegions("psixy")
        addFromRegions("dx")
        addFromRegions("dy")
        addFromRegions("poloidal_distance")
        addFromRegionsXArray("total_poloidal_distance")
        addFromRegions("Brxy")
        addFromRegions("Bzxy")
        addFromRegions("Bpxy")
        addFromRegions("Btxy")
        addFromRegions("Bxy")
        addFromRegions("hy")
        # if not self.user_options.orthogonal:
        #    addFromRegions('beta')
        #    addFromRegions('eta')
        addFromRegions("dphidy")
        addFromRegions("ShiftTorsion")
        addFromRegions("zShift")
        addFromRegionsXArray("ShiftAngle")
        # I think IntShiftTorsion should be the same as sinty in Hypnotoad1.
        # IntShiftTorsion should never be used. It is only for some 'BOUT-06 style
        # differencing'. IntShiftTorsion is not written by Hypnotoad1, so don't write
        # here. /JTO 19/5/2019
        if not self.user_options.shiftedmetric:
            addFromRegions("sinty")
        addFromRegions("g11")
        addFromRegions("g22")
        addFromRegions("g33")
        addFromRegions("g12")
        addFromRegions("g13")
        addFromRegions("g23")
        addFromRegions("J")
        addFromRegions("g_11")
        addFromRegions("g_22")
        addFromRegions("g_33")
        addFromRegions("g_12")
        addFromRegions("g_13")
        addFromRegions("g_23")
        if self.user_options.curvature_type in (
            "curl(b/B) with x-y derivatives",
            "curl(b/B)",
        ):
            addFromRegions("curl_bOverB_x")
            addFromRegions("curl_bOverB_y")
            addFromRegions("curl_bOverB_z")
        addFromRegions("bxcvx")
        addFromRegions("bxcvy")
        addFromRegions("bxcvz")

        if hasattr(next(iter(self.equilibrium.regions.values())), "pressure"):
            addFromRegions("pressure")

    def writeArray(self, name, array, f):
        f.write(name, BoutArray(array.centre, attributes=array.attributes))
        f.write(
            name + "_xlow", BoutArray(array.xlow[:-1, :], attributes=array.attributes)
        )
        f.write(
            name + "_ylow", BoutArray(array.ylow[:, :-1], attributes=array.attributes)
        )

    def writeCorners(self, name, array, f):
        f.write(
            name + "_corners",
            BoutArray(array.corners[:-1, :-1], attributes=array.attributes),
        )

    def writeArrayXDirection(self, name, array, f):
        f.write(name, BoutArray(array.centre[:, 0], attributes=array.attributes))

    def writeGridfile(self, filename):
        from boututils.datafile import DataFile

        with DataFile(filename, create=True, format="NETCDF4") as f:
            # Save unique ID for grid file
            import uuid

            f.write_file_attribute("grid_id", str(uuid.uuid1()))
            f.write("nx", self.nx)
            # ny for BOUT++ excludes boundary guard cells
            f.write("ny", self.ny_noguards)
            f.write("y_boundary_guards", self.user_options.y_boundary_guards)
            f.write("curvature_type", self.user_options.curvature_type)
            f.write("Bt_axis", self.equilibrium.Bt_axis)

            if hasattr(self.equilibrium, "psi_axis"):
                f.write("psi_axis", self.equilibrium.psi_axis)
            if hasattr(self.equilibrium, "psi_bdry"):
                f.write("psi_bdry", self.equilibrium.psi_bdry)
            if hasattr(self.equilibrium, "psi_axis_gfile"):
                f.write("psi_axis_gfile", self.equilibrium.psi_axis_gfile)
            if hasattr(self.equilibrium, "psi_bdry_gfile"):
                f.write("psi_bdry_gfile", self.equilibrium.psi_bdry_gfile)

            # write the 2d fields
            for name in self.fields_to_output:
                self.writeArray(name, self.__dict__[name], f)

            # write corner positions, as these may be useful for plotting, etc.
            for name in ["Rxy", "Zxy"]:
                self.writeCorners(name, self.__dict__[name], f)

            if self.user_options.orthogonal:
                # Also write hy as "hthe" for backward compatibility
                self.writeArray("hthe", self.hy, f)

            # write the 1d fields
            for name in self.arrayXDirection_to_output:
                self.writeArrayXDirection(name, self.__dict__[name], f)

            # Write topology-setting indices for BoutMesh
            eq_region0 = next(iter(self.equilibrium.regions.values()))

            if len(self.x_startinds) == 2:
                # No separatrix in grid: self.x_startinds = [0, nx]
                if eq_region0.separatrix_radial_index == 0:
                    # SOL only
                    ixseps1 = -1
                    ixseps2 = -1
                else:
                    # core only
                    ixseps1 = self.nx
                    ixseps2 = self.nx
            elif len(self.x_startinds) == 3:
                # One separatrix: self.x_startinds = [0, ixseps, nx]
                ixseps1 = self.x_startinds[1]

                # note: ixseps2 may be changed below for cases where the two separatrices
                # are in the same radial location
                ixseps2 = self.nx
            elif len(self.x_startinds) == 4:
                # Two separatrices
                if self.equilibrium.double_null_type == "lower":
                    ixseps1 = self.x_startinds[1]
                    ixseps2 = self.x_startinds[2]
                elif self.equilibrium.double_null_type == "upper":
                    ixseps1 = self.x_startinds[2]
                    ixseps2 = self.x_startinds[1]
                else:
                    raise ValueError(
                        'Expected either double_null_type=="lower" or '
                        'double_null_type="upper" when there are two separatrices.'
                    )
            else:
                raise ValueError("More than two separatrices not supported by BoutMesh")

            if len(self.y_regions_noguards) == 1:
                # No X-points
                jyseps1_1 = -1
                jyseps2_1 = self.ny // 2
                ny_inner = self.ny // 2
                jyseps1_2 = self.ny // 2
                jyseps2_2 = self.ny
            elif len(self.y_regions_noguards) == 2:
                raise ValueError("Unrecognized topology with 2 y-regions")
            elif len(self.y_regions_noguards) == 3:
                # single-null
                jyseps1_1 = self.y_regions_noguards[0] - 1
                jyseps2_1 = self.ny // 2
                ny_inner = self.ny // 2
                jyseps1_2 = self.ny // 2
                jyseps2_2 = sum(self.y_regions_noguards[:2]) - 1
            elif len(self.y_regions_noguards) == 4:
                # single X-point with all 4 legs ending on walls
                jyseps1_1 = self.y_regions_noguards[0] - 1
                jyseps2_1 = jyseps1_1
                ny_inner = sum(self.y_regions_noguards[:2])
                jyseps2_2 = sum(self.y_regions_noguards[:3]) - 1
                jyseps1_2 = jyseps2_2

                # for BoutMesh topology, this is equivalent to 2 X-points on top of each
                # other, so there are 2 separatrices, in the same radial location
                ixseps2 = ixseps1
            elif len(self.y_regions_noguards) == 5:
                raise ValueError("Unrecognized topology with 5 y-regions")
            elif len(self.y_regions_noguards) == 6:
                # double-null
                jyseps1_1 = self.y_regions_noguards[0] - 1
                jyseps2_1 = sum(self.y_regions_noguards[:2]) - 1
                ny_inner = sum(self.y_regions_noguards[:3])
                jyseps1_2 = sum(self.y_regions_noguards[:4]) - 1
                jyseps2_2 = sum(self.y_regions_noguards[:5]) - 1

                if ixseps2 == self.nx:
                    # this is a connected-double-null configuration, with two
                    # separatrices in the same radial location
                    ixseps2 = ixseps1

            f.write("ixseps1", ixseps1)
            f.write("ixseps2", ixseps2)
            f.write("jyseps1_1", jyseps1_1)
            f.write("jyseps2_1", jyseps2_1)
            f.write("ny_inner", ny_inner)
            f.write("jyseps1_2", jyseps1_2)
            f.write("jyseps2_2", jyseps2_2)

            # Create poloidal coordinate (single-valued everywhere, includes y-boundary
            # cells)
            y = MultiLocationArray(self.nx, self.ny)
            y.centre[:, 1:] = numpy.cumsum(self.dy.centre, axis=1)[:, :-1]
            # Set xlow from x=0 entries of centre because xlow and centre have different
            # x-sizes, but y is constant in x so only actually need values from a single
            # x-index.
            y.xlow = y.centre[0, numpy.newaxis, :]
            y.ylow[:, :-1] = y.centre - 0.5 * self.dy.centre[:, :]
            y.ylow[:, -1] = y.centre[:, -1] + 0.5 * self.dy.centre[:, -1]
            y.attributes["bout_type"] = "Field2D"
            self.writeArray("y-coord", y, f)

            # Create poloidal coordinate which goes from 0 to 2pi in the core region
            theta = deepcopy(y)
            myg = self.user_options.y_boundary_guards
            for t in [theta.centre, theta.xlow, theta.ylow]:
                # Make zero of theta half a point before the start of the core region
                t -= theta.ylow[0, numpy.newaxis, jyseps1_1 + myg + 1, numpy.newaxis]
                if jyseps2_1 != jyseps1_2:
                    # Has second divertor, subtract y-increment in upper divertor legs
                    # from outer regions to make theta continuous in the core
                    # Set from x=0 entries of centre because xlow and ylow have different
                    # x-sizes, but y is constant in x so only actually need values from a
                    # single x-index.
                    t[:, ny_inner + 2 * myg :] -= (
                        theta.ylow[
                            0, numpy.newaxis, jyseps1_2 + 3 * myg + 1, numpy.newaxis
                        ]
                        - theta.ylow[
                            0, numpy.newaxis, jyseps2_1 + myg + 1, numpy.newaxis
                        ]
                    )
            theta.attributes["bout_type"] = "Field2D"
            self.writeArray("theta", theta, f)

            # Create straight-field-line poloidal coordinate which goes from 0 to 2pi in
            # the core region and is proportional to zShift
            chi = 2.0 * numpy.pi * self.zShift / self.ShiftAngle
            # need to add ylow values separately because ShiftAngle does not have a ylow
            # member
            chi.ylow = 2.0 * numpy.pi * self.zShift.ylow / self.ShiftAngle.centre
            # set to NaN in divertor leg regions where chi is not valid
            for c in [chi.centre, chi.xlow, chi.ylow]:
                c[:, : jyseps1_1 + 1] = float("nan")
                c[:, jyseps2_1 + 1 : jyseps1_2 + 1] = float("nan")
                c[:, jyseps2_2 + 1 :] = float("nan")
            chi.attributes["bout_type"] = "Field2D"
            self.writeArray("chi", chi, f)

            # BOUT++ ParallelTransform that metrics are compatible with
            if self.user_options.shiftedmetric:
                # Toroidal coordinates with shifts to calculate parallel derivatives
                f.write_file_attribute("parallel_transform", "shiftedmetric")
            else:
                # Field-aligned coordinates
                f.write_file_attribute("parallel_transform", "identity")

            # Save hypnotoad_inputs as a variable rather than an attribute because they
            # are long and attributes are printed by 'ncdump -h' or by ncdump when
            # looking at a different variable, which would be inconvenient. It is not
            # likely that we need to load the hypnotoad inputs in BOUT++, so no reason
            # to save as an attribute.
            #
            # Some options may be duplicated because of saving both equilibrium
            # and mesh options, but no other way to make sure to get everything
            # from subclasses (i.e. TokamakEquilibrium)
            inputs_string = (
                "Equilibrium options:\n"
                + self.equilibrium.user_options.as_table()
                + "\nNonorthogonal equilibrium options:\n"
                + self.equilibrium.nonorthogonal_options.as_table()
                + "\nMesh options:\n"
                + self.user_options.as_table()
            )
            f.write("hypnotoad_inputs", inputs_string)
            options_dict = dict(self.equilibrium.user_options)
            options_dict.update(self.equilibrium.nonorthogonal_options)
            options_dict.update(self.user_options)
            f.write("hypnotoad_inputs_yaml", yaml.dump(options_dict))

            f.write_file_attribute("hypnotoad_version", self.version)
            if self.git_hash is not None:
                f.write_file_attribute("hypnotoad_git_hash", self.git_hash)
                f.write_file_attribute(
                    "hypnotoad_git_diff",
                    self.git_diff if self.git_diff is not None else "",
                )

            if hasattr(self.equilibrium, "geqdsk_filename"):
                # If grid was created from a geqdsk file, save the file name
                f.write_file_attribute(
                    "hypnotoad_geqdsk_filename", self.equilibrium.geqdsk_filename
                )

            if hasattr(self.equilibrium, "geqdsk_input"):
                # If grid was created from a geqdsk file, save the file contents
                #
                # Write as string variable and not attribute because the string will be
                # long and attributes are printed by 'ncdump -h' or by ncdump when
                # looking at a different variable, which would be inconvenient. It is
                # not likely that we need to load the geqdsk file contents in BOUT++, so
                # no reason to save as an attribute.
                f.write(
                    "hypnotoad_input_geqdsk_file_contents",
                    self.equilibrium.geqdsk_input,
                )

            # save Python and module versions to enable reproducibility
            f.write("Python_version", sys.version)
            f.write("module_versions", module_versions_formatted())

    def plot2D(self, f, title=None):
        from matplotlib import pyplot

        try:
            vmin = f.min()
            vmax = f.max()
            if vmin == vmax:
                vmin -= 0.1
                vmax += 0.1

            for region, indices in zip(
                self.regions.values(), self.region_indices.values()
            ):
                pyplot.pcolor(
                    region.Rxy.corners,
                    region.Zxy.corners,
                    f[indices],
                    vmin=vmin,
                    vmax=vmax,
                )

            pyplot.colorbar()
        except NameError:
            raise NameError(
                "Some variable has not been defined yet: have you called "
                "Mesh.geometry()?"
            )
