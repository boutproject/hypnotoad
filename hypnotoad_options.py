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
Configuration information for grid generation

Default values set either here or in __init__ constructors of various classes.
"""

from options import Options

""" Options to be set by the user """
HypnotoadOptions = Options(
    ## General options for the mesh
        # Grid is orthogonal
        orthogonal = True,

        # Grid generated for paralleltransform=ShiftedMetric
        shiftedmetric = True,

        # Number of y-boundary guard cells
        y_boundary_guards = 0,

        # Expression to use to calculate curvature operator 'bxcv'
        # Possible values: 'Curl(b/B)' or 'bxkappa'
        curvature_type = 'bxkappa',

    ## Radial spacing options
        psi_spacing_separatrix_multiplier = None,

    ## Input parameters for poloidal spacing functions
        # Method to use for poloidal spacing function:
        #  - 'sqrt' for getSqrtPoloidalSpacingFunction
        #  - 'monotonic' for getMonotonicPoloidalDistanceFunc
        poloidal_spacing_method = 'sqrt',

        # spacing at the X-point end of a region
        xpoint_poloidal_spacing_length = None,

        # spacing at the wall end of a region (used for orthogonal grids)
        target_poloidal_spacing_length = None,

        # spacing at the X-point end of a region (used for non-orthogonal grids)
        nonorthogonal_xpoint_poloidal_spacing_length = None,

        # spacing at the wall end of a region (used for non-orthogonal grids)
        nonorthogonal_target_poloidal_spacing_length = None,

        # range near the X-point over which fixed poloidal position changes to
        # orthogonal position
        nonorthogonal_xpoint_poloidal_spacing_range = None,
        nonorthogonal_xpoint_poloidal_spacing_range_inner = None,
        nonorthogonal_xpoint_poloidal_spacing_range_outer = None,

        # range near the wall over which fixed poloidal position changes to
        # orthogonal position
        nonorthogonal_target_poloidal_spacing_range = None,
        nonorthogonal_target_poloidal_spacing_range_inner = None,
        nonorthogonal_target_poloidal_spacing_range_outer = None,

        nonorthogonal_radial_range_power = 1.,

        # method used to determine poloidal spacing of non-orthogonal grid
        nonorthogonal_spacing_method = 'combined',

        # Small increment in psi used to find vector along grad(psi) at end of
        # separatrix segment
        poloidal_spacing_delta_psi = None,

    ## Accuracy options for following Grad(psi)
        follow_perpendicular_rtol = None,
        follow_perpendicular_atol = None,

    ## Options for refining grids
        refine_width = None,
        refine_atol = None,

    ## Accuracy options for FineContour
        finecontour_Nfine = 1000,
        finecontour_atol = 1.e-12,

    ## Accuracy options for poloidal spacing functions
        sfunc_checktol = 1.e-13,

    ## Accuracy options for geometry checking
        geometry_rtol = 1.e-10,

    ## Switches for diagnostics to investigate when something is not converging
    # Info for FineContour.__init__
        finecontour_diagnose = False,
        poloidalfunction_diagnose = False,

    ## Switches for fudges to get rid of spikes near the X-point
        cap_Bp_ylow_xpoint = False,
    )

"""
Configuration to be used internally, values set by processing values in HypnotoadOptions
"""
HypnotoadInternalOptions = Options(
    # List of number of radial points in MeshRegions associated with this
    # EquilibriumRegion
    nx = None,

    # Number of y-points in this region
    ny = None,

    # Label of the kind of region:
    # - 'wall.X' starts at a wall and ends at an X-point
    # - 'X.wall' starts at an X-point and ends at a wall
    # - 'X.X' starts at an X-point and ends at an X-point
    kind = None,

    ### Parameters used by spacing functions
    ## Parameters for sqrt spacing function
    # Distance for polynomial part of spacing function at lower end
    sqrt_b_lower = None,

    # Distance for polynomial part of spacing function at upper end
    sqrt_b_upper = None,

    # Distance for sqrt part of spacing function (if used) at lower end
    sqrt_a_lower = None,

    # Distance for sqrt part of spacing function (if used) at upper end
    sqrt_a_upper = None,

    ## Parameters for monotonic spacing function
    # Distance for spacing function at lower end
    monotonic_d_lower = None,

    # Distance for spacing function at upper end
    monotonic_d_upper = None,

    # Distance for perpendicular spacing function at lower end
    perp_d_lower = None,

    # Distance for perpendicular spacing function at upper end
    perp_d_upper = None,

    # Normalization factor for number of points in contours, used to scale grid
    # spacing with total number of points, to keep functions consistent when
    # resolution is changed
    N_norm = None,

    # Distance for transition between fixed-poloidal-spacing grid and orthogonal grid
    # at the lower end. If 'None' then the value of monotonic_d_lower will be used instead.
    nonorthogonal_range_lower = None,
    nonorthogonal_range_lower_inner = None,
    nonorthogonal_range_lower_outer = None,

    # Distance for transition between fixed-poloidal-spacing grid and orthogonal grid
    # at the upper end. If 'None' then the value of monotonic_d_upper will be used instead.
    nonorthogonal_range_upper = None,
    nonorthogonal_range_upper_inner = None,
    nonorthogonal_range_upper_outer = None,

    )
