Tips and tricks for fixing problems
===================================

Initial strategy
----------------

The general advice when you are having problems generating a grid is to try and
start from something that builds successfully, even if it is far away from the
ial grid you want, with the aim of adjusting gradually so that it is possible
to solve one problem at a time.

While debugging, it is suggested to set ``poloidalfunction_diagnose = True`` to
get more information on errors.

Suggestions for 'simple' initial settings:

* Make an orthogonal grid ``orthogonal = true``
* Do not include y-boundary cells ``y_boundary_guards = 0``. It is recommended
  to set ``y_boundary_guards`` to the ``MYG`` value you will use with BOUT++,
  but setting it to a value greater than 0 makes the grid extend beyond the
  divertor targets, which can cause problems.
* Set ``psinorm_*`` settings so that the grid is as close as possible to the
  separatrix (or separatrices). This minimises the amount of the magnetic field
  that needs to be gridded.

For speed, it is useful to set ``nx_*`` to small numbers, e.g. 2, as this
reduces the number of ``PsiContour`` objects that are needed. (The ``ny_*``
resolution makes little difference to the grid generation time.) If
``finecontour_Nfine`` is greater than its default value of 100, it may also
help to reduce this value while debugging.

Once a grid builds with these 'simple' settings, change them to the desired
values in the order of the list (it may be useful to extend the range of
``psinorm_*`` in stages, to help isolate where the problem is coming from).
Finally increase the resolution to the desired values.

Nonorthogonal grids
+++++++++++++++++++

When building nonorthogonal grids, some extra simplifications may be useful (try these before increasing ``psinorm_*``):

* Set the range over which the 'fixed' spacing (as opposed to orthogonal
  spacing) is used to fairly large values, by increasing the
  ``nonorthogonal_*_poloidal_spacing_range*`` settings (maybe to values between
  1 and 10).

It may be useful to relax these simplifications first, or extend the range of
``psinorm_*`` first. Try in both orders if there are problems.

Useful parameters to experiment with
------------------------------------

With the parameters mentioned above fixed, probably the most useful ones to
experiment with are the ``xpoint_poloidal_spacing_length`` and the
``target_*_poloidal_spacing_length`` (which applies to the region in which the
problem is occuring). There are probably some 'easiest to grid' values which
are neither too large nor too small, at a guess these are likely to be around
those which keep the poloidal spacing of the grid uniform.

If 'refine point' failures are occuring, the following might also help:

* Set ``refine_timeout = None`` (or a suitably large value). At some places in
  some equilibria the 'refinement' may take a long time but still succeed.
* Including more methods for refinement by setting ``refine_methods =
  ['integrate+newton', 'integrate', 'line']`` may make refinement more robust
  (and should never hurt, as extra methods are only tried once previous ones
  have failed).
* When ``'line'`` is included in ``refine_methods``, increasing (or
  occasionally decreasing) ``refine_width`` may help.
* Changing ``refine_atol`` might make some difference.

For failures at walls, adjusting ``finecontour_extend_prefactor`` might help
[should link to where the role of ``finecontour_extend_prefactor`` is explained
here].

When setup of ``FineContour`` objects fails, decreasing
``finecontour_overdamping_factor`` should make the iteration more robust but
slower.

``poloidal_spacing_delta_psi`` and ``xpoint_offset`` might need to be increased
if there are problems around the X-point, particularly for lower resolution
input equilibria.

Nonorthogonal grids
+++++++++++++++++++

In addition to the parameters for an orthogonal grid, experiment with
``nonorthogonal_*_poloidal_spacing_length`` (again the ones which apply to the
region where the problem is).

Other tips for any grid
-----------------------

Double Null
+++++++++++

When generating double-null grids, if the input equilibrium is not connected
you may get strange results around the second X-point if the radial resolution
is high enough to see the fact that it is disconnected. The reason is that
different code branches are used to set up grids for connected and disconnected
equilibria. The 'connected' branch can tolerate small amounts of disconnection
and will still produce output for larger amounts, but the output will not be
good if the second X-point is more than ~1 grid spacing away from the first
separatrix. Instead, set ``nx_intersep`` to a non-zero value to use the
'disconnected' code branch.

x-boundary cells
++++++++++++++++

BOUT++'s ``nx`` includes boundary cells. For example if you want 32 radial grid
points in the core and 32 in the SOL with 2 boundary points (``MXG=2``, which
is the BOUT++ default), then you need to set ``nx_core = 34`` and ``nx_sol =
34``.

y-boundary cells
++++++++++++++++

BOUT++ can load grid files with or without y-boundary cells included . The
y-boundary cells are filled by extrapolation if not provided, which can lead to
problems (e.g. negative Jacobian), so it is recommended to provide boundary
cells in the grid file. This can be done by setting the option
``y_boundary_guards`` to a non-zero value (``MYG=2`` is the default number of
y-boundary points in BOUT++). However, creating a grid in the boundaries
requires following psi-contours past the wall. If your equilibrium file
provides psi on a grid that extends a fair way past the wall then this should
usually be fine, but often equilibrium files do not extend much past the wall.
Extrapolating psi past the provided input usually does not provide a good
result. Extending past the wall may also be a problem if there is a coil very
close to the target. In these cases, the problem is minimised when ``ny`` is
large, corresponding to a small poloidal grid spacing, because then the grid
does not extend as far past the target. So the suggested workflow for setting
non-zero y_boundary_guards would be:

    1. explore parameters at low resolution with ``y_boundary_guards = 0``
       (which is the default)
    2. when you are happy with other settings, increase the various ``ny_*``
       parameters, and see if the grid generates successfully with
       ``y_boundary_guards = 0`` - if the grid reaches very close to the edge
       of the psi-contour plot at this point then you may have trouble on the
       next step
    3. increase ``y_boundary_guards`` to the desired value - if this fails or
       produces a badly distorted grid in the boundaries, consider whether your
       model can run with fewer y-guard/y-boundary cells, which would make this
       gridding easier - if it is possible obtaining an equilibrium file which
       provides psi a bit further out is likely to be helpful
    4. if all else fails, edit the grid file to ensure the values of the metric
       components, etc. in the y-boundary cells are sensible (whatever that
       means to you, at least make values that should be positive be positive).

Errors when refining PsiContour or FineContour objects
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Creation of PsiContour or FineContour objects may fail when refining points. This often happens if the contour extends into a region where psi is not smooth enough. Some possible work-arounds:

* if a coil or the centre column is close to a target then when
  ``y_boundary_guards>0`` so that the contours have to extend past the target,
  the contours may extend to a region where psi is not smooth enough. One thing
  to try is increasing the relevant ``ny`` setting or decreasing
  ``target_*_poloidal_spacing_length`` to decrease the grid spacing, so that
  contours extend less far past the target. If this does not help, it may be
  the underlying ``FineContour`` object that is extending too far - in this
  case try decreasing ``finecontour_extend_prefactor``;
  ``finecontour_extend_prefactor`` should be reduced as little as possible
  because if it is too small, distances may be calculated by extrapolating past
  the end of the ``FineContour``, which will give inaccurate results.
* if warnings about ``FineContour: maximum iterations (200) exceeded...``` are
  produced, try adding ``'line'`` as the first entry of the refine_methods
  option.  The ``'line'`` option seems to be more robust for the small changes
  needed for the iteration in ``FineContour.equaliseSpacing()``. You could also
  try increasing ``finecontour_maxits`` and/or decreasing
  ``finecontour_overdamping_factor`` (this factor must be between 0 and 1 -
  closer to 1 is faster if it converges, but less robust, while closer to 0 is
  more robust).

Staggered grids
+++++++++++++++

If you use grids that are staggered in the y-direction (``CELL_YLOW``
quantities in BOUT++) and have problems with simulations - possibly showing up
as evolution being unphysically fast right from the beginning of a simulation -
a possible cause of problems is spikes in the ``*_ylow`` metric coefficients
due to ``Bp`` being very close to zero adjacent to the X-point. A hack to work
around this is to set ``cap_Bp_ylow_xpoint = True`` - this option limits the
minimum of ``Bp`` on the ylow points to the average of the adjacent cell-centre
points.

Nonorthogonal grids
+++++++++++++++++++

More tips for generating nonorthogonal grids are collected in their own
section: :ref:`nonorthogonal-tips:Nonorthogonal tips`.
