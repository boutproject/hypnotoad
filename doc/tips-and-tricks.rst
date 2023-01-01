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
