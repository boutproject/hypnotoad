Introduction and overview
=========================

Hypnotoad is the grid generator for `BOUT++
<https://bout-dev.readthedocs.io>`_. 'Grid generator' here means a tool that
takes a magnetic equilibrium (for a tokamak), generates a set of flux-surface
aligned grid points for a block-structured logical grid, and calculates
geometric quantities (metric components, etc.) required by BOUT++ at each
point.

Aims
----

* BOUT++ does not support arbitrary non-uniform grids - the spacing must vary
  slowly from one grid point to the next. This gives a constraint on grid
  generation, e.g. there should not be a jump in ``dx`` when crossing the
  separatrix.
* Grid generation should be robust: if a small change is made to the inputs,
  the changes to the grid should be small and smoothly varying.
* Grid generation should be reproducible: the inputs should be saved in a
  single file, from which the grid file can be re-generated without interactive
  input.

Overview
--------

The process for grid generation is:

* Analyse the equilibrium, identifying the X-points and central O-point.
  Determine whether equilibrium is single-null or double-null.
* Select the set of flux surfaces to be gridded.
* Divide the poloidal cross section into regions that can be covered by a
  logically rectangular grid, and are then subdivided into subregions. Each
  subregion joins to at most one other subregion at each edge.
* Determine the poloidal positions of the points along each flux surface in
  each subregion.
* Use the equilibrium magnetic field to calculate geometric quantities defined
  by the locally field-aligned coordinate system used by BOUT++ (see `this
  description
  <https://bout-dev.readthedocs.io/en/latest/user_docs/coordinates.html#field-aligned-coordinates>`_).

These stages are described in more detail in the following sections.

The main focus of hypnotoad and of this manual is on grids for diverted tokamak
configurations. For other toroidally symmetric configurations, see
:ref:`other-configurations:Other configurations`. For stellarator (or other 3d)
configurations, see instead the `Zoidberg grid generator
<https://github.com/boutproject/zoidberg>`_.

Orthogonal vs. non-orthogonal grids
-----------------------------------

Hypnotoad can generate both orthogonal and non-orthogonal grids. 'Orthogonal'
here means that the radial and poloidal coordinates are orthogonal (in the R-Z
plane).

.. note:: This is not the same as the simulation coordinate system being
   orthogonal. BOUT++ uses locally field aligned coordinates, and these are
   orthogonal only in an unsheared slab configuration.

On an orthogonal grid several metric components are exactly zero, which makes
numerical schemes simpler. An orthogonal grid also effectively shifts the
position of the wall at the divertor targets to make it orthogonal to the flux
surfaces, which prevents strongly distorted grids near the targest -- this is
not physically correct, but may make simulations run more easily. Due to flux
expansion around the X-point, the poloidal spacing is either very large near
the X-point, or very small at the radial boundaries of the grid adjacent to the
X-point -- usually a compromise is chosen where the poloidal grid spacing is
larger than would be ideal near the X-point but smaller than would be near the
radial boundaries.

Nonorthogonal grids have many more free parameters than orthogonal grids,
making it challenging to generically choose a good grid generation scheme.
Hypnotoad's :ref:`nonorthogonal grids <nonorthogonal-grid:Nonorthogonal grid>` are constructed to:

* follow the wall exactly at the divertor targets
* have 'sensible' poloidal grid spacing around the X-points, avoiding the
  compromise noted above that is forced by orthogonal grids
* be close to orthogonal away from the targets and X-points
