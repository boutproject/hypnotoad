Geometry
========

If you are using the GUI, the geometrical quantities discussed here are
calculated (and then written to the grid file) when the 'Write Grid' button is
pressed.

Geometrical quantities for the standard BOUT++ locally field-aligned coordinate
system are calculated following the definitions in `the BOUT++ manual
<https://bout-dev.readthedocs.io/en/latest/user_docs/coordinates.html#field-aligned-coordinates>`_. 

.. note:: Only the *locally* field aligned coordinate system is supported. The
   'integrated shear' :math:`I` that would be needed for a globally field
   aligned coordinate system is not computed -- it is anyway not well defined
   in diverted tokamak topology as the reference location cannot be chosen in a
   consistent way for all regions (core, SOL and PFR).

For orthogonal grids, all the geometrical quantities can be defined in terms of:

* magnetic field components from the equilibrium and their :math:`R` and
  :math:`Z` derivatives. These can all be calculated at the grid points from
  the interpolating functions defined from the input data.
* the poloidal distance along the contours, which is used to calculate
  :math:`h_y = |\nabla y|` using the poloidal distances between cell faces, see
  :meth:`MeshRegion.calcHy() <hypnotoad.core.mesh.MeshRegion.calcHy>`.
* ``zShift``, (see `definition here
  <https://bout-dev.readthedocs.io/en/latest/user_docs/coordinates.html#zshift>`_),
  which is calculated by integrating :math:`B_t/RB_p`, evaluated again using
  the interpolated input data, along :class:`FineContour
  <hypnotoad.core.equilibrium.FineContour>`\s (see
  :meth:`MeshRegion.calcZShift() <hypnotoad.core.mesh.MeshRegion.calcZShift>`).

None of these depend on the BOUT++ grid itself (i.e. we are not using any
numerical derivatives calculated on a 'coarse grid'), so the calculation of the
geometric quantities at each grid point depends only on the R-Z position of the
grid point and is robustly independent of the parameters of the grid, so will
be calculated consistently and with no loss of accuracy when the grid
parameters are varied. The accuracy of the evaluation of the geometric
quantities only depends on the accuracy of the magnetic equilibrium input data,
and the number of points on fine contours ``finecontour_Nfine``.

One exception to the comments above is ``ShiftTorsion`` which does use a radial
derivative on the BOUT++ grid (as noted below in the nonorthogonal grids
paragraph, this should be accurate anyway). However ``ShiftTorsion`` is only
used in BOUT++'s' ``Curl()`` which is rarely (if ever) used. The calculation of
``ShiftTorsion`` should be verified before being used if it is ever needed.

For nonorthogonal grids, one additional quantity is required -- the angle
between the covariant and contravariant basis vectors :math:`\mathbf{e}_x` and
:math:`\nabla x` (which is the same as the angle between :math:`\mathbf{e}_y`
and :math:`\nabla y`) is also required. This angle is calculated by
approximating the :math:`x`-coordinate direction :math:`\mathbf{e}_x` with the
vector between the two grid points radially on either side of the point where
it is needed. This quantity does therefore depend on the BOUT++ grid. However,
the radial resolution is usually very high, since radially the perpendicular
variation of plasma turbulence has to be resolved (unlike the poloidal
direction where only parallel gradients have to be resolved), so this should be
a good approximation for any reasonable BOUT++ grid. This calculation is
performed in :meth:`MeshRegion.calcHy()
<hypnotoad.core.mesh.MeshRegion.calcHy>`.

**TODO:** document the expressions for metric coefficients, etc. in the
nonorthogonal coordinate system somewhere -- either here or in the BOUT++
manual.

Technical details
-----------------

The geometric quantities are calculated by various :class:`MeshRegion
<hypnotoad.core.mesh.MeshRegion>` methods that are called from
:meth:`Mesh.geometry() <hypnotoad.core.mesh.Mesh.geometry>`.
