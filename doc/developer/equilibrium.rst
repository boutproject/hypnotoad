Equilibrium implementations
===========================

:class:`Equilibrium <hypnotoad.core.equilibrium.Equilibrium>` is an abstract
(or partially abstract) base class. It defines the interface to magnetic field
data and topology. An implementation (derived class) which reads magnetic
equilibrium input and analyses the topology is needed for each particular
configuration.
:class:`TokamakEquilibrium <hypnotoad.cases.tokamak.TokamakEquilibrium>`,
:class:`CircularEquilibrium <hypnotoad.cases.circular.CircularEquilibrium>`,
and
:class:`TORPEXMagneticField <hypnotoad.cases.torpex.TORPEXMagneticField>` exist
already.

Any implementation of :class:`Equilibrium
<hypnotoad.core.equilibrium.Equilibrium>` must provide:

* ``self.psi``: function which takes two arguments, {R,Z}, and returns the
  value of psi at that position.
* ``self.f_R``: function which takes two arguments, {R,Z}, and returns the R
  component of the vector :math:`\nabla\psi/|\nabla\psi|^2`.
* ``self.f_Z``: function which takes two arguments, {R,Z}, and returns the Z
  component of the vector :math:`\nabla\psi/|\nabla\psi|^2`.
* ``self.Bp_R``: function which takes two arguments, {R,Z}, and returns the R
  component of the poloidal magnetic field.
* ``self.Bp_Z``: function which takes two arguments, {R,Z}, and returns the Z
  component of the poloidal magnetic field.
* ``self.x_points``: list of :class:`Point2D
  <hypnotoad.core.equilibrium.Point2D>` objects giving the position of the
  X-points ordered from primary X-point (nearest the core) outward.
* ``self.psi_sep``: values of psi on the separatrices, ordered the same as
  ``self.x_points``.
* ``self.fpol``: poloidal current function; takes one argument, psi, and returns fpol
  (function such that B_toroidal = fpol/R).
* ``self.fpolprime``: psi-derivative of fpol.
* ``self.Rmin``, ``self.Rmax``, ``self.Zmin``, ``self.Zmax``: positions of the
  corners of a bounding box for the gridding.
* ``self.regions``: ``OrderedDict`` of
  :class:`EquilibriumRegion <hypnotoad.core.equilibrium.EquilibriumRegion>`
  objects that specify the topology of the equilibrium.
* ``self.wall``: list of :class:`Point2D <hypnotoad.core.equilibrium.Point2D>`
  objects giving vertices of polygon representing the wall, in anti-clockwise
  order; assumed to be closed so last element and first are taken to be
  connected.
