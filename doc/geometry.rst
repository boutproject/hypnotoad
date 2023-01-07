Geometry
========

Comments on how geometrical quantities of the standard BOUT++ locally field
aligned coordinate system are calculated.

For orthogonal grids nothing depends on the grid as everything is calculated
from the position of a grid point (and poloidal distances using the
:class:`FineContour <hypnotoad.core.equilibrium.FineContour` objects).

For nonorthogonal grids, the calculation of the angle between poloidal and
radial directions uses the radial vector defined by the grid. Everything else
independent of the grid as for orthogonal grid. Inaccuracy should be small as
radial resolution is generally high compared to the length scale of the
magnetic equilibrium as the radial direction needs to capture the perpendicular
scales of the turbulence.
