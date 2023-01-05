MeshRegion notes
================

:class:`MeshRegion <hypnotoad.core.mesh.MeshRegion>` contains methods that
calculate geometric quantities for the standard BOUT++ locally field aligned
coordinate system. This was convenient, but ideally :class:`MeshRegion
<hypnotoad.core.mesh.MeshRegion>` should be generic, and the BOUT++ specific
functionality would be provided by a derived class, e.g. ``BoutMeshRegion``
that would be used by :class:`BoutMesh <hypnotoad.core.mesh.BoutMesh>`. In
practice, ``hypnotoad`` is only used with BOUT++ at present, so this is not an
issue at the moment.
