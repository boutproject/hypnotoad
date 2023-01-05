Grid file
=========

Scalar parameters
-----------------

Grid size and topology
++++++++++++++++++++++

``nx``

``ny``

``y_boundary_guards``

``ixseps1``, ``ixseps2``

``jyseps1_1``, ``jyseps2_1``, ``jyseps1_2``, ``jyseps2_2``

``ny_inner``

Equilibrium parameters
++++++++++++++++++++++

``Bt_axis``

``psi_axis``

``psi_bdry``

``psi_axis_gfile``

``psi_bdry_gfile``

Other options
+++++++++++++

``curvature_type``

1D arrays
---------

Poloidal coordinates
++++++++++++++++++++

BoutMesh writes three poloidal coordinates to the grid file:

.. note:: These coordinates are defined/created in BoutMesh because they
   require a global mesh, which is not required in Mesh where everything is
   defined only in terms of MeshRegions.

``y-coord``
        increments by ``dy`` between points and starts from zero at the
        beginning of the global grid. ``y`` includes boundary cells and is
        single-valued (at a given radial position) everywhere on the global
        grid. ``y`` has branch cuts adjacent to both X-points in the core, and
        adjacent to the X-point in the PFRs.

``theta``
        increments by ``dy`` between points and goes from 0 to 2pi in the core
        region. The lower inner divertor leg has negative values. The lower
        outer divertor leg has values >2pi. The upper inner leg (if it exists)
        has values increasing continuously from those in the inner SOL (these
        will overlap values in the outer core region). The outer upper leg (if
        it exists) has values continuous with those in the outer SOL (these
        will overlap values in the inner core region).

``chi``
        is a straight-field line poloidal coordinate proportional to the
        toroidal angle (i.e. to zShift). It goes from 0 to 2pi in the core, and
        is undefined on open field lines.

1D integral quantities
++++++++++++++++++++++

``total_poloidal_distance``

``ShiftAngle``

2D arrays
---------

Spatial positions
+++++++++++++++++

``Rxy``

``Zxy``

``Rxy_corners``, ``Zxy_corners``

Grid spacings
+++++++++++++

``dx``

``dy``

Magnetic field quantities
+++++++++++++++++++++++++

``psixy``

``Brxy``, ``Bzxy``

``Bpxy``, ``Btxy``

``Bxy``

Integral quantities
+++++++++++++++++++

``poloidal_distance``

``zShift``

``ShiftTorsion``

Coordinate related variables
++++++++++++++++++++++++++++

``hy``, ``hthe``

``dphidy``

Metric coefficients
+++++++++++++++++++

``g11``, ``g22``, ``g33``, ``g12``, ``g13``, ``g23``
        Note ``g12`` and ``g13`` vanish for orthogonal coordinates (although
        ``g13`` would be non-zero for globally field-aligned coordinates, which
        are not supported by hypnotoad).

``g_11``, ``g_22``, ``g_33``, ``g_12``, ``g_13``, ``g_23``
        Note ``g_12`` and ``g_13`` vanish for orthogonal coordinates (although
        they would both be non-zero for globally field-aligned coordinates,
        which are not supported by hypnotoad).

Jacobian
++++++++

``J``

Curvature
+++++++++

``curl_bOverB_x``, ``curl_bOverB_y``, ``curl_bOverB_z``

``bxcvx``, ``bxcvy``, ``bxcvz``

Equilibrium plama parameters
++++++++++++++++++++++++++++

``pressure``

Provenance tracking
-------------------

``hypnotoad_inputs``

``hypnotoad_inputs_yaml``

``hypnotoad_input_geqdsk_file_contents``

``Python_version``

``module_versions``
