Grid file
=========

If you are using the GUI, grid file is written when the 'Write Grid' button is
pressed.

Scalar parameters
-----------------

Grid size and topology
++++++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``nx``

     - Number of radial grid points in the global grid, including boundary cells.

   * - ``ny``

     - Number of poloidal grid points in the global grid, not including boundary
       cells.

   * - ``y_boundary_guards``

     - Number of poloidal boundary cells included in the grid at each divertor
       target.

   * - ``ixseps1``, ``ixseps2``

     - Radial grid indices of the separatrices. These give the index of the first
       grid point radially outside the separatrix. A single null grid has
       ``ixseps2 = nx``. A connected double null grid has ``ixseps1 = ixseps2``.
       For a disconnected double null, the primary X-point is the lower one if
       ``ixseps1 < ixseps2``, or the upper one if ``ixseps2 < ixseps1``.

   * - ``jyseps1_1``, ``jyseps2_1``, ``jyseps1_2``, ``jyseps2_2``

     - Poloidal grid indices of the X-points. These give the index of the last
       grid point poloidally before the poloidal position of the X-point is
       reached. For single null grids ``jyseps2_1 = jyseps1_2``.

   * - ``ny_inner``

     - For double null grids, gives the number of poloidal grid points before the
       upper target is reached, not including boundary cells.

Equilibrium parameters
++++++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``Bt_axis``

   * - ``psi_axis``

   * - ``psi_bdry``

   * - ``psi_axis_gfile``

   * - ``psi_bdry_gfile``

Other options
+++++++++++++

.. list-table::
   :widths: 30 70

   * - ``curvature_type``

1D arrays
---------

Poloidal coordinates
++++++++++++++++++++

BoutMesh writes three poloidal coordinates to the grid file:

.. note:: These coordinates are defined/created in BoutMesh because they
   require a global mesh, which is not required in Mesh where everything is
   defined only in terms of MeshRegions.

.. list-table::
   :widths: 30 70

   * - ``y-coord``

     -  Increments by ``dy`` between points and starts from zero at the
        beginning of the global grid. ``y`` includes boundary cells and is
        single-valued (at a given radial position) everywhere on the global
        grid. ``y`` has branch cuts adjacent to both X-points in the core, and
        adjacent to the X-point in the PFRs.

   * - ``theta``

     -  Increments by ``dy`` between points and goes from 0 to 2pi in the core
        region. The lower inner divertor leg has negative values. The lower
        outer divertor leg has values >2pi. The upper inner leg (if it exists)
        has values increasing continuously from those in the inner SOL (these
        will overlap values in the outer core region). The outer upper leg (if
        it exists) has values continuous with those in the outer SOL (these
        will overlap values in the inner core region).

   * - ``chi``

     -  Is a straight-field line poloidal coordinate proportional to the
        toroidal angle (i.e. to zShift). It goes from 0 to 2pi in the core, and
        is undefined on open field lines.

1D integral quantities
++++++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``total_poloidal_distance``

     - The total poloidal distance around a closed flux surface in the core.
       Not calculated on open flux surfaces.

   * - ``ShiftAngle``

     - The total toroidal angular displacement when following a field line one
       full poloidal turn around a closed flux surface. Not calculated on open
       flux surfaces.

2D arrays
---------

All 2D arrays are saved at the cell centre position, named with no suffix, e.g.
``Rxy``. For use with staggered grid codes, they are also saved at the 'lower'
cell face locations, with suffix ``_xlow``, e.g. ``Rxy_xlow`` for the
:math:`x`-direction cell faces and with suffix ``_ylow``, e.g. ``Rxy_ylow``,
for the :math:`y`-direction cell faces.

Spatial positions
+++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``Rxy``

     - Major radius.

   * - ``Zxy``

     - Height.

   * - ``Rxy_corners``, ``Zxy_corners``

     - Major radius and height of the lower-left corner of each grid cell. Not
       needed by BOUT++, but may be useful for post-processing.

   * - ``Rxy_lower_right_corners``, ``Zxy_lower_right_corners``,
       ``Rxy_upper_right_corners``, ``Zxy_upper_right_corners``,
       ``Rxy_upper_left_corners``, ``Zxy_upper_left_corners``

     - Major radius and height of the other three corners of each grid cell.
       Mostly redundant information with ``Rxy_corners`` and ``Zxy_corners``,
       but may make handling branch cuts and upper/outer boundaries more
       convenient. Not needed by BOUT++, but may be useful for post-processing.

Grid spacings
+++++++++++++

.. list-table::
   :widths: 30 70

   * - ``dx``

     - Coordinate spacing in the radial :math:`x` direction.

   * - ``dy``

     - Coordinate spacing in the poloidal :math:`y` direction.

Magnetic field quantities
+++++++++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``psixy``

     - Poloidal magnetic flux function, which is the poloidal magnetic flux divided by :math:`2\pi`.

   * - ``Brxy``, ``Bzxy``

     - Components of the magnetic field in the major-radial and vertical directions.

   * - ``Bpxy``, ``Btxy``

     - Components of the magnetic field in the poloidal and toroidal directions.

   * - ``Bxy``

     - Total magnetic field.

Integral quantities
+++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``poloidal_distance``

     - Poloidal distance (in metres) from the lower divertor target of each flux
       surface to the grid point (on open field lines), or from the poloidal
       location of the lower X-point (on closed field lines).

   * - ``zShift``

     - Toroidal displacement of a field line followed from some reference
       position to the poloidal location of the grid point.

   * - ``ShiftTorsion``

     - :math:`d^2\zeta/dxdy`, where :math:`zeta` is the toroidal angle. Only
       used in BOUT++ for the ``Curl()`` operator, which is rarely used. Note
       the calculation of this quantity has not been checked carefully, and
       should be verified if it is ever needed.

Coordinate related variables
++++++++++++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``hy``, ``hthe``

   * - ``dphidy``

Metric coefficients
+++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``g11``, ``g22``, ``g33``, ``g12``, ``g13``, ``g23``
        Note ``g12`` and ``g13`` vanish for orthogonal coordinates (although
        ``g13`` would be non-zero for globally field-aligned coordinates, which
        are not supported by hypnotoad).

   * - ``g_11``, ``g_22``, ``g_33``, ``g_12``, ``g_13``, ``g_23``
        Note ``g_12`` and ``g_13`` vanish for orthogonal coordinates (although
        they would both be non-zero for globally field-aligned coordinates,
        which are not supported by hypnotoad).

Jacobian
++++++++

.. list-table::
   :widths: 30 70

   * - ``J``

     - The Jacobian of the locally field aligned BOUT++ coordinate system.

Curvature
+++++++++

.. list-table::
   :widths: 30 70

   * - ``curl_bOverB_x``, ``curl_bOverB_y``, ``curl_bOverB_z``

     - Contravariant components (despite the slightly misleading variable
       names) of :math:`\nabla\times(\mathbf{b}/B)`, i.e.
       :math:`\nabla\times(\mathbf{b}/B)^x`,
       :math:`\nabla\times(\mathbf{b}/B)^y`, and
       :math:`\nabla\times(\mathbf{b}/B)^z`.

   * - ``bxcvx``, ``bxcvy``, ``bxcvz``

     - Contravariant components of the vector
       :math:`\frac{B}{2}\nabla\times\left(\frac{\mathbf{b}}{B}\right)`. Other
       forms (e.g. :math:`\mathbf{b}\times\mathbf{\kappa}`) could be
       implemented, for different settings of ``curvature_type``, but have not
       been implemented yet.

Equilibrium plama parameters
++++++++++++++++++++++++++++

.. list-table::
   :widths: 30 70

   * - ``pressure``

     - Pressure profile read from the geqdsk input file (if there was one).

Provenance tracking
-------------------

See :ref:`provenance-tracking:Provenance tracking`.
