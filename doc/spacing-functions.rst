Spacing functions
=================

'Spacing functions' map a grid index to a coordinate value, and are used to
position the points radially, as a function of :math:`psi`, and poloidally, as
a function of poloidal distance along the contours. The actual grid spacing is
(in the limit of large numbers of grid points) the derivative of the spacing
function with respect to grid index.

The spacing functions are defined in terms of the grid index normalised by the
total number of points (in the radial or poloidal direction) in the global
grid. For the poloidal direction the normalisation factor can be scaled by
changing the ``N_norm_prefactor`` setting, but it is not clear why this would
ever be useful (the same effect should be achievable by scaling the
``*_poloidal_spacing_length`` settings). This normalisation makes the spacing
consistent when the grid resolution is changed. For example, if the resolution
in one direction is doubled then the cell face positions on the original grid
are also cell face positions on the new grid, so that each original cell is
divided in two.

The spacing functions used are described in more detail in the following
sections: :ref:`radial-grid:Radial grid`, :ref:`orthogonal-grid:Orthogonal
grid`, and :ref:`nonorthogonal-grid:Nonorthogonal grid`.

Some general comments can be made here.

.. note:: At the level of :class:`PsiContour
   <hypnotoad.core.equilibrium.PsiContour>`, no distinction is made between
   cell centre and cell face (or cell face and cell corner) grid points, so in
   a region of the grid where there are :math:`n_\mathrm{c}` cell centre grid
   points in a certain direction, there will be :math:`n=(2n_\mathrm{c}+1)`
   points belonging to the :class:`PsiContour
   <hypnotoad.core.equilibrium.PsiContour>` on and inside the boundaries.

Let the spacing function be denoted :math:`s(i)` in a region where the index
:math:`i` is in the range :math:`0\leq i\leq i_\mathrm{max}`, and the
coordinate :math:`s` is in the range :math:`0\leq s\leq L`. All spacing
functions must satisfy

.. math::
   \begin{eqnarray}
   s(0) &=& 0 \\
   s(i_\mathrm{max}) &=& L
   \end{eqnarray}

Some of the spacing functions also have parameters that control the gradient at
one or both ends. Where possible, the gradient should match between the spacing
functions in adjoining regions, so that the grid spacing is continuous between
them. Therefore the gradient is controlled by parameters that set

.. math::
   \frac{ds}{di_\mathrm{norm}} = N_\mathrm{norm}\frac{ds}{di}

(where :math:`i_\mathrm{norm} = i/N_\mathrm{norm}`, and :math:`N_\mathrm{norm}`
is a global parameter which is the same in all regions,
:math:`N_\mathrm{norm}=n_{x,\mathrm{global}}` or
:math:`N_\mathrm{norm}=n_{y,\mathrm{global}}`) so that if the same parameter is
used for both regions at the boundary where they join, the grid spacing is the
same at the boundary (in the limit of large numbers of grid points).

The radial spacing functions also make the second derivative vanish
:math:`d^2s/di^2=0` at boundaries that are connected to another regions, so
that the grid spacing varies smoothly across the boundary, see
:ref:`radial-grid:Radial grid`.

.. _monotonic-requirement:

The spacing function must be monotonic within the range of indices where it is
used. Otherwise the grid would double back on itself (to put it another way,
some of the grid spacings would be negative), which makes no sense. Some of the
spacing functions are written to guarantee monotonicity by construction, as
noted below. Otherwise, monotonicity is checked and the spacing parameters need
to be adjusted if the condition is not satisfied.
