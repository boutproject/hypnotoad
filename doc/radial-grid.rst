Radial grid
===========

After the :class:`EquilibriumRegion
<hypnotoad.core.equilibrium.EquilibriumRegion>` objects are created, the radial
grids are constructed and stored in them.

Constructing the 'radial grid' means choosing the flux surfaces that are
gridded in the core, SOL, PFR, and (for disconnected double null
configurations) the region between the separatrices.

The parameters controlling the radial grid are:

* the number of radial grid points in the core (``nx_core``), SOL (``nx_sol``),
  and for disconnected double null inter-separatrix (``nx_inter_sep``)
* the 'inner' values of :math:`\psi` in the core (``psi_core`` or
  ``psinorm_core``) and PFRs (``psi_pf_lower``, ``psi_pf_upper`` or
  ``psinorm_pf_lower``, ``psinorm_pf_upper``)
* the 'outer' values of :math:`\psi` in the SOL (``psi_sol``, ``psi_sol_inner``
  or ``psinorm_sol``, ``psinorm_sol_inner``)
* ``psi_spacing_separatrix_multiplier`` which can be used to decrease (or
  increase) the grid spacing at the separatrices. Useful to increase resolution
  near separatrices so that the radial resolution is adequate near the
  X-points.

Note some of these settings are not used for single null configurations.

As the settings for the boundary value of :math:`\psi` are generally different
between core and PFR, and can also be different for inner and outer SOL regions
in a double null configuration, the grid values of :math:`\psi` are different
in different regions -- :math:`\psi` is not a 1d radial coordinate globally, for
the whole grid. This means that, for example, the SOL region connects at the
separatrix to the core and PFR(s), in different poloidal locations, which have
different :math:`\psi` grids. To ensure that the grid spacing is continuous and
smooth across the separatrix everywhere, the radial derivative of the grid
spacing is forced to vanish at the separatrices, so the radial gradient of the
grid spacing approaches zero on both sides of the boundary, ensuring the
spacing is smooth.

Technical details
-----------------

To make the radial grid spacing continuous and
smooth across the boundaries between regions (the separatrices) a gradient for
the :ref:`spacing function <spacing-functions:Spacing functions>` is imposed at
these boundaries, and the second derivative of the spacing function is
constrained to vanish. At the radial boundaries (where there is no connection
to another region) neither the gradient nor the second derivative of the
spacing function are constrained. All these requirements would be simple to
implement if the spacing function was a polynomial (they would reduce to some
simultaneous equations that need to be solved for the coefficents). However, it
was found that using a polynomial spacing function was not robust, as the
function could become :ref:`non-monotonic <monotonic-requirement>`. A more
complex algorithm has been implemented in :meth:`getSmoothMonotonicGridFunc()
<hypnotoad.core.equilibrium.Equilibrium.getSmoothMonotonicGridFunc>` that
guarantees the function is always monotonic, with the spacing function
involving the error function or trigonometric functions in some cases.
