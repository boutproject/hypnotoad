Construction of :class:`FineContour <hypnotoad.core.equilibrium.FineContour>` objects
=====================================================================================

:class:`FineContour <hypnotoad.core.equilibrium.FineContour>` objects are
created initially from :class:`PsiContour
<hypnotoad.core.equilibrium.PsiContour>` objects; this is largely for
historical reasons, as :class:`PsiContour
<hypnotoad.core.equilibrium.PsiContour>` was introduced first, and
:class:`FineContour <hypnotoad.core.equilibrium.FineContour>` added later to
improve accuracy and robustness.

After the :class:`FineContour <hypnotoad.core.equilibrium.FineContour>` is
created, or when it is modified (for example adding extra points beyond a
boundary, or decreasing the length when a wall intersects the initial contour),
it is then adjusted to ensure constant poloidal spacing of its points, using:

.. automethod:: hypnotoad.core.equilibrium.FineContour.equaliseSpacing
   :noindex:

The :class:`FineContour <hypnotoad.core.equilibrium.FineContour>` objects only
depend on the :math:`\psi`-value of the contour and number of points (set by
``finecontour_Nfine``), not on the poloidal spacing settings for the grid, etc.
(the initial guess for the iteration may depend on these, but this dependence
is eliminated by the iteration, up to the ``finecontour_atol`` tolerance). This
independence should ensure the consistency of the calculated poloidal
distances, metric coefficients, etc. when input parameters, e.g. resolution,
are varied.
