Poloidal distance calculations
==============================

Poloidal distance along flux surfaces is calculated on :class:`FineContour
<hypnotoad.core.equilibrium.FineContour>` objects using:

.. automethod:: hypnotoad.core.equilibrium.FineContour.calcDistance
   :noindex:

When the distance is needed at an arbitrary point on the contour (e.g. at the
positions of the :class:`PsiContour <hypnotoad.core.equilibrium.PsiContour>`
points), it is calculated based on the :class:`FineContour
<hypnotoad.core.equilibrium.FineContour>` distance using:

.. automethod:: hypnotoad.core.equilibrium.FineContour.getDistance
   :noindex:
