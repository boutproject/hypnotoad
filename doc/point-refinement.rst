Refinement of points to correct ψ
=================================

:class:`PsiContour <hypnotoad.core.equilibrium.PsiContour>` and
:class:`FineContour <hypnotoad.core.equilibrium.FineContour>` objects represent
flux surfaces -- surfaces of constant :math:`\psi` -- so all points beloning to
them should be at positions where :math:`\psi(R,Z)=\psi_\mathrm{C}`, where
:math:`\psi_\mathrm{C}` is the nominal :math:`\psi`-value of the contour and
:math:`\psi(R,Z)` is defined by the ``psi()`` method provided by
:class:`Equilibrum <hypnotoad.core.equilibrium.Equilibrium>`. However, the
initial placement of the points at each stage of the grid generation process is
given by the result of either an ODE integration or an interpolation or
extrapolation from the positions of some other points, so the point will not
(initially) lie exactly on the flux surface -- there will be some small error
in the :math:`\psi`-value.

To reduce the difference between the initial :math:`\psi` for each point and
the nominal :math:`\psi_\mathrm{C}` below a tolerance (set by ``refine_atol``),
:meth:`PsiContour.refinePoint()
<hypnotoad.core.equilibrium.PsiContour.refinePoint>` uses one of several
methods to refine the position of the point:

* :meth:`refinePointNewton()
  <hypnotoad.core.equilibrium.PsiContour.refinePointNewton>` - uses a Newton
  iteration. Usually converges quickly if the initial guess is close to the
  final point.
* :meth:`refinePointIntegrate()
  <hypnotoad.core.equilibrium.PsiContour.refinePointIntegrate>` - integrates
  along a vector :math:`(dR/d\psi,dZ/d\psi)` using ``scipy.integrate.solve_ivp``

.. automethod:: hypnotoad.core.equilibrium.PsiContour.refinePointNewton
   :noindex:

.. automethod:: hypnotoad.core.equilibrium.PsiContour.refinePointIntegrate
   :noindex:

.. automethod:: hypnotoad.core.equilibrium.PsiContour.refinePointLinesearch
   :noindex:

   The initial width is set by the ``refine_width`` setting.

A final option is available, ``'integrate+newton'`` which first uses
:meth:`refinePointIntegrate()
<hypnotoad.core.equilibrium.PsiContour.refinePointIntegrate>` to quickly
improve the initial guess, and then uses :meth:`refinePointNewton()
<hypnotoad.core.equilibrium.PsiContour.refinePointNewton>` to find an accurate
final position.

The choice of which method to use is set by the ``refine_method`` option
(``'newton'``, ``'integrate'``, ``'line'``, or ``'integrate+newton'``).
Multiple options can be passed in a list, in which case they will be tried in
order until one of them converges without error. This can be useful in case
different methods are more robust in different parts of the grid.
