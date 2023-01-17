#!/usr/bin/env python3

"""
This script is a hack to make a psi that gives a disconnected double null configuration.
The result is not a Grad-Shafranov equilibrium like the results from
CerfonFreidbergSymmetric.
"""

from CerfonFreidbergGeometry import CerfonFreidbergSymmetric, geqdsk

eq = CerfonFreidbergSymmetric.CerfonFreidbergSymmetric()

eq.initByName("NSTX")

eq.calculatePsi0FromCurrent(1.0e6)

orig_psi_symbolic = eq._getpsi
new_psi_symbolic = (
    lambda x, y, A, c1, c2, c3, c4, c5, c6, c7: orig_psi_symbolic(
        x, y, A, c1, c2, c3, c4, c5, c6, c7
    )
    + 1.0e-3 * y
)
eq._getpsi = new_psi_symbolic

wall = [(0.1, -1.8), (1.67, -1.8), (1.67, 1.8), (0.1, 1.8)]

with open("test_disconnected-double-null.eqdsk", "w") as f:
    geqdsk.write(eq, f, wall=wall)
