#!/usr/bin/env python3

from CerfonFreidbergGeometry import CerfonFreidbergSymmetric, geqdsk

eq = CerfonFreidbergSymmetric.CerfonFreidbergSymmetric()

eq.initByName("NSTX")

eq.calculatePsi0FromCurrent(1.0e6)

wall = [(0.1, -1.8), (1.67, -1.8), (1.67, 1.8), (0.1, 1.8)]

with open("test_connected-double-null.eqdsk", "w") as f:
    geqdsk.write(eq, f, wall=wall)
