#!/usr/bin/env python
#

import sys

if len(sys.argv) != 2:
    raise ValueError("Usage is {} geqdsk_file".format(sys.argv[0]))

filename = sys.argv[1]

from hypnotoad2 import tokamak

with open(filename, 'rt') as fh:
    eq = tokamak.read_geqdsk(fh)

# Make the regions, identifying single/double null etc.
    
eq.makeRegions(psinorm_pf=0.95, psinorm_sol=1.1)

try:
    import matplotlib.pyplot as plt
    
    eq.plotPotential(ncontours=40)
    for region in eq.regions.values():
        plt.plot([p.R for p in region.points], [p.Z for p in region.points], '-o')
    print("Close window to continue...")
    plt.show()
except:
    pass

# Create the mesh

from hypnotoad2.mesh import BoutMesh

mesh = BoutMesh(eq)
mesh.geometry()

eq.plotPotential(ncontours=40)
plt.plot(*eq.x_points[0], 'rx')
mesh.plotPoints(xlow=True, ylow=True, corners=True)
plt.show()

mesh.writeGridfile('bout.grd.nc')
