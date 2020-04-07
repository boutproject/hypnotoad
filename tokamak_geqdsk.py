#!/usr/bin/env python
#

import sys

if len(sys.argv) < 2 or len(sys.argv) > 3:
    raise ValueError("Usage is {} geqdsk_file [options.yaml]".format(sys.argv[0]))

filename = sys.argv[1]
if len(sys.argv) == 3:
    # Options yaml file
    import yaml
    with open(sys.argv[2], 'r') as inputfile:
        options = yaml.safe_load(inputfile)
else:
    options = {}

from hypnotoad2 import tokamak

with open(filename, 'rt') as fh:
    eq = tokamak.read_geqdsk(fh, options=options)

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
mesh.plotPoints(xlow = options.get("plot_xlow", True),
                ylow = options.get("plot_ylow", True),
                corners = options.get("plot_corners", True))
plt.show()

mesh.writeGridfile('bout.grd.nc')
