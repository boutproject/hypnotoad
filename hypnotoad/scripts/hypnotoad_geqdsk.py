#!/usr/bin/env python
#
# This script generates a BOUT++ grid file from a GEQDSK equilibrium,
# and optionally a set of inputs in a YAML file
#
# For example:
#  $ ./tokamak_geqdsk file.geqdsk  geqdsk_cdn.yaml
#
#

import sys
import warnings


def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        raise ValueError("Usage is {} geqdsk_file [options.yaml]".format(sys.argv[0]))

    filename = sys.argv[1]
    if len(sys.argv) == 3:
        # Options yaml file
        import yaml

        with open(sys.argv[2], "r") as inputfile:
            options = yaml.safe_load(inputfile)
    else:
        options = {}

    from ..cases import tokamak

    with open(filename, "rt") as fh:
        eq = tokamak.read_geqdsk(fh, settings=options, nonorthogonal_settings=options)

    if options.get("plot_regions", False):
        try:
            import matplotlib.pyplot as plt

            eq.plotPotential(ncontours=40)
            for region in eq.regions.values():
                plt.plot(
                    [p.R for p in region.points], [p.Z for p in region.points], "-o"
                )
            print("Close window to continue...")
            plt.show()
        except Exception as err:
            warnings.warn(str(err))

    # Create the mesh

    from ..core.mesh import BoutMesh

    mesh = BoutMesh(eq, settings=options)
    mesh.geometry()

    if options.get("plot_mesh", False):
        try:
            import matplotlib.pyplot as plt

            eq.plotPotential(ncontours=40)
            plt.plot(*eq.x_points[0], "rx")
            mesh.plotPoints(
                xlow=options.get("plot_xlow", True),
                ylow=options.get("plot_ylow", True),
                corners=options.get("plot_corners", True),
            )
            plt.show()
        except Exception as err:
            warnings.warn(str(err))

    mesh.writeGridfile(options.get("grid_file", "bout.grd.nc"))


if __name__ == "__main__":
    main()
