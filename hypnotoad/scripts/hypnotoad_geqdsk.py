#!/usr/bin/env python
#
# This script generates a BOUT++ grid file from a GEQDSK equilibrium,
# and optionally a set of inputs in a YAML file
#
# For example:
#  $ ./tokamak_geqdsk file.geqdsk  geqdsk_cdn.yaml
#
#

import warnings


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("inputfile", nargs="?", default=None)
    parser.add_argument("--pdb", action="store_true", default=False)
    args = parser.parse_args()

    filename = args.filename
    if args.inputfile is not None:
        # Options yaml file
        import yaml

        with open(args.inputfile, "r") as inputfile:
            options = yaml.safe_load(inputfile)
    else:
        options = {}

    from ..cases import tokamak
    from ..core.mesh import BoutMesh

    possible_options = (
        [opt for opt in tokamak.TokamakEquilibrium.user_options_factory.defaults]
        + [
            opt
            for opt in tokamak.TokamakEquilibrium.nonorthogonal_options_factory.defaults
        ]
        + [opt for opt in BoutMesh.user_options_factory.defaults]
    )
    unused_options = [opt for opt in options if opt not in possible_options]
    if unused_options != []:
        raise ValueError(
            f"There were options in the input file that are not used: {unused_options}"
        )

    if args.pdb:
        import pdb

        pdb.set_trace()

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

    mesh = BoutMesh(eq, options)
    mesh.calculateRZ()

    if options.get("plot_mesh", False):
        try:
            import matplotlib.pyplot as plt

            ax = eq.plotPotential(ncontours=40)
            ax.plot(*eq.x_points[0], "rx")
            mesh.plotPoints(
                xlow=options.get("plot_xlow", True),
                ylow=options.get("plot_ylow", True),
                corners=options.get("plot_corners", True),
                ax=ax,
            )
            plt.show()
        except Exception as err:
            warnings.warn(str(err))

    mesh.geometry()

    mesh.writeGridfile(options.get("grid_file", "bout.grd.nc"))


if __name__ == "__main__":
    main()
