#!/usr/bin/env python
#
# This script generates a BOUT++ grid file for a circular, concentric flux surface
# geometry, optionally using a set of inputs in a YAML file
#
# For example:
#  $ hypnotoad_circular settings.yml
#

import gc
import warnings


def main(*, add_noise=None):
    """
    Read (optional) input file, and write a grid file

    Parameters
    ----------
    add_noise : int, optional
        Only intended for use in tests. If an int is passed, used as the seed for a
        random number generator adding small amounts of noise to the EquilibriumRegion
        points before generating the grid.
    """

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("inputfile", nargs="?", default=None)
    parser.add_argument("--pdb", action="store_true", default=False)
    parser.add_argument("--plot-regions", action="store_true", default=False)
    parser.add_argument("--plot-mesh", action="store_true", default=False)
    args = parser.parse_args()

    if args.inputfile is not None:
        # Options yaml file
        import yaml

        with open(args.inputfile, "r") as inputfile:
            options = yaml.safe_load(inputfile)
    else:
        options = {}

    from ..cases.circular import CircularEquilibrium
    from ..core.mesh import BoutMesh

    possible_options = (
        [opt for opt in CircularEquilibrium.user_options_factory.defaults]
        + [opt for opt in CircularEquilibrium.nonorthogonal_options_factory.defaults]
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

    eq = CircularEquilibrium(settings=options, nonorthogonal_settings=options)

    if add_noise is not None:
        # Add machine-precision level noise for testing robustness of grid generation
        import numpy as np
        from ..core.equilibrium import Point2D

        rng = np.random.default_rng(add_noise)
        for region in eq.regions.values():
            region.points = [
                Point2D(
                    p.R + rng.normal(scale=1.0e-16), p.Z + rng.normal(scale=5.0e-17)
                )
                for p in region.points
            ]

    if args.plot_regions:
        try:
            import matplotlib.pyplot as plt

            eq.plotPotential(ncontours=40)
            for region in eq.regions.values():
                plt.plot(
                    [p.R for p in region.points],
                    [p.Z for p in region.points],
                    marker="o",
                )
            print("Close window to continue...")
            plt.show()
        except Exception as err:
            warnings.warn(str(err))

    # Create the mesh

    mesh = BoutMesh(eq, options)
    mesh.calculateRZ()

    if args.plot_mesh:
        try:
            import matplotlib.pyplot as plt

            ax = eq.plotPotential(ncontours=40)
            mesh.plotPoints(
                xlow=options.get("plot_xlow", True),
                ylow=options.get("plot_ylow", True),
                corners=options.get("plot_corners", True),
                ax=ax,
            )
            print("Close window to continue...")
            plt.show()
        except Exception as err:
            warnings.warn(str(err))

    mesh.geometry()

    mesh.writeGridfile(options.get("grid_file", "bout.grd.nc"))

    # Delete and garbage-collect hypnotoad objects here so that any ParallelMap
    # instances get deleted if they exists. ParallelMap.__del__() calls
    # terminate on the worker processes. Needs to happen before program exits
    # or program will hang waiting for parallel workers to finish (which they
    # never do because they are waiting for another job in
    # ParallelMap.worker_run()).
    del eq, mesh
    gc.collect()


if __name__ == "__main__":
    main()
