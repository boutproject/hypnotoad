#!/usr/bin/env python

import gc
import numpy as np


def create_tokamak(geometry="lsn", nx=65, ny=65):
    """
    Create an example, based on a simple analytic form for the poloidal flux.

    Inputs
    ------

    geometry  string    lsn, usn, cdn, udn, ldn, udn2
    nx        int       Number of points in major radius
    ny        int       Number of points in height

    Returns
    -------

    r1d[nx]       1D array of major radius [m]
    z1d[ny]       1D array of height [m]
    psi2d[nx,ny]  2D array of poloidal flux [Wb]
    """

    r1d = np.linspace(1.2, 1.8, nx)
    z1d = np.linspace(-0.5, 0.5, ny)
    r2d, z2d = np.meshgrid(r1d, z1d, indexing="ij")

    r0 = 1.5
    z0 = 0.3

    psi_functions = {
        "lsn": lambda R, Z: (
            np.exp(-((R - r0) ** 2 + (Z + z0 - 0.3) ** 2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z + z0 + 0.3) ** 2) / 0.3**2)
        ),
        "usn": lambda R, Z: (
            np.exp(-((R - r0) ** 2 + (Z - z0 - 0.3) ** 2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z - z0 + 0.3) ** 2) / 0.3**2)
        ),
        "cdn": lambda R, Z: (
            np.exp(-((R - r0) ** 2 + Z**2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z + 2 * z0) ** 2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z - 2 * z0) ** 2) / 0.3**2)
        ),
        "udn": lambda R, Z: (
            np.exp(-((R - r0) ** 2 + Z**2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z + 2 * z0 + 0.002) ** 2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z - 2 * z0) ** 2) / 0.3**2)
        ),
        "ldn": lambda R, Z: (
            -np.exp(-((R - r0) ** 2 + Z**2) / 0.3**2)
            - np.exp(-((R - r0) ** 2 + (Z + 2 * z0) ** 2) / 0.3**2)
            - np.exp(-((R - r0) ** 2 + (Z - 2 * z0 - 0.003) ** 2) / 0.3**2)
        ),
        # Double null, but with the secondary far from the plasma edge
        "udn2": lambda R, Z: (
            np.exp(-((R - r0) ** 2 + Z**2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z + 2 * z0 + 0.02) ** 2) / 0.3**2)
            + np.exp(-((R - r0) ** 2 + (Z - 2 * z0) ** 2) / 0.3**2)
        ),
    }

    if geometry not in psi_functions:
        raise ValueError(
            "geometry not recognised. Choices are {}".format(psi_functions.keys())
        )

    psi_func = psi_functions[geometry]

    return r1d, z1d, psi_func(r2d, z2d), psi_func(np.linspace(r0, 1.2 * r0, nx), 0.0)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "geometry",
        type=str,
        nargs="?",
        default="lsn",
        choices=["lsn", "usn", "cdn", "ldn", "udn", "udn2"],
    )
    parser.add_argument("--nx", type=int, default=65)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--np", "--number-of-processors", type=int, default=-1)
    parser.add_argument("--no-plot", action="store_true", default=False)
    args = parser.parse_args()

    if "sn" in args.geometry:
        filename = "single-null.yaml"
    elif "cdn" == args.geometry:
        filename = "connected-double-null.yaml"
    else:
        filename = "disconnected-double-null.yaml"

    # Read input options

    import yaml

    with open(filename, "r") as inputfile:
        options = yaml.safe_load(inputfile)

    if args.np >= 0:
        options.update(number_of_processors=args.np)

    # Generate an artificial poloidal flux function
    r1d, z1d, psi2d, psi1d = create_tokamak(
        geometry=args.geometry,
        nx=args.nx,
        ny=args.ny,
    )

    from hypnotoad import tokamak

    eq = tokamak.TokamakEquilibrium(
        r1d, z1d, psi2d, psi1d, [], settings=options  # psi1d, fpol
    )

    from hypnotoad.core.mesh import BoutMesh

    mesh = BoutMesh(eq, options)
    mesh.geometry()

    if not args.no_plot:
        import matplotlib.pyplot as plt

        eq.plotPotential(ncontours=40)

        plt.plot(*eq.x_points[0], "rx")

        mesh.plotPoints(xlow=True, ylow=True, corners=True)

        plt.show()

    del eq, mesh
    gc.collect()
