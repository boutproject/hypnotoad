#!/usr/bin/env python3


def get_arg_parser():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="""
        Plot the equilibrium stored in a geqdsk file
        """
    )
    parser.add_argument(
        "equilibrium_file", help="Path to equilibrium file in geqdsk format"
    )
    parser.add_argument(
        "--wall",
        action="store_false",
        default=True,
        help="Plot the wall contour if present",
    )
    parser.add_argument(
        "--wall-width",
        default=2.0,
        type=float,
        help="Width of the wall contour and separatrix contour",
    )
    parser.add_argument(
        "--separatrix",
        action="store_false",
        default=True,
        help="Plot the separatrix.  Note that for a disconnected double-null "
        "equilibrium, the thing plotted by default is not really the separatrix, as "
        "around the core the contour is a weighted average of the two separatrices "
        "constructed so that it joins the two X-points, while in the divertor legs it "
        "is the separatrix connected to the nearest X-point.",
    )
    parser.add_argument(
        "--separate-separatrices",
        action="store_true",
        default=False,
        help="Plot the two separatrices separately for a disconnected double-null "
        "equilibrium.",
    )
    parser.add_argument(
        "--separatrix-color",
        default="blue",
        nargs="*",
        help="Color (or list of colors) for the separatrix contour(s)",
    )
    parser.add_argument(
        "--psi-labels", action="store_true", default=False, help="Label psi contours"
    )
    parser.add_argument(
        "--n-points",
        default=100,
        type=int,
        help="Number points to use for interpolated array of psi values",
    )
    parser.add_argument(
        "--n-contours", default=40, type=int, help="Number of psi contours to plot"
    )
    parser.add_argument(
        "--color-contours",
        action="store_true",
        default=False,
        help="Color psi contours",
    )
    parser.add_argument(
        "--highlight-region",
        nargs=2,
        type=float,
        default=None,
        help="Highlight a region between the two values of normalised psi given by the "
        "arguments to this flag.",
    )
    parser.add_argument(
        "--highlight-color",
        default="orange",
        help="Color to use for region highlighted by `--hilight-region`",
    )
    parser.add_argument(
        "--show",
        action="store_false",
        default=True,
        help="Show plot in interactive window",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Name for output file. Suffix determines file format",
    )

    return parser


def main():
    args = get_arg_parser().parse_args()

    from ..cases import tokamak
    from matplotlib import pyplot as plt

    with open(args.equilibrium_file, "rt") as fh:
        eq = tokamak.read_geqdsk(fh)

    # Work out sensible aspect ratio for figure
    figwidth = 4.0
    figheight = figwidth * (eq.Zmax - eq.Zmin) / (eq.Rmax - eq.Rmin)
    plt.figure(figsize=(figwidth, figheight))

    if args.color_contours:
        colors = None
    else:
        colors = "grey"
    eq.plotPotential(
        npoints=args.n_points,
        ncontours=args.n_contours,
        labels=args.psi_labels,
        colors=colors,
        linestyles="-",
    )

    if args.highlight_region is not None:
        eq.plotHighlightRegion(args.highlight_region, color=args.highlight_color)

    if args.wall:
        eq.plotWall(linewidth=args.wall_width)

    if args.separatrix:
        eq.plotSeparatrix(
            scatter=False,
            separate_contours=args.separate_separatrices,
            linewidth=args.wall_width,
            linewidths=args.wall_width,
            color=args.separatrix_color,
            colors=args.separatrix_color,
        )

    if args.output is not None:
        plt.savefig(args.output, bbox_inches="tight")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
