#!/usr/bin/env python3


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description="""
        Script to plot grid cells for a BOUT++ grid, using the corner positions saved by
        hypnotoad.

        The CELL_CENTRE positions are shown in the plot as black dots, the cells are
        shown by joining the corners with black lines. The branch cuts (if shown) are
        thick red lines.
        """
    )
    parser.add_argument(
        "gridfile", help="Path to the grid file to plot grid cells from"
    )
    parser.add_argument(
        "--mxg",
        type=int,
        nargs="?",
        default=2,
        help="Number of boundary cells at x-boundaries (default 2)",
    )
    parser.add_argument(
        "--branch-cuts",
        action="store_true",
        default=False,
        help="Highlight the branch cuts?",
    )
    parser.add_argument(
        "--separatrix",
        action="store_true",
        default=False,
        help="Highlight the separatrix/separatrices?",
    )
    parser.add_argument(
        "--save-as",
        default=None,
        help="Filename to save plot to, suffix gives type of file",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        default=False,
        help="Skip showing the plot in a window?",
    )
    args = parser.parse_args()
    gridfile = args.gridfile
    mxg = args.mxg
    branch_cuts = args.branch_cuts
    separatrix = args.separatrix
    save_as = args.save_as
    no_show = args.no_show
    if mxg < 1:
        raise ValueError(f"mxg must be at least 1, got {mxg}")

    try:
        from xbout import open_boutdataset
    except ImportError:
        raise ImportError("xbout is required for this plotting script")

    from matplotlib import pyplot as plt

    ds = open_boutdataset(
        gridfile,
        geometry="toroidal",
        keep_xboundaries=True,
        keep_yboundaries=True,
        info=False,
        drop_variables="theta",
    )

    y_boundary_guards = ds.metadata.get("y_boundary_guards", 0)
    if y_boundary_guards < 1:
        raise ValueError(
            "Grid file does not include y-boundary cells. These are required for grid "
            "plotting"
        )

    plt.figure(constrained_layout=True)

    for r in ds.regions:
        ds_region = ds.bout.from_region(r, with_guards={"x": 1, "theta": 1})

        if ds_region.regions[r].connection_inner_x is None:
            xin = mxg
            xin_pol = mxg
        else:
            xin = None
            xin_pol = 1

        if ds_region.regions[r].connection_outer_x is None:
            if mxg == 1:
                xout_corner_rad = None
                xout_corner_pol = xout_corner_rad
            else:
                xout_corner_rad = -mxg + 1
                xout_corner_pol = xout_corner_rad
            xout_centre = -mxg
        else:
            xout_corner_rad = None
            xout_corner_pol = -1
            xout_centre = None

        if ds_region.regions[r].connection_lower_y is None:
            ylow = y_boundary_guards
        else:
            ylow = None

        if ds_region.regions[r].connection_upper_y is None:
            if y_boundary_guards == 1:
                yup_corner = None
            else:
                yup_corner = -y_boundary_guards + 1
            yup_centre = -y_boundary_guards
        else:
            yup_corner = None
            yup_centre = None

        # plot grid points at cell centres
        # Note: have to use `s` to set marker size, because xarray has used a `size`
        # kwarg for something else.
        ds_region.isel(
            x=slice(xin, xout_centre), theta=slice(ylow, yup_centre)
        ).plot.scatter("R", "Z", marker=".", color="k", s=1)

        # plot radial grid lines
        plt.plot(
            ds_region["Rxy_corners"].isel(
                x=slice(xin, xout_corner_rad), theta=slice(ylow, yup_corner)
            ),
            ds_region["Zxy_corners"].isel(
                x=slice(xin, xout_corner_rad), theta=slice(ylow, yup_corner)
            ),
            color="k",
        )

        # plot poloidal grid lines
        plt.plot(
            ds_region["Rxy_corners"].isel(
                x=slice(xin_pol, xout_corner_pol), theta=slice(ylow, yup_corner)
            ).T,
            ds_region["Zxy_corners"].isel(
                x=slice(xin_pol, xout_corner_pol), theta=slice(ylow, yup_corner)
            ).T,
            color="k",
        )

        if branch_cuts and ds_region.regions[r].connection_upper_y is not None:
            # Highlight the branch cuts
            # By arbitrary choice, plot from the region(s) below the branch cut, so
            # highlight the upper edge.
            plt.plot(
                ds_region["Rxy_corners"].isel(
                    x=slice(xin, xout_corner_rad), theta=-1
                ),
                ds_region["Zxy_corners"].isel(
                    x=slice(xin, xout_corner_rad), theta=-1
                ),
                color="r",
                linewidth=3,
                zorder=1000,
            )

        if separatrix and ds_region.regions[r].connection_outer_x is not None:
            # Highlight the branch cuts
            # By arbitrary choice, plot from the region(s) inside the separatrix, so
            # highlight the outer edge.
            plt.plot(
                ds_region["Rxy_corners"].isel(
                    x=-1, theta=slice(ylow, yup_corner)
                ),
                ds_region["Zxy_corners"].isel(
                    x=-1, theta=slice(ylow, yup_corner)
                ),
                color="b",
                linewidth=3,
                zorder=999,
            )

    plt.gca().set_aspect("equal")

    if save_as is not None:
        plt.savefig(save_as)

    if not no_show:
        plt.show()
