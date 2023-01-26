#!/usr/bin/env python3

from hypnotoad import tokamak, Mesh
from matplotlib import pyplot as plt
from netCDF4 import Dataset as NCFile
import numpy as np
import yaml


def get_arg_parser():
    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    settings_docs = tokamak.TokamakEquilibrium.user_options_factory.doc

    parser = ArgumentParser(
        description=f"""
        Script to get flux surfaces from the equilibrium in a geqdsk file, using
        hypnotoad.

        Flux surfaces are saved as a dimension ``(N,2)`` array where ``N`` is the number
        of points used to represent the flux surface. ``N`` will vary depending which
        region (core, SOL, PFR, etc.) the surface belongs to. All flux surfaces are
        saved into the `flux surfaces` group in the output file. The wall contour is
        also saved, as the variable `wall` in the root group.

        Options that can be set in the input file to control outputs of this script:

        * ``finecontour_Nfine`` - number of points on each flux surface will be 1, 2 or
          3 times this.
        * ``number_of_processors`` - this many processors will be used to parallelise
          some computations
        * ``nx_core`` - {settings_docs["nx_core"]}
        * ``nx_inter_sep`` - {settings_docs["nx_inter_sep"]}
        * ``nx_pf`` - {settings_docs["nx_pf"]}
        * ``nx_sol`` - {settings_docs["nx_sol"]}
        * ``nx_sol_inner`` - {settings_docs["nx_sol_inner"]}
        * ``nx_sol_outer`` - {settings_docs["nx_sol_outer"]}
        * ``psinorm_core`` - {settings_docs["psinorm_core"]}
        * ``psinorm_sol`` - {settings_docs["psinorm_sol"]}
        * ``psinorm_sol_inner`` - {settings_docs["psinorm_sol_inner"]}
        * ``psinorm_pf`` - {settings_docs["psinorm_pf"]}
        * ``psinorm_pf_lower`` - {settings_docs["psinorm_pf_lower"]}
        * ``psinorm_pf_upper`` - {settings_docs["psinorm_pf_upper"]}
        * ``psi_core`` - {settings_docs["psi_core"]}
        * ``psi_sol`` - {settings_docs["psi_sol"]}
        * ``psi_sol_inner`` - {settings_docs["psi_sol_inner"]}
        * ``psi_pf_lower`` - {settings_docs["psi_pf_lower"]}
        * ``psi_pf_upper`` - {settings_docs["psi_pf_upper"]}
        * ``psi_spacing_separatrix_multiplier`` -
          {settings_docs["psi_spacing_separatrix_multiplier"]}
        """,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "geqdsk", help="geqdsk file - contains magnetic equilibrium data"
    )
    parser.add_argument(
        "settings",
        default=None,
        help="YAML file with settings controlling generated flux surfaces",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="flux_surfaces.nc",
        help="Filename for flux surfaces output.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        default=False,
        help="Make a plot of the flux surfaces that are being written.",
    )
    parser.add_argument(
        "-s",
        "--save-plot",
        default=None,
        help="File name for a plot of the flux surfaces that are being written.",
    )

    return parser


def main():
    args = get_arg_parser().parse_args()
    geqdsk = args.geqdsk
    settings_filename = args.settings
    output_filename = args.output
    plot = args.plot
    plot_filename = args.save_plot

    # Load settings from YAML file
    with open(settings_filename, "r") as settingsfile:
        settings = yaml.safe_load(settingsfile)
    if settings is None:
        settings = {}

    # Always use nonorthogonal grids so that flux surfaces align with wall
    settings["orthogonal"] = False

    # Load equilibrium data from geqdsk and create hypnotoad.Equilibrium object
    with open(geqdsk, "rt") as fh:
        eq = tokamak.read_geqdsk(fh, settings=settings, nonorthogonal_settings=settings)

    mesh = Mesh(eq, settings)

    n_regions = len(mesh.regions)

    def get_connected_region_name(start_region_name, segment_ind, n_connected_regions):
        if n_regions == 6:
            # Single null
            if start_region_name == "inner_lower_divertor":
                if segment_ind == 0:
                    return "PFR"
                elif segment_ind == 1:
                    return "SOL"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for single-null."
                    )
            elif start_region_name == "upper_outer_divertor":
                if segment_ind == 0:
                    return "PFR"
                elif segment_ind == 1:
                    return "SOL"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for single-null."
                    )
            elif start_region_name == "core":
                if segment_ind == 0:
                    return "core"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for single-null."
                    )
            else:
                raise ValueError(
                    f"Unexpected start_region_name={start_region_name} for single-null."
                )
        elif n_regions == 12:
            # Connected double null
            if start_region_name == "inner_lower_divertor":
                if segment_ind == 0:
                    return "lower PFR"
                elif segment_ind == 1:
                    return "inner SOL"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for connected-double-null."
                    )
            elif start_region_name == "outer_upper_divertor":
                if segment_ind == 0:
                    return "upper PFR"
                elif segment_ind == 1:
                    return "outer SOL"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for connected-double-null."
                    )
            elif start_region_name == "inner_core" or start_region_name == "outer_core":
                if segment_ind == 0:
                    return "core"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for connected-double-null."
                    )
            else:
                raise ValueError(
                    f"Unexpected start_region_name={start_region_name} for "
                    "connected-double-null."
                )
        elif n_regions == 18:
            # Disconnected double null
            if start_region_name == "inner_lower_divertor":
                if segment_ind == 0:
                    return "lower PFR"
                elif segment_ind == 1:
                    if n_connected_regions == 2:
                        return "lower PFR intersep"
                    elif n_connected_regions == 4:
                        return "intersep"
                    else:
                        raise ValueError(
                            f"Unexpected n_connected_regions={n_connected_regions} in "
                            f"{start_region_name} intersep for connected-double-null."
                        )
                elif segment_ind == 2:
                    return "inner SOL"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for disconnected-double-null."
                    )
            elif start_region_name == "outer_upper_divertor":
                if segment_ind == 0:
                    return "upper PFR"
                elif segment_ind == 1:
                    if n_connected_regions == 2:
                        return "upper PFR intersep"
                    elif n_connected_regions == 4:
                        return "intersep"
                    else:
                        raise ValueError(
                            f"Unexpected n_connected_regions={n_connected_regions} in "
                            f"{start_region_name} intersep for connected-double-null."
                        )
                elif segment_ind == 2:
                    return "outer SOL"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for disconnected-double-null."
                    )
            elif start_region_name == "inner_core" or start_region_name == "outer_core":
                if segment_ind == 0:
                    return "core"
                else:
                    raise ValueError(
                        f"Unexpected segment_ind={segment_ind} in {start_region_name} "
                        "for disconnected-double-null."
                    )
            else:
                raise ValueError(
                    f"Unexpected start_region_name={start_region_name} for "
                    "disconnected-double-null."
                )
        else:
            raise ValueError(
                f"Unexpected number of regions ({n_regions}), cannot indentify "
                "topology."
            )

    output_file = NCFile(output_filename, "x")
    # Create 'dimension' whose first entry represents the R-coordinate of a point, and
    # second dimension the Z-coordinate
    output_file.createDimension("coordinate", 2)
    flux_surfaces_group = output_file.createGroup("flux surfaces")

    if plot or plot_filename is not None:
        fig, ax = plt.subplots()
        plot_properties = iter(plt.rcParams["axes.prop_cycle"])

    for (region_name, _), region_index in mesh.region_lookup.items():
        region = mesh.regions[region_index]
        if region.yGroupIndex != 0:
            # Will follow flux surfaces from region to region, so only want to start
            # from the first MeshRegion (poloidally) in a connected set
            continue

        # Get the MeshRegions that are poloidally connected to this one
        connected_regions = [region]
        next_region_ind = region.connections["upper"]
        while next_region_ind is not None and next_region_ind != region_index:
            next_region = mesh.regions[next_region_ind]
            connected_regions.append(next_region)
            next_region_ind = next_region.connections["upper"]

        # numpy array for R,Z points for each flux surface in this set
        Nfine = mesh.user_options.finecontour_Nfine
        n_connected_regions = len(connected_regions)
        ny = Nfine
        ny_total = n_connected_regions * ny - (n_connected_regions - 1)

        nx = len(region.contours)
        if region.connections["outer"] is not None:
            # Include separatrices in the region 'outside' them (i.e. the
            # inter-separatrix and SOL regions)
            nx = nx - 1

        positions = np.zeros([nx, ny_total, 2])
        for i_region, r in enumerate(connected_regions):
            ystart = max(i_region * (ny - 1), 0)
            yend = (i_region + 1) * (ny - 1) + 1
            for ix in range(nx):
                contour = r.contours[ix]
                fc = contour.get_fine_contour(psi=eq.psi)
                positions[ix, ystart:yend, :] = fc.positions[
                    fc.startInd : fc.endInd + 1, :
                ]

        name = get_connected_region_name(
            region_name, region.radialIndex, len(connected_regions)
        )
        this_x = f"{name} radial"
        this_y = f"{name} poloidal"
        output_file.createDimension(this_x, nx)
        output_file.createDimension(this_y, ny_total)
        flux_surfaces_group.createVariable(name, float, (this_x, this_y, "coordinate"))
        flux_surfaces_group[name][...] = positions

        if plot or plot_filename is not None:
            props = next(plot_properties)
            lines = ax.plot(positions[:, :, 0].T, positions[:, :, 1].T, **props)
            # Only label the first line so the legend only has one entry per region
            lines[0].set_label(name)

    # Add wall contour to output file
    n_wall = len(eq.wall)
    output_file.createDimension("wall index", n_wall)
    output_file.createVariable("wall", float, ("wall index", "coordinate"))
    output_file["wall"][...] = [[p.R, p.Z] for p in eq.wall]

    output_file.close()

    if plot or plot_filename is not None:
        ax.set_xlabel("R")
        ax.set_ylabel("Z")
        ax.legend()
    if plot_filename is not None:
        fig.savefig(plot_filename)
    if plot:
        plt.show()


if __name__ == "__main__":
    main()
