#!/usr/bin/env python3

"""
Compare two BOUT++ grid files

Grid files must have same nx and ny, although different y_boundary_guards values are
handled.

Absolute and relative differences are calculated on the grid, even though the points may
be at slightly different positions. An option to interpolate quantities from one grid
onto the grid points of the other might be useful but is not implemented yet. When
differences are plotted on the poloidal plane, they are (for simplicity) plotted at the
grid point positions of the first grid file.

Additional requirements: xbout
Can be installed with
```
$ pip install --user xbout
```
or
```
$ conda install xbout
```
"""

from matplotlib import pyplot as plt
import xarray as xr
from xbout import open_boutdataset
from xbout.region import _create_regions_toroidal
from xbout.utils import _set_attrs_on_all_vars


def check_missing_variables(ds1, ds2, *, ignore_ylow=False):

    variables = list(set(v for v in ds1).union(v for v in ds2))

    variables.sort()

    common_variables = []
    only_on_1 = []
    only_on_2 = []
    for v in variables:
        if v not in ds1:
            if not (ignore_ylow and "ylow" in v):
                only_on_2.append(v)
        elif v not in ds2:
            if not (ignore_ylow and "ylow" in v):
                only_on_1.append(v)
        else:
            common_variables.append(v)

    scalar_variables = set(v for v in ds1.metadata).union(v for v in ds2.metadata)

    common_scalar_variables = []
    for v in scalar_variables:
        if v not in ds1.metadata:
            if not (ignore_ylow and "ylow" in v):
                only_on_2.append(v)
        elif v not in ds2.metadata:
            if not (ignore_ylow and "ylow" in v):
                only_on_1.append(v)
        else:
            common_scalar_variables.append(v)

    if only_on_1 != []:
        print(f"{only_on_1} are present in {ds1.name} but not {ds2.name}")
    if only_on_2 != []:
        print(f"{only_on_2} are present in {ds2.name} but not {ds1.name}")

    return common_variables, common_scalar_variables


def check_scalars(ds1, ds2, variables):
    """
    Requires all variables to be present in both ds1 and ds2
    """

    for v in variables:
        if not ds1.metadata[v] == ds2.metadata[v]:
            print(
                f"{v} has different values: {ds1.metadata[v]} in {ds1.name} and "
                f"{ds2.metadata[v]} in {ds2.name}"
            )


def plot_grid_points(ds1, ds2, *, poloidal_plot, show_all):
    if poloidal_plot:
        # Rxy and Zxy get renamed when applying toroidal geometry to Datasets
        Rxy = "R"
        Zxy = "Z"
    else:
        Rxy = "Rxy"
        Zxy = "Zxy"

    plt.figure()
    plt.scatter(ds1[Rxy], ds1[Zxy], marker="x", label=ds1.name)
    plt.scatter(ds2[Rxy], ds2[Zxy], marker="+", label=ds2.name)
    l = plt.legend()
    l.set_draggable(True)

    if not show_all:
        plt.show()


def trim_yboundaries(ds):
    myg = ds.metadata["y_boundary_guards"]

    if myg == 0:
        return ds

    ycoord = ds1.metadata.get("bout_ydim", "y")

    ds = ds.isel({ycoord: slice(myg, -myg)})

    if ds.metadata["jyseps2_1"] != ds.metadata["jyseps1_2"]:
        # Has a second (upper) divertor that needs to be removed
        ny_inner = ds.metadata["ny_inner"]

        lower_part = ds.isel({ycoord: slice(ny_inner)})
        upper_part = ds.isel({ycoord: slice(ny_inner + 2 * myg, None)})

        ds = xr.concat([lower_part, upper_part], dim=ycoord)

    new_metadata = ds.metadata.copy()
    new_metadata["keep_yboundaries"] = 0
    ds = _set_attrs_on_all_vars(ds, "metadata", new_metadata)

    if hasattr(ds, "geometry") and ds.geometry == "toroidal":
        ds = _create_regions_toroidal(ds)

    return ds


def plot_arrays(ds1, ds2, variables, *, atol, poloidal_plot, show_all):

    if poloidal_plot:
        ds1 = ds1.drop(("x", "theta_coord"))
        ds2 = ds2.drop(("x", "theta_coord"))
    else:
        ds1 = ds1.drop(("x", "y"))
        ds2 = ds2.drop(("x", "y"))

    if not ds1.metadata["y_boundary_guards"] == ds2.metadata["y_boundary_guards"]:
        ds1_trimmed = trim_yboundaries(ds1)
        ds2_trimmed = trim_yboundaries(ds2)
    else:
        # No need to trim y-boundaries
        ds1_trimmed = ds1
        ds2_trimmed = ds2

    for v in variables:
        if poloidal_plot:
            fig, axes = plt.subplots(1, 4)
        else:
            fig, axes = plt.subplots(2, 2)
        axes = axes.flatten()
        da1 = ds1[v]
        da2 = ds2[v]
        da1_trimmed = ds1_trimmed[v]
        da2_trimmed = ds2_trimmed[v]

        # Special handling because IDL hypnotoad and Python hypnotoad define 'y'
        # differently (dy is different by a constant factor)
        if v == "hthe":
            da1 = da1 * ds1["dy"]
            da2 = da2 * ds2["dy"]
            da1_trimmed = da1_trimmed * ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed * ds2_trimmed["dy"]
            v = "hthe*dy"
        elif v == "J":
            da1 = da1 * ds1["dy"]
            da2 = da2 * ds2["dy"]
            da1_trimmed = da1_trimmed * ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed * ds2_trimmed["dy"]
            v = "J*dy"
        elif v == "bxcvy":
            da1 = da1 / ds1["dy"]
            da2 = da2 / ds2["dy"]
            da1_trimmed = da1_trimmed / ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed / ds2_trimmed["dy"]
            v = "bxcvy/dy"
        elif v == "g12":
            da1 = da1 / ds1["dy"]
            da2 = da2 / ds2["dy"]
            da1_trimmed = da1_trimmed / ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed / ds2_trimmed["dy"]
            v = "g12/dy"
        elif v == "g22":
            da1 = da1 / ds1["dy"] ** 2
            da2 = da2 / ds2["dy"] ** 2
            da1_trimmed = da1_trimmed / ds1_trimmed["dy"] ** 2
            da2_trimmed = da2_trimmed / ds2_trimmed["dy"] ** 2
            v = "g22/dy**2"
        elif v == "g23":
            da1 = da1 / ds1["dy"]
            da2 = da2 / ds2["dy"]
            da1_trimmed = da1_trimmed / ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed / ds2_trimmed["dy"]
            v = "g23/dy"
        elif v == "g_12":
            da1 = da1 * ds1["dy"]
            da2 = da2 * ds2["dy"]
            da1_trimmed = da1_trimmed * ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed * ds2_trimmed["dy"]
            v = "g_12*dy"
        elif v == "g_22":
            da1 = da1 * ds1["dy"] ** 2
            da2 = da2 * ds2["dy"] ** 2
            da1_trimmed = da1_trimmed * ds1_trimmed["dy"] ** 2
            da2_trimmed = da2_trimmed * ds2_trimmed["dy"] ** 2
            v = "g_22*dy**2"
        elif v == "g_23":
            da1 = da1 * ds1["dy"]
            da2 = da2 * ds2["dy"]
            da1_trimmed = da1_trimmed * ds1_trimmed["dy"]
            da2_trimmed = da2_trimmed * ds2_trimmed["dy"]
            v = "g_23*dy"

        vmin = min(da1.min().values, da2.min().values)
        vmax = max(da1.max().values, da2.max().values)

        absolute_difference = abs(da1_trimmed - da2_trimmed)
        abs_da1 = abs(da1_trimmed)
        relative_difference = absolute_difference / abs_da1.where(abs_da1 > atol, atol)

        fig.suptitle(v)
        if len(da1.shape) > 1:
            if poloidal_plot:
                absolute_difference.attrs["metadata"] = da1_trimmed.metadata
                absolute_difference.attrs["regions"] = da1_trimmed.regions
                absolute_difference = absolute_difference.assign_coords(
                    R=da1_trimmed["R"]
                )
                absolute_difference = absolute_difference.assign_coords(
                    Z=da1_trimmed["Z"]
                )
                relative_difference.attrs["metadata"] = da1_trimmed.metadata
                relative_difference.attrs["regions"] = da1_trimmed.regions
                relative_difference = relative_difference.assign_coords(
                    R=da1_trimmed["R"]
                )
                relative_difference = relative_difference.assign_coords(
                    Z=da1_trimmed["Z"]
                )

                da1.bout.pcolormesh(ax=axes[0], vmin=vmin, vmax=vmax)
                da2.bout.pcolormesh(ax=axes[1], vmin=vmin, vmax=vmax)
                absolute_difference.bout.pcolormesh(ax=axes[2])
                relative_difference.bout.pcolormesh(ax=axes[3])
            else:
                da1.plot(ax=axes[0], vmin=vmin, vmax=vmax)
                da2.plot(ax=axes[1], vmin=vmin, vmax=vmax)
                absolute_difference.plot(ax=axes[2])
                relative_difference.plot(ax=axes[3])
            axes[0].set_title(ds1.name)
            axes[1].set_title(ds2.name)
        else:
            da1.plot(ax=axes[0], label=ds1.name)
            da2.plot(ax=axes[0], label=ds2.name)
            l = axes[0].legend()
            l.set_draggable(True)
            absolute_difference.plot(ax=axes[2])
            relative_difference.plot(ax=axes[3])
        axes[2].set_title("absolute difference")
        axes[3].set_title("relative difference")

        if not show_all:
            plt.show()


if __name__ == "__main__":
    import argparse
    from sys import exit

    parser = argparse.ArgumentParser()
    parser.add_argument("ds1")
    parser.add_argument("ds2")
    parser.add_argument("--ignore-ylow", action="store_true", default=False)
    parser.add_argument("--atol", type=float, default=1.0e-8)
    parser.add_argument("--poloidal-plot", action="store_true", default=False)
    parser.add_argument("--show-all", action="store_true", default=False)
    args = parser.parse_args()

    ds1 = open_boutdataset(
        args.ds1, keep_xboundaries=True, keep_yboundaries=True, info=False
    ).load()
    ds1.attrs["name"] = args.ds1

    ds2 = open_boutdataset(
        args.ds2, keep_xboundaries=True, keep_yboundaries=True, info=False
    ).load()
    ds2.attrs["name"] = args.ds2

    if args.poloidal_plot:
        from xbout.geometries import apply_geometry

        ds1.metadata["MXG"] = 2
        ds1.metadata["MYG"] = ds1.metadata["y_boundary_guards"]
        ds2.metadata["MXG"] = 2
        ds2.metadata["MYG"] = ds2.metadata["y_boundary_guards"]

        coordinates = {"x": "psi_poloidal", "y": "theta_coord"}
        ds1 = apply_geometry(ds1, "toroidal", coordinates=coordinates)
        ds2 = apply_geometry(ds2, "toroidal", coordinates=coordinates)

    common_variables, common_scalar_variables = check_missing_variables(
        ds1, ds2, ignore_ylow=args.ignore_ylow
    )
    check_scalars(ds1, ds2, common_scalar_variables)
    plot_grid_points(ds1, ds2, poloidal_plot=args.poloidal_plot, show_all=args.show_all)
    plot_arrays(
        ds1,
        ds2,
        common_variables,
        atol=args.atol,
        poloidal_plot=args.poloidal_plot,
        show_all=args.show_all,
    )

    if args.show_all:
        plt.show()

    exit(0)
