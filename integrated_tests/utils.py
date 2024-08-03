import numpy as np
import os
from sys import argv
from xarray import open_dataset
import xarray.testing as xrt

from hypnotoad.scripts.hypnotoad_geqdsk import main as hyp_geqdsk

expected_different_attrs = [
    "grid_id",
    "hypnotoad_version",
    "hypnotoad_git_hash",
    "hypnotoad_git_diff",
]
expected_different_vars = [
    "hypnotoad_inputs",
    "hypnotoad_inputs_yaml",
    "Python_version",
    "module_versions",
    # Variables that have been added. These entries can be removed if/when the
    # expected output is re-generated.
    "penalty_mask",
    "closed_wall_R",
    "closed_wall_Z",
    "Jpar0",
    "Jpar0_xlow",
    "Jpar0_ylow",
]


def check_errors(ds1, ds2, *, rtol, atol):
    print("\nError check 1")
    print(
        "{0:>20}{1:>23}{2:>23}{3:>5}{4:>5}".format(
            "variable",
            "significant error",
            "value",
            "xind",
            "yind",
        )
    )
    for v in ds1:
        if len(ds1[v].dims) < 2:
            continue
        try:
            diff = ds1[v] - ds2[v]
            tol = atol + ds2[v] * rtol
            test = diff / tol
            try:
                loc = {k: v.values for k, v in abs(test).argmax(...).items()}
            except ValueError:
                pos = {"x": "error", "y": "error"}
            print(
                "{0:>20}{1:>23}{2:>23}{3:>5}{4:>5}".format(
                    v,
                    diff.isel(loc).values,
                    ds2[v].isel(loc).values,
                    loc["x"],
                    loc["y"],
                ),
                flush=True,
            )
        except np.core._exceptions.UFuncTypeError:
            pass

    print("\nError check 2")
    print(
        "{0:>20}{1:>23}{2:>23}{3:>23}{4:>23}".format(
            "variable",
            "relative error",
            "absolute error",
            "xind",
            "yind",
        ),
        flush=True,
    )
    for v in ds1:
        if len(ds1[v].dims) < 2:
            continue
        try:
            diff = ds1[v] - ds2[v]
            rdiff = diff / ds2[v]
            try:
                pos = {k: v.values for k, v in abs(rdiff).argmax(...).items()}
            except ValueError:
                pos = {"x": "error", "y": "error"}
            print(
                "{0:>20}{1:>23}{2:>23}{3:>23}{4:>23}".format(
                    v,
                    abs(rdiff).max().values,
                    abs(diff).max().values,
                    pos["x"],
                    pos["y"],
                ),
                flush=True,
            )
        except np.core._exceptions.UFuncTypeError:
            pass


def run_case(name, inputfile, expectedfile, *, rtol, atol, diagnose, add_noise=None):
    argv[2] = inputfile

    if os.path.isfile("bout.grd.nc"):
        os.remove("bout.grd.nc")

    hyp_geqdsk(add_noise=add_noise)

    expected = (
        open_dataset(expectedfile).load().drop(expected_different_vars, errors="ignore")
    )
    actual = (
        open_dataset("bout.grd.nc")
        .load()
        .drop(expected_different_vars, errors="ignore")
    )

    if diagnose:
        check_errors(expected, actual, rtol=rtol, atol=atol)

    xrt.assert_allclose(expected, actual, rtol=rtol, atol=atol)

    for attrname in expected_different_attrs:
        del expected.attrs[attrname]
        del actual.attrs[attrname]
    assert expected.attrs == actual.attrs

    variables = set(expected).union(set(actual))
    for v in variables:
        assert (
            expected[v].attrs == actual[v].attrs
        ), f"expect identical attributes for {v}"

    expected.close()
    actual.close()

    print(name, "case passed!", flush=True)

    return actual
