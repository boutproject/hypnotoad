#!/usr/bin/env python3

import os
from pathlib import Path
import sys
import xarray.testing as xrt

# Put the integrated_tests directory into sys.path so we can import from it
sys.path.append(str(Path(__file__).joinpath("..", "..", "..").resolve()))
from integrated_tests.utils import run_case

diagnose = False

rtol = 1.0e-9
atol = 5.0e-10

# make sure we are in the test directory
os.chdir(Path(__file__).parent)

while len(sys.argv) < 3:
    sys.argv.append(None)
sys.argv[1] = "../grid_files/test_connected-double-null.eqdsk"

grid_np1 = run_case(
    "nonorthogonal",
    "test_nonorthogonal.yml",
    "expected_nonorthogonal.grd.nc",
    rtol=rtol,
    atol=atol,
    diagnose=diagnose,
)

grid_np2 = run_case(
    "nonorthogonal",
    "test_nonorthogonal_np2.yml",
    "expected_nonorthogonal.grd.nc",
    rtol=rtol,
    atol=atol,
    diagnose=diagnose,
)

xrt.assert_identical(grid_np1, grid_np2)

sys.exit(0)
