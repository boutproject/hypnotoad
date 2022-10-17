#!/usr/bin/env python3

import os
from pathlib import Path
import sys

# Put the integrated_tests directory into sys.path so we can import from it
sys.path.append(str(Path(__file__).joinpath("..", "..", "..").resolve()))
from integrated_tests.utils import run_case

diagnose = False

rtol = 2.0e-9
atol = 5.0e-10

# make sure we are in the test directory
os.chdir(Path(__file__).parent)

while len(sys.argv) < 3:
    sys.argv.append(None)
sys.argv[1] = "../grid_files/test_connected-double-null.eqdsk"

run_case(
    "orthogonal",
    "../connected_doublenull_orthogonal/test_orthogonal.yml",
    "../connected_doublenull_orthogonal/expected_orthogonal.grd.nc",
    rtol=rtol,
    atol=atol,
    diagnose=diagnose,
    add_noise=246,
)

sys.exit(0)
