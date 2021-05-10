#!/usr/bin/env python3

import os
from pathlib import Path
import sys

diagnose = False

rtol = 1.0e-9
atol = 5.0e-10

# make sure we are in the test directory
os.chdir(Path(__file__).parent)

sys.path.insert(0, "..")
from utils import run_case  # noqa: E402

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
