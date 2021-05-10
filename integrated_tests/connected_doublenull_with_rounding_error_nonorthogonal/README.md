Integrated test with noise added, for nonorthogonal grid case.

Slightly modified version of the connected\_doublenull\_nonorthogonal test.
Purpose is to check robustness of grid generation. When expected output is
generated on the same system as the test is run, the results will match
regardless of tolerances, etc. This test introduces machine-precision level
errors to the initial state (after creating the TokamakEquilibrium but before
creating the BoutMesh) to check the sensitivity to rounding errors that might
be expected when running on a different system.

To run the test, execute `./runtest.py`.

Expected data requires git-lfs https://git-lfs.github.com/

Test requires xarray http://xarray.pydata.org/en/stable/

Test uses an artificial equilibrium, see ../grid\_files.
