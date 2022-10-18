> **Warning:** This example does not use/generate a geqdsk file, so does not
support the 'standard' workflow with `hypnotoad-gui` or `hypnotoad_geqdsk`.

To generate a tokamak grid for a simple poloidal flux, psi, given by an
analytic function:

    $ ./tokamak_example.py <type>

Where `<type>` is one of:
* `lsn`   Lower Single Null
* `usn`   Upper Single Null
* `cdn`   Connected Double Null
* `udn`   Upper Double Null
* `ldn`   Lower Double Null
* `udn2`  Upper Double Null, with a larger gap between separatrices.

The input files `single-null.yaml`, `connected-double-null.yaml` and
`disconnected-double-null.yaml` can be modified to change the grids for the
corresponding cases.

`tokamak_example.py` takes the additional arguments:
* `--np`, `--number-of-processors` sets the number of processors to use
    (default is to run in serial)
* `--nx` number of points (in major radius direction) to use when discretising
    psi on a Cartesian grid. Higher numbers give a more accurate psi. Does not
    affect the number of points in the generated grid.
* `--ny` number of points (in vertical direction) to use when discretising psi
    on a Cartesian grid. Higher numbers give a more accurate psi. Does not
    affect the number of points in the generated grid.
* `--no-plot` Do not produce plots. Means the script produces no output - only
    useful for checking if it can run without error in continuous integration
    testing.
