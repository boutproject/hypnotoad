#!/usr/bin/env python3

# Copyright 2020 J.T. Omotani
#
# Contact John Omotani john.omotani@ukaea.uk
#
# This file is part of Hypnotoad 2.
#
# Hypnotoad 2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Hypnotoad 2 is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# Hypnotoad 2.  If not, see <http://www.gnu.org/licenses/>.

"""
Create a BOUT++ grid for TORPEX from an input file giving coil currents and positions

Input file should contain coil parameters, for each coil:
    R: major radius in metres
    Z: major radius in metres
    I: clockwise current in Amps

Note: positions of cell corners are generated first, grid points are then put in the
centre of the cell.
"""


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--noplot", action="store_true")
    args = parser.parse_args()

    if not args.noplot:
        from matplotlib import pyplot

    from ..cases.torpex import createMesh

    filename = args.filename
    gridname = "torpex.grd.nc"

    mesh = createMesh(filename)

    try:
        mesh.geometry()
    except Exception as e:
        import traceback

        print("There was an exception in mesh.geometry:", str(e))
        print("****************************************")
        traceback.print_tb(e.__traceback__)
        print("****************************************")

    if not args.noplot:
        pyplot.figure()
        mesh.equilibrium.plotPotential()
        mesh.equilibrium.addWallToPlot()
        pyplot.plot(*mesh.equilibrium.x_points[0], "rx")
        mesh.plotPoints(xlow=True, ylow=True, corners=True)
        pyplot.show()

    mesh.writeGridfile(gridname)


if __name__ == "__main__":
    main()
