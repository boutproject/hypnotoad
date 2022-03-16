#!/usr/bin/env python3
#
# Inputs for hypnotoad are saved into grid files. This scripts reconstructs the original
# input files from a grid file

from argparse import ArgumentParser as AP
from netCDF4 import Dataset
from pathlib import Path


def main():
    parser = AP("Recreate input files from a hypnotoad grid file")
    parser.add_argument("grid_file", type=str)
    parser.add_argument(
        "-g",
        "--gfile-out",
        type=str,
        default=None,
        help="Name for recreated gfile. Defaults to name of original gfile.",
    )
    parser.add_argument("-y", "--yaml-out", type=str, default="hypnotoad.yaml")
    args = parser.parse_args()

    with Dataset(args.grid_file, "r") as gridfile:
        gfile_name = args.gfile_out
        if gfile_name is None:
            try:
                orig_path = Path(gridfile.getncattr("hypnotoad_geqdsk_filename"))
            except AttributeError:
                gfile_name = "hypnotoad.eqdsk"
            else:
                gfile_name = orig_path.name
        with open(gfile_name, "x") as gfile, open(args.yaml_out, "x") as yamlfile:
            gfile.write(gridfile["hypnotoad_input_geqdsk_file_contents"][...])
            yamlfile.write(gridfile["hypnotoad_inputs_yaml"][...])


if __name__ == "__main__":
    main()
