#!/usr/bin/env python3

import argparse
from Qt.QtWidgets import QApplication
import sys

from ..__init__ import __version__
from .gui import HypnotoadGui


def main():
    """
    BOUT++ grid generation
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    parser.parse_args()

    app = QApplication(sys.argv)
    window = HypnotoadGui()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
