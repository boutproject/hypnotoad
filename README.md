Installation
------------

#### From conda

``hypnotoad`` is available from the ``conda-forge`` channel. Install with

    $ conda install -c conda-forge hypnotoad

To use the GUI, one of PySide2 or PyQt5 is needed (PySide and PyQt4 may also
work, but have not been tested). Install with

    $ conda install -c conda-forge pyside2

or

    $ conda install -c conda-forge pyqt

#### From PyPi

``hypnotoad`` can be installed using ``pip`` by running

    $ pip install --user hypnotoad

To use the GUI, one of PySide2 or PyQt5 is needed (PySide and PyQt4 may also
work, but have not been tested). These can be fetched by choosing a variant.
For PySide2 use

    $ pip install --user hypnotoad[gui-pyside2]

or for PyQt5 use

    $ pip install --user hypnotoad[gui-PyQt5]

#### git repo

If you need to modify the hypnotoad code, or get development versions, clone
from github

    $ git clone git@github.com:boutproject/hypnotoad.git

You can install from the git repo with ``pip``, this is useful to get the
executables added to your path. If you use ``conda`` you may wish to first
install the dependencies using

    $ conda install boututils matplotlib netcdf4 numpy optionsfactory pyparsing pyqt pyyaml qt.py scipy

(replacing ``pyqt`` with ``pyside2`` if you prefer PySide2 to PyQt5) to ensure
they are not ``pip``-installed. Make sure to do an 'editable' install using
``-e`` or ``--editable`` option like

    $ cd hypnotoad
    $ pip install --user -e .

(you may also need to ``pip``-install ``PySide2`` or ``PyQt5`` to use the GUI).
If installing in a conda environment you do not need the ``--user`` argument.
This installs executables which use the code that's currently in the git repo,
so if you edit or update it you will see the updates. If you install with ``pip
install .`` (without the ``-e``) then ``pip`` can get confused because it can't
tell which version number is newer, as the git repo versions have a version
number based on the git hash, not a simple x.y.z; then pip may for example not
uninstall hypnotoad correctly.


Usage
-----

Options are read and set up in the Equilibrium (child-)class object, and passed
from there to the Mesh (child-)class object.

User-settable options, with their current values, are printed when an
Equilibrium object is created.  Internal options should not need to be set by
the user, but can be overridden with keyword arguments to the Equilibrium
constructor.

Hypnotoad can be run either as an executable (``hypnotoad_geqdsk``), which just
reads from an input file, using the gui (``hypnotoad-gui``) or interactively
from a Python shell. To ensure reproducibility, it is suggested to create your
final grid non-interactively. The interactive mode is intended to make it
easier to prototype the grid and find a good set of input parameters. Once you
have found a configuration you are happy with, you can save the current input
parameters using the save dialog in the gui, or with
``Equilibrium.saveOptions(filename='hypnotoad_options.yaml')`` from the Python
shell; this may be especially useful if you have changed some options from the
Python shell with keyword-arguments.

Grid generation can take a while with the default options, which are set for
high accuracy. When prototyping, it is suggested to temporarily use lower
accuracy. The following may be a good starting point:
- finecontour\_Nfine=100. This speeds up the creation of the internal,
  high-resolution, fixed-spacing representation of contours, and also
  calculations of distance along contours and some interpolation functions.
- gradPsiRtol=2.e-6 and gradPsiAtol=1.e-6. These control the maximum error on
  the integration along grad(psi) used to trace grid lines orthogonal to the
  flux surfaces. They do not usually make a huge difference, but affect the
  time spent in 'Following perpendicular'.
- If your wall is given by a large number of points (say more than 20) it might
  be worth creating a simpler one with fewer points for prototyping. This will
  speed up the 'finding wall intersections' stage. Note that the wall only
  matters where it intersects the grid.
- Decreasing the resolution of the grid will also help. The grid points will
  probably not be in exactly the same place, but the algorithms are intended to
  produce grid spacings that are inversely proportional to the total number of
  points, so the structure should be very similar.


Developing
----------

Developer documentation is [here](doc/developer/developer.md).
