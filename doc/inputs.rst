Inputs
======

Magnetic equilibrium
--------------------

The magnetic equilibrium should be provided in the form of a `geqdsk file
<https://fusion.gat.com/theory/Efitgeqdsk>`_ (note the link requires a GA
account to access). The data in the file should use `COCOS-1 conventions
<https://doi.org/10.1016/j.cpc.2012.09.010>`_.

The geqdsk file provides the poloidal magnetic flux function :math:`\psi`
(which is the poloidal magnetic flux divided by :math:`2\pi`), and
:math:`f_\mathrm{pol}` that determines the toroidal magnetic field as
:math:`B_\mathrm{toroidal} = f_\mathrm{pol}(\psi)\nabla\zeta` where
:math:`\zeta` is the toroidal angle (:math:`f_\mathrm{pol}(psi)` is sometimes
written as :math:`I(\psi)`).

The geqdsk file may also contain a set of points that describe the wall of the
tokamak. Only a contour in the poloidal plane can be represented, so 3d
features of the wall cannot be included, consistent with hypnotoad generating a
toroidally-symmetric grid. The wall contour is only used by hypnotoad at the
divertor targets, so the detail of the wall shape elsewhere is not important.

The geqdsk file is passed as the first argument to the command line program
``hypnotoad-geqdsk``, or when using ``hypnotoad-gui`` is entered in the 'geqdsk
file' box near the bottom of the window.

The geqdsk file used to create a grid file is saved in the grid file for
provenance tracking. It can be conveniently extracted from the grid file with
the ``hypnotoad-recreate-inputs`` program.

Settings
--------

Grid generation is controlled by settings, saved in a YAML file. The settings
can be conveniently edited by using ``hypnotoad-gui`` which can load (using the
'Options file' box near the bottom of the window) and save (using the 'File'
dialog menu or the buttons near the top of the screen) the settings to/from
YAML files. ``hypnotoad-gui`` provides an interactive interface where default
values are displayed (and defaults that depend on other settings are updated),
with tool-tips showing help text appearing if you hover the cursor over each
setting.

When generating a grid on the command line, the settings file is passed as the
second argument to ``hypnotoad-geqdsk``.

Relevant options will be described in the following sections. A complete list
of the options for both orthogonal and nonorthogonal grids is here:
:ref:`_temp/options:Tokamak options`. Options that only apply to nonorthogonal
grids are listed here: :ref:`_temp/nonorthogonal-options:Nonorthogonal
options`. Some other options used to set tolerances, etc. for the grid
generation process are here: :ref:`_temp/mesh-options:Mesh options`.

All the settings used for grid generation are saved (with all defaults
evaluated) in YAML format in the grid file for provenance tracking. They can be
conveniently extracted as a YAML file from the grid file with the
``hypnotoad-recreate-inputs`` program.
