Quickstart
==========

This quickstart describes how to use ``hypnotoad`` to create a grid for a
single-null or double-null diverted tokamak configuration.

Installation
------------

``hypnotoad`` can be installed using ``pip``::

    $ pip install --user hypnotoad

or conda (from the conda-forge channel)::

    $ conda install hypnotoad

Usage
-----

The magnetic equilibrium should be provided in the form of a `geqdsk file
<https://fusion.gat.com/theory/Efitgeqdsk>`_ (note the link requires a GA
account to access), using `COCOS-1 conventions
<https://doi.org/10.1016/j.cpc.2012.09.010>`_.

The settings controlling grid generation are provided in a YAML input file. It
is recommended to prepare this input file using the GUI, as described below.
For reproducibility, produce the final grid file by running::

    $ hypnotoad-geqdsk <path/to/geqdsk> <path/to/input.yml>

to ensure that the interactive nature of the GUI does not affect the grid file
(it should not, but bugs are possible!).

To prepare input:

1. Start the hypnotoad GUI::

   $ hypnotoad-gui

2. Load a geqdsk file using the 'Browse' button beside the 'geqdsk file' box
   near the bottom of the window.

3. (Optional) load settings from an existing input file, using the 'Browse'
   button beside the 'Options file' box near the bottom of the window.

4. Adjust the settings in the table in the left side of the window. Help text
   for each setting is shown in a pop-up if you hover the cursor over the
   setting.

5. Click 'Run' to create the set of grid points. Re-adjust the settings and
   repeat until you are satisfied.

6. Save the settings using the 'File' menu, or the buttons near the top of the
   window (hover over the buttons to see descriptions). Note that 'Save' in the
   'File' menu and the 'Save options to file' button will overwrite the
   existing input file (if one was loaded). Use 'Save as' or the 'Save options
   to file with new filename' button if you do not want to overwrite.

7. You can click 'Write Grid' to create a grid file. This step will perform
   calculations of the geometric quantities required by BOUT++ and may expose
   problems with the grid generation, even when the process was successful up
   to this step. Note that as mentioned above it is recommended to create the
   final grid on the command line using a settings file for reproducibility.

Unfortunately, many problems can occur during grid generation. For some hints
on how to solve or work around them, see :ref:`tips-and-tricks:Tips and tricks
for fixing problems`
