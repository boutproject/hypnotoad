What's new
==========

0.4.5 (unreleased)
------------------

### Bug fixes

- Fix `psi_axis` check when `reverse_current=True` (#131)\
  By [John Omotani](https://github.com/johnomotani)

0.4.4
-----

### Bug fixes

- Fix for circular geometries with non-constant q (#120)\
  By [John Omotani](https://github.com/johnomotani)

### New features

- Command to recreate input file and gfile from a grid file (#119)\
  By [John Omotani](https://github.com/johnomotani)
- Save all variables at `xlow` positions (as `*_xlow`), as these can be used by
  BOUT++ to create `CELL_XLOW` `Coordinates` objects. Also save the positions
  of the cell corners (`Rxy_corners` and `Zxy_corners`) as these may be useful
  for plotting (#118)\
  By [John Omotani](https://github.com/johnomotani)
- Add executable `hypnotoad_plot_grid_cells` that takes a grid file as input and uses
  `xbout` to plot the grid cells from the newly added corner positions (#118)\
  By [John Omotani](https://github.com/johnomotani)

0.4.3
-----

### Bug fixes

- Small fixes and minor optimizations for nonorthogonal gridding (#113)\
  By [John Omotani](https://github.com/johnomotani)
- Pass `psi` to keyword argument where this update had been missed previously
  (#113)\
  By [John Omotani](https://github.com/johnomotani)
- Fix plotting of poloidal lines in Mesh.plotPoints() (#113)\
  By [John Omotani](https://github.com/johnomotani)
- Minor fixes updating TORPEX case (#106)\
  By [John Omotani](https://github.com/johnomotani)
- Fix crash in hypnotoad-gui if input file did not contain an 'orthogonal'
  setting (#106)\
  By [John Omotani](https://github.com/johnomotani)

### New features

- Save poloidal distances along psi contours to grid file (#116)\
  By [John Omotani](https://github.com/johnomotani)
- Option to reverse sign of toroidal magnetic field Bt for TokamakEquilibrium (#114)\
  By [John Omotani](https://github.com/johnomotani)
- Enable DCT interpolation scheme in TokamakEquilibrum, can be selected with new option
  psi_interpolation_method (#113)\
  By [John Omotani](https://github.com/johnomotani)
- Add xpoint_offset option for TORPEX case (#113)\
  By [John Omotani](https://github.com/johnomotani)
- Clean up of examples. Input files from examples/tokamak can now be used with
  hypnotoad-gui (#106)\
  By [John Omotani](https://github.com/johnomotani)
- psi_axis and psi_bdry are added to the grid file for post-processing.
  The calculated values should be consistent with the EFIT data (#112)\
  By [Qian Xia](https://github.com/Vandoo)

0.4.1 (26th May 2021)
---------------------

### Bug fixes

- Check if an equilibrium can be gridded as a connected double-null before gridding.
  Prevents creation of invalid grids where second X-point is outside the first
  flux-surface in the SOL (#104)\
  By [John Omotani](https://github.com/johnomotani)
- Correct setting of ixseps1 and ixseps2 for upper disconnected double null
  case - should have ixseps1>ixseps2 but previously had ixseps1<ixseps2 (#109)\
  By [John Omotani](https://github.com/johnomotani)

0.4.0 (26th May 2021)
---------------------

### Breaking changes

- Method for calculating zShift made more accurate. This will change the calculated
  zShift for the same input file (#101)
- Default value for 'curvature_type' changed from "curl(b/B) with x-y derivatives" to
  "curl(b/B)". New default should be more accurate and more consistent when changing
  grid sizes, but will produce slightly different output from the same input file. If
  the old behaviour is needed, set 'curvature_type = "curl(b/B) with x-y derivatives"'
  explicitly (#100)

### New features

- Accuracy of calculation of zShift improved by integrating on FineContours
  rather than PsiContours (#101)\
  By [John Omotani](https://github.com/johnomotani)
- Circular, concentric flux surfaces magnetic geometry (#100)\
  By [John Omotani](https://github.com/johnomotani)
- Non-grid based method for calculating curvature (#100)\
  By [John Omotani](https://github.com/johnomotani)
- Parallelise most expensive loops (#99)\
  By [John Omotani](https://github.com/johnomotani)
- Fix position of start and end points of contours when refining and add more tolerance
  settings, to enable more reproducible grid generation (#95)\
  By [John Omotani](https://github.com/johnomotani)
- Options are saved as a YAML string in "hypnotoad_inputs_yaml" to make them
  easier to read in code later. (#98)\
  By [John Omotani](https://github.com/johnomotani)
- When exceptions are caught by the GUI, print the traceback as well as the
  exception message (#95)\
  By [John Omotani](https://github.com/johnomotani)

### Testing

- Integrated tests, based on an analytic, connected double-null equilibrium (#97, fixes
  #50)\
  By [John Omotani](https://github.com/johnomotani)

### Bug fixes

- Save all options to grid files. Previously only Equilibrium options were
  saved. Now also Mesh and nonorthogonal options (#98)\
  By [John Omotani](https://github.com/johnomotani)
- Setting to adjust extension of FineContours past targets, may help to avoid crashes on
  problematic equilibria (#96)\
  By [John Omotani](https://github.com/johnomotani)
- Diagnostic plots produced when some errors occur had invalid linestyles - use
  markers instead (#95, fixes #94)\
  By [John Omotani](https://github.com/johnomotani)

0.3.1 (11th February 2021)
--------------------------

### New features

- More robust calculation of distances in FineCountour.getDistance(), using
  closest approach to line segments. Can be important for grids with sharp angles (#87)\
  By [Ben Dudson](https://github.com/bendudson)

### Bug fixes

- Ensure FineContours always extend to the end of their parent PsiContour (#86, fixes
  #84)\
  By [Ben Dudson](https://github.com/bendudson)

0.3.0 (25th January 2021)
-------------------------

### Breaking changes

- Changed function used for determining radial positioning of grid points. Function now
  guaranteed to be monotonic, so is more robust. However this does change the output
  slightly compared to previous versions (#64)\
  By [John Omotani](https://github.com/johnomotani)
- Rename target_poloidal_spacing_length, nonorthogonal_target_poloidal_spacing_length,
  nonorthogonal_target_poloidal_spacing_range,
  nonorthogonal_target_poloidal_spacing_range_inner,
  nonorthogonal_target_poloidal_spacing_range_outer renamed to *_target_all_*, because
  extra settings were added to modify each of these parameters individually for each
  target (#75)\
  By [John Omotani](https://github.com/johnomotani)

### New features

- Python script to compare two grid files. Script uses xBOUT. Added to utils/
  subdirectory of repo, and not installed with hypnotoad package (to avoid adding
  dependency on xBOUT) (#83)\
  By [John Omotani](https://github.com/johnomotani)
- Option to start grid at upper-outer divertor instead of lower-inner (#80)\
  By [John Omotani](https://github.com/johnomotani)
- Smoothing copied from IDL hypnotoad for components of curvature vector (#79)\
  By [John Omotani](https://github.com/johnomotani)
- Check for unrecognised options in input files and raise an error if any are found
  (#76)\
  By [John Omotani](https://github.com/johnomotani)
- Extra settings added so spacings can be controlled separately at each target (#74)\
  By [John Omotani](https://github.com/johnomotani)
- EquilibriumRegion.getSqrtPoloidalDistanceFunc() upgraded to ensure that when
  it extrapolates the distance function is always monotonic. This is used when
  y_boundary_guards is greater than 0 (#73)\
  By [John Omotani](https://github.com/johnomotani)
- Command line argument for hypnotoad_geqdsk to call pdb.set_trace() to make it
  easier to debug exceptions with pdb (#72)\
  By [John Omotani](https://github.com/johnomotani)
- When grid file is created from a geqdsk input, save the filename, and the
  contents of the geqdsk file to the grid file (#71, closes #70)\
  By [John Omotani](https://github.com/johnomotani)
- UUID unique identifier saved into each grid file (#67, closes #66)\
  By [John Omotani](https://github.com/johnomotani)

### Bug fixes

- String outputs written as file attributes rather than variables (#69, fixes #68)\
  By [John Omotani](https://github.com/johnomotani)
- Failure when target_poloidal_spacing_length set to number (rather than the
  default None) when y_boundary_guards is non-zero (#64)\
  By [John Omotani](https://github.com/johnomotani)
- BoutMesh options now settable in GUI (#63)
  By [John Omotani](https://github.com/johnomotani)
- Changing settings in File->Preferences caused GUI to crash (#62, fixes #61)\
  By [John Omotani](https://github.com/johnomotani)


0.2.1 (12th January 2021)
-------------------------

### New features

- ``y-coord`` and ``theta`` poloidal coordinates written out by ``BoutMesh`` (#51,
  fixes #49)\
  By [John Omotani](https://github.com/johnomotani)

### Bug fixes

- Timeout if FineContour.refine() takes too long. Length of timeout set by
  refine_timeout option (#58)\
  By [John Omotani](https://github.com/johnomotani)


0.2.0 (25th June 2020)
----------------------

### New features

- Button allowing grid to be regenerated in the gui after nonorthogonal spacing
  options are changed (#26)\
  By [John Omotani](https://github.com/johnomotani)

### Bug fixes

- More robust generation of non-orthogonal grids for tokamak cases (#26)\
  By [John Omotani](https://github.com/johnomotani)


0.1.4 (21st June 2020)
----------------------

### Internal changes

- Make compatible with v5.13 of pyside2 (#45)\
  By [John Omotani](https://github.com/johnomotani)


0.1.3 (20th June 2020)
----------------------

### Internal changes

- Warn instead of failing in case of ImportError when setting "Qt5Agg" backend
  (#44)\
  By [John Omotani](https://github.com/johnomotani)


0.1.2 (20th June 2020)
----------------------

### Internal changes

- Use argparse in command line scripts (#42)\
  By [John Omotani](https://github.com/johnomotani)


0.1.1 (9th June 2020)
---------------------

### New features

- For orthogonal grids, save hy as 'hthe'. Allows backward compatibility with
  codes that compute metric coefficients for themselves. hy and hthe
  definitions are the same for orthogonal grids. They differ for non-orthogonal
  grids, so hy is *not* written as 'hthe' for non-orthogonal grids, to prevent
  silent errors (#39)\
  By [John Omotani](https://github.com/johnomotani)

- Options have associated `doc` attributes, visible as tool-tips in the GUI
  (#33)\
  By [John Omotani](https://github.com/johnomotani)

- Add preferences dialog to GUI (#25)\
  By [Peter Hill](https://github.com/ZedThree) and [John
  Omotani](https://github.com/johnomotani)

- Catch errors in HypnotoadGui.run(), allows changing settings and pressing Run
  button again if there was an error in grid generation (#24)\
  By [John Omotani](https://github.com/johnomotani) and [Peter
  Hill](https://github.com/ZedThree)

- Plot wall in gui window (#23)\
  By [John Omotani](https://github.com/johnomotani)


### Bug fixes

- If an empty string is passed to a value in the options table of the GUI,
  resets the option to its default value. Previously this caused a crash (#25)\
  By [John Omotani](https://github.com/johnomotani)


### Internal changes

- `options` package dependency removed, replaced by
  hypnotoad.utils.OptionsFactory, which creates immutable Options objects.
  WithMeta class used to store the value of each option along with (optionally)
  a `doc` attribute, required type, list of allowed values, and list of checks
  that the value must pass (else raise an exception) (#33)\
  By [John Omotani](https://github.com/johnomotani)

- Make options table keys in HypnotoadGui immutable. Also makes the use of
  Options object in Equilibrium/TokamakEquilibrium more consistent by ensuring
  self.user_options has only been delegated once from
  TokamakEquilibrium.default_options (i.e. the push() method has only been used
  once, otherwise update() or set() are called) (#25)\
  By [John Omotani](https://github.com/johnomotani)


0.1.0 (24th April 2020)
-----------------------

### New features

- Github Action to upload hypnotoad to PyPi on release (#19)\
  By [John Omotani](https://github.com/johnomotani)

- Version number auto-detected from git tags using versioneer (#18)\
  By [John Omotani](https://github.com/johnomotani)

- Graphical user interface, using Qt.py (#17)\
  By [Peter Hill](https://github.com/ZedThree)

- Github Actions automatically run pytest, flake8 and black (#10)\
  By [John Omotani](https://github.com/johnomotani)

- Support for tokamak grids (#2)\
  By [Ben Dudson](https://github.com/bendudson)


### Bug fixes

- Set wall=[] instead of wall=None when there is no wall (#12 #9)
  By [John Omotani](https://github.com/johnomotani)


### Documentation

- doc/whats-new.md documents changes (#22)\
  By [John Omotani](https://github.com/johnomotani)

- doc/developer/RELEASE\_HOWTO.md specifies release process
  By [John Omotani](https://github.com/johnomotani)


### Internal changes

- Rename package and repo from 'hypnotoad2' to 'hypnotoad' (#21)\
  By [John Omotani](https://github.com/johnomotani)


0.0.0 (24th March 2020)
-------------------

Python grid generator for BOUT++, supporting orthogonal or non-orthogonal
grids. Working for TORPEX X-point configuration.
By [John Omotani](https://github.com/johnomotani)
