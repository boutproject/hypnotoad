What's new
==========


0.1.1 (unreleased)
------------------

### New Features

- Add preferences dialog to GUI (#25)\
  By [Peter Hill](https://github.com/ZedThree) and [John
  Omotani](https://github.com/johnomotani)

- Catch errors in HypnotoadGui.run(), allows changing settings and pressing Run
  button again if there was an error in grid generation (#24)\
  By [John Omotani](https://github.com/johnomotani)

- Plot wall in gui window (#23)\
  By [John Omotani](https://github.com/johnomotani)


### Bug fixes

- If an empty string is passed to a value in the options table of the GUI,
  resets the option to its default value. Previously this caused a crash (#25)\
  By [John Omotani](https://github.com/johnomotani)


### Internal changes

- Make options table keys in HypnotoadGui immutable. Also makes the use of
  Options object in Equilibrium/TokamakEquilibrium more consistent by ensuring
  self.user_options has only been delegated once from
  TokamakEquilibrium.default_options (i.e. the push() method has only been used
  once, otherwise update() or set() are called) (#25)\
  By [John Omotani](https://github.com/johnomotani)


0.1.0 (24th April 2020)
-----------------------

### New Features

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
