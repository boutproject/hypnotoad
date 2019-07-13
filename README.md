Dependencies
------------

- options module ('pip3 install --user options')
- yaml module ('pip3 install --user PyYAML')
- scipy (recent enough version, tested with 1.3.0 'pip3 install --user --upgrade scipy')


Usage
-----

Options are read and set up in the Equilibrium (child-)class object, and passed
from there to the Mesh (child-)class object.

User-settable options, with their current values, are printed when an
Equilibrium object is created.  Internal options should not need to be set by
the user, but can be overridden with keyword arguments to the Equilibrium
constructor.

Hypnotoad2 can be run either as an executable, which just reads from an input
file, or interactively from a Python shell. To ensure reproducibility, it is
suggested to create your final grid non-interactively. The interactive mode is
intended to make it easier to prototype the grid and find a good set of input
parameters. Once you have found a configuration you are happy with, you can
save the current input parameters with
Equilibrium.saveOptions(filename='hypnotoad\_options.yaml'); this may be
especially useful if you have changed some options from the Python shell with
keyword-arguments.

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
