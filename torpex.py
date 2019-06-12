#!/usr/bin/env python3
"""
Create a BOUT++ grid for TORPEX from an input file giving coil currents and positions

Input file should contain coil parameters, for each coil:
    R: major radius in metres
    Z: major radius in metres
    I: clockwise current in Amps

Note: positions of cell corners are generated first, grid points are then put in the
centre of the cell.
"""

plotStuff = True

from collections import OrderedDict
import warnings

if plotStuff:
    from matplotlib import pyplot
import numpy

from hypnotoad2.mesh import BoutMesh
from hypnotoad2.equilibrium import setDefault, Equilibrium, Point2D, EquilibriumRegion, SolutionError
from hypnotoad2.hypnotoad_options import HypnotoadOptions, HypnotoadInternalOptions

# type for manipulating inforation about magnetic field coils
from collections import namedtuple
Coil = namedtuple('Coil', 'R, Z, I')

class TORPEXMagneticField(Equilibrium):
    """
    Magnetic configuration defined by coil positions and currents for the TORPEX device
    """

    # TORPEX wall is a circle radius 0.2 m around (1 m, 0 m)
    awall = 0.2
    Rcentre = 1.
    Zcentre = 0.

    # Bounding box
    Rmin = 0.8
    Rmax = 1.2
    Zmin = -.2
    Zmax = .2

    # Add TORPEX-specific options and default values
    user_options = HypnotoadOptions.add(
            nx_core = None,
            nx_sol = None,
            ny_inner_lower_divertor = None,
            ny_inner_upper_divertor = None,
            ny_outer_upper_divertor = None,
            ny_outer_lower_divertor = None,
            psi_core = None,
            psi_sol = None,
            psi_inner_sol = None,
            psi_pf = None,
            psi_lower_pf = None,
            psi_upper_pf = None,
            )

    def __init__(self, equilibOptions, meshOptions, **kwargs):

        self.equilibOptions = equilibOptions

        if 'Coils' in self.equilibOptions:
            self.coils = [Coil(**c) for c in self.equilibOptions['Coils']]

            self.magneticFunctionsFromCoils()

            Bt_axis = self.equilibOptions['Bt_axis']
        elif 'gfile' in self.equilibOptions:
            from hypnotoad2.dct_interpolation import DCT_2D

            # load a g-file
            try:
                from pyEquilibrium.geqdsk import Geqdsk
                gfile = Geqdsk(self.equilibOptions['gfile'])
            except AttributeError:
                from boututils.geqdsk import Geqdsk
                gfile = Geqdsk()
                gfile.openFile(self.equilibOptions['gfile'])

            R = numpy.linspace(gfile['rleft'], gfile['rleft'] + gfile['rdim'], gfile['nw'])
            Z = numpy.linspace(gfile['zmid'] - 0.5*gfile['zdim'], gfile['zmid'] + 0.5*gfile['zdim'], gfile['nh'])
            self.magneticFunctionsFromGrid(R, Z, gfile['psirz'])

            Bt_axis = gfile['bcentr']
        elif 'matfile' in self.equilibOptions:
            # Loading directly from the TORPEX-provided matlab file should be slightly
            # more accurate than going via a g-file because g-files don't save full
            # double-precision
            from hypnotoad2.dct_interpolation import DCT_2D

            from scipy.io import loadmat
            eqfile = loadmat(self.equilibOptions['matfile'])['eq']

            R = eqfile['R'][0, 0]
            Z = eqfile['Z'][0, 0]
            # TORPEX psi uses different sign convention from us
            psi = -eqfile['psi'][0, 0]

            Rinds = (R[0, :] >= self.Rmin - 0.02) * (R[0, :] <= self.Rmax + 0.02)
            Zinds = (Z[:, 0] >= self.Zmin - 0.02) * (Z[:, 0] <= self.Zmax + 0.02)

            R = R[:, Rinds]
            R = R[Zinds, :]
            Z = Z[:, Rinds]
            Z = Z[Zinds, :]
            psi = psi[:, Rinds]
            psi = psi[Zinds, :]

            self.magneticFunctionsFromGrid(R[0, :], Z[:, 0], psi)

            Bt = eqfile['Bphi'][0,0]
            RindMid = Bt.shape[1]//2
            ZindMid = Bt.shape[0]//2
            assert eqfile['R'][0,0][ZindMid, RindMid] == 1.
            assert eqfile['Z'][0,0][ZindMid, RindMid] == 0.
            Bt_axis = Bt[ZindMid, RindMid]
        else:
            raise ValueError('Failed to initialise psi function from inputs')

        # TORPEX plasma pressure so low fpol is constant
        self.fpol = lambda psi: Bt_axis / self.Rcentre

        # Set up options read from user input
        self.user_options = TORPEXMagneticField.user_options

        # Set sensible defaults for options
        self.user_options.set(
                xpoint_poloidal_spacing_length = 5.e-2,
                nonorthogonal_xpoint_poloidal_spacing_length = 5.e-2,
                follow_perpendicular_rtol = 2.e-8,
                follow_perpendicular_atol = 1.e-8,
                refine_width = 1.e-2,
                )

        default_options = self.user_options.copy()
        self.user_options.set(**meshOptions)
        self.user_options = self.user_options.push(kwargs)

        setDefault(self.user_options, 'psi_pf', self.user_options.psi_core)
        setDefault(self.user_options, 'psi_lower_pf', self.user_options.psi_pf)
        setDefault(self.user_options, 'psi_upper_pf', self.user_options.psi_pf)

        setDefault(self.user_options, 'psi_inner_sol', self.user_options.psi_sol)

        setDefault(self.user_options, 'poloidal_spacing_delta_psi',
                numpy.abs((self.user_options.psi_core - self.user_options.psi_sol)/20.))

        setDefault(self.user_options, 'nonorthogonal_xpoint_poloidal_spacing_range_inner',
                self.user_options.nonorthogonal_xpoint_poloidal_spacing_range)
        setDefault(self.user_options, 'nonorthogonal_xpoint_poloidal_spacing_range_outer',
                self.user_options.nonorthogonal_xpoint_poloidal_spacing_range)
        setDefault(self.user_options, 'nonorthogonal_target_poloidal_spacing_range_inner',
                self.user_options.nonorthogonal_target_poloidal_spacing_range)
        setDefault(self.user_options, 'nonorthogonal_target_poloidal_spacing_range_outer',
                self.user_options.nonorthogonal_target_poloidal_spacing_range)

        formatstring = '{:<50}|  {:<50}'
        print('\nOptions\n=======')
        print(formatstring.format('Name', 'Value'))
        print('----------------------------------------------------------------------------------------------------')
        for name, value in self.user_options.items():
            valuestring = str(value)
            if value == default_options[name]:
                valuestring += '\t(default)'
            print(formatstring.format(name, valuestring))
        print('')

        # Set up internal options
        # '.push(kwargs)' here lets the kwargs override any values (including for
        # 'internal' options that should not need to be set by the user) set as defaults
        # from HypnotoadOptions
        self.options = HypnotoadInternalOptions.push(kwargs)

        # Make a set of points representing the wall
        self.wall = [self.TORPEX_wall(theta) for theta in
                numpy.linspace(0., 2.*numpy.pi, 100, endpoint=False)]

        try:
            self.x_points = [self.findSaddlePoint(self.Rmin+0.05, self.Rmax-0.05, 0.8*self.Zmin,
                                                  0.8*self.Zmax)]
            self.psi_sep = [self.psi(*self.x_points[0])]
        except SolutionError:
            warnings.warn('Warning: failed to find X-point. Equilibrium generation will '
                    'fail')

    def TORPEX_wall(self, theta):
        """
        Return the location of the TORPEX wall parameterized by the angle theta
        anticlockwise around the centre of the vacuum vessel
        """
        return Point2D(self.Rcentre + self.awall*numpy.cos(theta),
                       self.Zcentre + self.awall*numpy.sin(theta))

    def addWallToPlot(self, npoints=None):
        if npoints is not None:
            theta = numpy.linspace(0., 2.*numpy.pi, npoints+1)
            pyplot.plot(*self.TORPEX_wall(theta))
        else:
            R = [p.R for p in self.wall]
            Z = [p.Z for p in self.wall]
            R.append(R[0])
            Z.append(Z[0])
            pyplot.plot(R, Z)

    def magneticFunctionsFromCoils(self):
        """
        Calculate the poloidal magnetic flux function psi = -R*A_phi, where A_phi is the
        toroidal (anti-clockwise) component of magnetic vector potential due to coils.
        See for example http://physics.usask.ca/~hirose/p812/notes/Ch3.pdf

        The currents in the coils are taken to be positive in the anti-clockwise direction
        here.

        Note e_R x e_phi = e_Z

        A radially increasing psi results in Bp going clockwise in the poloidal plane.
        """
        import sympy
        from sympy.functions.special.elliptic_integrals import elliptic_k, elliptic_e
        import scipy.special

        R,Z = sympy.symbols('R Z')
        mu0 = 4.e-7*sympy.pi

        A_phi = 0*R

        for coil in self.coils:
            # little-r is the vector position from the centre of the coil to (R,Z)
            # sinTheta is the angle between r and the axis through the centre of the coil
            rSquared = R**2 + (Z - coil.Z)**2
            r = sympy.sqrt(rSquared)
            sinTheta = R / r
            kSquared = 4*coil.R*r*sinTheta / (rSquared + coil.R**2 + 2*coil.R*r*sinTheta)
            A_phi += (
              coil.I*coil.R / sympy.sqrt(r**2 + coil.R**2 + 2*coil.R*r*sinTheta) / kSquared
              * ( (2-kSquared)*elliptic_k(kSquared) - 2*elliptic_e(kSquared) )
              )

        # multiply by costant pre-factor
        A_phi *= mu0/sympy.pi

        psi = -R*A_phi
        dpsidR = sympy.diff(psi, R)
        dpsidZ = sympy.diff(psi, Z)
        modGradpsiSquared = dpsidR**2 + dpsidZ**2
        B_R = dpsidZ/R
        B_Z = -dpsidR/R

        self.psi = sympy.lambdify([R,Z], psi, modules=['numpy',
            {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
        self.f_R = sympy.lambdify([R,Z], dpsidR/modGradpsiSquared, modules=['numpy',
            {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
        self.f_Z = sympy.lambdify([R,Z], dpsidZ/modGradpsiSquared, modules=['numpy',
            {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
        self.Bp_R = sympy.lambdify([R,Z], B_R, modules=['numpy',
            {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
        self.Bp_Z = sympy.lambdify([R,Z], B_Z, modules=['numpy',
            {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])

    def makeRegions(self, atol = 2.e-8, npoints=100):
        """
        Find the separatrix and create the regions to grid.

        For TORPEX, follow 4 legs away from the x-point, starting with a rough guess and
        then refining to the separatrix value of A_toroidal.
        """
        assert len(self.x_points) == 1, 'should be one X-point for TORPEX configuration'
        xpoint = self.x_points[0]

        boundary = self.findRoots_1d(
                lambda s: self.psi(*self.wallPosition(s)) - self.psi_sep[0], 4, 0., 1.)

        # put lower left leg first in list, go clockwise
        boundary = boundary[2::-1] + [boundary[3]]

        legnames = ['inner_lower_divertor', 'inner_upper_divertor',
                'outer_upper_divertor', 'outer_lower_divertor']
        kinds = ['wall.X', 'X.wall', 'wall.X', 'X.wall']

        # create input options for EquilibriumRegions
        legoptions = {}
        for i,name in enumerate(legnames):
            options = {}
            options['nx'] = [self.user_options.nx_core, self.user_options.nx_sol]
            options['ny'] = self.user_options['ny_'+name]
            options['kind'] = kinds[i]
            legoptions[name] = options

        # set hard-wired poloidal grid spacing options
        ny_total = sum(opt['ny'] for opt in legoptions.values())

        setDefault(self.options, 'N_norm', ny_total)
        self.regions = OrderedDict()
        wall_vectors = OrderedDict()
        s = numpy.linspace(10.*atol, 1., npoints)
        for i,boundary_position in enumerate(boundary):
            name = legnames[i]
            boundaryPoint = self.wallPosition(boundary_position)
            legR = xpoint.R + s*(boundaryPoint.R - xpoint.R)
            legZ = xpoint.Z + s*(boundaryPoint.Z - xpoint.Z)
            leg = EquilibriumRegion(self, legnames[i], 2, self.user_options,
                    self.options.push(legoptions[name]),
                    [Point2D(R,Z) for R,Z in zip(legR, legZ)], self.psi, self.psi_sep[0])
            self.regions[name] = leg.getRefined(atol=atol, width=0.02)
            wall_vectors[name] = self.wallVector(boundary_position)

        # Make the SeparatrixContours go around clockwise
        # Record the x-point position
        # Record the psi-values of segment boundaries
        # Record the desired radial grid spacing dpsidi at internal boundaries

        dpsidi_sep_inner = (self.user_options.psi_inner_sol - self.psi_sep[0]) / self.user_options.nx_sol
        dpsidi_sep_outer = (self.user_options.psi_sol - self.psi_sep[0]) / self.user_options.nx_sol
        dpsidi_sep_lower = (self.psi_sep[0] - self.user_options.psi_lower_pf) / self.user_options.nx_core
        dpsidi_sep_upper = (self.psi_sep[0] - self.user_options.psi_upper_pf) / self.user_options.nx_core
        if self.user_options.psi_lower_pf < self.user_options.psi_sol:
            dpsidi_sep = min(dpsidi_sep_inner, dpsidi_sep_outer, dpsidi_sep_lower,
                    dpsidi_sep_upper)
        else:
            dpsidi_sep = max(dpsidi_sep_inner, dpsidi_sep_outer, dpsidi_sep_lower,
                    dpsidi_sep_upper)

        # decrease (assuming the factor is <1) the spacing around the separatrix by the
        # factor psi_spacing_separatrix_multiplier
        if self.user_options.psi_spacing_separatrix_multiplier is not None:
            dpsidi_sep = self.user_options.psi_spacing_separatrix_multiplier * dpsidi_sep

        # lower PF
        lower_psi_func = self.getPolynomialGridFunc(self.user_options.nx_core,
                self.user_options.psi_lower_pf, self.psi_sep[0], grad_upper=dpsidi_sep)
        lower_psi_vals = self.make1dGrid(self.user_options.nx_core, lower_psi_func)

        # upper PF
        upper_psi_func = self.getPolynomialGridFunc(self.user_options.nx_core,
                self.user_options.psi_upper_pf, self.psi_sep[0], grad_upper=dpsidi_sep)
        upper_psi_vals = self.make1dGrid(self.user_options.nx_core, upper_psi_func)

        # inner SOL
        inner_psi_func = self.getPolynomialGridFunc(self.user_options.nx_sol,
                self.psi_sep[0], self.user_options.psi_inner_sol, grad_lower=dpsidi_sep)
        inner_psi_vals = self.make1dGrid(self.user_options.nx_sol, inner_psi_func)

        # outer SOL
        outer_psi_func = self.getPolynomialGridFunc(self.user_options.nx_sol,
                self.psi_sep[0], self.user_options.psi_sol, grad_lower=dpsidi_sep)
        outer_psi_vals = self.make1dGrid(self.user_options.nx_sol, outer_psi_func)

        def setupRegion(name, psi_vals1, psi_vals2, reverse):
            r = self.regions[name]
            r.psi_vals = [psi_vals1, psi_vals2]
            r.separatrix_radial_index = 1
            if reverse:
                r.reverse()
                r.xPointsAtEnd[1] = xpoint
                r.wallSurfaceAtStart = wall_vectors[name]
            else:
                r.xPointsAtStart[1] = xpoint
                r.wallSurfaceAtEnd = wall_vectors[name]

        setupRegion('inner_lower_divertor', lower_psi_vals, inner_psi_vals, True)
        setupRegion('inner_upper_divertor', upper_psi_vals, inner_psi_vals, False)
        setupRegion('outer_upper_divertor', upper_psi_vals, outer_psi_vals, True)
        setupRegion('outer_lower_divertor', lower_psi_vals, outer_psi_vals, False)

        # inner lower PF -> outer lower PF
        self.makeConnection('inner_lower_divertor', 0, 'outer_lower_divertor', 0)

        # inner lower SOL -> inner upper SOL
        self.makeConnection('inner_lower_divertor', 1, 'inner_upper_divertor', 1)

        # outer upper PF -> inner upper PF
        self.makeConnection('outer_upper_divertor', 0, 'inner_upper_divertor', 0)

        # outer upper SOL -> outer lower SOL
        self.makeConnection('outer_upper_divertor', 1, 'outer_lower_divertor', 1)

def parseInput(filename):
    import yaml

    with open(filename, 'r') as inputfile:
        inputs = yaml.safe_load(inputfile)
    mesh_inputs = inputs['Mesh']
    del inputs['Mesh']
    equilib_inputs = inputs
    if 'Coils' in equilib_inputs:
        print('Coils:',equilib_inputs['Coils'])
    elif 'gfile' in equilib_inputs:
        print('gfile:', equilib_inputs['gfile'])
    
    return equilib_inputs, mesh_inputs

def createMesh(filename, **kwargs):
    # parse input file
    equilibOptions, meshOptions = parseInput(filename)

    equilibrium = TORPEXMagneticField(equilibOptions, meshOptions, **kwargs)

    print('X-point',equilibrium.x_points[0],'with psi='+str(equilibrium.psi_sep[0]))

    equilibrium.makeRegions()

    return BoutMesh(equilibrium)

def createEqdsk(equilib, *, nR=None, Rmin=None, Rmax=None, nZ=None, Zmin=None, Zmax=None,
        filename='torpex_test.g'):
    from pyEquilibrium.geqdsk import Geqdsk

    R = numpy.linspace(Rmin, Rmax, nR)[numpy.newaxis, :]
    Z = numpy.linspace(Zmin, Zmax, nZ)[:, numpy.newaxis]

    gout = Geqdsk()
    gout.set('nw', nR)
    gout.set('nh', nZ)
    gout.set('rdim', Rmax - Rmin)
    gout.set('zdim', Zmax - Zmin)
    gout.set('rcentr', 0.5*(Rmax - Rmin))
    gout.set('rleft', Rmin)
    gout.set('zmid', 0.5*(Zmax + Zmin))
    gout.set('rmaxis', 1.)
    gout.set('zmaxis', 0.)
    # these values very arbitrary as don't have a magnetic axis
    gout.set('simag', equilib.psi(1., Zmax))
    gout.set('sibry', equilib.psi_sep[0])
    gout.set('bcentr', equilib.fpol(0.)/1.)
    gout.set('current', 0.)
    gout.set('xdum', 0.)

    gout.set('fpol', equilib.fpol(0.) * numpy.ones(nR)) # works for TORPEX because we assume fpol is constant - plasma response neglected
    gout.set('pres', numpy.zeros(nR))
    gout.set('ffprime', numpy.zeros(nR))
    gout.set('pprime', numpy.zeros(nR))
    gout.set('psirz', equilib.psi(R, Z))

    gout.set('rbbbs', [Rmin, Rmax, Rmax, Rmin])
    gout.set('zbbbs', [Zmin, Zmin, Zmax, Zmax])

    theta = numpy.linspace(0., 2.*numpy.pi, 100, endpoint=False)
    gout.set('rlim', [equilib.TORPEX_wall(t).R for t in theta])
    gout.set('zlim', [equilib.TORPEX_wall(t).Z for t in theta])

    gout.dump(filename)

if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]
    gridname = 'torpex.grd.nc'

    mesh = createMesh(filename)

    try:
        mesh.geometry()
    except Exception as e:
        import traceback
        print('There was an exception in mesh.geometry:', str(e))
        print('****************************************')
        traceback.print_tb(e.__traceback__)
        print('****************************************')

    if plotStuff:
        pyplot.figure()
        mesh.equilibrium.plotPotential()
        mesh.equilibrium.addWallToPlot()
        pyplot.plot(*mesh.equilibrium.x_points[0], 'rx')
        mesh.plotPoints(xlow=True, ylow=True, corners=True)
        pyplot.show()

    mesh.writeGridfile(gridname)

    exit(0)
