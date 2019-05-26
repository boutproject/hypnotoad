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

import numpy
from collections import OrderedDict
import warnings
from mesh import BoutMesh
from equilibrium import Equilibrium, Point2D, EquilibriumRegion, SolutionError
if plotStuff:
    from matplotlib import pyplot

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

    def __init__(self, equilibOptions, meshOptions):
        if 'Coils' in equilibOptions:
            from collections import namedtuple

            Coil = namedtuple('Coil', 'R, Z, I')
            self.coils = [Coil(**c) for c in equilibOptions['Coils']]

            self.magneticFunctionsFromCoils()

        elif 'gfile' in equilibOptions:
            from geqdsk import Geqdsk
            from dct_interpolation import DCT_2D

            # load a g-file
            gfile = Geqdsk
            try:
                gfile.read(equilibOptions['gfile'])
            except AttributeError:
                gfile.openFile(equilibOptions['gfile'])

            R = numpy.linspace(gfile['rleft'], gfile['rdim'], gfile['nw'])
            Z = numpy.linspace(gfile['zmid'] - 0.5*gfile['zdim'], gfile['zmid'] + 0.5*gfile['zdim'], gfile['nh'])
            self.magneticFunctionsFromGrid(R, Z, gfile['psirz'])

        else:
            raise ValueError('Failed to initialise psi function from inputs')

        Bt_axis = equilibOptions['Bt_axis']

        # TORPEX plasma pressure so low fpol is constant
        self.fpol = lambda psi: Bt_axis / self.Rcentre

        self.options = meshOptions

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

    def addWallToPlot(self, npoints=100):
        theta = numpy.linspace(0., 2.*numpy.pi, npoints+1)
        pyplot.plot(*self.TORPEX_wall(theta))

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

    def magneticFunctionsFromGrid(self, R, Z, psiRZ):
        from dct_interpolation import DCT_2D

        self._dct = DCT_2D(R, Z, psiRZ)

        self.psi = lambda R, Z: self._dct(R, Z)
        modGradpsiSquared = lambda R, Z: numpy.sqrt(self._dct.ddR(R, Z)**2
                                                         + self._dct.ddZ(R, Z)**2)
        self.f_R = lambda R, Z: self._dct.ddR(R, Z) / modGradpsiSquared(R, Z)
        self.f_Z = lambda R, Z: self._dct.ddZ(R, Z) / modGradpsiSquared(R, Z)
        self.Bp_R = lambda R, Z: self._dct.ddZ(R, Z) / R
        self.Bp_Z = lambda R, Z: -self._dct.ddR(R, Z) / R

    def makeRegions(self, atol = 2.e-8, npoints=100):
        """
        Find the separatrix and create the regions to grid.

        For TORPEX, follow 4 legs away from the x-point, starting with a rough guess and
        then refining to the separatrix value of A_toroidal.
        """
        wall_position = lambda s: self.TORPEX_wall(s*2.*numpy.pi)

        assert len(self.x_points) == 1, 'should be one X-point for TORPEX configuration'
        xpoint = self.x_points[0]

        boundary = self.findRoots_1d(
                lambda s: self.psi(*wall_position(s)) - self.psi_sep[0], 4, 0., 1.)

        # put lower left leg first in list, go clockwise
        boundary = boundary[2::-1] + [boundary[3]]

        boundaryPoints = tuple(wall_position(s) for s in boundary)

        legnames = ['inner_lower_divertor', 'inner_upper_divertor',
                'outer_upper_divertor', 'outer_lower_divertor']

        # create input options for EquilibriumRegions
        legoptions = {}
        nx_core = self.readOption('nx_core')
        nx_sol = self.readOption('nx_sol')
        self.y_boundary_guards = self.readOption('y_boundary_guards')
        for name in legnames:
            options = {}
            options['nx'] = [nx_core, nx_sol]
            options['ny'] = self.readOption('ny_'+name)
            options['y_boundary_guards'] = self.y_boundary_guards
            legoptions[name] = options

        psi_core = self.readOption('psi_core')
        psi_pf = self.readOption('psi_pf', psi_core)
        psi_lower_pf = self.readOption('psi_lower_pf', psi_pf)
        psi_upper_pf = self.readOption('psi_upper_pf', psi_pf)

        psi_sol = self.readOption('psi_sol')
        psi_inner_sol = self.readOption('psi_inner_sol', psi_sol)

        self.regions = OrderedDict()
        s = numpy.linspace(10.*atol, 1., npoints)
        for i,point in enumerate(boundaryPoints):
            name = legnames[i]
            legR = xpoint.R + s*(point.R - xpoint.R)
            legZ = xpoint.Z + s*(point.Z - xpoint.Z)
            leg = EquilibriumRegion(name, 2, legoptions[name],
                    [Point2D(R,Z) for R,Z in zip(legR, legZ)], self.psi, self.psi_sep[0])
            self.regions[name] = leg.getRefined(atol=atol, width=0.02)

        # Make the SeparatrixContours go around clockwise
        # Record the x-point position
        # Record the psi-values of segment boundaries
        # Record the desired radial grid spacing dpsidi at internal boundaries

        self.psi_spacing_separatrix_multiplier = self.readOption('psi_spacing_separatrix_multiplier', None)
        dpsidi_sep_inner = (psi_inner_sol - self.psi_sep[0]) / nx_sol
        dpsidi_sep_outer = (psi_sol - self.psi_sep[0]) / nx_sol
        dpsidi_sep_lower = (self.psi_sep[0] - psi_lower_pf) / nx_core
        dpsidi_sep_upper = (self.psi_sep[0] - psi_upper_pf) / nx_core
        if psi_lower_pf < psi_sol:
            dpsidi_sep = min(dpsidi_sep_inner, dpsidi_sep_outer, dpsidi_sep_lower,
                    dpsidi_sep_upper)
        else:
            dpsidi_sep = max(dpsidi_sep_inner, dpsidi_sep_outer, dpsidi_sep_lower,
                    dpsidi_sep_upper)

        # decrease (assuming the factor is <1) the spacing around the separatrix by the
        # factor psi_spacing_separatrix_multiplier
        if self.psi_spacing_separatrix_multiplier is not None:
            dpsidi_sep = self.psi_spacing_separatrix_multiplier * dpsidi_sep

        # lower PF
        lower_psi_func = self.getPolynomialGridFunc(nx_core, psi_lower_pf,
                self.psi_sep[0], grad_upper=dpsidi_sep)
        lower_psi_vals = self.make1dGrid(nx_core, lower_psi_func)

        # upper PF
        upper_psi_func = self.getPolynomialGridFunc(nx_core, psi_upper_pf,
                self.psi_sep[0], grad_upper=dpsidi_sep)
        upper_psi_vals = self.make1dGrid(nx_core, upper_psi_func)

        # inner SOL
        inner_psi_func = self.getPolynomialGridFunc(nx_sol, self.psi_sep[0],
                psi_inner_sol, grad_lower=dpsidi_sep)
        inner_psi_vals = self.make1dGrid(nx_sol, inner_psi_func)

        # outer SOL
        outer_psi_func = self.getPolynomialGridFunc(nx_sol, self.psi_sep[0],
                psi_sol, grad_lower=dpsidi_sep)
        outer_psi_vals = self.make1dGrid(nx_sol, outer_psi_func)


        # set poloidal grid spacing
        d_xpoint = self.readOption('xpoint_poloidal_spacing_length', 5.e-2)
        ny_total = sum(r.ny_noguards for r in self.regions.values())

        # inner lower
        r = self.regions['inner_lower_divertor']
        r.reverse()
        r.xPointsAtEnd[1] = xpoint
        r.psi_vals = [lower_psi_vals, inner_psi_vals]
        r.sfunc = self.getSqrtPoloidalDistanceFunc(r.distance[-1], 2*r.ny_noguards,
                ny_total, d_upper=0., d_sqrt_upper=d_xpoint)

        # inner upper
        r = self.regions['inner_upper_divertor']
        r.xPointsAtStart[1] = xpoint
        r.psi_vals = [upper_psi_vals, inner_psi_vals]
        r.sfunc = self.getSqrtPoloidalDistanceFunc(r.distance[-1], 2*r.ny_noguards,
                ny_total, d_lower=0., d_sqrt_lower=d_xpoint)

        # outer upper
        r = self.regions['outer_upper_divertor']
        r.reverse()
        r.xPointsAtEnd[1] = xpoint
        r.psi_vals = [upper_psi_vals, outer_psi_vals]
        r.sfunc = self.getSqrtPoloidalDistanceFunc(r.distance[-1], 2*r.ny_noguards,
                ny_total, d_upper=0., d_sqrt_upper=d_xpoint)

        # outer lower
        r = self.regions['outer_lower_divertor']
        r.xPointsAtStart[1] = xpoint
        r.psi_vals = [lower_psi_vals, outer_psi_vals]
        r.sfunc = self.getSqrtPoloidalDistanceFunc(r.distance[-1], 2*r.ny_noguards,
                ny_total, d_lower=0., d_sqrt_lower=d_xpoint)

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

def createMesh(filename):
    # parse input file
    equilibOptions, meshOptions = parseInput(filename)

    equilibrium = TORPEXMagneticField(equilibOptions, meshOptions)

    print('X-point',equilibrium.x_points[0],'with psi='+str(equilibrium.psi_sep[0]))

    equilibrium.makeRegions()

    return BoutMesh(equilibrium, meshOptions)

if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]
    gridname = 'torpex.grd.nc'

    mesh = createMesh(filename)

    mesh.geometry()

    if plotStuff:
        mesh.equilibrium.plotPotential()
        mesh.equilibrium.addWallToPlot()
        pyplot.plot(*mesh.equilibrium.x_points[0], 'rx')
        mesh.plotPoints(xlow=True, ylow=True, corners=True)
        pyplot.show()

    mesh.writeGridfile(gridname)

    exit(0)
