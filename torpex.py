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

# Torpex parameters
Rmin = 0.8
Rmax = 1.2
Zmin = -.2
Zmax = .2

# Could try to calculate this from range of A_toroidal in the domain or something, but
# easier just to eyeball it from contour plots.
# Only used to calculate some relative tolerances
typical_A = 1.e-7

import numpy
from scipy.optimize import minimize_scalar, brentq, root
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
if plotStuff:
    from matplotlib import pyplot

def TORPEX_wall(theta):
    """
    Return the location of the TORPEX wall parameterized by the angle theta
    anticlockwise around the centre of the vacuum vessel
    """
    # TORPEX wall is a circle radius 0.2 m around (1 m, 0 m)
    awall = 0.2
    Rcentre = 1.
    Zcentre = 0.
    return Point2D(Rcentre + awall*numpy.cos(theta), Zcentre + awall*numpy.sin(theta))

def addWallToPlot(npoints=100):
    theta = numpy.linspace(0., 2.*numpy.pi, npoints+1)
    pyplot.plot(*TORPEX_wall(theta))

class Point2D:
    """
    A point in 2d space.
    Can be added, subtracted, multiplied by scalar
    """
    def __init__(self, R, Z):
        self.R = R
        self.Z = Z

    def __add__(self, other):
        return Point2D(self.R+other.R, self.Z+other.Z)

    def __sub__(self, other):
        return Point2D(self.R-other.R, self.Z-other.Z)

    def __mul__(self, other):
        return Point2D(self.R*other, self.Z*other)

    def __rmul__(self, other):
        return Point2D(self.R*other, self.Z*other)

    def __truediv__(self, other):
        return Point2D(self.R/other, self.Z/other)

    def __iter__(self):
        """
        Along with __next__() allows Point2D class to be treated like a tuple, e.g.
        p = Point2D(1., 0.)
        val = f(*p)
        where f is a function that takes two arguments
        """
        self.iterStep = 0
        return self

    def __next__(self):
        if self.iterStep == 0:
            self.iterStep = 1
            return self.R
        elif self.iterStep == 1:
            self.iterStep = 2
            return self.Z
        else:
            raise StopIteration

    def __repr__(self):
        """
        Allow Point2D to be printed
        """
        return 'Point2D('+str(self.R)+','+str(self.Z)+')'

class MeshContour:
    """
    Represents a contour as a collection of points.
    Includes methods for interpolation.
    Mostly behaves like a list
    """
    def __init__(self, points, A_toroidal, Aval):
        self.points = points

        self.distance = [0.]
        for i in range(1,len(self.points)):
            self.distance.append(
                    self.distance[-1] + calc_distance(self.points[i-1], self.points[i]))

        # Function that evaluates the vector potential at R,Z
        self.A_toroidal = A_toroidal

        # Value of vector potential on this contour
        self.Aval = Aval

    def __iter__(self):
        return self.points.__iter__()

    def __str__(self):
        return self.points.__str__()

    def __getitem__(self, key):
        return self.points.__getitem__(key)

    def __setitem__(self, key, value):
        return self.points.__setitem__(key, value)

    def append(self, point):
        self.points.append(point)
        self.distance.append(
                self.distance[-1] + calc_distance(self.points[-2], self.points[-1]))

    def refine(self, *args, **kwargs):
        self = self.getRefined(*args, **kwargs)

    def getRefined(self, width=.2, atol=2.e-8):
        f = lambda R,Z: self.A_toroidal(R, Z) - self.Aval

        def perpLine(p, tangent, w):
            # p - point through which to draw perpLine
            # tangent - vector tangent to original curve, result will be perpendicular to this
            # w - width on either side of p to draw the perpLine to
            modTangent = numpy.sqrt(tangent.R**2 + tangent.Z**2)
            perpIdentityVector = Point2D(tangent.Z/modTangent, -tangent.R/modTangent)
            return lambda s: p + 2.*(s-0.5)*w*perpIdentityVector

        def refinePoint(p, tangent):
            converged = False
            w = width
            sp = []
            ep = []
            while not converged:
                try:
                    pline = perpLine(p, tangent, w)
                    sp.append(pline(0.))
                    ep.append(pline(1.))
                    snew, info = brentq(lambda s: f(*pline(s)), 0., 1., xtol=atol, full_output=True)
                    converged = info.converged
                except ValueError:
                    pass
                w /= 2.
                if w < atol:
                    if numpy.abs(f(*pline(0.))) < atol*typical_A and numpy.abs(f(*pline(1.))) < atol*typical_A:
                        # f is already so close to 0 that the point does not need refining
                        snew = 0.5
                        pass
                    else:
                        raise ValueError("Could not find interval to refine point")

            return pline(snew)

        newpoints = []
        newpoints.append(refinePoint(self.points[0], self.points[1] - self.points[0]))
        for i,p in enumerate(self.points[1:-1]):
            newpoints.append(refinePoint(p, self.points[i+1] - self.points[i-1]))
        newpoints.append(refinePoint(self.points[-1], self.points[-1] - self.points[-2]))

        return MeshContour(newpoints, self.A_toroidal, self.Aval)

    def interpFunction(self):
        distance = numpy.array(numpy.float64(self.distance))
        R = numpy.array(numpy.float64([p.R for p in self.points]))
        Z = numpy.array(numpy.float64([p.Z for p in self.points]))
        interpR = interp1d(distance, R, kind='cubic',
                           assume_sorted=True, fill_value='extrapolate')
        interpZ = interp1d(distance, Z, kind='cubic',
                           assume_sorted=True, fill_value='extrapolate')
        total_distance = distance[-1]
        return lambda s: Point2D(interpR(s*total_distance), interpZ(s*total_distance))

    def getRegridded(self, npoints, width=1.e-4, atol=2.e-8, sfunc=None, extend=0):
        """
        Interpolate onto set of npoints points, then refine positions.
        By default points are uniformly spaced, this can be changed by passing 'sfunc'
        which replaces the uniform interval 's' with 's=sfunc(s)'.
        'extend' extends the contour past its existing end by a number of points
        Returns a new MeshContour.
        """
        s = numpy.linspace(0., (npoints-1+extend)/(npoints-1), npoints+extend)
        if sfunc is not None:
            s = sfunc(s)
        interp = self.interpFunction()
        new_contour = MeshContour([interp(x) for x in s], self.A_toroidal, self.Aval)
        return new_contour.getRefined(width, atol)

    def plot(self, *args, **kwargs):
        pyplot.plot([x.R for x in self], [x.Z for x in self], *args, **kwargs)

class Mesh:
    """
    Mesh quantities to be written to a grid file for BOUT++
    """
    def __init__(self, meshOptions, A_toroidal, f_R, f_Z, Bp_R, Bp_Z, fpol, xpoint, A_xpoint,
                 separatrixLegs):
        self.orthogonal = meshOptions['orthogonal']
        self.nx = meshOptions['nx']
        self.ny = meshOptions['ny']
        self.ixseps = meshOptions['ixseps']
        self.jyseps1 = meshOptions['jyseps1']
        self.jyseps2 = meshOptions['jyseps2']
        self.ny_inner = meshOptions['ny_inner']
        self.psi_inner = meshOptions['psi_inner']
        self.psi_outer = meshOptions['psi_outer']
        self.y_boundary_guards = meshOptions['y_boundary_guards']

        self.A_toroidal = A_toroidal
        self.f_R = f_R
        self.f_Z = f_Z
        self.Bp_R = Bp_R
        self.Bp_Z = Bp_Z

        # poloidal current function, gives B_toroidal
        self.fpol = fpol

        self.xpoint = xpoint
        self.A_xpoint = A_xpoint

        assert self.orthogonal # non-orthogonal not implelemented yet

        # number of points along each leg
        self.npol_leg = []
        self.npol_leg.append(self.jyseps1+1)
        self.npol_leg.append(self.ny_inner - (self.jyseps1+1))
        self.npol_leg.append(self.jyseps2+1 - self.ny_inner)
        self.npol_leg.append(self.ny - (self.jyseps2+1))

        # number of radial points 'inside' separatrix for each leg
        self.nrad_pf = self.ixseps

        # number of radial points 'outside' separatrix for each leg
        self.nrad_sol = self.nx - self.ixseps

        # generate points for cell centres and faces
        # wider poloidal spacing along separatrix near X-point, so orthogonal grid does
        # not get too squashed
        sfunc = lambda s: s**0.5
        self.separatrixLegs = [leg.getRegridded(2*np+1, sfunc=sfunc,
                                                extend=2*self.y_boundary_guards)
                               for leg,np in zip(separatrixLegs, self.npol_leg)]

        dpsi_inner = numpy.abs(A_xpoint - self.psi_inner)/float(self.nrad_pf)
        dpsi_outer = numpy.abs(self.psi_outer - A_xpoint)/float(self.nrad_sol)
        # make grid spacing (dpsi_inner+dpsi_outer)/2. at separatrix
        # in index space for indices of cell faces, psi needs to go through psi_inner at
        # 0, A_xpoint at ixseps and psi_outer at nx+1
        # for now use quadratic fit, leave grid refinement for later...
        # psi(i) = c + b*i + a*i^2
        # psi(0) = c = psi_inner
        # psi(ixseps) = A_xpoint
        #   A_xpoint - psi_inner = ixseps(b + a*ixseps)
        #   (A_xpoint - psi_inner)/ixseps = b + a*ixseps
        # psi(nx+1) = psi_outer
        #   psi_outer - psi_inner = (nx+1)*(b + a*(nx+1))
        #   (psi_outer - psi_inner)/(nx+1) = b + a*(nx+1)
        # a*(nx+1-ixseps) = (psi_outer-psi_inner)/(nx+1) - (A_xpoint-psi_inner)/ixseps
        # a = ((psi_outer-psi_inner)/(nx+1) - (A_xpoint-psi_inner)/ixseps) / (nx+1-ixseps)
        # b = ( (psi_outer-psi_inner)/(nx+1)**2 - (A_xpoint-psi_inner)/ixseps**2 ) / (1/(nx+1) - 1/ixseps)
        a = ((self.psi_outer-self.psi_inner)/(self.nx+1) -
                (self.A_xpoint-self.psi_inner)/self.ixseps) / (self.nx+1-self.ixseps)
        b = ( (self.psi_outer-self.psi_inner)/(self.nx+1)**2 -
                (self.A_xpoint-self.psi_inner)/self.ixseps**2 ) / (1./(self.nx+1) - 1./self.ixseps)
        c = self.psi_inner
        psi_index = lambda i: a*i**2 + b*i + c
        psi_face_vals_inner = [psi_index(i) for i in range(self.nrad_pf+1)]
        psi_face_vals_outer = [psi_index(i) for i in range(self.ixseps, self.nx+2)]
        self.psi_vals_inner = []
        self.psi_vals_outer = []
        self.dx_list = numpy.zeros(2*self.nx+1)
        for i in range(self.nrad_pf):
            psi_m = psi_face_vals_inner[i]
            psi_p = psi_face_vals_inner[i+1]
            self.psi_vals_inner.append(psi_m)
            self.psi_vals_inner.append(0.5*(psi_m + psi_p))
            if i > 0:
                self.dx_list[2*i] = numpy.abs(0.5*(psi_p - psi_face_vals_inner[i-1]))
            else:
                self.dx_list[2*i] = numpy.abs(psi_p - psi_m)
            self.dx_list[2*i+1] = numpy.abs(psi_p - psi_m)
        self.psi_vals_inner.append(psi_face_vals_inner[-1])
        self.dx_list[2*self.nrad_pf+2] = numpy.abs(0.5*(psi_face_vals_outer[1] -
                                    psi_face_vals_inner[-2]))
        for i in range(self.nrad_sol):
            psi_m = psi_face_vals_outer[i]
            psi_p = psi_face_vals_outer[i+1]
            self.psi_vals_outer.append(psi_m)
            self.psi_vals_outer.append(0.5*(psi_m + psi_p))
            if i > 0:
                self.dx_list[2*self.nrad_pf+2*i] = numpy.abs(0.5*(psi_p - psi_face_vals_outer[i-1]))
            self.dx_list[2*self.nrad_pf+2*i+1] = numpy.abs(psi_p - psi_m)
        self.psi_vals_outer.append(psi_face_vals_outer[-1])
        self.dx_list[2*self.nx] = numpy.abs(psi_face_vals_outer[-1] - psi_face_vals_outer[-2])

        print('Mesh get points')
        self.contours_pf = [[] for i in range(4)]
        self.contours_sol = [[] for i in range(4)]

        for i,leg in enumerate(self.separatrixLegs):
            print('leg',i,'...')
            sep_points = list(leg)
            if i%2 == 0:
                sep_points.reverse()

            perp_points_inner = followPerpendicular(self.f_R, self.f_Z, sep_points[0],
                    self.A_xpoint, self.psi_vals_inner[-2::-1])
            perp_points_inner.reverse()
            for j,point in enumerate(perp_points_inner):
                self.contours_pf[i].append(MeshContour([point], self.A_toroidal, self.psi_vals_inner[j]))

            perp_points_outer = followPerpendicular(self.f_R, self.f_Z, sep_points[0],
                    self.A_xpoint, self.psi_vals_outer[1:])
            for j,point in enumerate(perp_points_outer):
                self.contours_sol[i].append(MeshContour([point], self.A_toroidal, self.psi_vals_outer[j+1]))

            for count,p in enumerate(sep_points[1:]):
                print('point',count)
                perp_points_inner = followPerpendicular(self.f_R, self.f_Z, p,
                        self.A_xpoint, self.psi_vals_inner[-2::-1])
                perp_points_inner.reverse()
                for j,point in enumerate(perp_points_inner):
                    self.contours_pf[i][j].append(point)

                perp_points_outer = followPerpendicular(self.f_R, self.f_Z, p,
                        self.A_xpoint, self.psi_vals_outer[1:])
                for j,point in enumerate(perp_points_outer):
                    self.contours_sol[i][j].append(point)

            leg.refine(width=1.e-5)

    def geometry(self):
        """
        Calculate geometrical quantities for BOUT++
        """

        self.x_regions = [0, self.ixseps, self.nx]
        self.y_regions = [0,
                     self.jyseps1 + 1 + self.y_boundary_guards,
                     self.ny_inner + 2*self.y_boundary_guards,
                     self.jyseps2 + 1 + 3*self.y_boundary_guards,
                     self.ny + 4*self.y_boundary_guards
                    ]
        self.Rxy = numpy.zeros([self.nx, self.ny + 4*self.y_boundary_guards])
        self.Rxy_ylow = numpy.zeros([self.nx, self.ny + 4*self.y_boundary_guards])
        self.Zxy = numpy.zeros([self.nx, self.ny + 4*self.y_boundary_guards])
        self.Zxy_ylow = numpy.zeros([self.nx, self.ny + 4*self.y_boundary_guards])

        self.Rcorners = numpy.zeros([self.nx + 1, self.ny + 4*self.y_boundary_guards + 1])
        self.Zcorners = numpy.zeros([self.nx + 1, self.ny + 4*self.y_boundary_guards + 1])

        # these store the corner positions at the end of the upper left divertor guard
        # cells which would otherwise be missed because they would need to be in the same
        # place in the arrays as the upper right divertor positions
        self.Rcorners_extra = numpy.zeros(self.nx+1)
        self.Zcorners_extra = numpy.zeros(self.nx+1)

        for i,contours in enumerate(self.contours_pf):
            if i == 3:
                y_extra = 1
                yup = None
            else:
                y_extra = 0
                yup = -1

            self.Rxy[self.x_regions[0]:self.x_regions[1], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.R for p in contour[1::2]] for contour in contours[1::2]]

            self.Rxy_ylow[self.x_regions[0]:self.x_regions[1], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.R for p in contour[0:-1:2]] for contour in contours[1::2]]

            self.Zxy[self.x_regions[0]:self.x_regions[1], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.Z for p in contour[1::2]] for contour in contours[1::2]]

            self.Zxy_ylow[self.x_regions[0]:self.x_regions[1], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.Z for p in contour[0:-1:2]] for contour in contours[1::2]]

            sep = list(self.separatrixLegs[i])
            if i%2 == 0:
                sep.reverse()
            self.Rcorners[self.x_regions[0]:self.x_regions[1], self.y_regions[i]:self.y_regions[i+1]+y_extra] = \
                    [[p.R for p in contour[0:yup:2]] for contour in contours[0::2]]
            self.Zcorners[self.x_regions[0]:self.x_regions[1], self.y_regions[i]:self.y_regions[i+1]+y_extra] = \
                    [[p.Z for p in contour[0:yup:2]] for contour in contours[0::2]]

            if i==1:
                self.Rcorners_extra[self.x_regions[0]:self.x_regions[1]+1] = \
                        [contour[-1].R for contour in contours[0::2]+[sep]]
                self.Zcorners_extra[self.x_regions[0]:self.x_regions[1]+1] = \
                        [contour[-1].Z for contour in contours[0::2]+[sep]]

        for i,contours in enumerate(self.contours_sol):
            if i == 3:
                y_extra = 1
                yup = None
            else:
                y_extra = 0
                yup = -1

            self.Rxy[self.x_regions[1]:self.x_regions[2], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.R for p in contour[1::2]] for contour in contours[1::2]]

            self.Rxy_ylow[self.x_regions[1]:self.x_regions[2], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.R for p in contour[0:-1:2]] for contour in contours[1::2]]

            self.Zxy[self.x_regions[1]:self.x_regions[2], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.Z for p in contour[1::2]] for contour in contours[1::2]]

            self.Zxy_ylow[self.x_regions[1]:self.x_regions[2], self.y_regions[i]:self.y_regions[i+1]] = \
                    [[p.Z for p in contour[0:-1:2]] for contour in contours[1::2]]

            sep = list(self.separatrixLegs[i])
            if i%2 == 0:
                sep.reverse()
            self.Rcorners[self.x_regions[1]:self.x_regions[2]+1, self.y_regions[i]:self.y_regions[i+1]+y_extra] = \
                    [[p.R for p in contour[0:yup:2]] for contour in [sep]+contours[0::2]]
            self.Zcorners[self.x_regions[1]:self.x_regions[2]+1, self.y_regions[i]:self.y_regions[i+1]+y_extra] = \
                    [[p.Z for p in contour[0:yup:2]] for contour in [sep]+contours[0::2]]

            if i==1:
                self.Rcorners_extra[self.x_regions[1]:self.x_regions[2]+1] = \
                        [contour[-1].R for contour in [sep]+contours[0::2]]
                self.Zcorners_extra[self.x_regions[1]:self.x_regions[2]+1] = \
                        [contour[-1].Z for contour in [sep]+contours[0::2]]

        self.psixy = self.A_toroidal(self.Rxy, self.Zxy)
        self.psixy_ylow = self.A_toroidal(self.Rxy_ylow, self.Zxy_ylow)
        self.dx = numpy.zeros([self.nx, self.ny + 2*self.y_boundary_guards])
        self.dx[:] = numpy.array(self.dx_list[1::2])[:, numpy.newaxis]
        self.dx_ylow = numpy.zeros([self.nx, self.ny + 2*self.y_boundary_guards])
        self.dx_ylow[:] = numpy.array(self.dx_list[0:-1:2])[:, numpy.newaxis]
        if self.psi_inner > self.psi_outer:
            # x-coordinate is -psixy so x always increases radially across grid
            self.bpsign = -1.
            self.xcoord = -self.psixy
        else:
            self.bpsign = 1.
            self.xcoord = self.psixy

        self.dy = 2.*numpy.pi/float(self.ny) * numpy.ones([self.nx, self.ny + 2*self.y_boundary_guards])
        self.dy_ylow = 2.*numpy.pi/float(self.ny) * numpy.ones([self.nx, self.ny + 2*self.y_boundary_guards])

        self.Brxy = self.Bp_R(self.Rxy, self.Zxy)
        self.Brxy_ylow = self.Bp_R(self.Rxy_ylow, self.Zxy_ylow)
        self.Bzxy = self.Bp_Z(self.Rxy, self.Zxy)
        self.Bzxy_ylow = self.Bp_Z(self.Rxy_ylow, self.Zxy_ylow)
        self.Bpxy = numpy.sqrt(self.Brxy**2 + self.Bzxy**2)
        self.Bpxy_ylow = numpy.sqrt(self.Brxy_ylow**2 + self.Bzxy_ylow**2)
        # determine direction - dot Bp with Grad(y) vector
        # evaluate in 'sol' at outer radial boundary
        Bp_dot_grady = self.Brxy[-1, self.jyseps2+1]*(self.Rxy[-1, self.jyseps2+2] - self.Rxy[-1, self.jyseps2]) +\
                       self.Bzxy[-1, self.jyseps2+1]*(self.Zxy[-1, self.jyseps2+2] - self.Zxy[-1, self.jyseps2])
        if Bp_dot_grady < 0.:
            print("Poloidal field is in opposite direction to Grad(theta) -> Bp negative")
            self.Bpxy = -self.Bpxy
            self.Bpxy_ylow = -self.Bpxy_ylow
            if self.bpsign > 0.:
                raise ValueError("Sign of Bp should be negative?")
        else:
            if self.bpsign < 0.:
                raise ValueError("Sign of Bp should be positive?")

        # Get toroidal field from poloidal current function fpol
        self.Btxy = self.fpol / self.Rxy
        self.Btxy_ylow = self.fpol / self.Rxy_ylow

        self.Bxy = numpy.sqrt(self.Bpxy**2 + self.Btxy**2)
        self.Bxy_ylow = numpy.sqrt(self.Bpxy_ylow**2 + self.Btxy_ylow**2)

        # poloidal arc-length between points
        self.hthe, self.hthe_ylow = self.get_hthe()

    def DDX(self, f):
        result = numpy.zeros([self.nx, self.ny + 2*self.y_boundary_guards])
        result[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.*self.dx[1:-1])
        result[0, :] = (-1.5*f[0, :] + 2.*f[1,:] - 0.5*f[2, :]) / self.dx[0]
        result[-1, :] = (1.5*f[-1, :] - 2.*f[-2,:] + 0.5*f[-3, :]) / self.dx[-1]
        return result

    def writeGridfile(self, filename):
        from boututils.datafile import DataFile

        with DataFile(filename, create=True) as f:
            f.write('nx', self.nx)
            f.write('ny', self.ny)
            f.write('y_boundary_guards', self.y_boundary_guards)
            f.write('Rxy', self.Rxy)
            f.write('Rxy_ylow', self.Rxy_ylow)
            f.write('Zxy', self.Zxy)
            f.write('Zxy_ylow', self.Zxy_ylow)
            f.write('psixy', self.psixy)
            f.write('psixy_ylow', self.psixy_ylow)
            f.write('dx', self.dx)
            f.write('dx_ylow', self.dx_ylow)
            f.write('dy', self.dy)
            f.write('dy_ylow', self.dy_ylow)
            f.write('Bpxy', self.Bpxy)
            f.write('Bpxy_ylow', self.Bpxy_ylow)
            f.write('Btxy', self.Btxy)
            f.write('Btxy_ylow', self.Btxy_ylow)
            f.write('Bxy', self.Bxy)
            f.write('Bxy_ylow', self.hthe_ylow)
            f.write('hthe', self.Bxy)
            f.write('hthe_ylow', self.hthe_ylow)

    def plot2D(self, f, title=None, ylow=False):
        try:
            vmin = f.min()
            vmax = f.max()

            # leg 0
            R = self.Rcorners[:, self.y_regions[0]:self.y_regions[1]+1].copy()
            Z = self.Zcorners[:, self.y_regions[0]:self.y_regions[1]+1].copy()
            # fix upper PF region corners
            R[self.x_regions[0]:self.x_regions[1], -1] = self.Rcorners[self.x_regions[0]:self.x_regions[1], self.y_regions[3]]
            Z[self.x_regions[0]:self.x_regions[1], -1] = self.Zcorners[self.x_regions[0]:self.x_regions[1], self.y_regions[3]]
            pyplot.pcolor(R, Z, f[:, self.y_regions[0]:self.y_regions[1]], vmin=vmin, vmax=vmax)

            # leg 1
            R = self.Rcorners[:, self.y_regions[1]:self.y_regions[2]+1].copy()
            Z = self.Zcorners[:, self.y_regions[1]:self.y_regions[2]+1].copy()
            # fix upper corners
            R[:, -1] = self.Rcorners_extra
            Z[:, -1] = self.Zcorners_extra
            pyplot.pcolor(R, Z, f[:, self.y_regions[1]:self.y_regions[2]], vmin=vmin, vmax=vmax)

            # leg 2
            R = self.Rcorners[:, self.y_regions[2]:self.y_regions[3]+1].copy()
            Z = self.Zcorners[:, self.y_regions[2]:self.y_regions[3]+1].copy()
            # fix upper PF region corners
            R[self.x_regions[0]:self.x_regions[1], -1] = self.Rcorners[self.x_regions[0]:self.x_regions[1], self.y_regions[1]]
            Z[self.x_regions[0]:self.x_regions[1], -1] = self.Zcorners[self.x_regions[0]:self.x_regions[1], self.y_regions[1]]
            pyplot.pcolor(R, Z, f[:, self.y_regions[2]:self.y_regions[3]], vmin=vmin, vmax=vmax)

            # leg 3
            R = self.Rcorners[:, self.y_regions[3]:self.y_regions[4]+1].copy()
            Z = self.Zcorners[:, self.y_regions[3]:self.y_regions[4]+1].copy()
            pyplot.pcolor(R, Z, f[:, self.y_regions[3]:self.y_regions[4]], vmin=vmin, vmax=vmax)

            pyplot.colorbar()
        except NameError:
            raise NameError('Some variable has not been defined yet: have you called Mesh.geometry()?')

    def get_hthe(self):
        """
        Get poloidal arc length on grid.
        Similar logic to plot2D in terms of R,Z on grid needed
        """
        hthe = numpy.zeros([self.nx, self.ny + 4*self.y_boundary_guards])
        hthe_ylow = numpy.zeros([self.nx, self.ny + 4*self.y_boundary_guards])

        # leg 0
        R = self.Rxy_ylow[:, self.y_regions[0]:self.y_regions[1]+1].copy()
        Z = self.Zxy_ylow[:, self.y_regions[0]:self.y_regions[1]+1].copy()
        # fix upper PF region points
        R[self.x_regions[0]:self.x_regions[1], -1] = self.Rxy_ylow[self.x_regions[0]:self.x_regions[1], self.y_regions[3]]
        Z[self.x_regions[0]:self.x_regions[1], -1] = self.Zxy_ylow[self.x_regions[0]:self.x_regions[1], self.y_regions[3]]
        hthe[:, self.y_regions[0]:self.y_regions[1]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)
        # now for ylow
        R = numpy.zeros([self.nx, self.y_regions[1]+1])
        Z = numpy.zeros([self.nx, self.y_regions[1]+1])
        R[:, 1:] = self.Rxy[:, self.y_regions[0]:self.y_regions[1]]
        Z[:, 1:] = self.Zxy[:, self.y_regions[0]:self.y_regions[1]]
        # 'fix' lower points
        R[:, 0] = self.Rxy_ylow[:, self.y_regions[0]]
        Z[:, 0] = self.Zxy_ylow[:, self.y_regions[0]]
        hthe_ylow[:, self.y_regions[0]:self.y_regions[1]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)
        hthe_ylow[:, self.y_regions[0]] *= 2. # double because only used distance from ylow point to next cell centre

        # leg 1
        R = self.Rxy_ylow[:, self.y_regions[1]:self.y_regions[2]+1].copy()
        Z = self.Zxy_ylow[:, self.y_regions[1]:self.y_regions[2]+1].copy()
        # fix upper points
        contours = self.contours_pf[1]
        R[self.x_regions[0]:self.x_regions[1], -1] = [contour[-1].R for contour in contours[0::2]]
        Z[self.x_regions[0]:self.x_regions[1], -1] = [contour[-1].Z for contour in contours[0::2]]
        contours = self.contours_sol[1]
        R[self.x_regions[1]:self.x_regions[2], -1] = [contour[-1].R for contour in contours[0::2]]
        Z[self.x_regions[1]:self.x_regions[2], -1] = [contour[-1].Z for contour in contours[0::2]]
        hthe[:, self.y_regions[1]:self.y_regions[2]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)
        # now for ylow
        R = self.Rxy[:, self.y_regions[1]-1:self.y_regions[2]].copy()
        Z = self.Zxy[:, self.y_regions[1]-1:self.y_regions[2]].copy()
        # fix PF lower points
        R[self.x_regions[0]:self.x_regions[1],0] = \
                self.Rxy[self.x_regions[0]:self.x_regions[1],self.y_regions[3]-1]
        Z[self.x_regions[0]:self.x_regions[1],0] = \
                self.Zxy[self.x_regions[0]:self.x_regions[1],self.y_regions[3]-1]
        hthe_ylow[:, self.y_regions[1]:self.y_regions[2]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)

        # leg 2
        R = self.Rxy_ylow[:, self.y_regions[2]:self.y_regions[3]+1].copy()
        Z = self.Zxy_ylow[:, self.y_regions[2]:self.y_regions[3]+1].copy()
        # fix upper PF region points
        R[self.x_regions[0]:self.x_regions[1], -1] = self.Rxy_ylow[self.x_regions[0]:self.x_regions[1], self.y_regions[1]]
        Z[self.x_regions[0]:self.x_regions[1], -1] = self.Zxy_ylow[self.x_regions[0]:self.x_regions[1], self.y_regions[1]]
        hthe[:, self.y_regions[2]:self.y_regions[3]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)
        # now for ylow
        R = self.Rxy[:, self.y_regions[2]-1:self.y_regions[3]].copy()
        Z = self.Zxy[:, self.y_regions[2]-1:self.y_regions[3]].copy()
        # 'fix' lower points
        R[:, 0] = self.Rxy_ylow[:, self.y_regions[2]]
        Z[:, 0] = self.Zxy_ylow[:, self.y_regions[2]]
        hthe_ylow[:, self.y_regions[2]:self.y_regions[3]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)
        hthe_ylow[:, self.y_regions[2]] *= 2. # double because only used distance from ylow point to next cell centre

        # leg 3
        R = numpy.zeros([self.nx, self.y_regions[4]-self.y_regions[3]+1])
        Z = numpy.zeros([self.nx, self.y_regions[4]-self.y_regions[3]+1])
        R[:,:-1] = self.Rxy_ylow[:, self.y_regions[3]:self.y_regions[4]]
        Z[:,:-1] = self.Zxy_ylow[:, self.y_regions[3]:self.y_regions[4]]
        # fix upper points
        contours = self.contours_pf[3]
        R[self.x_regions[0]:self.x_regions[1], -1] = [contour[-1].R for contour in contours[0::2]]
        Z[self.x_regions[0]:self.x_regions[1], -1] = [contour[-1].Z for contour in contours[0::2]]
        contours = self.contours_sol[3]
        R[self.x_regions[1]:self.x_regions[2], -1] = [contour[-1].R for contour in contours[0::2]]
        Z[self.x_regions[1]:self.x_regions[2], -1] = [contour[-1].Z for contour in contours[0::2]]
        hthe[:, self.y_regions[3]:self.y_regions[4]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)
        # now for ylow
        R = self.Rxy[:, self.y_regions[3]-1:self.y_regions[4]].copy()
        Z = self.Zxy[:, self.y_regions[3]-1:self.y_regions[4]].copy()
        # fix PF lower points
        R[self.x_regions[0]:self.x_regions[1],0] = \
                self.Rxy[self.x_regions[0]:self.x_regions[1],self.y_regions[1]-1]
        Z[self.x_regions[0]:self.x_regions[1],0] = \
                self.Zxy[self.x_regions[0]:self.x_regions[1],self.y_regions[1]-1]
        hthe_ylow[:, self.y_regions[3]:self.y_regions[4]] = \
                numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)

        return hthe, hthe_ylow

def parseInput(filename):
    import yaml
    from collections import namedtuple

    with open(filename, 'r') as inputfile:
        coil_inputs, mesh_inputs = yaml.safe_load_all(inputfile)
    print('Coils:',coil_inputs['Coils'])
    
    Coil = namedtuple('Coil', 'R, Z, I')
    return [Coil(**c) for c in coil_inputs['Coils']], coil_inputs['Bt_axis'], mesh_inputs['Mesh']

def magneticFunctions(coils):
    """
    Calculate toroidal (anticlockwise) component of magnetic vector potential due to coils
    See for example http://physics.usask.ca/~hirose/p812/notes/Ch3.pdf

    Note e_R x e_phi = e_Z
    """
    import sympy
    from sympy.functions.special.elliptic_integrals import elliptic_k, elliptic_e
    import scipy.special

    R,Z = sympy.symbols('R Z')
    mu0 = 4.e-7*sympy.pi

    potential = 0*R

    for coil in coils:
        # little-r is the vector position from the centre of the coil to (R,Z)
        # sinTheta is the angle between r and the axis through the centre of the coil
        rSquared = R**2 + (Z - coil.Z)**2
        r = sympy.sqrt(rSquared)
        sinTheta = R / r
        kSquared = 4*coil.R*r*sinTheta / (rSquared + coil.R**2 + 2*coil.R*r*sinTheta)
        potential += -(
          coil.I*coil.R / sympy.sqrt(r**2 + coil.R**2 + 2*coil.R*r*sinTheta) / kSquared
          * ( (2-kSquared)*elliptic_k(kSquared) - 2*elliptic_e(kSquared) )
          )

    # multiply by costant pre-factor
    potential *= mu0/sympy.pi

    dAdR = sympy.diff(potential, R)
    dAdZ = sympy.diff(potential, Z)
    modGradASquared = dAdR**2 + dAdZ**2

    A_func = sympy.lambdify([R,Z], potential, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    f_R_func = sympy.lambdify([R,Z], dAdR/modGradASquared, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    f_Z_func = sympy.lambdify([R,Z], dAdZ/modGradASquared, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    Bp_R_func = sympy.lambdify([R,Z], dAdZ, modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    Bp_Z_func = sympy.lambdify([R,Z], -1/R * sympy.diff(R*potential, R), modules=['numpy',
        {'elliptic_k':scipy.special.ellipk, 'elliptic_e':scipy.special.ellipe}])
    return (A_func, f_R_func, f_Z_func, Bp_R_func, Bp_Z_func)

def plotPotential(potential, npoints=100, ncontours=40):
    pyplot.figure()
    R = numpy.linspace(Rmin, Rmax, npoints)
    Z = numpy.linspace(Zmin, Zmax, npoints)
    contours = pyplot.contour(
            R, Z, potential(R[:,numpy.newaxis], Z[numpy.newaxis,:]).T, ncontours)
    pyplot.clabel(contours, inline=False, fmt='%1.3g')
    pyplot.axes().set_aspect('equal')

def calc_distance(p1, p2):
    d = p2 - p1
    return numpy.sqrt(d.R**2 + d.Z**2)

def findMinimum_1d(pos1, pos2, f, atol=1.e-14):
    coords = lambda s: pos1 + s*(pos2-pos1)
    result = minimize_scalar(lambda s: f(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
    if result.success:
        return coords(result.x)
    else:
        raise ValueError('findMinimum_1d failed')

def findMaximum_1d(pos1, pos2, f, atol=1.e-14):
    coords = lambda s: pos1 + s*(pos2-pos1)
    # minimize -f to find maximum
    result = minimize_scalar(lambda s: -f(*coords(s)), method='bounded', bounds=(0., 1.), options={'xatol':atol})
    if result.success:
        return coords(result.x)
    else:
        raise ValueError('findMaximum_1d failed')

def findExtremum_1d(pos1, pos2, f, rtol=1.e-5, atol=1.e-14):
    smallDistance = 10.*rtol*calc_distance(pos1, pos2)

    minpos = findMinimum_1d(pos1, pos2, f, atol)
    if calc_distance(pos1,minpos) > smallDistance and calc_distance(pos2,minpos) > smallDistance:
        # minimum is not at either end of the interval
        return minpos, True

    maxpos = findMaximum_1d(pos1, pos2, f, atol)
    if calc_distance(pos1,maxpos) > smallDistance and calc_distance(pos2,maxpos) > smallDistance:
        return maxpos, False

    raise ValueError("Neither minimum nor maximum found in interval")

def findSaddlePoint(f, atol=2.e-8):
    posTop, minTop = findExtremum_1d(Point2D(Rmin, Zmax), Point2D(Rmax, Zmax), f)
    posBottom, minBottom = findExtremum_1d(Point2D(Rmin, Zmin), Point2D(Rmax, Zmin), f)
    posLeft, minLeft = findExtremum_1d(Point2D(Rmin, .8*Zmin), Point2D(Rmin, .8*Zmax), f)
    posRight, minRight = findExtremum_1d(Point2D(Rmax, .8*Zmin), Point2D(Rmax, .8*Zmax), f)

    assert minTop == minBottom
    assert minLeft == minRight
    assert minTop != minLeft

    if minTop:
        vertSearch = findMaximum_1d
    else:
        vertSearch = findMinimum_1d

    if minLeft:
        horizSearch = findMaximum_1d
    else:
        horizSearch = findMinimum_1d

    extremumVert = Point2D(Rmin, Zmin)
    extremumHoriz = Point2D(Rmax, Zmax)

    count = 0
    while calc_distance(extremumVert, extremumHoriz) > atol:
        count = count+1

        extremumVert = vertSearch(posBottom, posTop, f, 0.5*atol)
        posLeft.Z = extremumVert.Z
        posRight.Z = extremumVert.Z

        extremumHoriz = horizSearch(posLeft, posRight, f, 0.5*atol)
        posBottom.R = extremumHoriz.R
        posTop.R = extremumHoriz.R

    print('findSaddlePoint took',count,'iterations to converge')

    return (extremumVert+extremumHoriz)/2.

def findRoots_1d(f, n, xmin, xmax, atol = 2.e-8, rtol = 1.e-5, maxintervals=1024):
    """
    Find n roots of a scalar function f(x) in the range xmin<=x<=xmax
    Assume they're not too close to each other - exclude a small region around each found
    root when searching for more.
    """
    smallDistance = rtol * (xmax - xmin)
    foundRoots = 0
    roots = []
    n_intervals = n
    while True:
        interval_points = numpy.linspace(xmin, xmax, n_intervals+1)
        interval_f = f(interval_points)
        lucky_roots = numpy.where(interval_f == 0.)
        if len(lucky_roots[0]) > 0:
            raise NotImplementedError("Don't handle interval points that happen to land "
                    "on a root yet!")
        intervals_with_roots = numpy.where(numpy.sign(interval_f[:-1]) !=
                                           numpy.sign(interval_f[1:]))[0]
        if len(intervals_with_roots) >= n:
            break
        n_intervals *= 2
        if n_intervals > maxintervals:
            raise ValueError("Could not find", n, "roots when checking", maxintervals,
                             "intervals")

    # find roots in the intervals
    for i in intervals_with_roots:
        root, info = brentq(f, interval_points[i], interval_points[i+1], xtol=atol,
                full_output=True)
        if not info.converged:
            raise ValueError("Root finding failed in {" + str(interval_points[i]) + "," +
                    str(interval_points[i+1]) + "} with end values {" + str(interval_f[i])
                    + "," + str(interval_f[i+1]))
        roots.append(root)
        foundRoots += 1

    return roots

def findSeparatrix(A_toroidal, xpoint, A_x, atol = 2.e-8, npoints=100):
    """
    Follow 4 legs away from the x-point, starting with a rough guess and then refining to
    the separatrix value of A_toroidal.
    """
    boundaryThetas = findRoots_1d(lambda theta: A_toroidal(*TORPEX_wall(theta)) - A_x, 4,
            0., 2.*numpy.pi)

    # put lower left leg first in list, go clockwise
    boundaryThetas = boundaryThetas[2::-1] + [boundaryThetas[3]]

    boundaryPoints = tuple(TORPEX_wall(theta) for theta in boundaryThetas)

    legs = []
    s = numpy.linspace(10.*atol, 1., npoints)
    for point in boundaryPoints:
        legR = xpoint.R + s*(point.R - xpoint.R)
        legZ = xpoint.Z + s*(point.Z - xpoint.Z)
        leg = MeshContour([Point2D(R,Z) for R,Z in zip(legR, legZ)], A_toroidal, A_x)
        leg = leg.getRefined(atol=atol)
        legs.append(leg)

    return legs

def followPerpendicular(f_R, f_Z, p0, A0, Avals, rtol=2.e-8, atol=1.e-8):
    """
    Follow a line perpendicular to Bp from point p0 until magnetic potential A_target is
    reached.
    """
    f = lambda A,x: (f_R(x[0], x[1]), f_Z(x[0], x[1]))
    Arange = (A0, Avals[-1])
    solution = solve_ivp(f, Arange, tuple(p0), t_eval=Avals, rtol=rtol, atol=atol,
            vectorized=True)

    return [Point2D(*p) for p in solution.y.T]

def createMesh(filename):
    # parse input file
    coils, Bt_axis, meshOptions = parseInput(filename)

    A_toroidal, f_R, f_Z, Bp_R, Bp_Z = magneticFunctions(coils)

    xpoint = findSaddlePoint(A_toroidal)
    A_xpoint = A_toroidal(*xpoint)
    print('X-point',xpoint,'with A_toroidal='+str(A_xpoint))

    # note legs are ordered in theta
    separatrixLegs = findSeparatrix(A_toroidal, xpoint, A_xpoint)

    fpol = Bt_axis / 1. # Major radius of TORPEX axis is 1m

    return Mesh(meshOptions, A_toroidal, f_R, f_Z, Bp_R, Bp_Z, fpol, xpoint, A_xpoint, separatrixLegs)


if __name__ == '__main__':
    from sys import argv, exit

    filename = argv[1]
    gridname = 'torpex.grd.nc'

    mesh = createMesh(filename)

    if plotStuff:
        plotPotential(mesh.A_toroidal)
        #plotPotential(lambda R,Z: A_toroidal(R,Z)-A_xpoint)
        addWallToPlot()
        pyplot.plot(*mesh.xpoint, 'rx')
        for l in mesh.separatrixLegs:
            l.plot('1')
        for contours in mesh.contours_pf:
            for contour in contours:
                contour.plot('x')
        for contours in mesh.contours_sol:
            for contour in contours:
                contour.plot('+')
        pyplot.show()

    mesh.geometry()
    mesh.writeGridFile(gridname)

    exit(0)
