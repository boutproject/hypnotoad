"""
Classes to handle Meshes and geometrical quantities for generating BOUT++ grids
"""

import numpy
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq

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

def calc_distance(p1, p2):
    d = p2 - p1
    return numpy.sqrt(d.R**2 + d.Z**2)

class MeshContour:
    """
    Represents a contour as a collection of points.
    Includes methods for interpolation.
    Mostly behaves like a list
    """
    def __init__(self, points, psi, Aval):
        self.points = points

        self.distance = [0.]
        for i in range(1,len(self.points)):
            self.distance.append(
                    self.distance[-1] + calc_distance(self.points[i-1], self.points[i]))

        # Function that evaluates the vector potential at R,Z
        self.psi = psi

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

    def reverse(self):
        self.points.reverse()
        self.distance.reverse()
        self.distance = [self.distance[0] - d for d in self.distance]

    def refine(self, *args, **kwargs):
        self = self.getRefined(*args, **kwargs)

    def getRefined(self, width=.2, atol=2.e-8):
        f = lambda R,Z: self.psi(R, Z) - self.Aval

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
                    if numpy.abs(f(*pline(0.))) < atol*self.A_xpoint and numpy.abs(f(*pline(1.))) < atol*self.A_xpoint:
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

        return MeshContour(newpoints, self.psi, self.Aval)

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

    def getRegridded(self, npoints, width=1.e-5, atol=2.e-8, sfunc=None, extend_lower=0,
            extend_upper=0):
        """
        Interpolate onto set of npoints points, then refine positions.
        By default points are uniformly spaced, this can be changed by passing 'sfunc'
        which replaces the uniform interval 's' with 's=sfunc(s)'.
        'extend_lower' and 'extend_upper' extend the contour past its existing ends by a
        number of points.
        Returns a new MeshContour.
        """
        s = numpy.linspace(-extend_lower/(npoints-1),
                (npoints-1+extend_upper)/(npoints-1), npoints+extend_lower+extend_upper)
        if sfunc is not None:
            s = sfunc(s)
        interp = self.interpFunction()
        new_contour = MeshContour([interp(x) for x in s], self.psi, self.Aval)
        return new_contour.getRefined(width, atol)

    def plot(self, *args, **kwargs):
        from matplotlib import pyplot
        pyplot.plot([x.R for x in self], [x.Z for x in self], *args, **kwargs)

class MeshRegion:
    """
    A simple rectangular region of a Mesh, that connects to one other region (or has a
    boundary) on each edge.
    Note that these regions include cell face and boundary points, so there are
    (2nx+1)*(2ny+1) points for an nx*ny grid.
    """
    def __init__(self, meshParent, myID, localnx, localny, separatrix, psi_vals,
            connections, isInner):
        print('creating region ', myID)

        # MeshContour objects containing the points in this region
        self.contours = contours

        # the Mesh object that owns this MeshRegion
        self.meshParent = meshParent

        # ID that Mesh uses to keep track of its MeshRegions
        self.myID = myID

        # sizes of the grid in this MeshRegion
        self.nx = localnx
        self.ny = localny

        # psi values for radial grid
        self.psi_vals = psi_vals

        # Dictionary that specifies whether a boundary is connected to another region or
        # is an actual boundary
        self.connections = connections

        # y-boundary guard cells needed if the region edge is a real boundary, i.e. not
        # connected to another region
        if connections['lower'] is None:
            self.y_guards_lower = meshParent.y_boundary_guards
        else:
            self.y_guards_lower = 0
        if connections['upper'] is None:
            self.y_guards_upper = meshParent.y_boundary_guards
        else:
            self.y_guards_upper = 0

        # include guard cells in contours (note, does not change passed-in separatrix
        # object)
        self.separatrix = separatrix.getRegridded(2*self.ny,
                extend_lower=self.y_guards_lower, extend_upper=self.y_guards_upper)

        # get points in this region
        self.contours = []
        if isInner:
            temp_psi_vals = self.psi_vals[::-1]
        else:
            temp_psi_vals = self.psi_vals
        perp_points = followPerpendicular(meshParent.f_R, meshParent.f_Z,
                self.separatrix[0], meshParent.psi_sep, temp_psi_vals)
        if isInner:
            perp_points.reverse()
        for i,point in enumerate(perp_points_inner):
            self.contours.append(MeshContour([point], mesh_parent.psi,
                self.psi_vals[i]))
        for p in sep[1:]:
            perp_points = followPerpendicular(meshParent.f_R, meshParent.f_Z, p,
                    meshParent.psi_sep, temp_psi_vals)
            if isInner:
                perp_points.reverse()
            for i,point in enumerate(perp_points_inner):
                self.contours_pf[i].append(point)

        # refine the contours to make sure they are at exactly the right psi-value
        for contour in contours:
            contour.refine()

    def geometry(self):
        """
        Calculate geometrical quantities for this region
        """

        self.Rcorners = numpy.zeros([self.nx + 1,
            self.ny + self.y_guards_lower + self.y_guards_upper + 1])
        self.Zcorners = numpy.zeros([self.nx + 1,
            self.ny + self.y_guards_lower + self.y_guards_upper + 1])

        self.Rxy = numpy.array([[p.R for p in contour[1::2]]
            for contour in contours[1::2]])

        self.Rxy_ylow = numpy.array([[p.R for p in contour[0:-1:2]]
            for contour in contours[1::2]])

        self.Zxy = numpy.array( [[p.Z for p in contour[1::2]]
            for contour in contours[1::2]])

        self.Zxy_ylow = numpy.array( [[p.Z for p in contour[0:-1:2]]
            for contour in contours[1::2]])

        self.Rcorners = numpy.array( [[p.R for p in contour[0::2]]
            for contour in contours[0::2]])
        self.Zcorners = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in contours[0::2]])

        self.psixy = self.psi(self.Rxy, self.Zxy)
        self.psixy_ylow = self.psi(self.Rxy_ylow, self.Zxy_ylow)

        self.dx = numpy.zeros([self.nx, self.ny + 2*self.y_boundary_guards])
        self.dx[:] = numpy.array(self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]
        self.dx_ylow = numpy.zeros([self.nx, self.ny + 2*self.y_boundary_guards])
        self.dx_ylow[:] = numpy.array(self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]

        if self.psi_vals[0] > self.psi_vals[-1]:
            # x-coordinate is -psixy so x always increases radially across grid
            self.bpsign = -1.
            self.xcoord = -self.psixy
        else:
            self.bpsign = 1.
            self.xcoord = self.psixy

        self.dy = meshParent.dy * numpy.ones([self.nx, self.ny + self.y_guards_lower
                                              + self.y_guards_upper])
        self.dy_ylow = meshParent.dy * numpy.ones([self.nx, self.ny + self.y_guards_lower
                                                   + self.y_guards_upper])

        self.Brxy = meshParent.Bp_R(self.Rxy, self.Zxy)
        self.Brxy_ylow = meshParent.Bp_R(self.Rxy_ylow, self.Zxy_ylow)
        self.Bzxy = meshParent.Bp_Z(self.Rxy, self.Zxy)
        self.Bzxy_ylow = meshParent.Bp_Z(self.Rxy_ylow, self.Zxy_ylow)
        self.Bpxy = numpy.sqrt(self.Brxy**2 + self.Bzxy**2)
        self.Bpxy_ylow = numpy.sqrt(self.Brxy_ylow**2 + self.Bzxy_ylow**2)
        # determine direction - dot Bp with Grad(y) vector
        # evaluate in 'sol' at outer radial boundary
        Bp_dot_grady = (
            self.Brxy[-1, self.ny/2 + self.y_guards_lower]
            *(self.Rxy[-1, self.ny/2 + self.y_guards_lower + 1]
                - self.Rxy[-1, self.ny/2 + self.y_guards_lower - 1])
            + self.Bzxy[-1, self.ny/2 + self.y_guards_lower]
              *(self.Zxy[-1, self.ny/2 + self.y_guards_lower + 1]
                  - self.Zxy[-1, self.ny/2 + self.y_guards_lower - 1]) )
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
        self.Btxy = meshParent.fpol / self.Rxy
        self.Btxy_ylow = meshParent.fpol / self.Rxy_ylow

        self.Bxy = numpy.sqrt(self.Bpxy**2 + self.Btxy**2)
        self.Bxy_ylow = numpy.sqrt(self.Bpxy_ylow**2 + self.Btxy_ylow**2)

        self.hthe = numpy.sqrt((Rxy_ylow[:,1:] - Rxy_ylow[:,:-1])**2
                               + (Zxy_ylow[:,1:] - Zxy_ylow[:,:-1])**2)
        # for hthe_ylow, need R, Z values from below the lower face of this region
        R = numpy.zeros([self.nx, self.ny + self.y_guards_lower + self.y_guards_upper + 1])
        R[:, 1:] = self.Rxy
        Z = numpy.zeros([self.nx, self.ny + self.y_guards_lower + self.y_guards_upper + 1])
        Z[:, 1:] = self.Zxy
        if self.connections['lower'] is not None:
            R[:, 0] = self.getNeighbour('lower').Rxy[:, -1]
            Z[:, 0] = self.getNeighbour('lower').Zxy[:, -1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:, 0] = 2.*self.Rxy_ylow[:, 0] - self.Rxy[:, 0]
            Z[:, 0] = 2.*self.Zxy_ylow[:, 0] - self.Zxy[:, 0]
        self.hthe_ylow = numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)

    def getNeighbour(self, face):
        return self.MeshParent.regions[self.connections[face]]

    def DDX(self, f):
        raise ValueError('not implemented for MeshRegion yet')
        result = numpy.zeros([self.nx, self.ny + 2*self.y_boundary_guards])
        result[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.*self.dx[1:-1])
        result[0, :] = (-1.5*f[0, :] + 2.*f[1,:] - 0.5*f[2, :]) / self.dx[0]
        result[-1, :] = (1.5*f[-1, :] - 2.*f[-2,:] + 0.5*f[-3, :]) / self.dx[-1]
        return result

class Mesh:
    """
    Mesh quantities to be written to a grid file for BOUT++
    """
    def __init__(self, meshOptions, psi, f_R, f_Z, Bp_R, Bp_Z, fpol, psi_sep,
                 separatrixLegs):
        self.orthogonal = meshOptions['orthogonal']
        self.nx_core = meshOptions['nx_core']
        self.nx_between = meshOptions['nx_between']
        self.nx_sol = meshOptions['nx_sol']
        self.ny_inner_lower_divertor = meshOptions['ny_inner_lower_divertor']
        self.ny_inner_core = meshOptions['ny_inner_core']
        self.ny_inner_upper_divertor = meshOptions['ny_inner_upper_divertor']
        self.ny_outer_upper_divertor = meshOptions['ny_outer_upper_divertor']
        self.ny_outer_core = meshOptions['ny_outer_core']
        self.ny_outer_lower_divertor = meshOptions['ny_outer_lower_divertor']
        self.psi_inner = meshOptions['psi_inner']
        self.psi_outer = meshOptions['psi_outer']
        self.y_boundary_guards = meshOptions['y_boundary_guards']

        self.psi = psi
        self.f_R = f_R
        self.f_Z = f_Z
        self.Bp_R = Bp_R
        self.Bp_Z = Bp_Z

        # poloidal current function, gives B_toroidal
        self.fpol = fpol

        self.psi_sep = psi_sep

        if self.nx_between > 0:
            raise ValueError("nx_between > 0 - there are 2 separatrices - need to find psi-value of second separatrix")
        assert self.orthogonal # non-orthogonal not implelemented yet

        # nx, ny both include boundary guard cells
        self.nx = self.nx_core + self.nx_between + self.nx_sol

        if ( (self.ny_inner_upper_divertor > 0 or self.ny_outer_upper_divertor > 0)
                and (self.ny_inner_lower_divertor > 0 or self.ny_inner_core > 0
                     or self.ny_outer_core > 0 or self.ny_outer_lower_divertor > 0) ):
            # there are two distinct divertor/limiter targets
            #  - the upper divertor legs are not empty, and also the other regions are not
            #    all empty
            self.upper_target_y_boundary_guards = self.y_boundary_guards

        self.ny = (self.ny_inner_lower_divertor + self.ny_inner_core
                   + self.ny_inner_upper_divertor + self.ny_outer_upper_divertor
                   + self.ny_outer_core + self.ny_outer_lower_divertor
                   + 2*self.y_boundary_guards + 2*self.upper_target_y_boundary_guards)

        # in index space for indices of cell faces, psi needs to go through psi_inner at
        # 0, psi_sep at nx_core and psi_outer at nx_core+nx_between+nx_sol+1
        # for now use quadratic fit, leave grid refinement for later...
        # psi(i) = c + b*i + a*i^2
        # psi(0) = c = psi_inner
        # psi(ixseps) = psi_sep
        #   psi_sep - psi_inner = ixseps(b + a*ixseps)
        #   (psi_sep - psi_inner)/ixseps = b + a*ixseps
        # psi(nx+1) = psi_outer
        #   psi_outer - psi_inner = (nx+1)*(b + a*(nx+1))
        #   (psi_outer - psi_inner)/(nx+1) = b + a*(nx+1)
        # a*(nx+1-ixseps) = (psi_outer-psi_inner)/(nx+1) - (psi_sep-psi_inner)/ixseps
        # a = ((psi_outer-psi_inner)/(nx+1) - (psi_sep-psi_inner)/ixseps) / (nx+1-ixseps)
        # b = ( (psi_outer-psi_inner)/(nx+1)**2 - (psi_sep-psi_inner)/ixseps**2 ) / (1/(nx+1) - 1/ixseps)
        nx = self.nx_core + self.nx_between + self.nx_sol
        a = ((self.psi_outer-self.psi_inner)/(nx+1) -
                (self.psi_sep-self.psi_inner)/self.nx_core) / (nx+1-self.nx_core)
        b = ( (self.psi_outer-self.psi_inner)/(nx+1)**2 -
                (self.psi_sep-self.psi_inner)/self.nx_core**2 ) / (1./(nx+1) -
                        1./self.nx_core)
        c = self.psi_inner
        psi_index = lambda i: a*i**2 + b*i + c
        psi_face_vals_inner = numpy.array(psi_index(i) for i in range(self.nx_core+1))
        psi_face_vals_between = numpy.array(psi_index(i) for i in
                range(self.nx_core+1, self.nx_core+self.nx_between+1))
        psi_face_vals_outer = numpy.array(psi_index(i) for i in
                range(self.nx_core+self.nx_between+1, self.nx_core+self.nx_between+self.nx_sol+2))
        self.psi_vals_inner = numpy.zeros(2*self.nx_core+1)
        self.psi_vals_between = numpy.zeros(2*self.nx_between+1)
        self.psi_vals_outer = numpy.zeros(2*self.nx_sol+1)
        self.psi_vals_inner[0::2] = psi_face_vals_inner
        self.psi_vals_inner[1::2] = 0.5*(psi_face_vals_inner[:-1] + psi_face_vals_inner[1:])
        self.psi_vals_between[0::2] = psi_face_vals_between
        self.psi_vals_between[1::2] = 0.5*(psi_face_vals_between[:-1] + psi_face_vals_between[1:])
        self.psi_vals_outer[0::2] = psi_face_vals_outer
        self.psi_vals_outer[1::2] = 0.5*(psi_face_vals_outer[:-1] + psi_face_vals_outer[1:])

        # Generate MeshRegion object for each section of the mesh
        # For region numbers see figure in 'BOUT++ Topology' section of BOUT++ manual
        self.regions = {}

        # Keep ranges of global indices for each region separately, because we don't want
        # MeshRegion objects to depend on global indices
        self.region_indices = {}
        x_regions = (slice(None, self.nx_core, None),
                     slice(self.nx_core, self.nx_core + self.nx_between, None),
                     slice(self.nx_core + self.nx_between, None, None))
        y_sizes = [0,
                   self.ny_inner_lower_divertor + self.y_boundary_guards,
                   self.ny_inner_core,
                   self.ny_inner_upper_divertor + self.upper_target_y_boundary_guards,
                   self.ny_outer_upper_divertor + self.upper_target_y_boundary_guards,
                   self.ny_outer_core,
                   self.ny_outer_lower_divertor + self.y_boundary_guards]
        y_startinds = numpy.cumsum(y_sizes)
        y_regions = (slice(y_startinds[i], y_startinds[i+1], None)
                     for i in range(len(y_startinds-1)))

        # Region 1 - inner lower PF
        # Region 2 - inner lower between separatrices
        # Region 3 - inner lower SOL
        if self.ny_lower_inner_divertor > 0:
            sep = separatrix['inner_lower'].getRegridded(2*self.ny_inner_lower_divertor+1,
                                                         sfunc=sfunc)
            sep.reverse()
            if self.nx_core > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_between > 0:
                    connections['outer'] = 2 # inner lower between separatrix region
                elif self.nx_sol > 0:
                    connections['outer'] = 3 # inner lower SOL
                else:
                    connections['outer'] = None
                connections['lower'] = None
                if self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 14 # outer lower PF
                else:
                    connections['upper'] = None
                self.regions[1] = MeshRegion(self, 1, self.nx_core,
                        self.ny_inner_lower_divertor, sep, self.psi_vals_inner,
                        connections, True)
                self.region_indices[1] = (numpy.index_exp[x_regions[0], y_regions[0]])
            if self.nx_between > 0:
                connections = {}
                if self.nx_core > 0:
                    connections['inner'] = 1 # inner lower PF
                else:
                    connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 3 # inner lower SOL
                else:
                    connections['outer'] = None
                connections['lower'] = None
                if self.ny_inner_core > 0:
                    connections['upper'] = 5 # inner between separatrix region
                elif self.ny_outer_core > 0:
                    connections['upper'] = 12 # outer between separatrix region
                elif self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 15 # outer lower between separatrix region
                else:
                    connections['upper'] = None
                self.regions[2] = MeshRegion(self, 2, self.nx_between,
                        self.ny_inner_lower_divertor, sep, self.psi_vals_between,
                        connections, False)
                self.region_indices[2] = (numpy.index_exp[x_regions[1], y_regions[0]])
            if self.nx_sol > 0:
                connections = {}
                if self.nx_between > 0:
                    connections['inner'] = 2 # inner lower between separatrix region
                elif self.nx_core > 0:
                    connections['inner'] = 1 # inner lower PF
                else:
                    connections['inner'] = None
                connections['outer'] = None
                connections['lower'] = None
                if self.ny_inner_core > 0:
                    connections['upper'] = 6 # inner SOL
                elif self.ny_inner_upper_divertor > 0:
                    connections['upper'] = 7 # inner upper SOL
                elif self.ny_outer_upper_divertor > 0:
                    # this probably shouldn't happen, but if there is no upper inner leg,
                    # but there is an upper outer leg, we can't connect to the outer SOL,
                    # so there must be a limiter/divertor plate going right up to the
                    # upper X-point
                    connections['upper'] = None
                elif self.ny_outer_core > 0:
                    connections['upper'] = 13 # outer SOL
                elif self.ny_lower_outer_divertor > 0:
                    connections['upper'] = 16 # outer lower SOL
                else:
                    connections['upper'] = None
                self.regions[3] = MeshRegion(self, 3, self.nx_sol,
                        self.ny_inner_lower_divertor, sep, self.psi_vals_outer,
                        connections, False)
                self.region_indices[3] = (numpy.index_exp[x_regions[2], y_regions[0]])

        # Region 4 - inner core
        # Region 5 - inner between separatrices
        # Region 6 - inner SOL
        if self.ny_inner_core > 0:
            sep = separatrix['inner_core'].getRegridded(2*self.ny_inner_core+1, sfunc=sfunc)
            if self.nx_core > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_between > 0:
                    connections['outer'] = 5 # inner between separatrix region
                elif self.nx_sol > 0:
                    connections['outer'] = 6 # inner SOL
                else:
                    connections['outer'] = None
                if self.ny_outer_core > 0:
                    connections['lower'] = 11 # outer core
                    connections['upper'] = 11 # outer core
                else:
                    connections['lower'] = 4 # periodic inner core
                    connections['upper'] = 4 # periodic inner core
                self.regions[4] = MeshRegion(self, 4, self.nx_core,
                        self.ny_inner_core, sep, self.psi_vals_inner, connections, True)
                self.region_indices[4] = (numpy.index_exp[x_regions[0], y_regions[1]])
            if self.nx_between > 0:
                connections = {}
                if self.nx_core > 0:
                    connections['inner'] = 4 # inner core
                else:
                    connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 6 # inner SOL
                else:
                    connections['outer'] = None
                if self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 2 # inner lower between separatrix region
                else:
                    connections['lower'] = None
                if self.ny_outer_core > 0:
                    connections['upper'] = 12 # outer between separatrix region
                elif self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 15 # outer lower between separatrix region
                else:
                    connections['upper'] = None
                self.regions[5] = MeshRegion(self, 5, self.nx_between,
                        self.ny_inner_core, sep, self.psi_vals_between, connections, False)
                self.region_indices[5] = (numpy.index_exp[x_regions[1], y_regions[1]])
            if self.nx_sol > 0:
                connections = {}
                if self.nx_between > 0:
                    connections['inner'] = 5 # inner between separatrix region
                elif self.nx_core > 0:
                    connections['inner'] = 4 # inner core
                else:
                    connections['inner'] = None
                connections['outer'] = None
                if self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 3 # inner lower SOL
                else:
                    connections['lower'] = None
                if self.ny_inner_upper_divertor > 0:
                    connections['upper'] = 7 # inner upper SOL
                elif self.ny_outer_upper_divertor > 0:
                    # this probably shouldn't happen, but if there is no upper inner leg,
                    # but there is an upper outer leg, we can't connect to the outer SOL,
                    # so there must be a limiter/divertor plate going right up to the
                    # upper X-point
                    connections['upper'] = None
                elif self.ny_outer_core > 0:
                    connections['upper'] = 13 # outer SOL
                elif self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 16 # outer lower SOL
                else:
                    connections['upper'] = None
                self.regions[6] = MeshRegion(self, 6, self.nx_sol,
                        self.ny_inner_core, sep, self.psi_vals_outer, connections, False)
                self.region_indices[6] = (numpy.index_exp[x_regions[2], y_regions[1]])

        # Region 7 - inner upper SOL
        # Region 8 - inner upper PF
        nx_upper_pf = self.nx_between + self.nx_core
        if self.ny_inner_upper_divertor > 0:
            sep = separatrix['inner_upper_divertor'].getRegridded(2*self.ny_inner_upper_divertor+1,
                                                                  sfunc=sfunc)
            if nx_upper_pf > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 7 # inner upper SOL
                else:
                    connections['outer'] = None
                if self.ny_outer_upper_divertor > 0:
                    connections['lower'] = 9 # outer upper PF
                else:
                    connections['lower'] = None
                connections['upper'] = None
                self.regions[7] = MeshRegion(self, 7, nx_upper_pf,
                        self.ny_inner_upper_divertor, sep, self.psi_vals_inner,
                        connections, True)
                self.region_indices[7] = (numpy.index_exp[:nx_upper_pf, y_regions[2]])
            if self.nx_sol > 0:
                connections = {}
                if nx_upper_pf > 0:
                    connections['inner'] = 7 # outer upper PF
                else:
                    connections['inner'] = None
                connections['outer'] = None
                if self.ny_inner_core > 0:
                    connections['lower'] = 6 # inner SOL
                elif self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 3 # inner lower SOL
                else:
                    connections['lower'] = None
                connections['upper'] = None
                self.regions[8] = MeshRegion(self, 8, self.nx_sol,
                        self.ny_inner_upper_divertor, sep, self.psi_vals_outer,
                        connections, False)
                self.region_indices[8] = (numpy.index_exp[nx_upper_pf:, y_regions[2]])

        # Region 9 - outer upper PF
        # Region 10 - outer upper SOL
        if self.ny_outer_upper_divertor > 0:
            sep = separatrix['outer_upper_divertor'].getRegridded(2*self.ny_inner_upper_divertor+1,
                                                                  sfunc=sfunc)
            sep.reverse()
            if nx_upper_pf > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 10 # outer upper SOL
                else:
                    connections['outer'] = None
                connections['lower'] = None
                if self.ny_outer_upper_divertor > 0:
                    connections['upper'] = 8 # inner upper PF
                else:
                    connections['upper'] = None
                self.regions[9] = MeshRegion(self, 9, nx_upper_pf,
                        self.ny_outer_upper_divertor, sep, self.psi_vals_inner,
                        connections, True)
                self.region_indices[9] = (numpy.index_exp[:nx_upper_pf, y_regions[3]])
            if self.nx_sol > 0:
                connections = {}
                connections['lower'] = None
                if self.ny_outer_core > 0:
                    connections['upper'] = 13 # outer SOL
                elif self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 16 # outer lower SOL
                else:
                    connections['upper'] = None
                if nx_upper_pf > 0:
                    connections['inner'] = 9 # outer upper PF
                else:
                    connections['inner'] = None
                connections['outer'] = None
                self.regions[10] = MeshRegion(self, 10, self.nx_sol,
                        self.ny_outer_upper_divertor, sep, self.psi_vals_outer,
                        connections, False)
                self.region_indices[10] = (numpy.index_exp[nx_upper_pf:, y_regions[3]])

        # Region 11 - outer core
        # Region 12 - outer between separatrices
        # Region 13 - outer SOL
        if self.ny_outer_core > 0:
            sep = separatrix['outer_core'].getRegridded(2*self.ny_outer_core+1,
                                                        sfunc=sfunc)
            if self.nx_core > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_between > 0:
                    connections['outer'] = 12 # outer between separatrix region
                elif self.nx_sol > 0:
                    connections['outer'] = 13 # outer SOL region
                else:
                    connections['outer'] = None
                if self.ny_inner_core > 0:
                    connections['lower'] = 4 # inner core
                    connections['upper'] = 4 # inner core
                else:
                    connections['lower'] = 11 # outer core
                    connections['upper'] = 11 # outer core
                self.regions[11] = MeshRegion(self, 11, self.nx_core,
                        self.ny_outer_core, sep, self.psi_vals_inner, connections, True)
                self.region_indices[11] = (numpy.index_exp[x_regions[0], y_regions[4]])
            if self.nx_between > 0:
                connections = {}
                if self.nx_core > 0:
                    connections['inner'] = 11 # outer core
                else:
                    connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 13 # outer SOL region
                else:
                    connections['outer'] = None
                if self.ny_inner_core > 0:
                    connections['lower'] = 5 # inner between separatrix region
                elif self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 2 # inner lower between separatrix region
                else:
                    connections['lower'] = None
                if self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 15 # outer lower between separatrix region
                else:
                    connections['upper'] = None
                self.regions[12] = MeshRegion(self, 12, self.nx_between,
                        self.ny_outer_core, sep, self.psi_vals_between, connections, False)
                self.region_indices[12] = (numpy.index_exp[x_regions[1], y_regions[4]])
            if self.nx_sol > 0:
                connections = {}
                if self.nx_between > 0:
                    connections['inner'] = 12 # outer between separatrix region
                elif self.nx_core > 0:
                    connections['inner'] = 11 # outer core
                else:
                    connections['inner'] = None
                connections['outer'] = None
                if self.ny_outer_upper_divertor > 0:
                    connections['lower'] = 10
                elif self.ny_inner_upper_divertor > 0:
                    # this probably shouldn't happen, but if there is no upper outer leg,
                    # but there is an upper inner leg, we can't connect to the inner SOL,
                    # so there must be a limiter/divertor plate going right up to the
                    # upper X-point
                    connections['lower'] = None
                elif self.ny_inner_core > 0:
                    connections['lower'] = 6 # inner SOL
                elif self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 3 # inner lower SOL
                else:
                    connections['lower'] = None
                if self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 16 # outer lower SOL
                else:
                    connections['upper'] = None
                self.regions[13] = MeshRegion(self, 13, self.nx_sol,
                        self.ny_outer_core, sep, self.psi_vals_outer, connections, False)
                self.region_indices[13] = (numpy.index_exp[x_regions[2], y_regions[4]])

        # Region 14 - outer lower PF
        # Region 15 - outer lower between separatrices
        # Region 16 - outer lower SOL
        if self.ny_outer_lower_divertor > 0:
            sep = separatrix['outer_lower'].getRegridded(2*self.ny_outer_lower_divertor+1,
                                                         sfunc=sfunc)
            if self.nx_core > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_between > 0:
                    connections['outer'] = 15 # outer lower between separatrix region
                elif self.nx_sol > 0:
                    connections['outer'] = 16 # outer lower SOL
                else:
                    connections['outer'] = None
                if self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 1 # inner lower PF
                else:
                    connections['lower'] = None
                self.regions[14] = MeshRegion(self, 14, self.nx_core,
                        self.ny_outer_lower_divertor, sep, self.psi_vals_inner,
                        connections, True)
                self.region_indices[14] = (numpy.index_exp[x_regions[0], y_regions[5]])
            if self.nx_between > 0:
                connections = {}
                if self.nx_core > 0:
                    connections['inner'] = 14 # outer lower PF
                else:
                    connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 16 # outer lower SOL
                else:
                    connections['outer'] = None
                if self.ny_outer_core > 0:
                    connections['lower'] = 12 # outer between separatrix region
                elif self.ny_inner_core > 0:
                    connections['lower'] = 5 # inner between separatrix region
                elif self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 2 # inner lower between separatrix region
                else:
                    connections['lower'] = None
                self.regions[15] = MeshRegion(self, 15, self.nx_between,
                        self.ny_outer_lower_divertor, sep, self.psi_vals_between,
                        connections, False)
                self.region_indices[15] = (numpy.index_exp[x_regions[1], y_regions[5]])
            if self.nx_sol > 0:
                connections = {}
                if self.nx_between > 0:
                    connections['inner'] = 15 # outer lower between separatrix region
                elif self.nx_core > 0:
                    connections['inner'] = 14 # outer lower PF
                else:
                    connections['inner'] = None
                connections['outer'] = None
                if self.ny_outer_core > 0:
                    connections['lower'] = 13 # outer SOL
                elif self.ny_outer_upper_divertor > 0:
                    connections['lower'] = 10 # outer upper SOL
                elif self.ny_inner_upper_divertor > 0:
                    # this probably shouldn't happen, but if there is no upper outer leg,
                    # but there is an upper inner leg, we can't connect to the inner SOL,
                    # so there must be a limiter/divertor plate going right up to the
                    # upper X-point
                    connections['lower'] = None
                elif self.ny_inner_core > 0:
                    connections['lower'] = 6 # inner SOL
                elif self.ny_inner_lower_divertor > 0:
                    connections['lower'] = 3 # inner lower SOL
                else:
                    connections['lower'] = None
                connections['upper'] = None
                self.regions[16] = MeshRegion(self, 16, self.nx_sol,
                        self.ny_outer_lower_divertor, sep, self.psi_vals_outer,
                        connections, False)
                self.region_indices[16] = (numpy.index_exp[x_regions[2], y_regions[5]])

        # constant spacing in y for now
        self.dy = 2.*numpy.pi / (self.ny_inner_lower_divertor + self.ny_inner_core +
                self.ny_inner_upper_divertor + self.ny_outer_upper_divertor +
                self.ny_outer_core + self.ny_outer_lower_divertor)

    def geometry(self):
        """
        Calculate geometrical quantities for BOUT++
        """

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
        from matplotlib import pyplot
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
