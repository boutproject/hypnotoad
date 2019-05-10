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
    def __init__(self, points, psi, psival):
        self.points = points

        self.distance = [0.]
        for i in range(1,len(self.points)):
            self.distance.append(
                    self.distance[-1] + calc_distance(self.points[i-1], self.points[i]))

        # Function that evaluates the vector potential at R,Z
        self.psi = psi

        # Value of vector potential on this contour
        self.psival = psival

    def __iter__(self):
        return self.points.__iter__()

    def __str__(self):
        return self.points.__str__()

    def __getitem__(self, key):
        return self.points.__getitem__(key)

    def __len__(self):
        return self.points.__len__()

    def append(self, point):
        self.points.append(point)
        self.distance.append(
                self.distance[-1] + calc_distance(self.points[-2], self.points[-1]))

    def reverse(self):
        self.points.reverse()
        self.distance.reverse()
        self.distance = [self.distance[0] - d for d in self.distance]

    def refine(self, *args, **kwargs):
        new = self.getRefined(*args, **kwargs)
        self.points = new.points
        self.distance = new.distance

    def getRefined(self, width=1.e-5, atol=2.e-8):
        f = lambda R,Z: self.psi(R, Z) - self.psival

        def perpLine(p, tangent, w):
            # p - point through which to draw perpLine
            # tangent - vector tangent to original curve, result will be perpendicular to this
            # w - width on either side of p to draw the perpLine to
            modTangent = numpy.sqrt(tangent.R**2 + tangent.Z**2)
            perpIdentityVector = Point2D(tangent.Z/modTangent, -tangent.R/modTangent)
            return lambda s: p + 2.*(s-0.5)*w*perpIdentityVector

        def refinePoint(p, tangent):
            if numpy.abs(f(*p)) < atol*numpy.abs(self.psival):
                # don't need to refine
                return p
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
                    raise ValueError("Could not find interval to refine point")

            return pline(snew)

        newpoints = []
        newpoints.append(refinePoint(self.points[0], self.points[1] - self.points[0]))
        for i,p in enumerate(self.points[1:-1]):
            newpoints.append(refinePoint(p, self.points[i+1] - self.points[i-1]))
        newpoints.append(refinePoint(self.points[-1], self.points[-1] - self.points[-2]))

        return MeshContour(newpoints, self.psi, self.psival)

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
        new_contour = MeshContour([interp(x) for x in s], self.psi, self.psival)
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
    def __init__(self, meshParent, myID, localNx, localNy, separatrix, psi_vals,
            connections, isInner):
        print('creating region ', myID)

        # the Mesh object that owns this MeshRegion
        self.meshParent = meshParent

        # ID that Mesh uses to keep track of its MeshRegions
        self.myID = myID

        # sizes of the grid in this MeshRegion, include boundary guard cells
        self.nx = localNx
        self.ny = localNy

        # psi values for radial grid
        self.psi_vals = psi_vals

        # MeshContour representing the separatrix segment associated with this region
        self.separatrix = separatrix

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

        # get points in this region
        self.contours = []
        if isInner:
            temp_psi_vals = self.psi_vals[::-1]
        else:
            temp_psi_vals = self.psi_vals
        perp_points = followPerpendicular(meshParent.equilibrium.f_R,
                meshParent.equilibrium.f_Z, self.separatrix[0],
                meshParent.equilibrium.psi_sep[0], temp_psi_vals)
        if isInner:
            perp_points.reverse()
        for i,point in enumerate(perp_points):
            self.contours.append(MeshContour([point], meshParent.equilibrium.psi,
                self.psi_vals[i]))
        for p in self.separatrix[1:]:
            perp_points = followPerpendicular(meshParent.equilibrium.f_R,
                    meshParent.equilibrium.f_Z, p, meshParent.equilibrium.psi_sep[0],
                    temp_psi_vals)
            if isInner:
                perp_points.reverse()
            for i,point in enumerate(perp_points):
                self.contours[i].append(point)

        # refine the contours to make sure they are at exactly the right psi-value
        for contour in self.contours:
            contour.refine()

    def fillRZ(self):
        """
        Fill the Rxy, Rxy_ylow and Zxy, Zxy_ylow arrays for this region

        xlow values include the outer point, after the final cell-centre grid point
        ylow values include the upper point, above the final cell-centre grid point
        """

        self.Rcorners = numpy.zeros([self.nx + 1, self.ny + 1])
        self.Zcorners = numpy.zeros([self.nx + 1, self.ny + 1])

        self.Rxy = numpy.array([[p.R for p in contour[1::2]]
            for contour in self.contours[1::2]])

        self.Rxy_ylow = numpy.array([[p.R for p in contour[0::2]]
            for contour in self.contours[1::2]])

        self.Rxy_extra_upper = numpy.array([contour[-1].R
            for contour in self.contours[1::2]])

        self.Rxy_xlow = numpy.array([[p.R for p in contour[1::2]]
            for contour in self.contours[0::2]])

        self.Rxy_extra_outer = numpy.array([p.R for p in self.contours[-1][1::2]])

        self.Zxy = numpy.array( [[p.Z for p in contour[1::2]]
            for contour in self.contours[1::2]])

        self.Zxy_ylow = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in self.contours[1::2]])

        self.Zxy_extra_upper = numpy.array([contour[-1].Z
            for contour in self.contours[1::2]])

        self.Zxy_xlow = numpy.array([[p.Z for p in contour[1::2]]
            for contour in self.contours[0::2]])

        self.Zxy_extra_outer = numpy.array([p.Z for p in self.contours[-1][1::2]])

        self.Rcorners = numpy.array( [[p.R for p in contour[0::2]]
            for contour in self.contours[0::2]])
        self.Zcorners = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in self.contours[0::2]])

    def geometry(self):
        """
        Calculate geometrical quantities for this region
        """

        self.psixy = self.meshParent.equilibrium.psi(self.Rxy, self.Zxy)
        self.psixy_ylow = self.meshParent.equilibrium.psi(self.Rxy_ylow, self.Zxy_ylow)
        self.psixy_xlow = self.meshParent.equilibrium.psi(self.Rxy_xlow, self.Zxy_xlow)
        self.psicorners = self.meshParent.equilibrium.psi(self.Rcorners, self.Zcorners)

        self.dx = numpy.zeros([self.nx, self.ny])
        self.dx[:] = numpy.array(self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]
        self.dx_ylow = numpy.zeros([self.nx, self.ny+1])
        self.dx_ylow[:] = numpy.array(self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]

        if self.psi_vals[0] > self.psi_vals[-1]:
            # x-coordinate is -psixy so x always increases radially across grid
            self.bpsign = -1.
            self.xcoord = -self.psixy
        else:
            self.bpsign = 1.
            self.xcoord = self.psixy

        self.dy = self.meshParent.dy_scalar * numpy.ones([self.nx, self.ny])
        self.dy_ylow = self.meshParent.dy_scalar * numpy.ones([self.nx, self.ny+1])
        self.dy_xlow = self.meshParent.dy_scalar * numpy.ones([self.nx+1, self.ny])
        self.dycorners = self.meshParent.dy_scalar * numpy.ones([self.nx+1, self.ny+1])

        self.Brxy = self.meshParent.equilibrium.Bp_R(self.Rxy, self.Zxy)
        self.Brxy_ylow = self.meshParent.equilibrium.Bp_R(self.Rxy_ylow, self.Zxy_ylow)
        self.Brxy_xlow = self.meshParent.equilibrium.Bp_R(self.Rxy_xlow, self.Zxy_xlow)
        self.Brcorners = self.meshParent.equilibrium.Bp_R(self.Rcorners, self.Zcorners)
        self.Bzxy = self.meshParent.equilibrium.Bp_Z(self.Rxy, self.Zxy)
        self.Bzxy_ylow = self.meshParent.equilibrium.Bp_Z(self.Rxy_ylow, self.Zxy_ylow)
        self.Bzxy_xlow = self.meshParent.equilibrium.Bp_Z(self.Rxy_xlow, self.Zxy_xlow)
        self.Bzcorners = self.meshParent.equilibrium.Bp_Z(self.Rcorners, self.Zcorners)
        self.Bpxy = numpy.sqrt(self.Brxy**2 + self.Bzxy**2)
        self.Bpxy_ylow = numpy.sqrt(self.Brxy_ylow**2 + self.Bzxy_ylow**2)
        self.Bpxy_xlow = numpy.sqrt(self.Brxy_xlow**2 + self.Bzxy_xlow**2)
        self.Bpcorners = numpy.sqrt(self.Brcorners**2 + self.Bzcorners**2)
        # determine direction - dot Bp with Grad(y) vector
        # evaluate in 'sol' at outer radial boundary
        Bp_dot_grady = (
            self.Brxy[-1, self.ny//2]
            *(self.Rxy[-1, self.ny//2 + 1] - self.Rxy[-1, self.ny//2 - 1])
            + self.Bzxy[-1, self.ny//2]
              *(self.Zxy[-1, self.ny//2 + 1] - self.Zxy[-1, self.ny//2 - 1]) )
        if Bp_dot_grady < 0.:
            print("Poloidal field is in opposite direction to Grad(theta) -> Bp negative")
            self.Bpxy = -self.Bpxy
            self.Bpxy_ylow = -self.Bpxy_ylow
            self.Bpxy_xlow = -self.Bpxy_xlow
            self.Bpcorners = -self.Bpcorners
            if self.bpsign > 0.:
                raise ValueError("Sign of Bp should be negative?")
        else:
            if self.bpsign < 0.:
                raise ValueError("Sign of Bp should be positive?")

        # Get toroidal field from poloidal current function fpol
        self.Btxy = self.meshParent.equilibrium.fpol(self.psixy) / self.Rxy
        self.Btxy_ylow = self.meshParent.equilibrium.fpol(self.psixy_ylow) / self.Rxy_ylow
        self.Btxy_xlow = self.meshParent.equilibrium.fpol(self.psixy_xlow) / self.Rxy_xlow
        self.Btcorners = self.meshParent.equilibrium.fpol(self.psicorners) / self.Rcorners

        self.Bxy = numpy.sqrt(self.Bpxy**2 + self.Btxy**2)
        self.Bxy_ylow = numpy.sqrt(self.Bpxy_ylow**2 + self.Btxy_ylow**2)

        self.hthe, self.hthe_ylow, self.hthe_xlow, self.hthecorners = self.calcHthe()

        #if not self.meshParent.orthogonal:
        #    # Calculate beta (angle between x and y coordinates), used for non-orthogonal grid
        #    # Also calculate radial grid spacing
        #    self.beta, self.hrad = self.calcBeta()
        #    self.beta_ylow, self.hrad_ylow = self.calcBeta(ylow=True)

        #    # eta is the polodial non-orthogonality parameter
        #    self.eta = numpy.sin(self.beta)
        #    self.eta_ylow = numpy.sin(self.beta_ylow)
        #else:
        #    self.beta = 0.
        #    self.eta = 0.

        # field line pitch
        self.pitch = self.hthe * self.Btxy / (self.Bpxy * self.Rxy)
        self.pitch_ylow = self.hthe_ylow * self.Btxy_ylow / (self.Bpxy_ylow
                                                             * self.Rxy_ylow)
        self.pitch_xlow = self.hthe_xlow * self.Btxy_xlow / (self.Bpxy_xlow
                                                             * self.Rxy_xlow)
        self.pitchcorners = self.hthecorners * self.Btcorners / (self.Bpcorners
                                                                 * self.Rcorners)

        self.dqdpsi = self.DDX_L2C(self.pitch_xlow)
        self.dqdpsi_ylow = self.DDX_L2C(self.pitchcorners, ylow=True)

    def calcHthe(self, ylow=False):
        # hthe = |Grad(theta)|
        # hthe = dtheta/ds at constant psi, phi when psi and theta are orthogonal
        # approx dtheta/sqrt((R(j+1/2)-R(j-1/2))**2 + (Z(j+1/2)-Z(j-1/2)**2)
        assert self.meshParent.orthogonal

        # get positions at j+/-0.5
        R = self.Rxy_ylow
        Z = self.Zxy_ylow

        hthe= self.dy/numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)

        # for hthe_ylow, need R, Z values from below the lower face of this region and
        # above the upper face
        R = numpy.zeros([self.nx, self.ny + 2])
        R[:,1:-1] = self.Rxy
        Z = numpy.zeros([self.nx, self.ny + 2])
        Z[:,1:-1] = self.Zxy
        if self.connections['lower'] is not None:
            R[:,0] = self.getNeighbour('lower').Rxy[:, -1]
            Z[:,0] = self.getNeighbour('lower').Zxy[:, -1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:,0] = 2.*self.Rxy_ylow[:,0] - self.Rxy[:,0]
            Z[:,0] = 2.*self.Zxy_ylow[:,0] - self.Zxy[:,0]
        if self.connections['upper'] is not None:
            R[:,-1] = self.getNeighbour('upper').Rxy[:,0]
            Z[:,-1] = self.getNeighbour('upper').Zxy[:,0]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # value will never even be passed to Mesh, since we stagger to YLOW)
            R[:,-1] = 2.*self.Rxy_ylow[:,-1] - self.Rxy[:,-1]
            Z[:,-1] = 2.*self.Zxy_ylow[:,-1] - self.Zxy[:,-1]

        hthe_ylow =  self.dy_ylow/numpy.sqrt((R[:,1:] - R[:,:-1])**2
                                             + (Z[:,1:] - Z[:,:-1])**2)

        # for hthe_xlow, need R, Z values from the cell corners
        R = self.Rcorners
        Z = self.Zcorners

        hthe_xlow =  self.dy_xlow/numpy.sqrt((R[:,1:] - R[:,:-1])**2
                                             + (Z[:,1:] - Z[:,:-1])**2)

        # for hthecorners, need R, Z values from xlow
        R = numpy.zeros([self.nx+1, self.ny+2])
        Z = numpy.zeros([self.nx+1, self.ny+2])
        R[:,1:-1] = self.Rxy_xlow
        Z[:,1:-1] = self.Zxy_xlow
        if self.connections['lower'] is not None:
            R[:,0] = self.getNeighbour('lower').Rxy_xlow[:,-1]
            Z[:,0] = self.getNeighbour('lower').Zxy_xlow[:,-1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:,0] = 2.*self.Rcorners[:,0] - self.Rxy_xlow[:,0]
            Z[:,0] = 2.*self.Zcorners[:,0] - self.Zxy_xlow[:,0]
        if self.connections['upper'] is not None:
            R[:,-1] = self.getNeighbour('upper').Rxy_xlow[:,0]
            Z[:,-1] = self.getNeighbour('upper').Zxy_xlow[:,0]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # will not even be stored in Mesh, since we stagger to YLOW)
            R[:,-1] = 2.*self.Rcorners[:,-1] - self.Rxy_xlow[:,-1]
            Z[:,-1] = 2.*self.Zcorners[:,-1] - self.Zxy_xlow[:,-1]

        hthecorners =  self.dycorners/numpy.sqrt((R[:,1:] - R[:,:-1])**2
                                             + (Z[:,1:] - Z[:,:-1])**2)

        return hthe, hthe_ylow, hthe_xlow, hthecorners

    def calcBeta(self, ylow=False):
        """
        beta is the angle between x and y coordinates, used for non-orthogonal grid.
        Also calculate radial grid spacing, hrad
        """
        raise ValueError("non-orthogonal grids not calculated yet")

        if not ylow:
            # need to multiply f_R and f_Z by bpsign because we want the radially-outward
            # vector perpendicular to psi contours, and if bpsign is negative then psi
            # increases inward instead of outward so (f_R,f_Z) would be in the opposite
            # direction
            # Actually want the angle of the vector in the y-direction, i.e. (f_Z,-f_R)
            angle_grad_psi = numpy.arctan2(
                    self.bpsign*self.meshParent.equilibrium.f_Z(self.Rxy, self.Zxy),
                    -self.bpsign*self.meshParent.equilibrium.f_R(self.Rxy, self.Zxy))

            R = numpy.zeros([self.nx + 1, self.ny])
            R[:-1, :] = self.Rxy_xlow
            R[-1, :] = self.Rxy_extra_outer
            Z = numpy.zeros([self.nx + 1, self.ny])
            Z[:-1 :] = self.Zxy_xlow
            Z[-1, :] = self.Zxy_extra_outer
            # could calculate radial grid spacing - is it ever needed?
            hrad = numpy.sqrt((R[1:,:] - R[:-1,:])**2 + (Z[1:,:] - Z[:-1,:])**2)

            dR = R[1:,:] - R[:-1,:]
            dZ = Z[1:,:] - Z[:-1,:]
            angle_dr = numpy.arctan2(dR, dZ)
        else:
            # need to multiply f_R and f_Z by bpsign because we want the radially-outward
            # vector perpendicular to psi contours, and if bpsign is negative then psi
            # increases inward instead of outward so (f_R,f_Z) would be in the opposite
            # direction
            # Actually want the angle of the vector in the y-direction, i.e. (f_Z,-f_R)
            angle_grad_psi = numpy.arctan2(
                    self.bpsign*self.meshParent.equilibrium.f_Z(self.Rxy_ylow, self.Zxy_ylow),
                    -self.bpsign*self.meshParent.equilibrium.f_R(self.Rxy_ylow, self.Zxy_ylow))

            # could calculate radial grid spacing - is it ever needed?
            ## for hrad at ylow, can use Rcorners and Zcorners
            hrad = numpy.sqrt((self.Rcorners[1:,:-1] - self.Rcorners[:-1,:-1])**2 +
                              (self.Zcorners[1:,:-1] - self.Zcorners[:-1,:-1])**2)

            dR = self.Rcorners[1:,:-1] - self.Rcorners[:-1,:-1]
            dZ = self.Zcorners[1:,:-1] - self.Zcorners[:-1,:-1]
            angle_dr = numpy.arctan2(dR, dZ)

        return (angle_grad_psi - angle_dr - numpy.pi/2.), hrad

    def getNeighbour(self, face):
        if self.connections[face] is None:
            return None
        else:
            return self.meshParent.regions[self.connections[face]]

    def DDX_L2C(self, f, ylow=False):
        # assume the 'xlow' quantity f has nx+1 values and includes the outer point after
        # the last cell-centre grid point.
        assert f.shape[0] == self.nx + 1

        if not ylow:
            dx = self.dx
        else:
            dx = self.dx_ylow
        result = (f[1:, :] - f[:-1, :]) / dx
        return result

class Mesh:
    """
    Mesh quantities to be written to a grid file for BOUT++
    """
    def __init__(self, equilibrium, meshOptions):
        self.meshOptions = meshOptions
        self.orthogonal = self.readOption('orthogonal')
        self.nx_core = self.readOption('nx_core')
        self.nx_between = self.readOption('nx_between')
        self.nx_sol = self.readOption('nx_sol')
        self.ny_inner_lower_divertor = self.readOption('ny_inner_lower_divertor')
        self.ny_inner_core = self.readOption('ny_inner_core')
        self.ny_inner_upper_divertor = self.readOption('ny_inner_upper_divertor')
        self.ny_outer_upper_divertor = self.readOption('ny_outer_upper_divertor')
        self.ny_outer_core = self.readOption('ny_outer_core')
        self.ny_outer_lower_divertor = self.readOption('ny_outer_lower_divertor')

        self.psi_core = self.readOption('psi_core')
        self.psi_lower_pf = self.readOption('psi_lower_pf', self.psi_core)
        self.psi_upper_pf = self.readOption('psi_upper_pf', self.psi_core)

        self.psi_sol = self.readOption('psi_sol')
        self.psi_inner_sol = self.readOption('psi_inner_sol', self.psi_sol)
        # this option can only be set different from psi_sol in a double-null
        # configuration (i.e. if there are upper divertor legs)
        assert self.ny_inner_upper_divertor > 0 or self.ny_outer_upper_divertor > 0

        self.psi_spacing_separatrix_multiplier = self.readOption('psi_spacing_separatrix_multiplier', None)
        if self.psi_spacing_separatrix_multiplier is not None:
            if self.nx_between > 0:
                raise ValueError("Cannot use psi_spacing_separatrix_multiplier when "
                                 "there are points between two separatrices - to get "
                                 "the same effect, increase the number of points in the "
                                 "between-separatrix region.")

        self.y_boundary_guards = self.readOption('y_boundary_guards')

        self.equilibrium = equilibrium

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

        if self.nx_between > 0:
            # Use uniform spacing of psi in index space in the region between the two
            # separatrices
            assert self.psi_spacing_separatrix_multiplier is None
            assert len(self.psi_sep) == 2
            dpsidi_sep0 = (self.psi_sep[1] - self.psi_sep[0]) / self.nx_between
        else:
            # In index space for indices of cell faces, psi needs to go through psi_inner at
            # 0, psi_sep at nx_core and psi_outer at nx_core+nx_between+nx_sol+1.
            # Estimate base value of dpsi/di as the lower of the average gradients on either
            # side
            if self.psi_core < self.psi_sol:
                dpsidi_sep0 = min((self.equilibrium.psi_sep[0] - self.psi_core) / self.nx_core,
                                  (self.psi_sol - self.equilibrium.psi_sep[0]) / self.nx_sol)
            else:
                dpsidi_sep0 = max((self.equilibrium.psi_sep[0] - self.psi_core) / self.nx_core,
                                  (self.psi_sol - self.equilibrium.psi_sep[0]) / self.nx_sol)

        # decrease (presumably) the spacing around the separatrix by the factor
        # psi_spacing_separatrix_multiplier
        if self.psi_spacing_separatrix_multiplier is not None:
            dpsidi_sep = self.psi_spacing_separatrix_multiplier * dpsidi_sep0
        else:
            dpsidi_sep = dpsidi_sep0

        # fit quadratics on both sides, with gradient dpsidi_sep and values psi_sep at the
        # separatrices and values psi_inner at the inner boundary and psi_outer at the
        # outer boundary

        def getPsiFuncInner(psival):
            # for core region:
            # psi(i) = a*i^2 + b*i + c
            # psi(0) = psival = c
            # psi(nx_core) = psi_sep[0] = a*nx_core**2 + b*nx_core + c
            # dpsidi(nx_core) = dpsidi_sep = 2*a*nx_core + b
            # nx_core*a = dpsidi_sep - psi_sep[0]/nx_core + psival/nx_core
            # b = dpsidi_sep - 2*(dpsidi_sep - psi_sep[0]/nx_core - psival/nx_core)
            a = ( (dpsidi_sep*self.nx_core - self.equilibrium.psi_sep[0] + psival)
                    / self.nx_core**2)
            b = dpsidi_sep - 2.*a*self.nx_core
            c = psival
            return lambda i: a*i**2 + b*i + c

        def getPsiFuncBetweenAndInner(psival):
            if self.nx_between > 0:
                # disconnected double-null

                nx_upper_pf = self.nx_core + self.nx_between
                # for upper PF region:
                # psi(i) = a*i^2 + b*i + c
                # psi(0) = psival = c
                # psi(nx_upper_pf) = psi_sep[1] = a*nx_upper_pf**2 + b*nx_upper_pf + c
                # dpsidi(nx_upper_pf) = dpsidi_sep = 2*a*nx_upper_pf + b
                # nx_upper_pf*a = dpsidi_sep - psi_sep[1]/nx_upper_pf + psival/nx_upper_pf
                # b = dpsidi_sep - 2*(dpsidi_sep - psi_sep[1]/nx_upper_pf - psival/nx_upper_pf)
                a = ( (dpsidi_sep*nx_upper_pf - self.equilibrium.psi_sep[1] + psival)
                        / nx_upper_pf**2)
                b = dpsidi_sep - 2.*a*nx_upper_pf
                c = psival
                return lambda i: a*i**2 + b*i + c

            else:
                # connected double-null
                assert self.nx_between == 0

                nx_upper_pf = self.nx_core
                # for upper PF region:
                # psi(i) = a*i^2 + b*i + c
                # psi(0) = psival = c
                # psi(nx_upper_pf) = psi_sep[0] = a*nx_upper_pf**2 + b*nx_upper_pf + c
                # dpsidi(nx_upper_pf) = dpsidi_sep = 2*a*nx_upper_pf + b
                # nx_upper_pf*a = dpsidi_sep - psi_sep[0]/nx_upper_pf + psival/nx_upper_pf
                # b = dpsidi_sep - 2*(dpsidi_sep - psi_sep[0]/nx_upper_pf - psival/nx_upper_pf)
                a = ( (dpsidi_sep*nx_upper_pf - self.equilibrium.psi_sep[0] + psival)
                        / nx_upper_pf**2)
                b = dpsidi_sep - 2.*a*nx_upper_pf
                c = psival
                return lambda i: a*i**2 + b*i + c

        def getPsiFuncBetween():
            # Constant spacing between the separatrices
            # psi(0) = psi_sep[0]
            # psi(nx_between) = psi_sep[1]
            return lambda i: dpsidi_sep*i + self.equilibrium.psi_sep[0]

        def getPsiFuncOuter(psival):
            if self.nx_between > 0:
                psi_sep = self.equilibrium.psi_sep[1]
            else:
                assert self.nx_between == 0
                psi_sep = self.equilibrium.psi_sep[0]
            # for between separatrix and sol regions:
            # psi(i) = a*i^2 + b*i + c
            # psi(0) = psi_sep = c
            # psi(nx_sol) = psival = a*nx_sol**2 + b*nx_sol + c
            # dpsidi(0) = dpsidi_sep = b
            # a = (psival - dpsidi_sep*nx_sol - psi_inner)/nx_sol**2
            c2 = psi_sep
            b2 = dpsidi_sep
            a2 = (psival - b2*self.nx_sol - c2) / self.nx_sol**2
            return lambda i: a2*i**2 + b2*i + c2

        psi_face_vals_core = numpy.array(
                [getPsiFuncInner(self.psi_core)(i) for i in range(self.nx_core+1)])
        psi_face_vals_lower_pf = numpy.array(
                [getPsiFuncInner(self.psi_lower_pf)(i) for i in range(self.nx_core+1)])
        psi_face_vals_upper_pf = numpy.array(
                [getPsiFuncBetweenAndInner(self.psi_upper_pf)(i)
                    for i in range(self.nx_core+self.nx_between+1)])
        psi_face_vals_between = numpy.array(
                [getPsiFuncBetween()(i) for i in range(0, self.nx_between+1)])
        psi_face_vals_outer_sol = numpy.array(
                [getPsiFuncOuter(self.psi_sol)(i)
                    for i in range(self.nx_between, self.nx_sol+1)])
        psi_face_vals_inner_sol = numpy.array(
                [getPsiFuncOuter(self.psi_inner_sol)(i)
                    for i in range(self.nx_between, self.nx_sol+1)])
        if self.nx_core > 0:
            self.psi_vals_core = numpy.zeros(2*self.nx_core+1)
            self.psi_vals_lower_pf = numpy.zeros(2*self.nx_core+1)
            self.psi_vals_upper_pf = numpy.zeros(2*self.nx_core+1)
        else:
            self.psi_vals_core = numpy.zeros(0)
            self.psi_vals_lower_pf = numpy.zeros(0)
            self.psi_vals_upper_pf = numpy.zeros(0)
        if self.nx_between > 0:
            self.psi_valsbetween = numpy.zeros(2*self.nx_between+1)
        else:
            self.psi_vals_between = numpy.zeros(0)
        if self.nx_sol > 0:
            self.psi_vals_outer_sol = numpy.zeros(2*self.nx_sol+1)
            self.psi_vals_inner_sol = numpy.zeros(2*self.nx_sol+1)
        else:
            self.psi_vals_outer_sol = numpy.zeros(0)
            self.psi_vals_inner_sol = numpy.zeros(0)

        self.psi_vals_core[0::2] = psi_face_vals_core
        self.psi_vals_core[1::2] = 0.5*(psi_face_vals_core[:-1] + psi_face_vals_core[1:])
        self.psi_vals_lower_pf[0::2] = psi_face_vals_lower_pf
        self.psi_vals_lower_pf[1::2] = 0.5*(psi_face_vals_lower_pf[:-1] +
                                            psi_face_vals_lower_pf[1:])
        self.psi_vals_upper_pf[0::2] = psi_face_vals_upper_pf
        self.psi_vals_upper_pf[1::2] = 0.5*(psi_face_vals_upper_pf[:-1] +
                                            psi_face_vals_upper_pf[1:])
        self.psi_vals_between[0::2] = psi_face_vals_between
        self.psi_vals_between[1::2] = 0.5*(psi_face_vals_between[:-1] +
                                                 psi_face_vals_between[1:])
        self.psi_vals_outer_sol[0::2] = psi_face_vals_outer_sol
        self.psi_vals_outer_sol[1::2] = 0.5*(psi_face_vals_outer_sol[:-1] +
                                             psi_face_vals_outer_sol[1:])
        self.psi_vals_inner_sol[0::2] = psi_face_vals_inner_sol
        self.psi_vals_inner_sol[1::2] = 0.5*(psi_face_vals_inner_sol[:-1] +
                                             psi_face_vals_inner_sol[1:])

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
        y_regions = tuple(slice(y_startinds[i], y_startinds[i+1], None)
                     for i in range(len(y_startinds)-1))

        # functions that set poloidal grid spacing:
        # - to use in divertor legs - sqrt of arc length in poloidal plane
        sfunc_leg = lambda s: s**0.5
        # - to use in core regions
        sfunc_core = lambda s: 0.5*(s**0.5 + 1.-(1.-s)**0.5)

        # Region 1 - inner lower PF
        # Region 2 - inner lower between separatrices
        # Region 3 - inner lower SOL
        def get_sep(sepname, ny, connections, sfunc, reverse):
            """
            Utility function to wrap up adding guard cells to separatrix when necessary.
            """
            if reverse:
                # going to reverse the separatrix section after extending, so need to pass
                # lower/upper guards 'backwards'
                if connections['lower'] is None:
                    upper_guards = self.y_boundary_guards
                else:
                    upper_guards = 0
                if connections['upper'] is None:
                    lower_guards = self.y_boundary_guards
                else:
                    lower_guards = 0
            else:
                if connections['lower'] is None:
                    lower_guards = self.y_boundary_guards
                else:
                    lower_guards = 0
                if connections['upper'] is None:
                    upper_guards = self.y_boundary_guards
                else:
                    upper_guards = 0
            sep = self.equilibrium.separatrix[sepname].getRegridded(2*ny+1,
                    extend_lower=2*lower_guards, extend_upper=2*upper_guards, sfunc=sfunc)
            if reverse:
                sep.reverse()
            return sep, ny + lower_guards + upper_guards

        if self.ny_inner_lower_divertor > 0:
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
                sep, localNy = get_sep('inner_lower_divertor',
                        self.ny_inner_lower_divertor, connections, sfunc_leg, True)
                self.regions[1] = MeshRegion(self, 1, self.nx_core, localNy, sep,
                        self.psi_vals_lower_pf, connections, True)
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
                sep, localNy = get_sep('inner_lower_divertor',
                        self.ny_inner_lower_divertor, connections, sfunc_leg, True)
                self.regions[2] = MeshRegion(self, 2, self.nx_between, localNy, sep,
                        self.psi_vals_between, connections, False)
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
                sep, localNy = get_sep('inner_lower_divertor',
                        self.ny_inner_lower_divertor, connections, sfunc_leg, True)
                self.regions[3] = MeshRegion(self, 3, self.nx_sol, localNy, sep,
                        self.psi_vals_inner_sol, connections, False)
                self.region_indices[3] = (numpy.index_exp[x_regions[2], y_regions[0]])

        # Region 4 - inner core
        # Region 5 - inner between separatrices
        # Region 6 - inner SOL
        if self.ny_inner_core > 0:
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
                sep, localNy = get_sep('inner_core', self.ny_inner_core, connections,
                        sfunc_core, False)
                self.regions[4] = MeshRegion(self, 4, self.nx_core, localNy, sep,
                        self.psi_vals_core, connections, True)
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
                sep, localNy = get_sep('inner_core', self.ny_inner_core, connections,
                        sfunc_core, False)
                self.regions[5] = MeshRegion(self, 5, self.nx_between, localNy, sep,
                        self.psi_vals_between, connections, False)
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
                sep, localNy = get_sep('inner_core', self.ny_inner_core, connections,
                        sfunc_core, False)
                self.regions[6] = MeshRegion(self, 6, self.nx_sol, localNy, sep,
                        self.psi_vals_inner_sol, connections, False)
                self.region_indices[6] = (numpy.index_exp[x_regions[2], y_regions[1]])

        # Region 7 - inner upper SOL
        # Region 8 - inner upper PF
        nx_upper_pf = self.nx_between + self.nx_core
        if self.ny_inner_upper_divertor > 0:
            if nx_upper_pf > 0:
                connections = {}
                connections['inner'] = None
                if self.nx_sol > 0:
                    connections['outer'] = 8 # inner upper SOL
                else:
                    connections['outer'] = None
                if self.ny_outer_upper_divertor > 0:
                    connections['lower'] = 9 # outer upper PF
                else:
                    connections['lower'] = None
                connections['upper'] = None
                sep, localNy = get_sep('inner_upper_divertor',
                        self.ny_inner_upper_divertor, connections, sfunc_leg, False)
                self.regions[7] = MeshRegion(self, 7, nx_upper_pf, localNy, sep,
                        self.psi_vals_upper_pf, connections, True)
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
                sep, localNy = get_sep('inner_upper_divertor',
                        self.ny_inner_upper_divertor, connections, sfunc_leg, False)
                self.regions[8] = MeshRegion(self, 8, self.nx_sol, localNy, sep,
                        self.psi_vals_inner_sol, connections, False)
                self.region_indices[8] = (numpy.index_exp[nx_upper_pf:, y_regions[2]])

        # Region 9 - outer upper PF
        # Region 10 - outer upper SOL
        if self.ny_outer_upper_divertor > 0:
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
                sep, localNy = get_sep('outer_upper_divertor',
                        self.ny_outer_upper_divertor, connections, sfunc_leg, True)
                self.regions[9] = MeshRegion(self, 9, nx_upper_pf, localNy, sep,
                        self.psi_vals_upper_pf, connections, True)
                self.region_indices[9] = (numpy.index_exp[:nx_upper_pf, y_regions[3]])
            if self.nx_sol > 0:
                connections = {}
                if nx_upper_pf > 0:
                    connections['inner'] = 9 # outer upper PF
                else:
                    connections['inner'] = None
                connections['outer'] = None
                connections['lower'] = None
                if self.ny_outer_core > 0:
                    connections['upper'] = 13 # outer SOL
                elif self.ny_outer_lower_divertor > 0:
                    connections['upper'] = 16 # outer lower SOL
                else:
                    connections['upper'] = None
                sep, localNy = get_sep('outer_upper_divertor',
                        self.ny_outer_upper_divertor, connections, sfunc_leg, True)
                self.regions[10] = MeshRegion(self, 10, self.nx_sol, localNy, sep,
                        self.psi_vals_outer_sol, connections, False)
                self.region_indices[10] = (numpy.index_exp[nx_upper_pf:, y_regions[3]])

        # Region 11 - outer core
        # Region 12 - outer between separatrices
        # Region 13 - outer SOL
        if self.ny_outer_core > 0:
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
                sep, localNy = get_sep('outer_core', self.ny_outer_core, connections,
                        sfunc_core, False)
                self.regions[11] = MeshRegion(self, 11, self.nx_core, localNy, sep,
                        self.psi_vals_core, connections, True)
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
                sep, localNy = get_sep('outer_core', self.ny_outer_core, connections,
                        sfunc_core, False)
                self.regions[12] = MeshRegion(self, 12, self.nx_between, localNy, sep,
                        self.psi_vals_between, connections, False)
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
                sep, localNy = get_sep('outer_core', self.ny_outer_core, connections,
                        sfunc_core, False)
                self.regions[13] = MeshRegion(self, 13, self.nx_sol, localNy, sep,
                        self.psi_vals_outer_sol, connections, False)
                self.region_indices[13] = (numpy.index_exp[x_regions[2], y_regions[4]])

        # Region 14 - outer lower PF
        # Region 15 - outer lower between separatrices
        # Region 16 - outer lower SOL
        if self.ny_outer_lower_divertor > 0:
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
                connections['upper'] = None
                sep, localNy = get_sep('outer_lower_divertor',
                        self.ny_outer_lower_divertor, connections, sfunc_leg, False)
                self.regions[14] = MeshRegion(self, 14, self.nx_core, localNy, sep,
                        self.psi_vals_lower_pf, connections, True)
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
                connections['upper'] = None
                sep, localNy = get_sep('outer_lower_divertor',
                        self.ny_outer_lower_divertor, connections, sfunc_leg, False)
                self.regions[15] = MeshRegion(self, 15, self.nx_between, localNy, sep,
                        self.psi_vals_between, connections, False)
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
                sep, localNy = get_sep('outer_lower_divertor',
                        self.ny_outer_lower_divertor, connections, sfunc_leg, False)
                self.regions[16] = MeshRegion(self, 16, self.nx_sol, localNy, sep,
                        self.psi_vals_outer_sol, connections, False)
                self.region_indices[16] = (numpy.index_exp[x_regions[2], y_regions[5]])

        # constant spacing in y for now
        self.dy_scalar = 2.*numpy.pi / (self.ny_inner_lower_divertor + self.ny_inner_core
                + self.ny_inner_upper_divertor + self.ny_outer_upper_divertor
                + self.ny_outer_core + self.ny_outer_lower_divertor)

        # create groups that connect in x
        self.x_groups = []
        region_set = set(self.regions.values())
        while region_set:
            for region in region_set:
                if region.connections['inner'] is None:
                    break
            group = []
            while True:
                group.append(region)
                region_set.remove(region)
                region = region.getNeighbour('outer')
                if region is None or group.count(region) > 0:
                    # reached boundary or have all regions in a periodic group
                    break
            self.x_groups.append(group)

        # create groups that connect in y
        self.y_groups = []
        region_set = set(self.regions.values())
        while region_set:
            for region in region_set:
                if region.connections['lower'] is None:
                    break
                # note, if no region with connections['lower']=None is found, then some
                # arbitrary region will be 'region' after this loop. This is OK, as this
                # region must be part of a periodic group, which we will handle.
            group = []
            while True:
                group.append(region)
                region_set.remove(region)
                region = region.getNeighbour('upper')
                if region is None or group.count(region) > 0:
                    # reached boundary or have all regions in a periodic group
                    break
            self.y_groups.append(group)

    def readOption(self, name, default=None):
        print('reading option', name, end='')
        try:
            value = self.meshOptions[name]
            print(':', value)
            return value
        except KeyError:
            print(' - not set - setting default:', default)
            return default

    def geometry(self):
        """
        Calculate geometrical quantities for BOUT++
        """
        for region in self.regions.values():
            region.fillRZ()

        def addFromRegion(f, f_region, regionID):
            f[self.region_indices[regionID]] = f_region

        self.Rxy = numpy.zeros([self.nx, self.ny])
        self.Rxy_ylow = numpy.zeros([self.nx, self.ny])
        self.Rxy_xlow = numpy.zeros([self.nx, self.ny])
        self.Zxy = numpy.zeros([self.nx, self.ny])
        self.Zxy_ylow = numpy.zeros([self.nx, self.ny])
        self.Zxy_xlow = numpy.zeros([self.nx, self.ny])
        self.psixy = numpy.zeros([self.nx, self.ny])
        self.psixy_ylow = numpy.zeros([self.nx, self.ny])
        self.dx = numpy.zeros([self.nx, self.ny])
        self.dx_ylow = numpy.zeros([self.nx, self.ny])
        self.dy = numpy.zeros([self.nx, self.ny])
        self.dy_ylow = numpy.zeros([self.nx, self.ny])
        self.Brxy = numpy.zeros([self.nx, self.ny])
        self.Brxy_ylow = numpy.zeros([self.nx, self.ny])
        self.Bzxy = numpy.zeros([self.nx, self.ny])
        self.Bzxy_ylow = numpy.zeros([self.nx, self.ny])
        self.Bpxy = numpy.zeros([self.nx, self.ny])
        self.Bpxy_ylow = numpy.zeros([self.nx, self.ny])
        self.Btxy = numpy.zeros([self.nx, self.ny])
        self.Btxy_ylow = numpy.zeros([self.nx, self.ny])
        self.Bxy = numpy.zeros([self.nx, self.ny])
        self.Bxy_ylow = numpy.zeros([self.nx, self.ny])
        self.hthe = numpy.zeros([self.nx, self.ny])
        self.hthe_ylow = numpy.zeros([self.nx, self.ny])
        #if not self.orthogonal:
        #    self.beta = numpy.zeros([self.nx, self.ny])
        #    self.beta_ylow = numpy.zeros([self.nx, self.ny])
        #    self.eta = numpy.zeros([self.nx, self.ny])
        #    self.eta_ylow = numpy.zeros([self.nx, self.ny])
        self.pitch = numpy.zeros([self.nx, self.ny])
        self.pitch_ylow = numpy.zeros([self.nx, self.ny])
        self.dqdpsi = numpy.zeros([self.nx, self.ny])
        self.dqdpsi_ylow = numpy.zeros([self.nx, self.ny])

        for region in self.regions.values():
            region.geometry()

            addFromRegion(self.Rxy, region.Rxy, region.myID)
            addFromRegion(self.Rxy_ylow, region.Rxy_ylow[:,:-1], region.myID)
            addFromRegion(self.Rxy_xlow, region.Rxy_xlow[:-1,:], region.myID)
            addFromRegion(self.Zxy, region.Zxy, region.myID)
            addFromRegion(self.Zxy_ylow, region.Zxy_ylow[:,:-1], region.myID)
            addFromRegion(self.Zxy_xlow, region.Zxy_xlow[:-1,:], region.myID)
            addFromRegion(self.psixy, region.psixy, region.myID)
            addFromRegion(self.psixy_ylow, region.psixy_ylow[:,:-1], region.myID)
            addFromRegion(self.dx, region.dx, region.myID)
            addFromRegion(self.dx_ylow, region.dx_ylow[:,:-1], region.myID)
            addFromRegion(self.dy, region.dy, region.myID)
            addFromRegion(self.dy_ylow, region.dy_ylow[:,:-1], region.myID)
            addFromRegion(self.Brxy, region.Brxy, region.myID)
            addFromRegion(self.Brxy_ylow, region.Brxy_ylow[:,:-1], region.myID)
            addFromRegion(self.Bzxy, region.Bzxy, region.myID)
            addFromRegion(self.Bzxy_ylow, region.Bzxy_ylow[:,:-1], region.myID)
            addFromRegion(self.Bpxy, region.Bpxy, region.myID)
            addFromRegion(self.Bpxy_ylow, region.Bpxy_ylow[:,:-1], region.myID)
            addFromRegion(self.Btxy, region.Btxy, region.myID)
            addFromRegion(self.Btxy_ylow, region.Btxy_ylow[:,:-1], region.myID)
            addFromRegion(self.Bxy, region.Bxy, region.myID)
            addFromRegion(self.Bxy_ylow, region.Bxy_ylow[:,:-1], region.myID)
            addFromRegion(self.hthe, region.hthe, region.myID)
            addFromRegion(self.hthe_ylow, region.hthe_ylow[:,:-1], region.myID)
            #if not self.orthogonal:
            #    addFromRegion(self.beta, region.beta, region.myID)
            #    addFromRegion(self.beta_ylow, region.beta_ylow[:,:-1], region.myID)
            #    addFromRegion(self.eta, region.eta, region.myID)
            #    addFromRegion(self.eta_ylow, region.eta_ylow[:,:-1], region.myID)
            addFromRegion(self.pitch, region.pitch, region.myID)
            addFromRegion(self.pitch_ylow, region.pitch_ylow[:,:-1], region.myID)
            addFromRegion(self.dqdpsi, region.dqdpsi, region.myID)
            addFromRegion(self.dqdpsi_ylow, region.dqdpsi_ylow[:,:-1], region.myID)

    def writeGridfile(self, filename):
        from boututils.datafile import DataFile

        with DataFile(filename, create=True) as f:
            f.write('nx', self.nx)
            # ny for BOUT++ excludes boundary guard cells
            f.write('ny', self.ny - 2*self.y_boundary_guards
                          - 2*self.upper_target_y_boundary_guards)
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
            f.write('hthe', self.hthe)
            f.write('hthe_ylow', self.hthe_ylow)
            #if not self.orthogonal:
            #    f.write('beta', self.beta)
            #    f.write('beta_ylow', self.beta_ylow)
            #    f.write('eta', self.eta)
            #    f.write('eta_ylow', self.eta_ylow)
            f.write('pitch', self.pitch)
            f.write('pitch_ylow', self.pitch_ylow)

    def plot2D(self, f, title=None):
        from matplotlib import pyplot

        try:
            vmin = f.min()
            vmax = f.max()

            for region, indices in zip(self.regions.values(), self.region_indices.values()):
                pyplot.pcolor(region.Rcorners, region.Zcorners, f[indices],
                              vmin=vmin, vmax=vmax)

            pyplot.colorbar()
        except NameError:
            raise NameError('Some variable has not been defined yet: have you called Mesh.geometry()?')

    def plotPoints(self, xlow=False, ylow=False, corners=False):
        from matplotlib import pyplot
        from cycler import cycle

        colors = cycle(pyplot.rcParams['axes.prop_cycle'].by_key()['color'])
        for region in self.regions.values():
            c = next(colors)
            pyplot.scatter(region.Rxy, region.Zxy, marker='x', c=c)
            if xlow:
                pyplot.scatter(region.Rxy_xlow, region.Zxy_xlow, marker='1', c=c)
            if ylow:
                pyplot.scatter(region.Rxy_ylow, region.Zxy_ylow, marker='2', c=c)
            if corners:
                pyplot.scatter(region.Rcorners, region.Zcorners, marker='+', c=c)

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
