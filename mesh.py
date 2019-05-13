"""
Classes to handle Meshes and geometrical quantities for generating BOUT++ grids
"""

import numpy
import numbers
from scipy.integrate import solve_ivp
import warnings
from equilibrium import Point2D, PsiContour, EquilibriumRegion

class MultiLocationArray(numpy.lib.mixins.NDArrayOperatorsMixin):
    """
    Container for arrays representing points at different cell locations
    Not all have to be filled.
    """
    _centre_array = None
    _xlow_array = None
    _ylow_array = None
    _corners_array = None

    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny

    @property
    def centre(self):
        if self._centre_array is None:
            self._centre_array = numpy.zeros([self.nx, self.ny])
        return self._centre_array

    @centre.setter
    def centre(self, value):
        if self._centre_array is None:
            self._centre_array = numpy.zeros([self.nx, self.ny])
        self._centre_array[...] = value

    @property
    def xlow(self):
        if self._xlow_array is None:
            self._xlow_array = numpy.zeros([self.nx + 1, self.ny])
        return self._xlow_array

    @xlow.setter
    def xlow(self, value):
        if self._xlow_array is None:
            self._xlow_array = numpy.zeros([self.nx + 1, self.ny])
        self._xlow_array[...] = value

    @property
    def ylow(self):
        if self._ylow_array is None:
            self._ylow_array = numpy.zeros([self.nx, self.ny + 1])
        return self._ylow_array

    @ylow.setter
    def ylow(self, value):
        if self._ylow_array is None:
            self._ylow_array = numpy.zeros([self.nx, self.ny + 1])
        self._ylow_array[...] = value

    @property
    def corners(self):
        if self._corners_array is None:
            self._corners_array = numpy.zeros([self.nx + 1, self.ny + 1])
        return self._corners_array

    @corners.setter
    def corners(self, value):
        if self._corners_array is None:
            self._corners_array = numpy.zeros([self.nx + 1, self.ny + 1])
        self._corners_array[...] = value

    # The following __array_ufunc__ implementation allows the MultiLocationArray class to
    # be handled by Numpy functions, and add, subtract, etc. like an ndarray.
    # The implementation is mostly copied from the example in
    # https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.mixins.NDArrayOperatorsMixin.html#numpy.lib.mixins.NDArrayOperatorsMixin

    # One might also consider adding the built-in list type to this
    # list, to support operations like np.add(array_like, list)
    _HANDLED_TYPES = (numpy.ndarray, numbers.Number)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.get('out', ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use MultiLocationArray instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle MultiLocationArray objects.
            if not isinstance(x, self._HANDLED_TYPES + (MultiLocationArray,)):
                return NotImplemented

        result = MultiLocationArray(self.nx, self.ny)

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._centre_array is not None:
            this_inputs = tuple(x._centre_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.centre if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].centre = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.centre = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._xlow_array is not None:
            this_inputs = tuple(x._xlow_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.xlow if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].xlow = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.xlow = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._ylow_array is not None:
            this_inputs = tuple(x._ylow_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.ylow if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].ylow = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.ylow = this_result

        # Defer to the implementation of the ufunc on unwrapped values.
        if self._corners_array is not None:
            this_inputs = tuple(x._corners_array if isinstance(x, MultiLocationArray)
                    else x for x in inputs)
            if out:
                kwargs['out'] = tuple(
                    x.corners if isinstance(x, MultiLocationArray) else x
                    for x in out)
            this_result = getattr(ufunc, method)(*this_inputs, **kwargs)

            if type(this_result) is tuple:
                # multiple return values
                if not type(result) is tuple:
                    result = tuple(MultiLocationArray(self.nx, self.ny)
                            for x in this_result)
                for i,x in enumerate(this_result):
                    result[i].corners = x

            elif method == 'at':
                # no return value
                result = None
            else:
                # one return value
                result.corners = this_result

        return result

class MeshRegion:
    """
    A simple rectangular region of a Mesh, that connects to one other region (or has a
    boundary) on each edge.
    Note that these regions include cell face and boundary points, so there are
    (2nx+1)*(2ny+1) points for an nx*ny grid.
    """
    def __init__(self, meshParent, myID, equilibriumRegion, psi_vals, connections, radialIndex):
        print('creating region', myID, '-',
                equilibriumRegion.name+'('+str(radialIndex)+')')

        # the Mesh object that owns this MeshRegion
        self.meshParent = meshParent

        # ID that Mesh uses to keep track of its MeshRegions
        self.myID = myID

        # psi values for radial grid
        self.psi_vals = psi_vals

        # EquilibriumRegion representing the segment associated with this region
        self.equilibriumRegion = equilibriumRegion

        # sizes of the grid in this MeshRegion, include boundary guard cells
        self.nx = self.equilibriumRegion.nx[radialIndex]
        self.ny = self.equilibriumRegion.ny(radialIndex)
        self.ny_noguards = self.equilibriumRegion.ny_noguards

        # Dictionary that specifies whether a boundary is connected to another region or
        # is an actual boundary
        self.connections = connections

        # Number of this region, counting radially outward
        self.radialIndex = radialIndex

        # # y-boundary guard cells needed if the region edge is a real boundary, i.e. not
        # # connected to another region
        # if self.connections['lower'] is None:
        #     self.y_guards_lower = self.equilibriumRegion.y_boundary_guards
        # else:
        #     self.y_guards_lower = 0
        # if self.connections['upper'] is None:
        #     self.y_guards_upper = self.equilibriumRegion.y_boundary_guards
        # else:
        #     self.y_guards_upper = 0

        # get points in this region
        self.contours = []
        if self.radialIndex == 0:
            temp_psi_vals = self.psi_vals[::-1]
        else:
            temp_psi_vals = self.psi_vals
        perp_points = followPerpendicular(meshParent.equilibrium.f_R,
                meshParent.equilibrium.f_Z, self.equilibriumRegion[0],
                meshParent.equilibrium.psi_sep[0], temp_psi_vals)
        if self.radialIndex == 0:
            perp_points.reverse()
        for i,point in enumerate(perp_points):
            self.contours.append(PsiContour([point], meshParent.equilibrium.psi,
                self.psi_vals[i]))
        for p in self.equilibriumRegion[1:]:
            perp_points = followPerpendicular(meshParent.equilibrium.f_R,
                    meshParent.equilibrium.f_Z, p, meshParent.equilibrium.psi_sep[0],
                    temp_psi_vals)
            if self.radialIndex == 0:
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

        self.Rxy = MultiLocationArray(self.nx, self.ny)
        self.Zxy = MultiLocationArray(self.nx, self.ny)

        self.Rxy.centre = numpy.array([[p.R for p in contour[1::2]]
            for contour in self.contours[1::2]])

        self.Rxy.ylow = numpy.array([[p.R for p in contour[0::2]]
            for contour in self.contours[1::2]])

        self.Rxy.xlow = numpy.array([[p.R for p in contour[1::2]]
            for contour in self.contours[0::2]])

        self.Zxy.centre = numpy.array( [[p.Z for p in contour[1::2]]
            for contour in self.contours[1::2]])

        self.Zxy.ylow = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in self.contours[1::2]])

        self.Zxy.xlow = numpy.array([[p.Z for p in contour[1::2]]
            for contour in self.contours[0::2]])

        self.Rxy.corners = numpy.array( [[p.R for p in contour[0::2]]
            for contour in self.contours[0::2]])
        self.Zxy.corners = numpy.array( [[p.Z for p in contour[0::2]]
            for contour in self.contours[0::2]])

        # Fix up the corner values at the X-points. Because the PsiContour have to start
        # slightly away from the X-point in order for the integrator to go in the right
        # direction, the points that should be at the X-point will be slighly displaced,
        # and will not be consistent between regions. So replace these points with the
        # X-point position instead.
        xpoint = self.equilibriumRegion.xPointsAtStart[self.radialIndex]
        if xpoint is not None:
            self.Rxy.corners[0,0] = xpoint.R
            self.Zxy.corners[0,0] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtStart[self.radialIndex+1]
        if xpoint is not None:
            self.Rxy.corners[-1,0] = xpoint.R
            self.Zxy.corners[-1,0] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtEnd[self.radialIndex]
        if xpoint is not None:
            self.Rxy.corners[0,-1] = xpoint.R
            self.Zxy.corners[0,-1] = xpoint.Z

        xpoint = self.equilibriumRegion.xPointsAtEnd[self.radialIndex+1]
        if xpoint is not None:
            self.Rxy.corners[-1,-1] = xpoint.R
            self.Zxy.corners[-1,-1] = xpoint.Z

    def getRZBoundary(self):
        # Upper value of ylow array logically overlaps with the lower value in the upper
        # neighbour. They should be close, but aren't guaranteed to be identical already
        # because they were taken from separate PsiContour obects. Use the value from the
        # upper neighbour to ensure consistency.
        # Also do similarly for the corner arrays.
        # Don't need to do this for the x-boundaries, because there the PsiContour
        # objects are shared between neighbouring regions.
        #
        # This needs to be a separate method from fillRZ() so that it can be called after
        # all regions have filled their Rxy and Zxy arrays.
        if self.connections['upper'] is not None:
            up = self.getNeighbour('upper')
            self.Rxy.ylow[:,-1] = up.Rxy.ylow[:,0]
            self.Zxy.ylow[:,-1] = up.Zxy.ylow[:,0]
            self.Rxy.corners[:,-1] = up.Rxy.corners[:,0]
            self.Zxy.corners[:,-1] = up.Zxy.corners[:,0]

    def geometry(self):
        """
        Calculate geometrical quantities for this region
        """

        self.psixy = self.meshParent.equilibrium.psi(self.Rxy, self.Zxy)

        self.dx = MultiLocationArray(self.nx, self.ny)
        self.dx.centre = numpy.array(self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]
        self.dx.ylow = numpy.array(self.psi_vals[2::2] - self.psi_vals[:-2:2])[:, numpy.newaxis]

        if self.psi_vals[0] > self.psi_vals[-1]:
            # x-coordinate is -psixy so x always increases radially across grid
            self.bpsign = -1.
            self.xcoord = -self.psixy
        else:
            self.bpsign = 1.
            self.xcoord = self.psixy

        self.dy = MultiLocationArray(self.nx, self.ny)
        self.dy.centre = self.meshParent.dy_scalar
        self.dy.ylow = self.meshParent.dy_scalar
        self.dy.xlow = self.meshParent.dy_scalar
        self.dy.corners = self.meshParent.dy_scalar

        self.Brxy = self.meshParent.equilibrium.Bp_R(self.Rxy, self.Zxy)
        self.Bzxy = self.meshParent.equilibrium.Bp_Z(self.Rxy, self.Zxy)
        self.Bpxy = numpy.sqrt(self.Brxy**2 + self.Bzxy**2)
        # determine direction - dot Bp with Grad(y) vector
        # evaluate in 'sol' at outer radial boundary
        Bp_dot_grady = (
            self.Brxy.centre[-1, self.ny//2]
            *(self.Rxy.centre[-1, self.ny//2 + 1] - self.Rxy.centre[-1, self.ny//2 - 1])
            + self.Bzxy.centre[-1, self.ny//2]
              *(self.Zxy.centre[-1, self.ny//2 + 1] - self.Zxy.centre[-1, self.ny//2 - 1]) )
        if Bp_dot_grady < 0.:
            print("Poloidal field is in opposite direction to Grad(theta) -> Bp negative")
            self.Bpxy = -self.Bpxy
            if self.bpsign > 0.:
                raise ValueError("Sign of Bp should be negative?")
        else:
            if self.bpsign < 0.:
                raise ValueError("Sign of Bp should be positive?")

        # Get toroidal field from poloidal current function fpol
        self.Btxy = self.meshParent.equilibrium.fpol(self.psixy) / self.Rxy

        self.Bxy = numpy.sqrt(self.Bpxy**2 + self.Btxy**2)

        self.hthe = self.calcHthe()

        #if not self.meshParent.orthogonal:
        #    # Calculate beta (angle between x and y coordinates), used for non-orthogonal grid
        #    # Also calculate radial grid spacing
        #    self.beta, self.hrad = self.calcBeta()

        #    # eta is the polodial non-orthogonality parameter
        #    self.eta = numpy.sin(self.beta)
        #else:
        #    self.beta.centre = 0.
        #    self.eta.centre = 0.

        # field line pitch
        self.pitch = self.hthe * self.Btxy / (self.Bpxy * self.Rxy)

        self.dqdpsi = MultiLocationArray(self.nx, self.ny)
        self.dqdpsi.centre = self.DDX_L2C(self.pitch.xlow)
        self.dqdpsi.ylow = self.DDX_L2C(self.pitch.corners, ylow=True)

    def calcHthe(self, ylow=False):
        # hthe = |Grad(theta)|
        # hthe = dtheta/ds at constant psi, phi when psi and theta are orthogonal
        # approx dtheta/sqrt((R(j+1/2)-R(j-1/2))**2 + (Z(j+1/2)-Z(j-1/2)**2)
        assert self.meshParent.orthogonal

        # get positions at j+/-0.5
        R = self.Rxy.ylow
        Z = self.Zxy.ylow

        hthe = MultiLocationArray(self.nx, self.ny)
        hthe.centre = self.dy.centre/numpy.sqrt((R[:,1:] - R[:,:-1])**2 + (Z[:,1:] - Z[:,:-1])**2)

        # for hthe_ylow, need R, Z values from below the lower face of this region and
        # above the upper face
        R = numpy.zeros([self.nx, self.ny + 2])
        R[:,1:-1] = self.Rxy.centre
        Z = numpy.zeros([self.nx, self.ny + 2])
        Z[:,1:-1] = self.Zxy.centre
        if self.connections['lower'] is not None:
            R[:,0] = self.getNeighbour('lower').Rxy.centre[:, -1]
            Z[:,0] = self.getNeighbour('lower').Zxy.centre[:, -1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:,0] = 2.*self.Rxy.ylow[:,0] - self.Rxy.centre[:,0]
            Z[:,0] = 2.*self.Zxy.ylow[:,0] - self.Zxy.centre[:,0]
        if self.connections['upper'] is not None:
            R[:,-1] = self.getNeighbour('upper').Rxy.centre[:,0]
            Z[:,-1] = self.getNeighbour('upper').Zxy.centre[:,0]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # value will never even be passed to Mesh, since we stagger to YLOW)
            R[:,-1] = 2.*self.Rxy.ylow[:,-1] - self.Rxy.centre[:,-1]
            Z[:,-1] = 2.*self.Zxy.ylow[:,-1] - self.Zxy.centre[:,-1]

        hthe.ylow =  self.dy.ylow/numpy.sqrt((R[:,1:] - R[:,:-1])**2
                                              + (Z[:,1:] - Z[:,:-1])**2)

        # for hthe_xlow, need R, Z values from the cell corners
        R = self.Rxy.corners
        Z = self.Zxy.corners

        hthe.xlow =  self.dy.xlow/numpy.sqrt((R[:,1:] - R[:,:-1])**2
                                              + (Z[:,1:] - Z[:,:-1])**2)

        # for hthecorners, need R, Z values from xlow
        R = numpy.zeros([self.nx+1, self.ny+2])
        Z = numpy.zeros([self.nx+1, self.ny+2])
        R[:,1:-1] = self.Rxy.xlow
        Z[:,1:-1] = self.Zxy.xlow
        if self.connections['lower'] is not None:
            R[:,0] = self.getNeighbour('lower').Rxy.xlow[:,-1]
            Z[:,0] = self.getNeighbour('lower').Zxy.xlow[:,-1]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # corresponding value at the upper boundary does not even exist, since we
            # stagger to YLOW)
            R[:,0] = 2.*self.Rxy.corners[:,0] - self.Rxy.xlow[:,0]
            Z[:,0] = 2.*self.Zxy.corners[:,0] - self.Zxy.xlow[:,0]
        if self.connections['upper'] is not None:
            R[:,-1] = self.getNeighbour('upper').Rxy.xlow[:,0]
            Z[:,-1] = self.getNeighbour('upper').Zxy.xlow[:,0]
        else:
            # dumb extrapolation, but should not need the affected guard cell value (the
            # will not even be stored in Mesh, since we stagger to YLOW)
            R[:,-1] = 2.*self.Rxy.corners[:,-1] - self.Rxy.xlow[:,-1]
            Z[:,-1] = 2.*self.Zxy.corners[:,-1] - self.Zxy.xlow[:,-1]

        hthe.corners =  self.dy.corners/numpy.sqrt((R[:,1:] - R[:,:-1])**2
                                                    + (Z[:,1:] - Z[:,:-1])**2)

        return hthe

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
            dx = self.dx.centre
        else:
            dx = self.dx.ylow
        result = (f[1:, :] - f[:-1, :]) / dx
        return result

class Mesh:
    """
    Mesh represented by a collection of connected MeshRegion objects
    """
    def __init__(self, equilibrium, meshOptions):
        self.meshOptions = meshOptions
        self.orthogonal = self.readOption('orthogonal')

        self.equilibrium = equilibrium

        assert self.orthogonal # non-orthogonal not implelemented yet

        # Generate MeshRegion object for each section of the mesh
        self.regions = {}

        # functions that set poloidal grid spacing:
        # - to use in divertor legs with X-point at start - sqrt of arc length in poloidal plane
        sfunc_leg_start = lambda s: s**0.5
        # - to use in divertor legs with X-point at end - sqrt of arc length in poloidal plane
        sfunc_leg_end = lambda s: 1.-(1.-s)**0.5
        # - to use in core regions
        sfunc_core = lambda s: 0.5*(s**0.5 + 1.-(1.-s)**0.5)
        # - to use in regions with no X-point (divertor targets at both ends)
        sfunc_noX = lambda s: s

        # Make consecutive numbering scheme for regions
        regionlist = []
        self.region_lookup = {}
        for reg_name,eq_reg in equilibrium.regions.items():
            for i in range(eq_reg.nSegments):
                region_number = len(regionlist)
                regionlist.append((reg_name, i))
                self.region_lookup[(reg_name, i)] = region_number
        # Get connections between regions
        connections = {}
        psi_vals = {}
        for region_id,(eq_reg,i) in enumerate(regionlist):
            connections[region_id] = {}
            region = equilibrium.regions[eq_reg]
            c = region.connections[i]
            for key, val in c.items():
                if val is not None:
                    connections[region_id][key] = self.region_lookup[val]
                else:
                    connections[region_id][key] = None
            warnings.warn("Using linear function for psi")
            psi_vals[region_id] = numpy.linspace(region.psi_boundaries[i],
                    region.psi_boundaries[i+1], 2*region.nx[i]+1)

        for eq_region in self.equilibrium.regions.values():
            for i in range(eq_region.nSegments):
                region_id = self.region_lookup[(eq_region.name,i)]
                if (eq_region.connections[i]['lower'] is None and
                        eq_region.connections[i]['upper'] is None):
                    sfunc = sfunc_noX
                elif eq_region.connections[i]['lower'] is None:
                    sfunc = sfunc_leg_end
                elif eq_region.connections[i]['upper'] is None:
                    sfunc = sfunc_leg_start
                else:
                    sfunc = sfunc_core
                eq_region_with_boundaries = eq_region.getRegridded(radialIndex=i,
                        sfunc=sfunc)
                self.regions[region_id] = MeshRegion(self, region_id,
                        eq_region_with_boundaries, psi_vals[region_id],
                        connections[region_id], i)

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
        for region in self.regions.values():
            region.getRZBoundary()
        for region in self.regions.values():
            region.geometry()

    def plotPoints(self, xlow=False, ylow=False, corners=False):
        from matplotlib import pyplot
        from cycler import cycle

        colors = cycle(pyplot.rcParams['axes.prop_cycle'].by_key()['color'])
        for region in self.regions.values():
            c = next(colors)
            pyplot.scatter(region.Rxy.centre, region.Zxy.centre, marker='x', c=c,
                    label=region.myID)
            if xlow:
                pyplot.scatter(region.Rxy.xlow, region.Zxy.xlow, marker='1', c=c)
            if ylow:
                pyplot.scatter(region.Rxy.ylow, region.Zxy.ylow, marker='2', c=c)
            if corners:
                pyplot.scatter(region.Rxy.corners, region.Zxy.corners, marker='+', c=c)
        pyplot.legend()

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

class BoutMesh(Mesh):
    """
    Mesh quantities to be written to a grid file for BOUT++

    Requires that the MeshRegion members fit together into a global logically-rectangular
    Mesh, with the topology assumed by BOUT++ (allowing complexity up to
    disconnected-double-null).

    For compatibility with BOUT++, the regions in the OrderedDict equilibrium.regions must
    be in the order: inner_lower_divertor, inner_core, inner_upper_divertor,
    outer_upper_divertor, outer_core, outer_lower_divertor. This ensures the correct
    positioning in the global logically rectangular grid. Regions are allowed to not be
    present (if they would have size 0).
    """
    def __init__(self, equilibrium, meshOptions):

        super().__init__(equilibrium, meshOptions)

        # nx, ny both include boundary guard cells
        eq_region0 = next(iter(self.equilibrium.regions.values()))
        self.nx = sum(eq_region0.nx)

        self.ny = sum(r.ny(0) for r in self.equilibrium.regions.values())

        self.ny_noguards = sum(r.ny_noguards for r in self.equilibrium.regions.values())

        # if self.nx_between > 0:
        #     # Use uniform spacing of psi in index space in the region between the two
        #     # separatrices
        #     assert self.psi_spacing_separatrix_multiplier is None
        #     assert len(self.psi_sep) == 2
        #     dpsidi_sep0 = (self.psi_sep[1] - self.psi_sep[0]) / self.nx_between
        # else:
        #     # In index space for indices of cell faces, psi needs to go through psi_inner at
        #     # 0, psi_sep at nx_core and psi_outer at nx_core+nx_between+nx_sol+1.
        #     # Estimate base value of dpsi/di as the lower of the average gradients on either
        #     # side
        #     if self.psi_core < self.psi_sol:
        #         dpsidi_sep0 = min((self.equilibrium.psi_sep[0] - self.psi_core) / self.nx_core,
        #                           (self.psi_sol - self.equilibrium.psi_sep[0]) / self.nx_sol)
        #     else:
        #         dpsidi_sep0 = max((self.equilibrium.psi_sep[0] - self.psi_core) / self.nx_core,
        #                           (self.psi_sol - self.equilibrium.psi_sep[0]) / self.nx_sol)

        # # decrease (presumably) the spacing around the separatrix by the factor
        # # psi_spacing_separatrix_multiplier
        # if self.psi_spacing_separatrix_multiplier is not None:
        #     dpsidi_sep = self.psi_spacing_separatrix_multiplier * dpsidi_sep0
        # else:
        #     dpsidi_sep = dpsidi_sep0

        # # fit quadratics on both sides, with gradient dpsidi_sep and values psi_sep at the
        # # separatrices and values psi_inner at the inner boundary and psi_outer at the
        # # outer boundary

        # def getPsiFuncInner(psival):
        #     # for core region:
        #     # psi(i) = a*i^2 + b*i + c
        #     # psi(0) = psival = c
        #     # psi(nx_core) = psi_sep[0] = a*nx_core**2 + b*nx_core + c
        #     # dpsidi(nx_core) = dpsidi_sep = 2*a*nx_core + b
        #     # nx_core*a = dpsidi_sep - psi_sep[0]/nx_core + psival/nx_core
        #     # b = dpsidi_sep - 2*(dpsidi_sep - psi_sep[0]/nx_core - psival/nx_core)
        #     a = ( (dpsidi_sep*self.nx_core - self.equilibrium.psi_sep[0] + psival)
        #             / self.nx_core**2)
        #     b = dpsidi_sep - 2.*a*self.nx_core
        #     c = psival
        #     return lambda i: a*i**2 + b*i + c

        # def getPsiFuncBetweenAndInner(psival):
        #     if self.nx_between > 0:
        #         # disconnected double-null

        #         nx_upper_pf = self.nx_core + self.nx_between
        #         # for upper PF region:
        #         # psi(i) = a*i^2 + b*i + c
        #         # psi(0) = psival = c
        #         # psi(nx_upper_pf) = psi_sep[1] = a*nx_upper_pf**2 + b*nx_upper_pf + c
        #         # dpsidi(nx_upper_pf) = dpsidi_sep = 2*a*nx_upper_pf + b
        #         # nx_upper_pf*a = dpsidi_sep - psi_sep[1]/nx_upper_pf + psival/nx_upper_pf
        #         # b = dpsidi_sep - 2*(dpsidi_sep - psi_sep[1]/nx_upper_pf - psival/nx_upper_pf)
        #         a = ( (dpsidi_sep*nx_upper_pf - self.equilibrium.psi_sep[1] + psival)
        #                 / nx_upper_pf**2)
        #         b = dpsidi_sep - 2.*a*nx_upper_pf
        #         c = psival
        #         return lambda i: a*i**2 + b*i + c

        #     else:
        #         # connected double-null
        #         assert self.nx_between == 0

        #         nx_upper_pf = self.nx_core
        #         # for upper PF region:
        #         # psi(i) = a*i^2 + b*i + c
        #         # psi(0) = psival = c
        #         # psi(nx_upper_pf) = psi_sep[0] = a*nx_upper_pf**2 + b*nx_upper_pf + c
        #         # dpsidi(nx_upper_pf) = dpsidi_sep = 2*a*nx_upper_pf + b
        #         # nx_upper_pf*a = dpsidi_sep - psi_sep[0]/nx_upper_pf + psival/nx_upper_pf
        #         # b = dpsidi_sep - 2*(dpsidi_sep - psi_sep[0]/nx_upper_pf - psival/nx_upper_pf)
        #         a = ( (dpsidi_sep*nx_upper_pf - self.equilibrium.psi_sep[0] + psival)
        #                 / nx_upper_pf**2)
        #         b = dpsidi_sep - 2.*a*nx_upper_pf
        #         c = psival
        #         return lambda i: a*i**2 + b*i + c

        # def getPsiFuncBetween():
        #     # Constant spacing between the separatrices
        #     # psi(0) = psi_sep[0]
        #     # psi(nx_between) = psi_sep[1]
        #     return lambda i: dpsidi_sep*i + self.equilibrium.psi_sep[0]

        # def getPsiFuncOuter(psival):
        #     if self.nx_between > 0:
        #         psi_sep = self.equilibrium.psi_sep[1]
        #     else:
        #         assert self.nx_between == 0
        #         psi_sep = self.equilibrium.psi_sep[0]
        #     # for between separatrix and sol regions:
        #     # psi(i) = a*i^2 + b*i + c
        #     # psi(0) = psi_sep = c
        #     # psi(nx_sol) = psival = a*nx_sol**2 + b*nx_sol + c
        #     # dpsidi(0) = dpsidi_sep = b
        #     # a = (psival - dpsidi_sep*nx_sol - psi_inner)/nx_sol**2
        #     c2 = psi_sep
        #     b2 = dpsidi_sep
        #     a2 = (psival - b2*self.nx_sol - c2) / self.nx_sol**2
        #     return lambda i: a2*i**2 + b2*i + c2

        # psi_face_vals_core = numpy.array(
        #         [getPsiFuncInner(self.psi_core)(i) for i in range(self.nx_core+1)])
        # psi_face_vals_lower_pf = numpy.array(
        #         [getPsiFuncInner(self.psi_lower_pf)(i) for i in range(self.nx_core+1)])
        # psi_face_vals_upper_pf = numpy.array(
        #         [getPsiFuncBetweenAndInner(self.psi_upper_pf)(i)
        #             for i in range(self.nx_core+self.nx_between+1)])
        # psi_face_vals_between = numpy.array(
        #         [getPsiFuncBetween()(i) for i in range(0, self.nx_between+1)])
        # psi_face_vals_outer_sol = numpy.array(
        #         [getPsiFuncOuter(self.psi_sol)(i)
        #             for i in range(self.nx_between, self.nx_sol+1)])
        # psi_face_vals_inner_sol = numpy.array(
        #         [getPsiFuncOuter(self.psi_inner_sol)(i)
        #             for i in range(self.nx_between, self.nx_sol+1)])
        # if self.nx_core > 0:
        #     self.psi_vals_core = numpy.zeros(2*self.nx_core+1)
        #     self.psi_vals_lower_pf = numpy.zeros(2*self.nx_core+1)
        #     self.psi_vals_upper_pf = numpy.zeros(2*self.nx_core+1)
        # else:
        #     self.psi_vals_core = numpy.zeros(0)
        #     self.psi_vals_lower_pf = numpy.zeros(0)
        #     self.psi_vals_upper_pf = numpy.zeros(0)
        # if self.nx_between > 0:
        #     self.psi_valsbetween = numpy.zeros(2*self.nx_between+1)
        # else:
        #     self.psi_vals_between = numpy.zeros(0)
        # if self.nx_sol > 0:
        #     self.psi_vals_outer_sol = numpy.zeros(2*self.nx_sol+1)
        #     self.psi_vals_inner_sol = numpy.zeros(2*self.nx_sol+1)
        # else:
        #     self.psi_vals_outer_sol = numpy.zeros(0)
        #     self.psi_vals_inner_sol = numpy.zeros(0)

        # self.psi_vals_core[0::2] = psi_face_vals_core
        # self.psi_vals_core[1::2] = 0.5*(psi_face_vals_core[:-1] + psi_face_vals_core[1:])
        # self.psi_vals_lower_pf[0::2] = psi_face_vals_lower_pf
        # self.psi_vals_lower_pf[1::2] = 0.5*(psi_face_vals_lower_pf[:-1] +
        #                                     psi_face_vals_lower_pf[1:])
        # self.psi_vals_upper_pf[0::2] = psi_face_vals_upper_pf
        # self.psi_vals_upper_pf[1::2] = 0.5*(psi_face_vals_upper_pf[:-1] +
        #                                     psi_face_vals_upper_pf[1:])
        # self.psi_vals_between[0::2] = psi_face_vals_between
        # self.psi_vals_between[1::2] = 0.5*(psi_face_vals_between[:-1] +
        #                                          psi_face_vals_between[1:])
        # self.psi_vals_outer_sol[0::2] = psi_face_vals_outer_sol
        # self.psi_vals_outer_sol[1::2] = 0.5*(psi_face_vals_outer_sol[:-1] +
        #                                      psi_face_vals_outer_sol[1:])
        # self.psi_vals_inner_sol[0::2] = psi_face_vals_inner_sol
        # self.psi_vals_inner_sol[1::2] = 0.5*(psi_face_vals_inner_sol[:-1] +
        #                                      psi_face_vals_inner_sol[1:])

        # For region numbers see figure in 'BOUT++ Topology' section of BOUT++ manual -
        # except here '7' is inner upper PF and '8' is inner upper SOL, which is the other
        # way around from that figure, but means all regions are labelled from inner-x to
        # outer-x

        # Keep ranges of global indices for each region separately, because we don't want
        # MeshRegion objects to depend on global indices
        assert all([r.nx == eq_region0.nx for r in self.equilibrium.regions.values()])
        x_sizes = [0] + list(eq_region0.nx)
        x_startinds = numpy.cumsum(x_sizes)
        x_regions = tuple(slice(x_startinds[i], x_startinds[i+1], None)
                     for i in range(len(x_startinds)-1))
        y_total = 0
        y_regions = {}
        for regname, region in self.equilibrium.regions.items():
            # all segments must have the same ny, i.e. same number of y-boundary guard
            # cells
            this_ny = region.ny(0)
            assert all(region.ny(i) == this_ny for i in range(region.nSegments))

            y_total_new = y_total + this_ny
            reg_slice = slice(y_total, y_total_new, None)
            y_total = y_total_new
            y_regions[regname] = reg_slice

        self.region_indices = {}
        for reg_name in self.equilibrium.regions:
            for i in range(len(x_regions)):
                self.region_indices[self.region_lookup[(reg_name, i)]] = numpy.index_exp[
                        x_regions[i], y_regions[reg_name]]

        # constant spacing in y for now
        self.dy_scalar = 2.*numpy.pi / self.ny_noguards

    def geometry(self):
        # Call geometry() method of base class
        super().geometry()

        def addFromRegion(f, f_region, regionID):
            if f_region._centre_array is not None:
                f.centre[self.region_indices[regionID]] = f_region.centre
            if f_region._xlow_array is not None:
                f.xlow[self.region_indices[regionID]] = f_region.xlow[:-1,:]
            if f_region._ylow_array is not None:
                f.ylow[self.region_indices[regionID]] = f_region.ylow[:,:-1]
            if f_region._corners_array is not None:
                f.corners[self.region_indices[regionID]] = f_region.corners[:-1,:-1]

        self.Rxy = MultiLocationArray(self.nx, self.ny)
        self.Zxy = MultiLocationArray(self.nx, self.ny)
        self.psixy = MultiLocationArray(self.nx, self.ny)
        self.dx = MultiLocationArray(self.nx, self.ny)
        self.dy = MultiLocationArray(self.nx, self.ny)
        self.Brxy = MultiLocationArray(self.nx, self.ny)
        self.Bzxy = MultiLocationArray(self.nx, self.ny)
        self.Bpxy = MultiLocationArray(self.nx, self.ny)
        self.Btxy = MultiLocationArray(self.nx, self.ny)
        self.Bxy = MultiLocationArray(self.nx, self.ny)
        self.hthe = MultiLocationArray(self.nx, self.ny)
        #if not self.orthogonal:
        #    self.beta = MultiLocationArray(self.nx, self.ny)
        #    self.eta = MultiLocationArray(self.nx, self.ny)
        self.pitch = MultiLocationArray(self.nx, self.ny)
        self.dqdpsi = MultiLocationArray(self.nx, self.ny)

        for region in self.regions.values():
            addFromRegion(self.Rxy, region.Rxy, region.myID)
            addFromRegion(self.Zxy, region.Zxy, region.myID)
            addFromRegion(self.psixy, region.psixy, region.myID)
            addFromRegion(self.dx, region.dx, region.myID)
            addFromRegion(self.dy, region.dy, region.myID)
            addFromRegion(self.Brxy, region.Brxy, region.myID)
            addFromRegion(self.Bzxy, region.Bzxy, region.myID)
            addFromRegion(self.Bpxy, region.Bpxy, region.myID)
            addFromRegion(self.Btxy, region.Btxy, region.myID)
            addFromRegion(self.Bxy, region.Bxy, region.myID)
            addFromRegion(self.hthe, region.hthe, region.myID)
            #if not self.orthogonal:
            #    addFromRegion(self.beta, region.beta, region.myID)
            #    addFromRegion(self.eta, region.eta, region.myID)
            addFromRegion(self.pitch, region.pitch, region.myID)
            addFromRegion(self.dqdpsi, region.dqdpsi, region.myID)

    def writeArray(self, name, array, f):
        f.write(name, array.centre)
        f.write(name+'_ylow', array.ylow[:, :-1])

    def writeGridfile(self, filename):
        from boututils.datafile import DataFile

        with DataFile(filename, create=True) as f:
            f.write('nx', self.nx)
            # ny for BOUT++ excludes boundary guard cells
            f.write('ny', self.ny_noguards)
            f.write('y_boundary_guards', self.equilibrium.y_boundary_guards)
            self.writeArray('Rxy', self.Rxy, f)
            self.writeArray('Zxy', self.Zxy, f)
            self.writeArray('psixy', self.psixy, f)
            self.writeArray('dx', self.dx, f)
            self.writeArray('dy', self.dy, f)
            self.writeArray('Bpxy', self.Bpxy, f)
            self.writeArray('Btxy', self.Btxy, f)
            self.writeArray('Bxy', self.Bxy, f)
            self.writeArray('hthe', self.hthe, f)
            #if not self.orthogonal:
            #    self.writeArray('beta', self.beta, f)
            #    self.writeArray('eta', self.eta, f)
            self.writeArray('pitch', self.pitch, f)

    def plot2D(self, f, title=None):
        from matplotlib import pyplot

        try:
            vmin = f.min()
            vmax = f.max()

            for region, indices in zip(self.regions.values(), self.region_indices.values()):
                pyplot.pcolor(region.Rxy.corners, region.Zxy.corners, f[indices],
                              vmin=vmin, vmax=vmax)

            pyplot.colorbar()
        except NameError:
            raise NameError('Some variable has not been defined yet: have you called Mesh.geometry()?')
