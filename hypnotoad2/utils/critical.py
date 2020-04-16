"""
Routines to find critical points (O- and X-points)

Copyright 2016 Ben Dudson, University of York. Email: benjamin.dudson@york.ac.uk

This file is part of FreeGS.

FreeGS is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FreeGS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with FreeGS.  If not, see <http://www.gnu.org/licenses/>.

"""


from scipy import interpolate
from numpy.linalg import inv
from numpy import (
    dot,
    linspace,
    argmin,
    abs,
    clip,
    amax,
    zeros,
)
import numpy as np


def find_critical(R, Z, psi, discard_xpoints=False):
    """
    Find critical points

    Inputs
    ------

    R - R(nr, nz) 2D array of major radii
    Z - Z(nr, nz) 2D array of heights
    psi - psi(nr, nz) 2D array of psi values

    Returns
    -------

    Two lists of critical points

    opoint, xpoint

    Each of these is a list of tuples with (R, Z, psi) points

    The first tuple is the primary O-point (magnetic axis)
    and primary X-point (separatrix)

    """

    # Get a spline interpolation function
    f = interpolate.RectBivariateSpline(R[:, 0], Z[0, :], psi)

    # Find candidate locations, based on minimising Bp^2
    Bp2 = (f(R, Z, dx=1, grid=False) ** 2 + f(R, Z, dy=1, grid=False) ** 2) / R ** 2

    # Get grid resolution, which determines a reasonable tolerance
    # for the Newton iteration search area
    dR = R[1, 0] - R[0, 0]
    dZ = Z[0, 1] - Z[0, 0]
    radius_sq = 9 * (dR ** 2 + dZ ** 2)

    # Find local minima

    J = zeros([2, 2])

    xpoint = []
    opoint = []

    nx, ny = Bp2.shape
    for i in range(2, nx - 2):
        for j in range(2, ny - 2):
            if (
                (Bp2[i, j] < Bp2[i + 1, j + 1])
                and (Bp2[i, j] < Bp2[i + 1, j])
                and (Bp2[i, j] < Bp2[i + 1, j - 1])
                and (Bp2[i, j] < Bp2[i - 1, j + 1])
                and (Bp2[i, j] < Bp2[i - 1, j])
                and (Bp2[i, j] < Bp2[i - 1, j - 1])
                and (Bp2[i, j] < Bp2[i, j + 1])
                and (Bp2[i, j] < Bp2[i, j - 1])
            ):

                # Found local minimum

                R0 = R[i, j]
                Z0 = Z[i, j]

                # Use Newton iterations to find where
                # both Br and Bz vanish
                R1 = R0
                Z1 = Z0

                count = 0
                while True:

                    Br = -f(R1, Z1, dy=1, grid=False) / R1
                    Bz = f(R1, Z1, dx=1, grid=False) / R1

                    if Br ** 2 + Bz ** 2 < 1e-6:
                        # Found a minimum. Classify as either
                        # O-point or X-point

                        # Evaluate D = fxx * fyy - (fxy)^2
                        if False:
                            D = (
                                f(R1, Z1, dx=2)[0][0] * f(R1, Z1, dy=2)[0][0]
                                - (f(R1, Z1, dx=1, dy=1)[0][0]) ** 2
                            )

                            # print("D0 = %e" % D)

                        if False:  # abs(D) < 1:
                            # D small, so need to use another method
                            # print("Small discriminant D (%e)" % (D,))

                            # Try second derivative in index space

                            dR = R[1, 0] - R[0, 0]
                            dZ = Z[0, 1] - Z[0, 0]
                            d2dr2 = (
                                psi[i + 1, j] - 2.0 * psi[i, j] + psi[i - 1, j]
                            ) / dR ** 2
                            d2dz2 = (
                                psi[i, j + 1] - 2.0 * psi[i, j] + psi[i, j - 1]
                            ) / dZ ** 2
                            d2drdz = (
                                (psi[i + 1, j + 1] - psi[i + 1, j - 1]) / (2.0 * dZ)
                                - (psi[i - 1, j + 1] - psi[i - 1, j - 1]) / (2.0 * dZ)
                            ) / (2.0 * dR)
                            D = d2dr2 * d2dz2 - d2drdz ** 2

                            # print("D1 = %e" % D)

                        if True:
                            dR = R[1, 0] - R[0, 0]
                            dZ = Z[0, 1] - Z[0, 0]
                            d2dr2 = (
                                psi[i + 2, j] - 2.0 * psi[i, j] + psi[i - 2, j]
                            ) / (2.0 * dR) ** 2
                            d2dz2 = (
                                psi[i, j + 2] - 2.0 * psi[i, j] + psi[i, j - 2]
                            ) / (2.0 * dZ) ** 2
                            d2drdz = (
                                (psi[i + 2, j + 2] - psi[i + 2, j - 2]) / (4.0 * dZ)
                                - (psi[i - 2, j + 2] - psi[i - 2, j - 2]) / (4.0 * dZ)
                            ) / (4.0 * dR)
                            D = d2dr2 * d2dz2 - d2drdz ** 2

                            # print("D2 = %e" % D)

                        if D < 0.0:
                            # Found X-point
                            # print("Found X-point at %e, %e (f=%e, D=%e)"
                            #       % (R1,Z1, f(R1,Z1)[0][0], D) )
                            xpoint.append((R1, Z1, f(R1, Z1)[0][0]))
                        else:
                            # Found O-point
                            # print("Found O-point at %e, %e (f=%e, D=%e)"
                            #       % (R1,Z1, f(R1,Z1)[0][0], D) )
                            opoint.append((R1, Z1, f(R1, Z1)[0][0]))
                        break

                    # Jacobian matrix
                    # J = ( dBr/dR, dBr/dZ )
                    #     ( dBz/dR, dBz/dZ )

                    J[0, 0] = -Br / R1 - f(R1, Z1, dy=1, dx=1)[0][0] / R1
                    J[0, 1] = -f(R1, Z1, dy=2)[0][0] / R1
                    J[1, 0] = -Bz / R1 + f(R1, Z1, dx=2) / R1
                    J[1, 1] = f(R1, Z1, dx=1, dy=1)[0][0] / R1

                    d = dot(inv(J), [Br, Bz])

                    R1 = R1 - d[0]
                    Z1 = Z1 - d[1]

                    count += 1
                    # If (R1,Z1) is too far from (R0,Z0) then discard
                    # or if we've taken too many iterations
                    if ((R1 - R0) ** 2 + (Z1 - Z0) ** 2 > radius_sq) or (count > 100):
                        # Discard this point
                        break

    # Remove duplicates
    def remove_dup(points):
        result = []
        for n, p in enumerate(points):
            dup = False
            for p2 in result:
                if (p[0] - p2[0]) ** 2 + (p[1] - p2[1]) ** 2 < 1e-5:
                    dup = True  # Duplicate
                    break
            if not dup:
                result.append(p)  # Add to the list
        return result

    xpoint = remove_dup(xpoint)
    opoint = remove_dup(opoint)

    if len(opoint) == 0:
        # Can't order primary O-point, X-point so return
        print("Warning: No O points found")
        return opoint, xpoint

    # Find primary O-point by sorting by distance from middle of domain
    Rmid = 0.5 * (R[-1, 0] + R[0, 0])
    Zmid = 0.5 * (Z[0, -1] + Z[0, 0])
    opoint.sort(key=lambda x: (x[0] - Rmid) ** 2 + (x[1] - Zmid) ** 2)

    # Draw a line from the O-point to each X-point. Psi should be
    # monotonic; discard those which are not

    if True:  # discard_xpoints:
        Ro, Zo, Po = opoint[0]  # The primary O-point
        xpt_keep = []
        for xpt in xpoint:
            Rx, Zx, Px = xpt

            rline = linspace(Ro, Rx, num=50)
            zline = linspace(Zo, Zx, num=50)

            pline = f(rline, zline, grid=False)

            if Px < Po:
                pline *= -1.0  # Reverse, so pline is maximum at X-point

            # Now check that pline is monotonic
            # Tried finding maximum (argmax) and testing
            # how far that is from the X-point. This can go
            # wrong because psi can be quite flat near the X-point
            # Instead here look for the difference in psi
            # rather than the distance in space

            maxp = amax(pline)
            if (maxp - pline[-1]) / (maxp - pline[0]) > 0.001:
                # More than 0.1% drop in psi from maximum to X-point
                # -> Discard
                continue

            ind = argmin(pline)  # Should be at O-point
            if (rline[ind] - Ro) ** 2 + (zline[ind] - Zo) ** 2 > 1e-4:
                # Too far, discard
                continue
            xpt_keep.append(xpt)
        xpoint = xpt_keep

    # Sort X-points by distance to primary O-point in psi space
    psi_axis = opoint[0][2]
    xpoint.sort(key=lambda x: (x[2] - psi_axis) ** 2)

    return opoint, xpoint


def find_psisurface(eq, r0, z0, r1, z1, psival=1.0, n=100, axis=None):
    """
    eq      - Equilibrium object
    (r0,z0) - Start location inside separatrix
    (r1,z1) - Location outside separatrix

    n - Number of starting points to use
    """
    # Clip (r1,z1) to be inside domain
    # Shorten the line so that the direction is unchanged
    if abs(r1 - r0) > 1e-6:
        rclip = clip(r1, eq.Rmin, eq.Rmax)
        z1 = z0 + (z1 - z0) * abs((rclip - r0) / (r1 - r0))
        r1 = rclip

    if abs(z1 - z0) > 1e-6:
        zclip = clip(z1, eq.Zmin, eq.Zmax)
        r1 = r0 + (r1 - r0) * abs((zclip - z0) / (z1 - z0))
        z1 = zclip

    r = linspace(r0, r1, n)
    z = linspace(z0, z1, n)

    if axis is not None:
        axis.plot(r, z)

    psidiff = eq.psi(r, z) - psival

    # Find the first index this crosses zero
    ind = np.argmax(psidiff[1:] * psidiff[0:-1] < 0.0)
    # between ind and ind-1

    f = psidiff[ind] / (psidiff[ind] - psidiff[ind - 1])

    r = (1.0 - f) * r[ind] + f * r[ind - 1]
    z = (1.0 - f) * z[ind] + f * z[ind - 1]

    if axis is not None:
        axis.plot(r, z, "bo")

    return r, z
