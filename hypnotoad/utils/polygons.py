"""
Routines for geometric calculations with polygons

Copyright 2019 Ben Dudson, University of York. Email: benjamin.dudson@york.ac.uk

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


def intersect(r1, z1, r2, z2, closed1=True, closed2=True):
    """Test if two polynomials intersect. The polynomials consist of
    (r1, z1) and (r2, z2) line segments. All inputs are expected to be lists.

    Returns True or False.
    """

    assert len(r1) == len(z1)
    assert len(r2) == len(z2)

    n1 = len(r1) if closed1 else len(r1) - 1
    n2 = len(r2) if closed2 else len(r2) - 1

    for i in range(n1):
        for j in range(n2):
            # Test for intersection between two line segments:
            # (r1[i],z1[i]) -- (r1[i+1],z1[i+1])
            # (r1[j],z1[j]) -- (r1[j+1],z1[j+1])
            # Note that since polynomials are closed the indices wrap around
            ip = (i + 1) % n1
            jp = (j + 1) % n2

            a = r1[ip] - r1[i]
            b = r2[jp] - r2[j]
            c = z1[ip] - z1[i]
            d = z2[jp] - z2[j]

            dr = r2[jp] - r1[i]
            dz = z2[jp] - z1[i]

            det = a * d - b * c

            if abs(det) < 1e-6:
                continue  # Almost certainly doesn't intersect

            alpha = (d * dr - b * dz) / det  # Location along line 1 [0,1]
            beta = (a * dz - c * dr) / det  # Location along line 2 [0,1]

            if (alpha > 0.0) & (alpha < 1.0) & (beta > 0.0) & (beta < 1.0):
                return True
    # Tested all combinations, none intersect
    return False


def area(polygon):
    """
    Calculate the area of a polygon. Can be positive (clockwise) or negative
    (anticlockwise)

    Input

    polygon   [ (r1, z1), (r2, z2), ... ]
    """
    nvert = len(polygon)  # Number of vertices

    # Integrate using trapezium rule. The sign of (r2-r1) ensures that
    # positive and negative areas leave only the area of the polygon.
    area = 0.0
    for i in range(nvert):
        r1, z1 = polygon[i]
        r2, z2 = polygon[(i + 1) % nvert]  # Next vertex in periodic list
        area += (r2 - r1) * (z1 + z2)  # 2*area
    return 0.5 * area


def clockwise(polygon):
    """
    Detect whether a polygon is clockwise or anti-clockwise
    True -> clockwise
    False -> anticlockwise

    Input

    polygon   [ (r1, z1), (r2, z2), ... ]
    """
    # Work out the winding direction by calculating the area
    return area(polygon) > 0
