"""
Routines for geometric calculations with polygons

License
-------

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
