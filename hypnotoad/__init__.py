# Copyright 2019 J.T. Omotani
#
# Contact John Omotani john.omotani@ukaea.uk
#
# This file is part of Hypnotoad 2.
#
# Hypnotoad 2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Hypnotoad 2 is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# Hypnotoad 2.  If not, see <http://www.gnu.org/licenses/>.

from .cases import tokamak, torpex
from .core.equilibrium import Point2D, EquilibriumRegion, Equilibrium, SolutionError
from .core.mesh import MultiLocationArray, MeshRegion, Mesh, BoutMesh
from .__version__ import get_versions

__version__ = get_versions()["version"]

__all__ = [
    "tokamak",
    "torpex",
    "Point2D",
    "EquilibriumRegion",
    "Equilibrium",
    "SolutionError",
    "MultiLocationArray",
    "MeshRegion",
    "Mesh",
    "BoutMesh",
    "__version__",
]
