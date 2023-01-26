#!/usr/bin/env python3

from hypnotoad.cases import tokamak
from hypnotoad.core.equilibrium import Point2D
from hypnotoad.core.mesh import BoutMesh
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path

geqdsk_path = Path(
    __file__,
    "..",
    "..",
    "..",
    "integrated_tests",
    "grid_files",
    "test_connected-double-null.eqdsk",
).resolve()

image_dir = Path(__file__, "..").resolve()

options = {
    "finecontour_atol": 1.0e-10,
    "ny_sol": 16,
    "orthogonal": True,
    "psinorm_core": 0.8,
    "psinorm_pf": 0.95,
    "psinorm_sol": 2.0,
    "psinorm_sol_inner": 1.03,
    "target_all_poloidal_spacing_length": 1.5,
    "xpoint_poloidal_spacing_length": 0.2,
    "y_boundary_guards": 2,
}

with open(geqdsk_path, "rt") as fh:
    eq = tokamak.read_geqdsk(fh, settings=options, nonorthogonal_settings=options)

mesh = BoutMesh(eq, options)
mesh.calculateRZ()


# Plot a PsiContour
plt.figure(figsize=(4, 3), constrained_layout=True)

mreg = mesh.regions[mesh.region_lookup[("outer_lower_divertor", 1)]]
contour = mreg.contours[7]

eq.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
contour.plot(psi=eq.psi, linestyle="", marker="x", color="red")
plt.xlim([0.0, 1.4])
plt.ylim([-2.0, -1.0])
plt.savefig(image_dir.joinpath("cdn-PsiContour.svg"))


# Plot a FineContour
plt.figure(figsize=(4, 8), constrained_layout=True)

finecontour = contour.get_fine_contour(psi=eq.psi)

eq.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
contour.plot(psi=eq.psi, linestyle="", marker="x", color="red")
finecontour.plot(psi=eq.psi, linestyle="", marker=".", markersize=1, color="blue")
plt.xlim([0.8, 1.0])
plt.ylim([-2.0, -1.5])
plt.savefig(image_dir.joinpath("cdn-FineContour.svg"))


# Plot an EquilibriumRegion
plt.figure(figsize=(4, 8), constrained_layout=True)

region_name = "outer_core"
eqreg = eq.regions[region_name]
reg_number = 0
for reg in eq.regions.keys():
    if reg == region_name:
        break
    else:
        reg_number = reg_number + 1

ax = eq.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
eqreg.plot(psi=eq.psi, linestyle="", marker="+", markersize=10, color=f"C{reg_number}")
plt.savefig(image_dir.joinpath("cdn-EquilibriumRegion.svg"))


# Plot the Equilibrium
plt.figure(figsize=(4, 8), constrained_layout=True)

ax = eq.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
for eqreg in eq.regions.values():
    eqreg.plot(psi=eq.psi, linestyle="", marker="+", markersize=10)
plt.savefig(image_dir.joinpath("cdn-Equilibrium.svg"))


# Plot a MeshRegion
plt.figure(figsize=(4, 8), constrained_layout=True)

mregion_name = ("inner_core", 0)
mreg_number = mesh.region_lookup[mregion_name]
mreg = mesh.regions[mreg_number]

ax = eq.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
ax.scatter(
    mreg.Rxy.centre,
    mreg.Zxy.centre,
    marker="x",
    color=f"C{mreg_number}",
)
plt.savefig(image_dir.joinpath("cdn-MeshRegion.svg"))


# Plot the Mesh
plt.figure(figsize=(4, 8), constrained_layout=True)
ax = mesh.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
mesh.plotPoints(ax=ax, legend=False)
plt.savefig(image_dir.joinpath("cdn-Mesh.svg"))

# Plot critical points
plt.figure(figsize=(4, 8), constrained_layout=True)

ax = eq.plotPotential(labels=False, colors="gray", linestyles="solid")
eq.plotWall()
for p in eq.x_points:
    plt.scatter(*p, marker="x", c="red", s=200, zorder=10, linewidth=3)
plt.scatter(*eq.o_point, marker="+", c="blue", s=200, zorder=10, linewidth=3)
plt.savefig(image_dir.joinpath("cdn-critical-points.svg"))


# Nonorthogonal figures
nonorth_options = {
    "finecontour_atol": 1.0e-10,
    "ny_sol": 16,
    "orthogonal": False,
    "psinorm_core": 0.8,
    "psinorm_pf": 0.95,
    "psinorm_sol": 2.0,
    "psinorm_sol_inner": 1.03,
    "target_all_poloidal_spacing_length": 1.5,
    "y_boundary_guards": 2,
    "nonorthogonal_radial_range_power": 4,
    "nonorthogonal_target_all_poloidal_spacing_range": 0.05,
    "nonorthogonal_target_all_poloidal_spacing_range_inner": 0.1,
    "nonorthogonal_target_all_poloidal_spacing_range_outer": 0.05,
    "nonorthogonal_xpoint_poloidal_spacing_range_inner": 0.05,
    "nonorthogonal_xpoint_poloidal_spacing_range_outer": 0.02,
}

with open(geqdsk_path, "rt") as fh:
    nonorth_eq = tokamak.read_geqdsk(
        fh, settings=nonorth_options, nonorthogonal_settings=nonorth_options
    )

nonorth_mesh = BoutMesh(nonorth_eq, nonorth_options)
nonorth_mesh.calculateRZ()

# Plot the Mesh
plt.figure(figsize=(4, 8), constrained_layout=True)
ax = nonorth_mesh.plotPotential(labels=False, colors="gray", linestyles="solid")
nonorth_eq.plotWall()
nonorth_mesh.plotPoints(ax=ax, legend=False)
plt.savefig(image_dir.joinpath("cdn-nonorth-Mesh.svg"))


# Nonorthogonal spacing function illustrations
# Standard, 'combined' poloidal spacing
nonorth_options = {
    "nx_core": 20,
    "nx_sol": 20,
    "ny_inner_divertor": 8,
    "ny_outer_divertor": 8,
    "ny_sol": 64,
    "orthogonal": False,
    "psinorm_core": 0.8,
    "psinorm_pf": 0.95,
    "psinorm_sol": 2.0,
    "psinorm_sol_inner": 1.03,
    "target_all_poloidal_spacing_length": 1.5,
    "nonorthogonal_radial_range_power": 4,
    "nonorthogonal_target_all_poloidal_spacing_range": 0.05,
    "nonorthogonal_target_all_poloidal_spacing_range_inner": 0.1,
    "nonorthogonal_target_all_poloidal_spacing_range_outer": 0.2,
    "nonorthogonal_xpoint_poloidal_spacing_range_inner": 0.05,
    "nonorthogonal_xpoint_poloidal_spacing_range_outer": 0.02,
}

with open(geqdsk_path, "rt") as fh:
    nonorth_eq = tokamak.read_geqdsk(
        fh, settings=nonorth_options, nonorthogonal_settings=nonorth_options
    )

# Hack to add a funkier wall to make the the examples nicer
nonorth_eq.wall = [
    Point2D(0.1, -1.8),
    Point2D(0.85, -1.8),
    Point2D(0.95, -1.7),
    Point2D(1.67, -1.7),
    Point2D(1.67, 1.7),
    Point2D(0.95, 1.7),
    Point2D(0.85, 1.8),
    Point2D(0.1, 1.8),
]
closed_wall = nonorth_eq.wall + [nonorth_eq.wall[0]]
nonorth_eq.closed_wallarray = np.array([(p.R, p.Z) for p in closed_wall])

nonorth_mesh = BoutMesh(nonorth_eq, nonorth_options)
nonorth_mesh.calculateRZ()

# Plot the Mesh
plt.figure(figsize=(4, 3), constrained_layout=True)
ax = nonorth_mesh.plotPotential(labels=False, colors="gray", linestyles="solid")
ax.set_xlim((0.5, 1.1))
ax.set_ylim((-1.85, -1.4))
ax.set_aspect("equal")
nonorth_eq.plotWall()
nonorth_mesh.plotPoints(ax=ax, legend=False, plot_types="radial")
plt.savefig(image_dir.joinpath("cdn-nonorth-combined.svg"))

# 'Target' fixed-poloidal spacing
nonorth_options = {
    "finecontour_atol": 1.0e-10,
    "nx_core": 20,
    "nx_sol": 20,
    "ny_inner_divertor": 8,
    "ny_outer_divertor": 8,
    "ny_sol": 64,
    "orthogonal": False,
    "psinorm_core": 0.8,
    "psinorm_pf": 0.95,
    "psinorm_sol": 2.0,
    "psinorm_sol_inner": 1.03,
    "target_all_poloidal_spacing_length": 1.5,
    "nonorthogonal_radial_range_power": 4,
    "nonorthogonal_target_all_poloidal_spacing_range": 1000,
    "nonorthogonal_xpoint_poloidal_spacing_range": 0.01,
    "nonorthogonal_xpoint_poloidal_spacing_range_inner": 0.025,
    "nonorthogonal_xpoint_poloidal_spacing_range_outer": 0.01,
}

with open(geqdsk_path, "rt") as fh:
    nonorth_eq = tokamak.read_geqdsk(
        fh, settings=nonorth_options, nonorthogonal_settings=nonorth_options
    )

# Hack to add a funkier wall to make the the examples nicer
nonorth_eq.wall = [
    Point2D(0.1, -1.8),
    Point2D(0.85, -1.8),
    Point2D(0.95, -1.7),
    Point2D(1.67, -1.7),
    Point2D(1.67, 1.7),
    Point2D(0.95, 1.7),
    Point2D(0.85, 1.8),
    Point2D(0.1, 1.8),
]
closed_wall = nonorth_eq.wall + [nonorth_eq.wall[0]]
nonorth_eq.closed_wallarray = np.array([(p.R, p.Z) for p in closed_wall])

nonorth_mesh = BoutMesh(nonorth_eq, nonorth_options)
nonorth_mesh.calculateRZ()

# Plot the Mesh
plt.figure(figsize=(4, 3), constrained_layout=True)
ax = nonorth_mesh.plotPotential(labels=False, colors="gray", linestyles="solid")
ax.set_xlim((0.5, 1.1))
ax.set_ylim((-1.85, -1.4))
ax.set_aspect("equal")
nonorth_eq.plotWall()
nonorth_mesh.plotPoints(ax=ax, legend=False, plot_types="radial")
plt.savefig(image_dir.joinpath("cdn-nonorth-target-spacing.svg"))

# 'X-point' fixed-perpendicular spacing
nonorth_options = {
    "finecontour_atol": 1.0e-10,
    "nx_core": 20,
    "nx_sol": 20,
    "ny_inner_divertor": 8,
    "ny_outer_divertor": 8,
    "ny_sol": 64,
    "orthogonal": False,
    "psinorm_core": 0.8,
    "psinorm_pf": 0.95,
    "psinorm_sol": 2.0,
    "psinorm_sol_inner": 1.03,
    "target_all_poloidal_spacing_length": 1.5,
    "nonorthogonal_radial_range_power": 4,
    "nonorthogonal_target_all_poloidal_spacing_range": 0.01,
    "nonorthogonal_xpoint_poloidal_spacing_range": 100.0,
    "nonorthogonal_xpoint_poloidal_spacing_range_inner": 100.0,
    "nonorthogonal_xpoint_poloidal_spacing_range_outer": 100.0,
}

with open(geqdsk_path, "rt") as fh:
    nonorth_eq = tokamak.read_geqdsk(
        fh, settings=nonorth_options, nonorthogonal_settings=nonorth_options
    )

# Hack to add a funkier wall to make the the examples nicer
nonorth_eq.wall = [
    Point2D(0.1, -1.8),
    Point2D(0.85, -1.8),
    Point2D(0.95, -1.7),
    Point2D(1.67, -1.7),
    Point2D(1.67, 1.7),
    Point2D(0.95, 1.7),
    Point2D(0.85, 1.8),
    Point2D(0.1, 1.8),
]
closed_wall = nonorth_eq.wall + [nonorth_eq.wall[0]]
nonorth_eq.closed_wallarray = np.array([(p.R, p.Z) for p in closed_wall])

nonorth_mesh = BoutMesh(nonorth_eq, nonorth_options)
nonorth_mesh.calculateRZ()

# Plot the Mesh
plt.figure(figsize=(4, 3), constrained_layout=True)
ax = nonorth_mesh.plotPotential(labels=False, colors="gray", linestyles="solid")
ax.set_xlim((0.5, 1.1))
ax.set_ylim((-1.85, -1.4))
ax.set_aspect("equal")
nonorth_eq.plotWall()
nonorth_mesh.plotPoints(ax=ax, legend=False, plot_types="radial")
plt.savefig(image_dir.joinpath("cdn-nonorth-xpoint-spacing.svg"))
