
#%%
# %load_ext autoreload
# %autoreload 2


import os
import hypnotoad
from hypnotoad.cases.dipole import DipoleEquilibrium
from hypnotoad.core.mesh import BoutMesh
#%%

options = {}



options["dpeqfile"] = os.path.join(hypnotoad.__path__[0],"../examples/dipole/dipole_eq.h5")
options["y_boundary_guards"] = 0
options["ny"] = 32
options["nx"] = 34
options["psi_spacing_separatrix_multiplier"] =  0.1   # Smaller -> pack near separatrix

# Poloidal grid spacing

options["target_all_poloidal_spacing_length"] = 0.3   # Smaller -> pack near targets
options["xpoint_poloidal_spacing_length"]= 0.25
#options["psi_inner"] = -2.0
# xpoint_poloidal_spacing_length
#
# finecontour_Nfine: 500
# finecontour_atol: 5.0e-16
# finecontour_diagnose: false
# finecontour_extend_prefactor: 2.0
# finecontour_maxits: 1000
# finecontour_overdamping_factor: 0.8
# follow_perpendicular_atol: 1.0e-20
# follow_perpendicular_rtol: 1.0e-16
# geometry_rtol: 1.0e-10
# leg_refine_atol: 2.0e-16
# leg_refine_maxits: 10000
options["nonorthogonal_radial_range_power"]: 2
options["nonorthogonal_spacing_method"]: combined
# nonorthogonal_target_all_poloidal_spacing_length: 1.0
# nonorthogonal_target_all_poloidal_spacing_range: 0.05
# nonorthogonal_target_all_poloidal_spacing_range_inner: 0.05
# nonorthogonal_target_all_poloidal_spacing_range_outer: 0.05
# nonorthogonal_target_inner_lower_poloidal_spacing_length: 1.0
# nonorthogonal_target_inner_lower_poloidal_spacing_range: 0.05
# nonorthogonal_target_inner_lower_poloidal_spacing_range_inner: 0.05
# nonorthogonal_target_inner_lower_poloidal_spacing_range_outer: 0.05
# nonorthogonal_target_inner_upper_poloidal_spacing_length: 1.0
# nonorthogonal_target_inner_upper_poloidal_spacing_range: 0.05
# nonorthogonal_target_inner_upper_poloidal_spacing_range_inner: 0.05
# nonorthogonal_target_inner_upper_poloidal_spacing_range_outer: 0.05
# nonorthogonal_target_outer_lower_poloidal_spacing_length: 1.0
# nonorthogonal_target_outer_lower_poloidal_spacing_range: 0.05
# nonorthogonal_target_outer_lower_poloidal_spacing_range_inner: 0.05
# nonorthogonal_target_outer_lower_poloidal_spacing_range_outer: 0.05
# nonorthogonal_target_outer_upper_poloidal_spacing_length: 1.0
# nonorthogonal_target_outer_upper_poloidal_spacing_range: 0.05
# nonorthogonal_target_outer_upper_poloidal_spacing_range_inner: 0.05
# nonorthogonal_target_outer_upper_poloidal_spacing_range_outer: 0.05
# nonorthogonal_xpoint_poloidal_spacing_length: 1.0
# nonorthogonal_xpoint_poloidal_spacing_range: 0.015
# nonorthogonal_xpoint_poloidal_spacing_range_inner: 0.05
# nonorthogonal_xpoint_poloidal_spacing_range_outer: 0.02

eq = DipoleEquilibrium(settings=options, nonorthogonal_settings=options)

import matplotlib.pyplot as plt

eq.plotPotential(ncontours=40)
for region in eq.regions.values():
    plt.plot(
        [p.R for p in region.points],
        [p.Z for p in region.points],
        marker="o",
    )
plt.show()

mesh = BoutMesh(eq, options)
mesh.calculateRZ()


import matplotlib.pyplot as plt

ax = eq.plotPotential(ncontours=40)
mesh.plotPoints(
    xlow=options.get("plot_xlow", True),
    ylow=options.get("plot_ylow", True),
    corners=options.get("plot_corners", True),
    ax=ax,
)
plt.show()


mesh.geometry()
mesh.writeGridfile(options.get("grid_file", "bout.dipole.grd.nc"))
