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
    "eqdsk_files",
    "SPR45_IXD1_AH-JTOmod-add-wall.geqdsk"
).resolve()

image_dir = Path(__file__, "..").resolve()

options = {"psinorm_core": 0.99,
    # "finecontour_atol": 1.0e-10,
    # "ny_sol": 16,
    # "orthogonal": True,
    # "psinorm_pf": 0.99,
    # "psinorm_pf_lower": 0.99,

    "psinorm_sol": 1.01,
#     "psinorm_sol_inner": 1.01,
#     "target_all_poloidal_spacing_length": 1.5,
#     "xpoint_poloidal_spacing_length": 0.2,
#     "y_boundary_guards": 2,
}

with open(geqdsk_path, "rt") as fh:
    eq = tokamak.read_geqdsk(fh, settings=options, nonorthogonal_settings=options)
print("eq object is: ")
print(eq)
print("making boutmesh...")
mesh = BoutMesh(eq, options)
mesh.calculateRZ()
print("made boutmesh")

# mreg_upper_coarse = mesh.regions[mesh.region_lookup[("outer_upper_divertor", 1)]]
# mreg_lower_coarse = mesh.regions[mesh.region_lookup[("outer_lower_divertor", 1)]]
# mreg_core_coarse  = mesh.regions[mesh.region_lookup[("outer_core", 1)]]
# contour_upper_coarse = mreg_upper_coarse.contours[1]
# contour_lower_coarse = mreg_lower_coarse.contours[1]
# contour_core_coarse  = mreg_core_coarse.contours[1]

# fc_upper_coarse = contour_upper_coarse.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=20)
# fc_lower_coarse = contour_lower_coarse.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=20)
# fc_core_coarse  = contour_core_coarse.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=20)
# print("total distance parallel upper COARSE: ", fc_upper_coarse.parallel_distance[-1])
# print("total distance parallel lower COARSE: ", fc_lower_coarse.parallel_distance[-1])
# print("total distance parallel core  COARSE: ", fc_core_coarse.parallel_distance[-1])
# length_ratio = fc_core_coarse.parallel_distance[-1] / fc_upper_coarse.parallel_distance[-1]
# Nfine_core = 500
# Nfine_upper = Nfine_core // length_ratio
# Nfine_lower = Nfine_upper
# print("Nfine upper: ", Nfine_upper)
# print("Nfine lower: ", Nfine_lower)
# print("Nfine core: ", Nfine_core)

mreg_upper = mesh.regions[mesh.region_lookup[("outer_upper_divertor", 1)]]
mreg_lower = mesh.regions[mesh.region_lookup[("outer_lower_divertor", 1)]]
mreg_core  = mesh.regions[mesh.region_lookup[("outer_core", 1)]]

print("NUMBER OF UPPER, CORE AND LOWER CONTOURS: ", len(mreg_upper.contours), len(mreg_core.contours), len(mreg_lower.contours))
contour_upper_0 = mreg_upper.contours[0]
contour_lower_0 = mreg_lower.contours[0]
contour_core_0  = mreg_core.contours[0]
contour_upper_1 = mreg_upper.contours[1]
contour_lower_1 = mreg_lower.contours[1]
contour_core_1  = mreg_core.contours[1]
contour_upper_2 = mreg_upper.contours[2]
contour_lower_2 = mreg_lower.contours[2]
contour_core_2  = mreg_core.contours[2]
contour_upper_3 = mreg_upper.contours[3]
contour_lower_3 = mreg_lower.contours[3]
contour_core_3  = mreg_core.contours[3]
contour_upper_end = mreg_upper.contours[10]
contour_lower_end = mreg_lower.contours[10]
contour_core_end = mreg_core.contours[10]


# THESE HAVE THE OLD Nfine VALUES HARDCODED TO MAKE PARALLEL DISTNACES MATCH BETWEEN REGIONS ROUGHLY
finecontour_upper_0 = contour_upper_0.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=100)
finecontour_lower_0 = contour_lower_0.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=100)
finecontour_core_0  = contour_core_0.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=139)
finecontour_upper_1 = contour_upper_1.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=100)
finecontour_lower_1 = contour_lower_1.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=100)
finecontour_core_1  = contour_core_1.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=231)
finecontour_upper_2 = contour_upper_2.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=75)
finecontour_lower_2 = contour_lower_2.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=75)
finecontour_core_2  = contour_core_2.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=180)
finecontour_upper_3 = contour_upper_3.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=75)
finecontour_lower_3 = contour_lower_3.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=75)
finecontour_core_3  = contour_core_3.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=184)

# finecontour_upper_0 = contour_upper_0.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1000)
# finecontour_lower_0 = contour_lower_0.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1000)
# finecontour_core_0  = contour_core_0.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1390)
# finecontour_upper_1 = contour_upper_1.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1000)
# finecontour_lower_1 = contour_lower_1.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1000)
# finecontour_core_1  = contour_core_1.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=2310)
# finecontour_upper_2 = contour_upper_2.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=750)
# finecontour_lower_2 = contour_lower_2.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=750)
# finecontour_core_2  = contour_core_2.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1800)
# finecontour_upper_3 = contour_upper_3.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=750)
# finecontour_lower_3 = contour_lower_3.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=750)
# finecontour_core_3  = contour_core_3.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=1840)
finecontour_upper_end = contour_upper_end.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=100)
finecontour_lower_end = contour_lower_end.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=100)
finecontour_core_end  = contour_core_end.get_fine_contour_parallel(psi=eq.psi, equilibrium=eq, Nfine=271)

print("finecontour_core_0.positions: ", finecontour_core_0.positions)
B_values_core_0  = np.sqrt(eq.B2(finecontour_core_0.positions[:,0],  finecontour_core_0.positions[:,1]))
print("B_values_core: ", B_values_core_0)
B_values_upper_0 = np.sqrt(eq.B2(finecontour_upper_0.positions[:,0], finecontour_upper_0.positions[:,1]))
print("B_values_upper: ", B_values_upper_0)
B_values_lower_0 = np.sqrt(eq.B2(finecontour_lower_0.positions[:,0], finecontour_lower_0.positions[:,1]))
print("B_values_lower: ", B_values_lower_0)
B_concatenated_0 = np.concatenate((B_values_upper_0, B_values_core_0, B_values_lower_0))
print("B_concatenated: ", B_concatenated_0)
print("finecontour_core_1.positions: ", finecontour_core_1.positions)
B_values_core_1  = np.sqrt(eq.B2(finecontour_core_1.positions[:,0],  finecontour_core_1.positions[:,1]))
print("B_values_core: ", B_values_core_1)
B_values_upper_1 = np.sqrt(eq.B2(finecontour_upper_1.positions[:,0], finecontour_upper_1.positions[:,1]))
print("B_values_upper: ", B_values_upper_1)
B_values_lower_1 = np.sqrt(eq.B2(finecontour_lower_1.positions[:,0], finecontour_lower_1.positions[:,1]))
print("B_values_lower: ", B_values_lower_1)
B_concatenated_1 = np.concatenate((B_values_upper_1, B_values_core_1, B_values_lower_1))
print("B_concatenated: ", B_concatenated_1)
print("finecontour_core_2.positions: ", finecontour_core_2.positions)
B_values_core_2  = np.sqrt(eq.B2(finecontour_core_2.positions[:,0],  finecontour_core_2.positions[:,1]))
print("B_values_core: ", B_values_core_2)
B_values_upper_2 = np.sqrt(eq.B2(finecontour_upper_2.positions[:,0], finecontour_upper_2.positions[:,1]))
print("B_values_upper: ", B_values_upper_2)
B_values_lower_2 = np.sqrt(eq.B2(finecontour_lower_2.positions[:,0], finecontour_lower_2.positions[:,1]))
print("B_values_lower: ", B_values_lower_2)
B_concatenated_2 = np.concatenate((B_values_upper_2, B_values_core_2, B_values_lower_2))
print("B_concatenated: ", B_concatenated_2)
print("finecontour_core_3.positions: ", finecontour_core_3.positions)
B_values_core_3  = np.sqrt(eq.B2(finecontour_core_3.positions[:,0],  finecontour_core_3.positions[:,1]))
print("B_values_core: ", B_values_core_3)
B_values_upper_3 = np.sqrt(eq.B2(finecontour_upper_3.positions[:,0], finecontour_upper_3.positions[:,1]))
print("B_values_upper: ", B_values_upper_3)
B_values_lower_3 = np.sqrt(eq.B2(finecontour_lower_3.positions[:,0], finecontour_lower_3.positions[:,1]))
print("B_values_lower: ", B_values_lower_3)
B_concatenated_3 = np.concatenate((B_values_upper_3, B_values_core_3, B_values_lower_3))
print("B_concatenated: ", B_concatenated_3)
print("finecontour_core_end.positions: ", finecontour_core_end.positions)
B_values_core_end  = np.sqrt(eq.B2(finecontour_core_end.positions[:,0],  finecontour_core_end.positions[:,1]))
print("B_values_core: ", B_values_core_end)
B_values_upper_end = np.sqrt(eq.B2(finecontour_upper_end.positions[:,0], finecontour_upper_end.positions[:,1]))
print("B_values_upper: ", B_values_upper_end)
B_values_lower_end = np.sqrt(eq.B2(finecontour_lower_end.positions[:,0], finecontour_lower_end.positions[:,1]))
print("B_values_lower: ", B_values_lower_end)
B_concatenated_end = np.concatenate((B_values_upper_end, B_values_core_end, B_values_lower_end))
print("B_concatenated: ", B_concatenated_end)


figwidth = 4.0
figheight = figwidth * (eq.Zmax - eq.Zmin) / (eq.Rmax - eq.Rmin)
fig, ax = plt.subplots(figsize=(figwidth, figheight), constrained_layout=True)


colors = "grey"
eq.plotPotential(
    npoints=50,
    ncontours=100,
    labels=True,
    colors=colors,
    linestyles="-",
)
eq.plotWall(axis=ax)

finecontour_upper_0.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="green")
finecontour_lower_0.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="green")
finecontour_core_0.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="green")
finecontour_upper_1.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="blue")
finecontour_lower_1.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="blue")
finecontour_core_1.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="blue")
finecontour_upper_2.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="red")
finecontour_lower_2.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="red")
finecontour_core_2.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="red")
finecontour_upper_3.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="orange")
finecontour_lower_3.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="orange")
finecontour_core_3.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="orange")
finecontour_upper_end.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="pink")
finecontour_lower_end.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="pink")
finecontour_core_end.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=0.001, color="pink")


plt.savefig(image_dir.joinpath("testing_wholecontourplot.pdf"), bbox_inches="tight")



end = len(finecontour_upper_0.parallel_distance)-1
parallel_gaps_upper = finecontour_upper_0.parallel_distance[1:end] - finecontour_upper_0.parallel_distance[0 : end - 1]
print(parallel_gaps_upper)
print("total distance poloidal upper 0: ", finecontour_upper_0.distance[-1])
print("total distance parallel upper 0: ", finecontour_upper_0.parallel_distance[-1])
end = len(finecontour_core_0.parallel_distance)-1
parallel_gaps_core = finecontour_core_0.parallel_distance[1:end] - finecontour_core_0.parallel_distance[0 : end - 1]
print(parallel_gaps_core)
print("total distance poloidal core 0: ", finecontour_core_0.distance[-1])
print("total distance parallel core 0: ", finecontour_core_0.parallel_distance[-1])
end = len(finecontour_lower_0.parallel_distance)-1
parallel_gaps_lower = finecontour_lower_0.parallel_distance[1:end] - finecontour_lower_0.parallel_distance[0 : end - 1]
print(parallel_gaps_lower)
print("total distance poloidal lower 0: ", finecontour_lower_0.distance[-1])
print("total distance parallel lower 0: ", finecontour_lower_0.parallel_distance[-1])
print("final length of B field values 0: ", len(B_concatenated_0))
# delete the first 8 and last 8 elements of B_concatenated
B_concatenated_trimmed = B_concatenated_0#[8:-8]
print("final B field values trimmed 0: ", B_concatenated_trimmed.tolist())
print("length of final B field values trimmed 0: ", len(B_concatenated_trimmed))




end = len(finecontour_upper_1.parallel_distance)-1
parallel_gaps_upper = finecontour_upper_1.parallel_distance[1:end] - finecontour_upper_1.parallel_distance[0 : end - 1]
print(parallel_gaps_upper)
print("total distance poloidal upper 1: ", finecontour_upper_1.distance[-1])
print("total distance parallel upper 1: ", finecontour_upper_1.parallel_distance[-1])
end = len(finecontour_core_1.parallel_distance)-1
parallel_gaps_core = finecontour_core_1.parallel_distance[1:end] - finecontour_core_1.parallel_distance[0 : end - 1]
print(parallel_gaps_core)
print("total distance poloidal core 1: ", finecontour_core_1.distance[-1])
print("total distance parallel core 1: ", finecontour_core_1.parallel_distance[-1])
end = len(finecontour_lower_1.parallel_distance)-1
parallel_gaps_lower = finecontour_lower_1.parallel_distance[1:end] - finecontour_lower_1.parallel_distance[0 : end - 1]
print(parallel_gaps_lower)
print("total distance poloidal lower 1: ", finecontour_lower_1.distance[-1])
print("total distance parallel lower 1: ", finecontour_lower_1.parallel_distance[-1])
print("final length of B field values 1: ", len(B_concatenated_1))
# delete the first 8 and last 8 elements of B_concatenated
B_concatenated_trimmed = B_concatenated_1#[8:-8]
print("final B field values trimmed 1: ", B_concatenated_trimmed.tolist())
print("length of final B field values trimmed 1: ", len(B_concatenated_trimmed))

end = len(finecontour_upper_2.parallel_distance)-1
parallel_gaps_upper = finecontour_upper_2.parallel_distance[1:end] - finecontour_upper_2.parallel_distance[0 : end - 1]
print(parallel_gaps_upper)
print("total distance poloidal upper 2: ", finecontour_upper_2.distance[-1])
print("total distance parallel upper 2: ", finecontour_upper_2.parallel_distance[-1])
end = len(finecontour_core_2.parallel_distance)-1
parallel_gaps_core = finecontour_core_2.parallel_distance[1:end] - finecontour_core_2.parallel_distance[0 : end - 1]
print(parallel_gaps_core)
print("total distance poloidal core 2: ", finecontour_core_2.distance[-1])
print("total distance parallel core 2: ", finecontour_core_2.parallel_distance[-1])
end = len(finecontour_lower_2.parallel_distance)-1
parallel_gaps_lower = finecontour_lower_2.parallel_distance[1:end] - finecontour_lower_2.parallel_distance[0 : end - 1]
print(parallel_gaps_lower)
print("total distance poloidal lower 2: ", finecontour_lower_2.distance[-1])
print("total distance parallel lower 2: ", finecontour_lower_2.parallel_distance[-1])
print("final length of B field values 2: ", len(B_concatenated_2))
# delete the first 8 and last 8 elements of B_concatenated
B_concatenated_trimmed = B_concatenated_2#[8:-8]
print("final B field values trimmed 2: ", B_concatenated_trimmed.tolist())
print("length of final B field values trimmed 2: ", len(B_concatenated_trimmed))

end = len(finecontour_upper_3.parallel_distance)-1
parallel_gaps_upper = finecontour_upper_3.parallel_distance[1:end] - finecontour_upper_3.parallel_distance[0 : end - 1]
print(parallel_gaps_upper)
print("total distance poloidal upper 3: ", finecontour_upper_3.distance[-1])
print("total distance parallel upper 3: ", finecontour_upper_3.parallel_distance[-1])
end = len(finecontour_core_3.parallel_distance)-1
parallel_gaps_core = finecontour_core_3.parallel_distance[1:end] - finecontour_core_3.parallel_distance[0 : end - 1]
print(parallel_gaps_core)
print("total distance poloidal core 3: ", finecontour_core_3.distance[-1])
print("total distance parallel core 3: ", finecontour_core_3.parallel_distance[-1])
end = len(finecontour_lower_3.parallel_distance)-1
parallel_gaps_lower = finecontour_lower_3.parallel_distance[1:end] - finecontour_lower_3.parallel_distance[0 : end - 1]
print(parallel_gaps_lower)
print("total distance poloidal lower 3: ", finecontour_lower_3.distance[-1])
print("total distance parallel lower 3: ", finecontour_lower_3.parallel_distance[-1])
print("final length of B field values 3: ", len(B_concatenated_3))
# delete the first 8 and last 8 elements of B_concatenated
B_concatenated_trimmed = B_concatenated_3#[8:-8]
print("final B field values trimmed 3: ", B_concatenated_trimmed.tolist())
print("length of final B field values trimmed 3: ", len(B_concatenated_trimmed))

end = len(finecontour_upper_end.parallel_distance)-1
parallel_gaps_upper = finecontour_upper_end.parallel_distance[1:end] - finecontour_upper_end.parallel_distance[0 : end - 1]
print(parallel_gaps_upper)
print("total distance poloidal upper end: ", finecontour_upper_end.distance[-1])
print("total distance parallel upper end: ", finecontour_upper_end.parallel_distance[-1])
end = len(finecontour_core_end.parallel_distance)-1
parallel_gaps_core = finecontour_core_end.parallel_distance[1:end] - finecontour_core_end.parallel_distance[0 : end - 1]
print(parallel_gaps_core)
print("total distance poloidal core end: ", finecontour_core_end.distance[-1])
print("total distance parallel core end: ", finecontour_core_end.parallel_distance[-1])
end = len(finecontour_lower_end.parallel_distance)-1
parallel_gaps_lower = finecontour_lower_end.parallel_distance[1:end] - finecontour_lower_end.parallel_distance[0 : end - 1]
print(parallel_gaps_lower)
print("total distance poloidal lower end: ", finecontour_lower_end.distance[-1])
print("total distance parallel lower end: ", finecontour_lower_end.parallel_distance[-1])
print("final length of B field values end: ", len(B_concatenated_3))
# delete the first 8 and last 8 elements of B_concatenated
B_concatenated_trimmed = B_concatenated_end#[8:-8]
print("final B field values trimmed end: ", B_concatenated_trimmed.tolist())
print("length of final B field values trimmed end: ", len(B_concatenated_trimmed))



"""
B FIELD VALUES!!!!!!! FOR ROUGH ANALYTIC EXPRESSION TO USE IN MOMENT KINETICS CODE FOR NOW.

[0.78972903 0.75695292 0.70006378 0.64873391 0.60296929 0.56265807
 0.52764883 0.49771051 0.47254802 0.45175239 0.4348647  0.42137203
 0.41075044 0.40250482 0.39616853 0.39134206 0.38768643 0.38491617
 0.3828036  0.38116832 0.37986752 0.37878817 0.37784351 0.37696688
 0.37610704 0.37522499 0.37429135 0.37328425 0.37218778 0.37099083
 0.37086713 0.36955515 0.36813187 0.36659682 0.36495118 0.36319782
 0.3613413  0.35938773 0.35734462 0.35522187 0.35303035 0.35078131
 0.34848701 0.3461607  0.3438166  0.34147184 0.33914271 0.33684509
 0.33459477 0.3324075  0.33030175 0.32829312 0.32639448 0.32461682
 0.32297202 0.32147022 0.3201153  0.31890779 0.31785036 0.3169406
 0.31617012 0.31553065 0.31501543 0.31461056 0.31430325 0.31408499
 0.31394226 0.31386616 0.31385165 0.31389187 0.31398776 0.31414021
 0.31435371 0.31463648 0.31499569 0.31544148 0.31598186 0.316624
 0.31737082 0.31822233 0.319173   0.32021182 0.32132355 0.32248658
 0.32367637 0.32486497 0.32602158 0.32711671 0.32811998 0.32900373
 0.32974257 0.33031596 0.33070766 0.33090632 0.33090632 0.33070766
 0.33031596 0.32974257 0.32900373 0.32811998 0.32711671 0.32602158
 0.32486497 0.32367637 0.32248658 0.32132355 0.32021182 0.319173
 0.31822233 0.31737082 0.316624   0.31598186 0.31544148 0.31499569
 0.31463648 0.31435371 0.31414021 0.31398776 0.31389187 0.31385165
 0.31386616 0.31394226 0.31408499 0.31430325 0.31461056 0.31501543
 0.31553065 0.31617012 0.3169406  0.31785036 0.31890779 0.3201153
 0.32147022 0.32297202 0.32461682 0.32639448 0.32829312 0.33030175
 0.3324075  0.33459477 0.33684509 0.33914271 0.34147184 0.3438166
 0.3461607  0.34848701 0.35078131 0.35303035 0.35522187 0.35734462
 0.35938773 0.3613413  0.36319782 0.36495118 0.36659682 0.36813187
 0.36955515 0.37086713 0.37099083 0.37218778 0.37328425 0.37429135
 0.37522499 0.37610704 0.37696688 0.37784351 0.37878817 0.37986752
 0.38116832 0.3828036  0.38491617 0.38768643 0.39134206 0.39616853
 0.40250482 0.41075044 0.42137203 0.4348647  0.45175239 0.47254802
 0.49771051 0.52764883 0.56265807 0.60296929 0.64873391 0.70006378
 0.75695292 0.78972903]
"""
"""





mreg_upper = mesh.regions[mesh.region_lookup[("outer_upper_divertor", 1)]]
mreg_lower = mesh.regions[mesh.region_lookup[("outer_lower_divertor", 1)]]
mreg_core  = mesh.regions[mesh.region_lookup[("outer_core", 1)]]
contour_upper = mreg_upper.contours[5]
contour_lower = mreg_lower.contours[5]
contour_core  = mreg_core.contours[5]
finecontour_upper = contour_upper.get_fine_contour(psi=eq.psi, equilibrium=eq)
finecontour_lower = contour_lower.get_fine_contour(psi=eq.psi, equilibrium=eq)
finecontour_core  = contour_core.get_fine_contour(psi=eq.psi, equilibrium=eq)




figwidth = 4.0
figheight = figwidth * (eq.Zmax - eq.Zmin) / (eq.Rmax - eq.Rmin)
fig, ax = plt.subplots(figsize=(figwidth, figheight), constrained_layout=True)


colors = "grey"
eq.plotPotential(
    npoints=50,
    ncontours=100,
    labels=True,
    colors=colors,
    linestyles="-",
)
eq.plotWall(axis=ax)
contour_upper.plot(psi=eq.psi, ax=ax, linestyle="", marker="x", color="red")
finecontour_upper.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=1, color="blue")
contour_lower.plot(psi=eq.psi, ax=ax, linestyle="", marker="x", color="red")
finecontour_lower.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=1, color="blue")
contour_core.plot(psi=eq.psi, ax=ax, linestyle="", marker="x", color="red")
finecontour_core.plot(psi=eq.psi, ax=ax, linestyle="", marker=".", markersize=1, color="blue")

plt.savefig(image_dir.joinpath("testing_wholecontourplot.pdf"), bbox_inches="tight")

end = len(finecontour_upper.parallel_distance)-1
parallel_gaps_upper = finecontour_upper.parallel_distance[1:end] - finecontour_upper.parallel_distance[0 : end - 1]
print(parallel_gaps_upper)

contour_upper_parallel_distance = np.array(contour_upper.get_parallel_distance(psi=eq.psi, equilibrium=eq))
end = len(contour_upper_parallel_distance)-1
print(contour_upper_parallel_distance)
parallel_gaps_contour_upper = contour_upper_parallel_distance[1:end] - contour_upper_parallel_distance[0 : end - 1]
print(parallel_gaps_contour_upper)
"""

"""
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
"""