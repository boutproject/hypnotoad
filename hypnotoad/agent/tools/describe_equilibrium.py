def describe_equilibrium(gridfile) -> dict:
    """
    Extract physics-relevant metrics from a GEQDSK file that
    inform good mesh settings choices.
    """
    from ...geqdsk import _geqdsk
    from ...utils import critical
    import numpy as np

    # Read header line, discard single characters
    with open(gridfile, "rt") as fh:
        header = fh.readline()
    header_tok = [tok for tok in header.split() if len(tok) > 1]

    with open(gridfile, "rt") as fh:
        data = _geqdsk.read(fh)

    # Range of psi normalises psi derivatives
    psi_bdry_gfile = data["sibdry"]
    psi_axis_gfile = data["simagx"]

    R1D = np.linspace(
        data["rleft"], data["rleft"] + data["rdim"], data["nx"], endpoint=True
    )

    Z1D = np.linspace(
        data["zmid"] - 0.5 * data["zdim"],
        data["zmid"] + 0.5 * data["zdim"],
        data["ny"],
        endpoint=True,
    )

    psi2D = data["psi"]
    # Find critical points (O- and X-points)
    R2D, Z2D = np.meshgrid(R1D, Z1D, indexing="ij")
    opoints, xpoints = critical.find_critical(R2D, Z2D, psi2D, 1.0e-6, 1000)

    warning_list = []
    if len(opoints) == 0:
        warning_list.append("No O-points found in input magnetic field.")
        magnetic_axis = None
    else:
        magnetic_axis = {"R": opoints[0][0], "Z": opoints[0][1], "psinorm": 0.0}

    if len(xpoints) == 0:
        warning_list.append("No X-points found in input magnetic field.")
        xpoints = None
    else:
        xpoints = [
            {
                "R": r,
                "Z": z,
                "psinorm": (psi - psi_axis_gfile) / (psi_bdry_gfile - psi_axis_gfile),
            }
            for r, z, psi in xpoints[:3]
        ]  # Maximum 3

    return {
        "header": header_tok,
        "psi_increasing": psi_bdry_gfile > psi_axis_gfile,
        "magnetic_axis": magnetic_axis,
        "xpoints": xpoints,
        "warnings": warning_list,
    }
