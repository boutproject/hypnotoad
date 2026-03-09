"""
Functions to inspect a BoutMesh object and return summary information
"""

import numpy as np
from typing import Optional


def _region_stats(arr: np.ndarray, name: str) -> dict:
    """Compute statistics for a 2D array of cell sizes or metric values."""
    if np.any(~np.isfinite(arr)):
        return {"error": f"{name} contains non-finite values"}
    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "uniformity_ratio": (
            float(np.max(arr) / np.min(arr)) if np.min(arr) > 0 else None
        ),
    }


def _max_adjacent_ratio(arr2d: np.ndarray) -> float:
    """
    Maximum ratio between adjacent cell sizes in either direction.
    Values >> 1 indicate sudden jumps that cause numerical diffusion.
    """
    ratios = []
    for axis in (0, 1):
        slc_a = [slice(None)] * 2
        slc_b = [slice(None)] * 2
        slc_a[axis] = slice(None, -1)
        slc_b[axis] = slice(1, None)
        a = arr2d[tuple(slc_a)]
        b = arr2d[tuple(slc_b)]
        # Avoid division by zero; take max(a/b, b/a)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(b > 0, a / b, np.nan)
            r = np.fmax(r, np.where(a > 0, b / a, np.nan))
        finite_r = r[np.isfinite(r)]
        if len(finite_r):
            ratios.append(float(np.max(finite_r)))
    return max(ratios) if ratios else float("nan")


def _get_array(region, name: str) -> Optional[np.ndarray]:
    arr = getattr(region, name, None)
    if arr is None:
        return None
    return np.asarray(arr.centre)


def _collect_region_arrays(mesh) -> dict[str, dict[str, np.ndarray]]:
    """
    BoutMesh stores data in regions accessible via mesh.regions (a dict of
    name -> MeshRegion). Each region has 2D arrays for coordinates and metrics.
    Falls back to whole-mesh arrays if regions are not available.
    """
    arrays_by_region = {}

    # Try structured regions first (preferred)
    regions = getattr(mesh, "regions", None)
    if regions:
        for rindx, region in regions.items():
            # Mesh uses indices rindx
            arrays_by_region[region.name] = {
                "neighbors": {
                    edge_name: regions[nindx].name if nindx is not None else None
                    for edge_name, nindx in mesh.connections[rindx].items()
                },
                "Rxy": _get_array(region, "Rxy"),
                "Zxy": _get_array(region, "Zxy"),
                "dx": _get_array(region, "dx"),
                "dy": _get_array(region, "dy"),
                "J": _get_array(region, "J"),
                "g11": _get_array(region, "g11"),
                "g22": _get_array(region, "g22"),
                "g33": _get_array(region, "g33"),
                "g_11": _get_array(region, "g_11"),
                "g_22": _get_array(region, "g_22"),
                "g_33": _get_array(region, "g_33"),
                "Bxy": _get_array(region, "Bxy"),
            }
    else:
        # Fall back to whole-mesh attributes
        arrays_by_region["global"] = {
            "Rxy": _get_array(mesh, "Rxy"),
            "Zxy": _get_array(mesh, "Zxy"),
            "dx": _get_array(mesh, "dx"),
            "dy": _get_array(mesh, "dy"),
            "J": _get_array(mesh, "J"),
            "g11": _get_array(mesh, "g11"),
            "g22": _get_array(mesh, "g22"),
            "g_11": _get_array(mesh, "g_11"),
            "g_22": _get_array(mesh, "g_22"),
            "Bxy": _get_array(mesh, "Bxy"),
        }

    return arrays_by_region


def inspect_mesh(mesh, detail="summary") -> dict:
    """
    Inspect a BoutMesh object after calculateRZ() and geometry() have been called.
    Returns a structured diagnostic dict suitable as a tool result.
    """
    arrays_by_region = _collect_region_arrays(mesh)

    region_diagnostics = {}
    all_J = []
    all_dx = []

    errors = []
    warnings = []

    for rname, arrs in arrays_by_region.items():
        dx = arrs.get("dx")
        dy = arrs.get("dy")
        J = arrs.get("J")
        g22 = arrs.get("g22")
        g_22 = arrs.get("g_22")

        rd = {
            "neighbors": arrs.get("neighbors", None),
            "nx": dx.shape[0],
            "ny": dx.shape[1],
            "valid": True,
        }

        if dx is not None:
            rd["dx"] = _region_stats(dx, f"{rname}.dx")
            rd["max_adjacent_dx_ratio"] = _max_adjacent_ratio(dx)
            all_dx.append(dx.ravel())

        if dy is not None:
            rd["dy"] = _region_stats(dy, f"{rname}.dy")
            rd["max_adjacent_dy_ratio"] = _max_adjacent_ratio(dy)

            if g22 is not None:
                dlpol = dy / np.sqrt(g22)
                rd["dlpol_poloidal_cell_size"] = _region_stats(dlpol, f"{rname}.dlpol")
                rd["max_adjacent_dlpol_ratio"] = _max_adjacent_ratio(dlpol)

            if g_22 is not None:
                dlpar = dy * np.sqrt(g_22)
                rd["dlpar_parallel_cell_size"] = _region_stats(dlpar, f"{rname}.dlpar")
                rd["max_adjacent_dlpar_ratio"] = _max_adjacent_ratio(dlpar)

        if J is not None:
            rd["J"] = _region_stats(J, f"{rname}.J")
            rd["n_negative_jacobian"] = int(np.sum(J <= 0))
            all_J.append(J.ravel())

        for varname in [
            "dx",
            "dy",
            "J",
            "dlpol_poloidal_cell_size",
            "dlpar_parallel_cell_size",
        ]:
            if "error" in rd[varname]:
                errors.append(rd[varname]["error"])
                rd["valid"] = False

        region_diagnostics[rname] = rd

    # Check for small dx
    for rname, rd in region_diagnostics.items():
        if ("min" in rd["dx"]) and rd["dx"]["min"] < 1e-8:
            errors.append(f"Too small dx in region {rname}: {rd['dx']['min']}")
            rd["valid"] = False

    # Large non-uniformity
    for rname, rd in region_diagnostics.items():
        ratio = rd["dx"].get("uniformity_ratio", None)
        if ratio is None:
            continue
        if ratio > 1e2:
            warnings.append(f"Large dx uniformity_ratio in region {rname}: {ratio}")

    # Check for dx in neighboring regions
    for rname, rd in region_diagnostics.items():
        dx_mean = rd["dx"].get("mean", None)
        if dx_mean is None:
            continue

        # Which setting affects the size of this region?
        if (rd["neighbors"].get("inner", None) is not None) and (
            rd["neighbors"].get("outer", None) is not None
        ):
            nx_name = "nx_inter_sep"
        elif rd["neighbors"].get("inner", None) is None:
            if "core" in rname:
                nx_name = "nx_core"
            else:
                nx_name = "nx_pf"
        else:
            nx_name = "nx_sol"

        for direction in ["inner", "outer"]:
            nname = rd["neighbors"].get(direction, None)
            if nname is None:
                continue
            ndx_mean = region_diagnostics[nname]["dx"].get("mean", None)
            if ndx_mean is None:
                continue
            if dx_mean > 10 * ndx_mean:
                warnings.append(
                    f"mean dx in region {rname} is {dx_mean / ndx_mean} times dx in neighbor region {nname}. Increase setting {nx_name}."
                )
            if dx_mean < 0.1 * ndx_mean:
                warnings.append(
                    f"mean dx in region {rname} is {dx_mean / ndx_mean} times dx in neighbor region {nname}. Decrease setting {nx_name}."
                )

    if detail == "summary":
        # One-liner per region: name, size, pass/fail
        region_out = {
            name: {
                "nx": rd.get("nx"),
                "ny": rd.get("ny"),
                "valid": rd["valid"],
            }
            for name, rd in region_diagnostics.items()
        }
        return {
            "detail": "summary",
            "valid": len(errors) == 0,
            "errors": errors,
            "warnings": warnings,
            "regions": region_out,
            "hint": (
                "Call inspect_mesh(detail='full') to investigate warnings"
                if len(warnings) > 0
                else None
            ),
        }

    # Full
    return {
        "detail": "full",
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "regions": region_diagnostics,
    }
