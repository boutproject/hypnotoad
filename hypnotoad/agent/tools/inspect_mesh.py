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
        "uniformity_ratio": float(np.max(arr) / np.min(arr))
        if np.min(arr) > 0
        else None,
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
    arr = getattr(region, "Rxy", None)
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

    Attribute names (Rxy, dx, J, etc.) are looked up via multiple candidate names
    to be robust against BoutMesh internals. If a quantity cannot be found,
    it is reported as null rather than raising an exception.
    """
    arrays_by_region = _collect_region_arrays(mesh)

    region_diagnostics = {}
    all_J = []
    all_dx = []

    for rname, arrs in arrays_by_region.items():
        dx = arrs.get("dx")
        dy = arrs.get("dy")
        J = arrs.get("J")

        rd = {"neighbors": arrs.get("neighbors", None)}

        if dx is not None:
            rd["nx"] = dx.shape[1] if dx.ndim == 2 else dx.shape[0]
            rd["dx_stats"] = _region_stats(dx, f"{rname}.dx")
            rd["max_adjacent_dx_ratio"] = _max_adjacent_ratio(dx)
            all_dx.append(dx.ravel())

        if dy is not None:
            rd["ny"] = dy.shape[0] if dy.ndim == 2 else dy.shape[0]
            rd["dy_stats"] = _region_stats(dy, f"{rname}.dy")
            rd["max_adjacent_dy_ratio"] = _max_adjacent_ratio(dy)

        if J is not None:
            rd["J_stats"] = _region_stats(J, f"{rname}.J")
            rd["n_negative_jacobian"] = int(np.sum(J <= 0))
            all_J.append(J.ravel())

        region_diagnostics[rname] = rd

    if detail == "summary":
        # One-liner per region: name, size, pass/fail
        region_out = {
            name: {
                "nx": rd.get("nx"),
                "ny": rd.get("ny"),
                "ok": rd.get("n_negative_jacobian", 0) == 0,
            }
            for name, rd in region_diagnostics.items()
        }
        return {
            "detail": "summary",
            "valid": True,
            "regions": region_out,
            "hint": None,  # "Call inspect_mesh(detail='standard') to investigate warnings"
            # if warnings else None,
        }

    # Full
    return {
        "status": "ok",
        "regions": region_diagnostics,
    }
