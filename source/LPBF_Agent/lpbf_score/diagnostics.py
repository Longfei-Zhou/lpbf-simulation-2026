"""Transparent diagnostic metrics integrated from the earlier pysimle evaluator.

The legacy evaluator contained useful ideas (G/V, G*V, cooling-rate and remelt
distributions), but also data-derived scoring windows and domain-size proxies that
could make the same dataset define its own pass range.  This module keeps only the
directly computed, unit-labelled diagnostics.  They are intentionally not folded
into PASS/FAIL until a material/machine-specific experimental window is supplied.
"""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import pandas as pd


def _finite(values: Iterable[float], *, positive: bool = False) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    mask = np.isfinite(array)
    if positive:
        mask &= array > 0
    return array[mask]


def distribution(values: Iterable[float], *, positive: bool = False) -> dict[str, Any]:
    """Return a JSON-safe robust distribution summary."""
    array = _finite(values, positive=positive)
    if not len(array):
        return {
            "count": 0,
            "min": None,
            "p10": None,
            "median": None,
            "mean": None,
            "p90": None,
            "max": None,
            "std": None,
        }
    return {
        "count": int(len(array)),
        "min": float(np.min(array)),
        "p10": float(np.quantile(array, 0.10)),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "p90": float(np.quantile(array, 0.90)),
        "max": float(np.max(array)),
        "std": float(np.std(array, ddof=1)) if len(array) > 1 else 0.0,
    }


def _column(frame: pd.DataFrame, name: str) -> np.ndarray:
    if name not in frame.columns:
        return np.array([], dtype=float)
    return pd.to_numeric(frame[name], errors="coerce").to_numpy(float)


def _active_mask(frame: pd.DataFrame) -> np.ndarray:
    """Select points that experienced a thermal/solidification event.

    Custom target grids deliberately retain numMelt=0 points.  Including those
    zeros in G, V and tSol distributions would make a coverage failure look like a
    microstructure observation.  numMelt is therefore the preferred activity mask.
    """
    if "numMelt" in frame.columns:
        values = _column(frame, "numMelt")
        return np.isfinite(values) & (values > 0)
    candidates = []
    for name in ("tSol", "G", "V", "dTdt"):
        values = _column(frame, name)
        if len(values):
            candidates.append(np.isfinite(values) & (np.abs(values) > 0))
    if not candidates:
        return np.ones(len(frame), dtype=bool)
    return np.logical_or.reduce(candidates)


def build_process_diagnostics(
    solid: pd.DataFrame,
    snapshots: pd.DataFrame,
) -> dict[str, Any]:
    """Build diagnostic-only thermal, solidification and melt-pool summaries."""
    active = _active_mask(solid)
    active_count = int(active.sum())
    solidification: dict[str, Any] = {
        "all_point_count": int(len(solid)),
        "thermally_active_point_count": active_count,
        "thermally_active_fraction": (
            float(active_count / len(solid)) if len(solid) else None
        ),
    }

    for column, key, positive, absolute in (
        ("G", "temperature_gradient_k_m", True, False),
        ("V", "solidification_velocity_m_s", True, False),
        ("dTdt", "cooling_rate_abs_k_s", True, True),
        # tSol is the absolute event time in the 3DThesis run, not a local
        # liquid-residence duration.  Keep the name explicit to prevent the
        # legacy report's misleading "solidification duration" interpretation.
        ("tSol", "solidification_event_time_s", True, False),
    ):
        values = _column(solid, column)
        if len(values):
            selected = values[active]
            if absolute:
                selected = np.abs(selected)
            solidification[key] = distribution(selected, positive=positive)

    g = _column(solid, "G")
    v = _column(solid, "V")
    cooling = _column(solid, "dTdt")
    if len(g) and len(v):
        valid = active & np.isfinite(g) & np.isfinite(v) & (g > 0) & (v > 0)
        ratio = g[valid] / v[valid]
        product = g[valid] * v[valid]
        solidification["g_over_v_k_s_m2"] = distribution(ratio, positive=True)
        solidification["g_times_v_k_s"] = distribution(product, positive=True)
        if len(cooling):
            observed = np.abs(cooling[valid])
            valid_compare = np.isfinite(observed) & (observed > 0)
            relative_error = np.abs(observed[valid_compare] - product[valid_compare]) / np.maximum(
                observed[valid_compare], 1e-12
            )
            solidification["dtdt_vs_gv_median_relative_error"] = (
                float(np.median(relative_error)) if len(relative_error) else None
            )

    remelt: dict[str, Any] = {"available": False}
    remelt_values = _column(solid, "numMelt")
    if len(remelt_values):
        finite = remelt_values[np.isfinite(remelt_values)]
        remelt = {
            "available": True,
            "distribution": distribution(finite),
            "count_0": int(np.sum(finite == 0)),
            "count_1": int(np.sum(finite == 1)),
            "count_2": int(np.sum(finite == 2)),
            "count_3_or_more": int(np.sum(finite >= 3)),
            "fraction_0": float(np.mean(finite == 0)) if len(finite) else None,
            "fraction_2_or_more": float(np.mean(finite >= 2)) if len(finite) else None,
            "fraction_3_or_more": float(np.mean(finite >= 3)) if len(finite) else None,
        }

    snapshot_summary: dict[str, Any] = {
        "snapshot_count": int(len(snapshots)),
    }
    for column, key in (
        ("width_m", "melt_pool_width_m"),
        ("depth_m", "melt_pool_depth_m"),
        ("length_m", "melt_pool_length_m"),
        ("molten_volume_m3", "melt_pool_volume_m3"),
        ("max_temperature_k_diagnostic_only", "peak_temperature_k_diagnostic_only"),
    ):
        values = _column(snapshots, column)
        if len(values):
            snapshot_summary[key] = distribution(values, positive=True)

    return {
        "status": "DIAGNOSTIC_ONLY_NOT_SCORED",
        "integration_source": "pysimle concepts, recomputed with LPBF_Agent definitions",
        "solidification": solidification,
        "remelting": remelt,
        "snapshots": snapshot_summary,
        "interpretation": {
            "g_over_v": "Relative solidification-morphology indicator; no universal higher-is-better direction.",
            "g_times_v": "Cooling-rate consistency proxy; compare against dTdt and calibrated material targets.",
            "scope": "These values support diagnosis and future calibration, but do not change release status by themselves.",
        },
        "legacy_methods_rejected": [
            "self-normalisation using the same evaluated dataset",
            "whole-domain mean temperature as a quality target",
            "temperature range divided by domain span as a gradient field",
            "unreferenced 60/40 stage weights and legacy 316L process windows",
        ],
    }
