#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Posterior spectral diagnostics for POSEIDON retrievals.

This helper regenerates posterior spectra from the saved MultiNest output
using ``POSEIDON.retrieval.retrieved_samples(...)`` so we do not depend only
on the median retrieved spectrum when propagating MSE and reduced chi-squared.

The module is intentionally usable in two ways:

1. Imported from notebooks:
   from posterior_diagnostics import run_case

2. Executed from the terminal:
   python posterior_diagnostics.py --mode uncontam --branch phoenix

The legacy metrics columns (``MSE``, ``chi2``, ``chi2_reduced``) remain tied to
the median retrieved spectrum for backwards compatibility with the existing
logs, while new ``posterior_*`` columns summarize the propagated distribution
computed from posterior sampled spectra.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
SPHINX_DIR = BASE_DIR / "sphinx_injection"

PLANET_NAME = "Trappist-1e"
INSTRUMENTS = ["JWST_NIRSpec_PRISM"]

WL_MIN_UM = 0.4
WL_MAX_UM = 6.0
R_GRID = 4000

P_MIN_BAR = 1.0e-2
P_MAX_BAR = 2.0
N_LAYERS = 100
P_REF_BAR = 1.0

OPACITY_TREATMENT = "opacity_sampling"
T_FINE = np.arange(200.0, 400.0 + 10.0, 10.0)
LOG_P_FINE = np.arange(-2.0, 2.0 + 0.2, 0.2)

SIGMA_Q16 = 0.15865525393145707
SIGMA_Q84 = 0.8413447460685429

DEFAULT_POSTERIOR_SAMPLES = 1000
DEFAULT_N_FREE_PARAMS = 11


@dataclass(frozen=True)
class CaseConfig:
    """Resolved configuration for one retrieval case."""

    mode_label: str
    branch: str
    source_label: str | None
    n_transits: int
    f_spot_case: float
    f_fac_case: float
    posterior_samples: int
    n_free_params: int
    clean_path: Path
    csv_path: Path
    npz_path: Path
    data_dir: Path
    output_dir: Path
    observation_file: str
    model_name: str


def _parse_bool(value: str) -> bool:
    """Parse a shell-style boolean flag."""
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y", "on"}:
        return True
    if normalized in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def _env_or_default(name: str, cast, default):
    """Read one environment variable with casting and a default value."""
    raw_value = os.getenv(name)
    if raw_value is None:
        return default

    try:
        return cast(raw_value)
    except ValueError as exc:
        raise ValueError(f"Invalid value for {name!r}: {raw_value!r}") from exc


def _resolve_path(raw_path: str | None, default_path: Path) -> Path:
    """Resolve user-provided or default paths relative to this directory."""
    if not raw_path:
        return default_path.resolve()

    path = Path(raw_path)
    if not path.is_absolute():
        path = (BASE_DIR / path).resolve()
    else:
        path = path.resolve()

    return path


def build_observation_filename(
    n_transits: int,
    f_spot_case: float,
    f_fac_case: float,
    *,
    reconstructed: bool = False,
    source_label: str | None = None,
) -> str:
    """Build one observation filename following the existing project naming."""
    stem = f"pandexo_output_{n_transits}transits_"
    if source_label:
        stem += f"{source_label}_"
    stem += f"fspot{f_spot_case:.2f}_ffac{f_fac_case:.2f}"
    if reconstructed:
        stem += "_recon"
    return f"{stem}.dat"


def build_model_name(
    mode_label: str,
    n_transits: int,
    f_spot_case: float,
    f_fac_case: float,
    *,
    source_label: str | None = None,
) -> str:
    """Build one POSEIDON model name following the existing project naming."""
    prefix = f"{source_label}_" if source_label else ""
    return (
        f"{prefix}{mode_label}_{n_transits}T_"
        f"{f_spot_case:.2f}spot-{f_fac_case:.2f}fac"
    )


def load_clean_two_cols(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a two-column clean spectrum file and sort it by wavelength."""
    arr = np.genfromtxt(path, comments="#", dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.shape[1] < 2:
        raise ValueError(
            f"Expected at least 2 columns (wl_um, depth) in clean spectrum: {path}"
        )

    wl_clean = arr[:, 0].astype(float)
    y_clean = arr[:, 1].astype(float)
    order = np.argsort(wl_clean)
    return wl_clean[order], y_clean[order]


def bin_average_with_halfbins(
    wl_src: np.ndarray,
    y_src: np.ndarray,
    centers: np.ndarray,
    halfwidths: np.ndarray,
    *,
    nsamp: int = 256,
) -> np.ndarray:
    """Band-average one source spectrum over observed bin half-widths."""
    wl_src = np.asarray(wl_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    centers = np.asarray(centers, dtype=float)
    halfwidths = np.asarray(halfwidths, dtype=float)

    order = np.argsort(wl_src)
    wl_sorted = wl_src[order]
    y_sorted = y_src[order]

    out = np.empty_like(centers, dtype=float)
    for idx, (center, halfwidth) in enumerate(zip(centers, halfwidths)):
        start = center - halfwidth
        stop = center + halfwidth
        x = np.linspace(start, stop, nsamp)
        y = np.interp(x, wl_sorted, y_sorted)
        out[idx] = np.trapz(y, x) / (stop - start)

    return out


def metric_summary(values: np.ndarray, prefix: str) -> dict[str, float]:
    """Summarize one posterior metric distribution."""
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]

    if finite.size == 0:
        return {
            f"{prefix}_mean": np.nan,
            f"{prefix}_std": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_q16": np.nan,
            f"{prefix}_q84": np.nan,
        }

    return {
        f"{prefix}_mean": float(np.mean(finite)),
        f"{prefix}_std": float(np.std(finite, ddof=0)),
        f"{prefix}_median": float(np.quantile(finite, 0.5)),
        f"{prefix}_q16": float(np.quantile(finite, SIGMA_Q16)),
        f"{prefix}_q84": float(np.quantile(finite, SIGMA_Q84)),
    }


def compute_scalar_metrics(
    y_model: np.ndarray,
    y_true: np.ndarray,
    err: np.ndarray,
    *,
    n_free_params: int,
) -> dict[str, float | int]:
    """Compute MSE, RMSE, chi-squared and reduced chi-squared."""
    y_model = np.asarray(y_model, dtype=float)
    y_true = np.asarray(y_true, dtype=float)
    err = np.asarray(err, dtype=float)

    if not (y_model.shape == y_true.shape == err.shape):
        raise ValueError(
            "Metric inputs must share the same shape: "
            f"model={y_model.shape}, truth={y_true.shape}, err={err.shape}"
        )
    if np.any(err <= 0.0):
        raise ValueError("Found non-positive uncertainties in err; chi2 is undefined.")

    residuals = y_model - y_true
    n_points = int(residuals.size)
    dof = int(max(n_points - int(n_free_params), 0))

    mse = float(np.mean(residuals**2))
    rmse = float(np.sqrt(mse))
    chi2 = float(np.sum((residuals / err) ** 2))
    chi2_reduced = float(chi2 / dof) if dof > 0 else np.nan

    return {
        "N": n_points,
        "p": int(n_free_params),
        "dof": dof,
        "MSE": mse,
        "rmse": rmse,
        "rmse_ppm": float(1.0e6 * rmse),
        "chi2": chi2,
        "chi2_reduced": chi2_reduced,
    }


def summarize_posterior_metrics(
    ymodel_samples: np.ndarray,
    clean_binned: np.ndarray,
    err_obs: np.ndarray,
    *,
    dof: int,
) -> dict[str, float | int]:
    """Propagate MSE and chi-squared through posterior sampled spectra."""
    ymodel_samples = np.asarray(ymodel_samples, dtype=float)
    clean_binned = np.asarray(clean_binned, dtype=float)
    err_obs = np.asarray(err_obs, dtype=float)

    if ymodel_samples.ndim != 2:
        raise ValueError(
            f"Expected ymodel_samples to be 2D, received shape {ymodel_samples.shape}"
        )
    if ymodel_samples.shape[1] != clean_binned.size:
        raise ValueError(
            "Observed-grid sample length does not match clean truth bins: "
            f"{ymodel_samples.shape[1]} vs {clean_binned.size}"
        )
    if err_obs.size != clean_binned.size:
        raise ValueError(
            f"Observed uncertainties length mismatch: {err_obs.size} vs {clean_binned.size}"
        )
    if np.any(err_obs <= 0.0):
        raise ValueError("Found non-positive uncertainties in err_obs.")

    residuals = ymodel_samples - clean_binned[None, :]
    mse_samples = np.mean(residuals**2, axis=1)
    rmse_samples = np.sqrt(mse_samples)
    chi2_samples = np.sum((residuals / err_obs[None, :]) ** 2, axis=1)

    if dof > 0:
        chi2r_samples = chi2_samples / float(dof)
    else:
        chi2r_samples = np.full_like(chi2_samples, np.nan, dtype=float)

    summary: dict[str, float | int] = {
        "posterior_samples_used": int(ymodel_samples.shape[0]),
    }
    summary.update(metric_summary(mse_samples, "posterior_MSE"))
    summary.update(metric_summary(rmse_samples, "posterior_RMSE"))
    summary.update(metric_summary(1.0e6 * rmse_samples, "posterior_RMSE_ppm"))
    summary.update(metric_summary(chi2_samples, "posterior_chi2"))
    summary.update(metric_summary(chi2r_samples, "posterior_chi2_reduced"))
    return summary


def find_multinest_basename(output_dir: Path, planet_name: str, model_name: str) -> str:
    """Resolve the absolute MultiNest basename needed by retrieved_samples."""
    candidates = [
        output_dir / "POSEIDON_output" / planet_name / "retrievals" / "MultiNest_raw" / model_name,
        output_dir / "POSEIDON_output" / planet_name / "retrievals" / model_name,
        BASE_DIR / "POSEIDON_output" / planet_name / "retrievals" / "MultiNest_raw" / model_name,
        BASE_DIR / "POSEIDON_output" / planet_name / "retrievals" / model_name,
    ]

    for candidate in candidates:
        eq_file = Path(str(candidate) + "-post_equal_weights.dat")
        if eq_file.exists():
            return str(candidate.resolve())

    searched = "\n".join(str(Path(str(candidate) + "-post_equal_weights.dat")) for candidate in candidates)
    raise FileNotFoundError(
        "Could not find the MultiNest equal-weights posterior file required by "
        "retrieved_samples(...).\nSearched:\n"
        f"{searched}"
    )


def build_case_config(args: argparse.Namespace) -> CaseConfig:
    """Convert CLI arguments into one fully resolved retrieval case."""
    source_label = "sphinx" if args.branch == "sphinx" else None
    output_dir = SPHINX_DIR if args.branch == "sphinx" else BASE_DIR
    data_dir = BASE_DIR / ("observations_sphinx" if args.branch == "sphinx" else "observations")

    reconstructed = args.mode == "recon"

    observation_file = (
        args.observation_file
        if args.observation_file
        else build_observation_filename(
            args.n_transits,
            args.f_spot,
            args.f_fac,
            reconstructed=reconstructed,
            source_label=source_label,
        )
    )

    model_name = (
        args.model_name
        if args.model_name
        else build_model_name(
            args.mode,
            args.n_transits,
            args.f_spot,
            args.f_fac,
            source_label=source_label,
        )
    )

    npz_default = output_dir / "posterior_diagnostics" / f"{model_name}_posterior_samples.npz"
    csv_default = output_dir / "chi2_log_posterior.csv"
    clean_default = BASE_DIR / "pandexo_spec.txt"

    return CaseConfig(
        mode_label=args.mode,
        branch=args.branch,
        source_label=source_label,
        n_transits=args.n_transits,
        f_spot_case=args.f_spot,
        f_fac_case=args.f_fac,
        posterior_samples=args.posterior_samples,
        n_free_params=args.n_free_params,
        clean_path=_resolve_path(args.clean_path, clean_default),
        csv_path=_resolve_path(args.csv_path, csv_default),
        npz_path=_resolve_path(args.npz_path, npz_default),
        data_dir=data_dir.resolve(),
        output_dir=output_dir.resolve(),
        observation_file=observation_file,
        model_name=model_name,
    )


def build_poseidon_context(case: CaseConfig) -> dict[str, Any]:
    """Create the POSEIDON objects needed to resample posterior spectra."""
    from POSEIDON.constants import M_E, R_E, R_Sun
    from POSEIDON.core import (
        create_planet,
        create_star,
        define_model,
        load_data,
        read_opacities,
        set_priors,
        wl_grid_constant_R,
    )

    r_star = 0.1192 * R_Sun
    t_star = 2566.0
    metallicity_star = 0.00
    log_g_star = 5.2396

    if case.mode_label == "contam":
        star = create_star(
            r_star,
            t_star,
            log_g_star,
            metallicity_star,
            stellar_grid="phoenix",
            stellar_contam="two_spots",
            f_spot=case.f_spot_case,
            T_spot=0.86 * t_star,
            f_fac=case.f_fac_case,
            T_fac=t_star + 100.0,
        )
    else:
        star = create_star(
            r_star,
            t_star,
            log_g_star,
            metallicity_star,
            stellar_grid="phoenix",
        )

    r_planet = 0.917985 * R_E
    m_planet = 0.6356 * M_E
    t_eq_planet = 255.0

    planet = create_planet(
        PLANET_NAME,
        r_planet,
        mass=m_planet,
        T_eq=t_eq_planet,
    )

    wl_model = wl_grid_constant_R(WL_MIN_UM, WL_MAX_UM, R_GRID)
    data = load_data(
        str(case.data_dir),
        datasets=[case.observation_file],
        instruments=INSTRUMENTS,
        wl_model=wl_model,
    )

    bulk_species = ["N2"]
    param_species = ["H2O", "CH4", "CO2", "O3"]

    model_kwargs: dict[str, Any] = {
        "PT_profile": "isotherm",
        "cloud_model": "cloud-free",
    }
    if case.mode_label == "contam":
        model_kwargs["stellar_contam"] = "two_spots"

    model = define_model(
        case.model_name,
        bulk_species,
        param_species,
        **model_kwargs,
    )

    prior_types: dict[str, str] = {
        "T": "uniform",
        "R_p_ref": "uniform",
        "log_H2O": "uniform",
        "log_CH4": "uniform",
        "log_CO2": "uniform",
        "log_O3": "uniform",
    }
    prior_ranges: dict[str, list[float]] = {
        "T": [200.0, 400.0],
        "R_p_ref": [0.85 * r_planet, 1.15 * r_planet],
        "log_H2O": [-8.0, -1.0],
        "log_CH4": [-8.0, -1.0],
        "log_CO2": [-5.0, -1.0],
        "log_O3": [-8.0, -1.0],
    }

    if case.mode_label == "contam":
        prior_types.update(
            {
                "f_spot": "uniform",
                "f_fac": "uniform",
                "T_phot": "uniform",
                "T_fac": "uniform",
                "T_spot": "uniform",
            }
        )
        prior_ranges.update(
            {
                "f_spot": [0.0, 0.26],
                "f_fac": [0.0, 0.70],
                "T_phot": [0.9 * t_star, 1.1 * t_star],
                "T_fac": [t_star, t_star + 150.0],
                "T_spot": [0.8 * t_star, 0.9 * t_star],
            }
        )

    priors = set_priors(planet, star, model, data, prior_types, prior_ranges)
    opac = read_opacities(model, wl_model, OPACITY_TREATMENT, T_FINE, LOG_P_FINE)

    p_grid = np.logspace(np.log10(P_MAX_BAR), np.log10(P_MIN_BAR), N_LAYERS)

    return {
        "planet": planet,
        "star": star,
        "model": model,
        "priors": priors,
        "opac": opac,
        "data": data,
        "wl_model": wl_model,
        "p_grid": p_grid,
        "p_ref": P_REF_BAR,
        "r_planet": r_planet,
    }


def sample_posterior_spectra(case: CaseConfig, context: dict[str, Any]) -> dict[str, np.ndarray | str]:
    """Regenerate posterior spectra from saved MultiNest samples."""
    from POSEIDON.retrieval import retrieved_samples

    retrieval_basename = find_multinest_basename(case.output_dir, PLANET_NAME, case.model_name)

    (
        _t_low2,
        _t_low1,
        _t_median,
        _t_high1,
        _t_high2,
        _log_x_low2,
        _log_x_low1,
        _log_x_median,
        _log_x_high1,
        _log_x_high2,
        spec_low2,
        spec_low1,
        spec_median,
        spec_high1,
        spec_high2,
        _t_best,
        spectrum_best,
        ymodel_best,
        ymodel_samples,
    ) = retrieved_samples(
        planet=context["planet"],
        star=context["star"],
        model=context["model"],
        opac=context["opac"],
        data=context["data"],
        retrieval_name=retrieval_basename,
        wl=context["wl_model"],
        P=context["p_grid"],
        P_ref_set=context["p_ref"],
        R_p_ref_set=context["r_planet"],
        P_param_set=context["p_ref"],
        He_fraction=0.17,
        N_slice_EM=2,
        N_slice_DN=4,
        spectrum_type="transmission",
        T_phot_grid=None,
        T_het_grid=None,
        log_g_phot_grid=None,
        log_g_het_grid=None,
        I_phot_grid=None,
        I_het_grid=None,
        y_p=np.array([0.0]),
        F_s_obs=None,
        constant_gravity=False,
        chemistry_grid=None,
        N_output_samples=case.posterior_samples,
    )

    return {
        "retrieval_basename": retrieval_basename,
        "wl_out": np.asarray(context["wl_model"], dtype=float),
        "spec_low2": np.asarray(spec_low2, dtype=float),
        "spec_low1": np.asarray(spec_low1, dtype=float),
        "spec_median": np.asarray(spec_median, dtype=float),
        "spec_high1": np.asarray(spec_high1, dtype=float),
        "spec_high2": np.asarray(spec_high2, dtype=float),
        "spectrum_best": np.asarray(spectrum_best, dtype=float),
        "ymodel_best": np.asarray(ymodel_best, dtype=float),
        "ymodel_samples": np.asarray(ymodel_samples, dtype=float),
    }


def build_metrics_row(
    case: CaseConfig,
    median_metrics: dict[str, float | int],
    posterior_metrics: dict[str, float | int],
) -> dict[str, Any]:
    """Build the CSV row with legacy and posterior summary columns."""
    row: dict[str, Any] = {
        "planet_name": PLANET_NAME,
        "model_name": case.model_name,
        "observation": case.observation_file,
        "branch": case.branch,
        "mode_label": case.mode_label,
        "source_label": case.source_label or "phoenix",
        "n_transits": case.n_transits,
        "spot_fraction": case.f_spot_case,
        "facula_fraction": case.f_fac_case,
        "metrics_basis": "median_retrieved_spectrum",
        "posterior_npz": str(case.npz_path),
    }
    row.update(median_metrics)
    row.update(posterior_metrics)
    return row


def upsert_metrics_csv(log_path: Path, row: dict[str, Any]) -> pd.DataFrame:
    """Insert or replace one row keyed by planet/model/observation/branch."""
    columns = [
        "planet_name",
        "model_name",
        "observation",
        "branch",
        "mode_label",
        "source_label",
        "n_transits",
        "spot_fraction",
        "facula_fraction",
        "metrics_basis",
        "N",
        "p",
        "dof",
        "MSE",
        "rmse",
        "rmse_ppm",
        "chi2",
        "chi2_reduced",
        "posterior_samples_used",
        "posterior_MSE_mean",
        "posterior_MSE_std",
        "posterior_MSE_median",
        "posterior_MSE_q16",
        "posterior_MSE_q84",
        "posterior_RMSE_mean",
        "posterior_RMSE_std",
        "posterior_RMSE_median",
        "posterior_RMSE_q16",
        "posterior_RMSE_q84",
        "posterior_RMSE_ppm_mean",
        "posterior_RMSE_ppm_std",
        "posterior_RMSE_ppm_median",
        "posterior_RMSE_ppm_q16",
        "posterior_RMSE_ppm_q84",
        "posterior_chi2_mean",
        "posterior_chi2_std",
        "posterior_chi2_median",
        "posterior_chi2_q16",
        "posterior_chi2_q84",
        "posterior_chi2_reduced_mean",
        "posterior_chi2_reduced_std",
        "posterior_chi2_reduced_median",
        "posterior_chi2_reduced_q16",
        "posterior_chi2_reduced_q84",
        "posterior_npz",
    ]

    log_path.parent.mkdir(parents=True, exist_ok=True)
    new_row = pd.DataFrame([{column: row.get(column, np.nan) for column in columns}], columns=columns)

    if log_path.exists():
        df_log = pd.read_csv(log_path)
        for column in columns:
            if column not in df_log.columns:
                df_log[column] = np.nan
        df_log = df_log[columns]

        key_mask = (
            (df_log["planet_name"].astype(str) == str(row["planet_name"]))
            & (df_log["model_name"].astype(str) == str(row["model_name"]))
            & (df_log["observation"].astype(str) == str(row["observation"]))
            & (df_log["branch"].astype(str) == str(row["branch"]))
        )

        if key_mask.any():
            df_log = df_log.loc[~key_mask].copy()

        df_log = pd.concat([df_log, new_row], ignore_index=True)
    else:
        df_log = new_row

    df_log.to_csv(log_path, index=False, float_format="%.10g")
    return df_log


def save_posterior_archive(
    case: CaseConfig,
    sampled: dict[str, np.ndarray | str],
    wl_data: np.ndarray,
    half_bin: np.ndarray,
    clean_binned: np.ndarray,
    err_obs: np.ndarray,
) -> None:
    """Write a compact archive with propagated spectra and grid metadata."""
    case.npz_path.parent.mkdir(parents=True, exist_ok=True)

    metadata = {
        "planet_name": PLANET_NAME,
        "model_name": case.model_name,
        "observation": case.observation_file,
        "branch": case.branch,
        "mode_label": case.mode_label,
        "n_transits": case.n_transits,
        "spot_fraction": case.f_spot_case,
        "facula_fraction": case.f_fac_case,
        "posterior_samples": case.posterior_samples,
        "retrieval_basename": sampled["retrieval_basename"],
    }

    np.savez_compressed(
        case.npz_path,
        metadata_json=json.dumps(metadata, indent=2),
        wl_model=np.asarray(sampled["wl_out"], dtype=float),
        wl_data=np.asarray(wl_data, dtype=float),
        half_bin=np.asarray(half_bin, dtype=float),
        clean_binned=np.asarray(clean_binned, dtype=float),
        err_obs=np.asarray(err_obs, dtype=float),
        spec_low2=np.asarray(sampled["spec_low2"], dtype=float),
        spec_low1=np.asarray(sampled["spec_low1"], dtype=float),
        spec_median=np.asarray(sampled["spec_median"], dtype=float),
        spec_high1=np.asarray(sampled["spec_high1"], dtype=float),
        spec_high2=np.asarray(sampled["spec_high2"], dtype=float),
        spectrum_best=np.asarray(sampled["spectrum_best"], dtype=float),
        ymodel_best=np.asarray(sampled["ymodel_best"], dtype=float),
        ymodel_samples=np.asarray(sampled["ymodel_samples"], dtype=float),
    )


def run_case(case: CaseConfig, *, verbose: bool = True) -> dict[str, Any]:
    """Run posterior propagation for one case and persist the diagnostics."""
    from POSEIDON.instrument import bin_spectrum_to_data
    from POSEIDON.utility import read_data

    if not case.clean_path.is_file():
        raise FileNotFoundError(f"Clean reference spectrum not found: {case.clean_path}")

    obs_path = case.data_dir / case.observation_file
    if not obs_path.is_file():
        raise FileNotFoundError(f"Observation file not found: {obs_path}")

    context = build_poseidon_context(case)
    sampled = sample_posterior_spectra(case, context)

    wl_clean, y_clean = load_clean_two_cols(case.clean_path)
    wl_data, half_bin, _y_obs, err_obs = read_data(str(case.data_dir), case.observation_file)
    clean_binned = bin_average_with_halfbins(wl_clean, y_clean, wl_data, half_bin)

    median_model_binned = bin_spectrum_to_data(
        np.asarray(sampled["spec_median"], dtype=float),
        np.asarray(sampled["wl_out"], dtype=float),
        context["data"],
    )
    median_metrics = compute_scalar_metrics(
        median_model_binned,
        clean_binned,
        err_obs,
        n_free_params=case.n_free_params,
    )

    posterior_metrics = summarize_posterior_metrics(
        np.asarray(sampled["ymodel_samples"], dtype=float),
        clean_binned,
        err_obs,
        dof=int(median_metrics["dof"]),
    )

    save_posterior_archive(case, sampled, wl_data, half_bin, clean_binned, err_obs)
    row = build_metrics_row(case, median_metrics, posterior_metrics)
    df_log = upsert_metrics_csv(case.csv_path, row)

    if verbose:
        print("---- Posterior diagnostics summary ----")
        print(f"Model name               : {case.model_name}")
        print(f"Observation              : {case.observation_file}")
        print(f"Branch                   : {case.branch}")
        print(f"Median-spectrum MSE      : {float(median_metrics['MSE']):.6e}")
        print(f"Median-spectrum chi2_red : {float(median_metrics['chi2_reduced']):.6f}")
        print(f"Posterior samples used   : {int(posterior_metrics['posterior_samples_used'])}")
        print(
            "Posterior MSE median     : "
            f"{float(posterior_metrics['posterior_MSE_median']):.6e}"
        )
        print(
            "Posterior chi2_red median: "
            f"{float(posterior_metrics['posterior_chi2_reduced_median']):.6f}"
        )
        print(f"Saved posterior archive  : {case.npz_path}")
        print(f"Updated metrics log      : {case.csv_path}")

    return {
        "case": case,
        "median_metrics": median_metrics,
        "posterior_metrics": posterior_metrics,
        "log_dataframe": df_log,
        "posterior_archive": case.npz_path,
    }


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser used by the diagnostics runner."""
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate posterior sampled spectra from a POSEIDON retrieval and "
            "propagate MSE / reduced chi-squared against the clean truth spectrum."
        )
    )
    parser.add_argument(
        "--mode",
        choices=("uncontam", "contam", "recon"),
        default=os.getenv("POSTERIOR_MODE", os.getenv("RETRIEVAL_MODE_LABEL", "uncontam")),
        help="Retrieval strategy label used in the existing project naming.",
    )
    parser.add_argument(
        "--branch",
        choices=("phoenix", "sphinx"),
        default=os.getenv("POSTERIOR_BRANCH", "phoenix"),
        help="Choose baseline PHOENIX outputs or the sphinx_injection branch.",
    )
    parser.add_argument(
        "--n-transits",
        dest="n_transits",
        type=int,
        default=_env_or_default("RETRIEVAL_N_TRANSITS", int, 10),
        help="Number of transits for the case naming.",
    )
    parser.add_argument(
        "--f-spot",
        dest="f_spot",
        type=float,
        default=_env_or_default("RETRIEVAL_F_SPOT_CASE", float, 0.26),
        help="Injected spot covering fraction for the case naming.",
    )
    parser.add_argument(
        "--f-fac",
        dest="f_fac",
        type=float,
        default=_env_or_default("RETRIEVAL_F_FAC_CASE", float, 0.70),
        help="Injected facula covering fraction for the case naming.",
    )
    parser.add_argument(
        "--posterior-samples",
        type=int,
        default=_env_or_default("POSTERIOR_SAMPLES", int, DEFAULT_POSTERIOR_SAMPLES),
        help="Number of posterior spectra to regenerate with retrieved_samples(...).",
    )
    parser.add_argument(
        "--n-free-params",
        type=int,
        default=_env_or_default("POSTERIOR_N_FREE_PARAMS", int, DEFAULT_N_FREE_PARAMS),
        help="Number of free parameters used when computing the degrees of freedom.",
    )
    parser.add_argument(
        "--clean-path",
        default=os.getenv("POSTERIOR_CLEAN_PATH"),
        help="Optional override for the clean truth spectrum path.",
    )
    parser.add_argument(
        "--csv-path",
        default=os.getenv("POSTERIOR_CSV_PATH"),
        help="Optional override for the metrics CSV output path.",
    )
    parser.add_argument(
        "--npz-path",
        default=os.getenv("POSTERIOR_NPZ_PATH"),
        help="Optional override for the posterior samples archive path.",
    )
    parser.add_argument(
        "--observation-file",
        default=os.getenv("POSTERIOR_OBSERVATION_FILE"),
        help="Optional explicit observation filename. Defaults to the existing naming convention.",
    )
    parser.add_argument(
        "--model-name",
        default=os.getenv("POSTERIOR_MODEL_NAME"),
        help="Optional explicit model name. Defaults to the existing naming convention.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=_env_or_default("POSTERIOR_QUIET", _parse_bool, False),
        help="Reduce console output.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    case = build_case_config(args)
    run_case(case, verbose=not args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
