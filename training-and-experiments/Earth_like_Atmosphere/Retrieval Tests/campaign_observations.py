#!/usr/bin/env python3
"""Generate PandExo observations and G-DAE reconstructions for the campaign.

The script starts from the clean spectrum in `pandexo_spec.txt`, applies the
configured stellar-contamination curve when needed, runs PandExo to produce a
noisy JWST/NIRSpec Prism observation, and then saves the matching G-DAE
reconstruction used by the `gdae` retrieval strategy.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from campaign_common import (
    CAMPAIGN_DIR,
    CONTAMINATION_DIR,
    EARTH_DIR,
    GDAE_MODEL_PATH,
    N_TRANSITS,
    CLEAN_SPECTRUM_PATH,
    CaseConfig,
    branch_dir,
    get_case,
    iter_cases,
    normalize_branch,
    normalize_test_id,
    observations_dir,
    pandexo_inputs_dir,
)

CONTAMINATION_PATTERN = re.compile(r"fspot(?P<f_spot>[0-9.]+)_ffac(?P<f_fac>[0-9.]+)\.txt$")
TRIM_IDX = 18  # Drop the short-wavelength bins PandExo returns outside the useful range.
DEFAULT_REF_FLAT: float | None = None
REF_FLAT_EPS = 1e-12
REF_FLAT_SOURCE_PATH = EARTH_DIR / 'spec_data' / 'airless_data.csv'

DEFAULT_UQ_PASSES = 100
DEFAULT_UQ_SEED = 12345

def load_two_column_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a wavelength/depth or wavelength/epsilon table sorted by wavelength."""
    arr = np.loadtxt(path)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}, got {arr.shape}.")
    wl = arr[:, 0].astype(np.float64)
    depth = arr[:, 1].astype(np.float64)
    order = np.argsort(wl)
    return wl[order], depth[order]


def load_epsilon_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    wl, epsilon = load_two_column_spectrum(path)
    return wl, epsilon


def build_contamination_index() -> dict[tuple[str, float, float], Path]:
    """Index PHOENIX and SPHINX epsilon curves by source and coverage fractions."""
    index: dict[tuple[str, float, float], Path] = {}
    for path in sorted(CONTAMINATION_DIR.glob("*TRAPPIST-1_contam_fspot*.txt")):
        match = CONTAMINATION_PATTERN.search(path.name)
        if not match:
            continue
        source = "sphinx" if path.name.startswith("sphinx_") else "original"
        f_spot = float(match.group("f_spot"))
        f_fac = float(match.group("f_fac"))
        index[(source, f_spot, f_fac)] = path
    return index


def find_contamination_curve(case: CaseConfig) -> Path:
    if case.contam_source == "clean":
        raise ValueError("Clean cases do not use a contamination curve.")
    index = build_contamination_index()
    key = (case.contam_source, case.f_spot, case.f_fac)
    if key not in index:
        raise KeyError(
            "No contamination curve found for "
            f"source={case.contam_source}, f_spot={case.f_spot}, f_fac={case.f_fac}."
        )
    return index[key]


def save_contaminated_input_spectrum(case: CaseConfig, test_id: str) -> Path:
    """Create the two-column spectrum PandExo will use for one campaign case."""
    output_dir = pandexo_inputs_dir(test_id, case.branch)
    output_dir.mkdir(parents=True, exist_ok=True)
    wl_clean, depth_clean = load_two_column_spectrum(CLEAN_SPECTRUM_PATH)

    if case.contam_source == "clean":
        depth_out = depth_clean
    else:
        wl_eps, epsilon = load_epsilon_curve(find_contamination_curve(case))
        epsilon_interp = np.interp(
            wl_clean,
            wl_eps,
            epsilon,
            left=epsilon[0],
            right=epsilon[-1],
        )
        depth_out = depth_clean * epsilon_interp

    out_path = output_dir / f"{case.observation_stem}_input.txt"
    np.savetxt(out_path, np.column_stack((wl_clean, depth_out)), fmt="%.10e")
    return out_path


def save_pandexo_spectrum_to_dat(input_spec_path: Path, output_dat_path: Path) -> None:
    """Run PandExo once and save a POSEIDON-compatible four-column observation.

    Output columns are wavelength, half-bin width, noisy transit depth, and
    transit-depth uncertainty.
    """
    import pandexo.engine.justdoit as jdi

    exo_dict = jdi.load_exo_dict()
    exo_dict["observation"].update(
        {
            "sat_level": 80,
            "sat_unit": "%",
            "baseline_unit": "total",
            "baseline": 0.9535 * 3 * 60 * 60,
            "noise_floor": 0,
            "noccultations": N_TRANSITS,
        }
    )
    exo_dict["star"].update(
        {
            "type": "phoenix",
            "mag": 11.354,
            "ref_wave": 1.25,
            "temp": 2566,
            "metal": 0.0,
            "logg": 5.2396,
        }
    )
    exo_dict["planet"].update(
        {
            "type": "user",
            "w_unit": "um",
            "f_unit": "rp^2/r*^2",
            "transit_duration": 0.9535 * 60 * 60,
            "td_unit": "s",
            "exopath": str(input_spec_path),
        }
    )

    inst_dict = jdi.load_mode_dict("NIRSpec Prism")
    inst_dict["configuration"]["detector"].update({"subarray": "sub512", "ngroup": 6})

    results = jdi.run_pandexo(exo_dict, inst_dict)
    final_spec = results["FinalSpectrum"]

    waves_trim = final_spec["wave"][TRIM_IDX:]
    spec_rand_trim = final_spec["spectrum_w_rand"][TRIM_IDX:]
    err_trim = final_spec["error_w_floor"][TRIM_IDX:]
    wave_err = np.gradient(waves_trim) / 2.0

    output_dat_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output_dat_path,
        np.column_stack((waves_trim, wave_err, spec_rand_trim, err_trim)),
        fmt="%.6e",
    )


def bin_average_with_halfbins(
    wl_src: np.ndarray,
    y_src: np.ndarray,
    centers: np.ndarray,
    halfwidths: np.ndarray,
    nsamp: int = 256,
) -> np.ndarray:
    """Average a high-resolution spectrum over observation bins."""
    wl_src = np.asarray(wl_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    centers = np.asarray(centers, dtype=np.float64)
    halfwidths = np.asarray(halfwidths, dtype=np.float64)

    order = np.argsort(wl_src)
    wl_sorted = wl_src[order]
    y_sorted = y_src[order]

    out = np.empty_like(centers, dtype=np.float64)
    for idx, (center, halfwidth) in enumerate(zip(centers, halfwidths)):
        left = center - halfwidth
        right = center + halfwidth
        sample_wl = np.linspace(left, right, nsamp)
        sample_flux = np.interp(sample_wl, wl_sorted, y_sorted)
        out[idx] = np.trapz(sample_flux, sample_wl) / (right - left)
    return out.astype(np.float32)


def infer_ref_flat_from_airless(
    path: Path = REF_FLAT_SOURCE_PATH,
    *,
    atol: float = 1e-12,
    rtol: float = 1e-8,
) -> float:
    """Infer a scalar flat reference depth from the experiment airless spectrum."""
    n_points = int(np.loadtxt(EARTH_DIR / 'waves.txt').shape[0])
    airless_df = pd.read_csv(path)
    airless_numeric = airless_df.select_dtypes(include=[np.number])
    depth = airless_numeric.iloc[0, -n_points:].to_numpy(dtype=np.float64)
    ref_flat = float(np.median(depth))
    if not np.allclose(depth, ref_flat, atol=atol, rtol=rtol):
        raise ValueError("Airless spectrum is not flat enough to define a scalar ref_flat.")
    if not np.isfinite(ref_flat) or ref_flat <= 0.0:
        raise ValueError(f"Invalid ref_flat inferred from {path}: {ref_flat}")
    return ref_flat


def resolve_ref_flat(ref_flat: float | None = DEFAULT_REF_FLAT) -> float:
    """Return an explicit scalar ref_flat or infer it from the experiment airless spectrum."""
    if ref_flat is None:
        return infer_ref_flat_from_airless()

    value = float(ref_flat)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"ref_flat must be a positive finite float, got {ref_flat}")
    return value


def to_log_ratio_1d(
    values: np.ndarray,
    ref_flat: float | None = DEFAULT_REF_FLAT,
    eps: float = REF_FLAT_EPS,
) -> np.ndarray:
    """Transform physical transit depths to log(depth / ref_flat)."""
    resolved_ref_flat = resolve_ref_flat(ref_flat)
    values = np.asarray(values, dtype=np.float32)
    values = np.clip(values, eps, None)
    return np.log(values / resolved_ref_flat).astype(np.float32)


def from_log_ratio_1d(
    values_log: np.ndarray,
    ref_flat: float | None = DEFAULT_REF_FLAT,
) -> np.ndarray:
    """Map log(depth / ref_flat) back to physical transit depth."""
    resolved_ref_flat = resolve_ref_flat(ref_flat)
    values_log = np.asarray(values_log, dtype=np.float32)
    return (np.exp(values_log) * resolved_ref_flat).astype(np.float32)


def load_autoencoder(model_path: Path = GDAE_MODEL_PATH) -> Any:
    from tensorflow import keras

    return keras.models.load_model(model_path)


def propagate_reconstruction_uncertainty(
    model: Any,
    y_obs: np.ndarray,
    y_err: np.ndarray,
    ref_flat: float,
    *,
    n_passes: int = DEFAULT_UQ_PASSES,
    seed: int | None = DEFAULT_UQ_SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate the notebook-consistent predictive reconstruction and error.

    MC dropout is evaluated repeatedly for the fixed observed spectrum. Every
    prediction is transformed back to physical transit-depth space before its
    mean and epistemic standard deviation are calculated. The total diagonal
    error follows the analysis notebook contract:

        sigma_total**2 = sigma_epi**2 + (0.5 * sigma_inst)**2

    The campaign files do not contain a wavelength covariance matrix, so only
    the diagonal of this predictive uncertainty is written to ``*_recon.dat``.
    """
    y_obs = np.asarray(y_obs, dtype=np.float32)
    y_err = np.asarray(y_err, dtype=np.float32)
    if y_obs.ndim != 1 or y_err.shape != y_obs.shape:
        raise ValueError(
            f"Expected 1D y_obs/y_err with matching shapes, got {y_obs.shape} and {y_err.shape}."
        )
    if not np.all(np.isfinite(y_obs)):
        raise ValueError("Observed transit depths must be finite.")
    if np.any(~np.isfinite(y_err)) or np.any(y_err <= 0.0):
        raise ValueError("Observed uncertainties must be finite and strictly positive.")
    if int(n_passes) < 2:
        raise ValueError(f"n_passes must be at least 2 to estimate a standard deviation, got {n_passes}.")

    resolved_ref_flat = resolve_ref_flat(ref_flat)
    if seed is not None:
        from tensorflow import keras

        keras.utils.set_random_seed(seed)

    observed_log = to_log_ratio_1d(y_obs, resolved_ref_flat).reshape(1, -1)
    repeated_log = np.repeat(observed_log, int(n_passes), axis=0)
    predictions_log = model(repeated_log, training=True).numpy()
    predictions = from_log_ratio_1d(predictions_log, resolved_ref_flat)

    # Average after returning every realization to physical transit-depth
    # space. In general, mean(exp(x)) is not exp(mean(x)).
    recon_mean = np.mean(predictions, axis=0, dtype=np.float64).astype(np.float32)
    sigma_epi = np.std(predictions, axis=0, ddof=1).astype(np.float32)
    sigma_ale = np.float32(0.5) * y_err
    recon_sigma = np.sqrt(sigma_epi**2 + sigma_ale**2).astype(np.float32)
    recon_sigma = np.maximum(recon_sigma, REF_FLAT_EPS)
    return recon_mean, recon_sigma


def reconstruct_observation_file(
    observation_path: Path,
    output_dir: Path,
    model_path: Path = GDAE_MODEL_PATH,
    *,
    ref_flat: float | None = DEFAULT_REF_FLAT,
    uq_passes: int = DEFAULT_UQ_PASSES,
    uq_seed: int | None = DEFAULT_UQ_SEED,
) -> Path:
    """Reconstruct one noisy observation and propagate its uncertainty."""
    observation = np.loadtxt(observation_path)
    if observation.ndim != 2 or observation.shape[1] < 4:
        raise ValueError(f"Expected four columns in {observation_path}, got {observation.shape}.")

    wl = observation[:, 0].astype(np.float32)
    d_wl = observation[:, 1].astype(np.float32)
    y_obs = observation[:, 2].astype(np.float32)
    y_err = observation[:, 3].astype(np.float32)

    wl_clean, y_clean = load_two_column_spectrum(CLEAN_SPECTRUM_PATH)
    y_clean_binned = bin_average_with_halfbins(wl_clean, y_clean, wl, d_wl)
    # Contract log_ratio_ref_flat_v1: the Earth-like G-DAE was trained with
    # wavelengths in strictly descending order. Enforce the contract instead
    # of relying on the input file's incidental ordering.
    descending = np.argsort(wl)[::-1]
    wl = wl[descending].copy()
    d_wl = d_wl[descending].copy()
    y_obs = y_obs[descending].copy()
    y_err = y_err[descending].copy()
    y_clean_binned = y_clean_binned[descending].copy()

    if uq_seed is not None:
        from tensorflow import keras

        # Keras creates the Dropout seed generators while deserializing the
        # model, so seed before loading it rather than only before inference.
        keras.utils.set_random_seed(uq_seed)

    model = load_autoencoder(model_path)
    resolved_ref_flat = resolve_ref_flat(ref_flat)
    y_recon, y_recon_err = propagate_reconstruction_uncertainty(
        model,
        y_obs,
        y_err,
        resolved_ref_flat,
        n_passes=uq_passes,
        seed=uq_seed,
    )

    recon_path = output_dir / f"{observation_path.stem}_recon.dat"
    np.savetxt(
        recon_path,
        np.column_stack((wl, d_wl, y_recon, y_recon_err)),
        fmt="%.10e",
    )
    return recon_path


def export_case_observation(
    case: CaseConfig,
    test_id: str,
    *,
    overwrite: bool = False,
    reconstruct_only: bool = False,
    uq_passes: int = DEFAULT_UQ_PASSES,
    uq_seed: int | None = DEFAULT_UQ_SEED,
) -> tuple[Path, Path]:
    """Generate the noisy observation and reconstruction for one configured case."""
    obs_dir = observations_dir(test_id, case.branch)
    obs_dir.mkdir(parents=True, exist_ok=True)
    observation_path = obs_dir / case.observation_file
    recon_path = obs_dir / case.reconstruction_file

    if observation_path.exists() and recon_path.exists() and not overwrite:
        print(f"Skipping existing observation and reconstruction: {test_id}/{case.branch}/{case.case_id}")
        return observation_path, recon_path

    if not reconstruct_only:
        input_path = save_contaminated_input_spectrum(case, test_id)
        save_pandexo_spectrum_to_dat(input_path, observation_path)
        print(f"Saved observation     : {observation_path}")
    elif not observation_path.exists():
        raise FileNotFoundError(f"Cannot reconstruct missing observation: {observation_path}")

    recon_path = reconstruct_observation_file(
        observation_path,
        obs_dir,
        uq_passes=uq_passes,
        uq_seed=uq_seed,
    )
    print(f"Saved reconstruction : {recon_path}")
    return observation_path, recon_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-id", required=True, help="Campaign test id, e.g. test_02.")
    parser.add_argument("--branch", choices=("phoenix", "sphinx", "all"), default="all")
    parser.add_argument("--f-spot", type=float, default=None)
    parser.add_argument("--f-fac", type=float, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--reconstruct-only",
        action="store_true",
        help="Do not run PandExo; rebuild *_recon.dat from existing observations.",
    )
    parser.add_argument(
        "--uq-passes",
        type=int,
        default=DEFAULT_UQ_PASSES,
        help=f"Monte Carlo passes for propagated reconstruction errors (default: {DEFAULT_UQ_PASSES}).",
    )
    parser.add_argument(
        "--uq-seed",
        type=int,
        default=DEFAULT_UQ_SEED,
        help=f"NumPy and Keras seed for reproducible uncertainty propagation (default: {DEFAULT_UQ_SEED}).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    test_id = normalize_test_id(args.test_id)
    (CAMPAIGN_DIR / "plots").mkdir(parents=True, exist_ok=True)

    if args.branch == "all":
        branches = ("phoenix", "sphinx")
    else:
        branches = (normalize_branch(args.branch),)

    selected_cases: list[CaseConfig] = []
    for branch in branches:
        if args.f_spot is not None or args.f_fac is not None:
            if args.f_spot is None or args.f_fac is None:
                raise ValueError("--f-spot and --f-fac must be supplied together.")
            selected_cases.append(get_case(branch, args.f_spot, args.f_fac))
        else:
            selected_cases.extend(iter_cases(branch))

    for case in selected_cases:
        branch_dir(test_id, case.branch).mkdir(parents=True, exist_ok=True)
        export_case_observation(
            case,
            test_id,
            overwrite=args.overwrite,
            reconstruct_only=args.reconstruct_only,
            uq_passes=args.uq_passes,
            uq_seed=args.uq_seed,
        )


if __name__ == "__main__":
    main()
