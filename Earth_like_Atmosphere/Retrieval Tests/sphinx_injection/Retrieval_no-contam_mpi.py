#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POSEIDON — Retrieval Pipeline for the SPHINX-injected No-contam strategy.
This script is configured for MPI parallelization and ignores stellar effects.

Usage example:
    mpirun -n 8 python -u Retrieval_no-contam_mpi.py
"""

import os
import time
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
from mpi4py import MPI

# Keep plotting and linear algebra backends MPI-friendly
os.environ.setdefault("MPLBACKEND", "Agg")  # Non-interactive backend for batch runs
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# Silence noisy pysynphot warnings during batch execution
warnings.filterwarnings("ignore", category=UserWarning, module=r"pysynphot")

# MPI state used for rank-aware logging and root-only post-processing
comm = MPI.COMM_WORLD
_rank = comm.Get_rank()

BASE_DIR = Path(__file__).resolve().parent
PARENT_DIR = BASE_DIR.parent
TIMES_PATH = BASE_DIR / "Times"

MODE_LABEL = "uncontam"
INJECTION_LABEL = "sphinx"


def _get_env_setting(name: str, cast, default):
    """Read one optional environment setting and cast it to the target type."""
    raw_value = os.getenv(name)
    if raw_value is None:
        return default

    try:
        return cast(raw_value)
    except ValueError as exc:
        raise ValueError(
            f"Invalid value for {name!r}: {raw_value!r}"
        ) from exc


N_TRANSITS = _get_env_setting("RETRIEVAL_N_TRANSITS", int, 10)
F_SPOT_CASE = _get_env_setting("RETRIEVAL_F_SPOT_CASE", float, 0.26)
F_FAC_CASE = _get_env_setting("RETRIEVAL_F_FAC_CASE", float, 0.70)


def build_observation_filename(
    n_transits: int,
    f_spot_case: float,
    f_fac_case: float,
    reconstructed: bool = False,
    source_label: Optional[str] = None,
) -> str:
    """Build the observation filename for one stellar-heterogeneity case."""
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
    source_label: Optional[str] = None,
) -> str:
    """Build a POSEIDON model name consistent with the selected case."""
    prefix = f"{source_label}_" if source_label else ""
    return (
        f"{prefix}{mode_label}_{n_transits}T_"
        f"{f_spot_case:.2f}spot-{f_fac_case:.2f}fac"
    )


def ensure_dataset_exists(data_dir: Path, obs_file: str) -> None:
    """Raise a clear error when the requested observation file is missing."""
    obs_path = data_dir / obs_file
    if not obs_path.is_file():
        raise FileNotFoundError(f"Observation not found: {obs_path}")


def append_timing_entry(
    times_path: Path,
    n_transits: int,
    f_spot_case: float,
    f_fac_case: float,
    elapsed_minutes: float,
    mode_label: str,
) -> None:
    """Append one retrieval runtime entry using the local Times-file format."""
    if not times_path.exists():
        times_path.write_text("N_T     ffacfspot  Time     Mode\n", encoding="utf-8")

    with times_path.open("a", encoding="utf-8") as stream:
        stream.write(
            f"{n_transits:<7d} "
            f"{f_fac_case:.2f}-{f_spot_case:.2f}   "
            f"{elapsed_minutes:>6.2f}   "
            f"\"{mode_label}\"\n"
        )

# Override POSEIDON shared-memory allocation to use float64 buffers
import POSEIDON.utility as _U


def _shared_memory_array_force64(node_rank, node_comm, shape):
    """
    Replaces POSEIDON.utility.shared_memory_array to force float64 allocation.
    Standard implementation for MPICH using Shared_query(MPI.PROC_NULL).
    """
    dtype = np.float64
    itemsize = np.dtype(dtype).itemsize

    # Only the first rank on each node reserves the shared buffer
    nbytes = int(np.prod(shape)) * itemsize if node_rank == 0 else 0

    # Expose the buffer through an MPI shared-memory window
    win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=node_comm)

    # All ranks attach to the segment created by the first rank
    buf, _ = win.Shared_query(MPI.PROC_NULL)

    # Wrap the raw buffer as a NumPy array
    arr = np.ndarray(buffer=buf, dtype=dtype, shape=shape)

    return arr, win


# Install the allocator override before importing retrieval helpers
_U.shared_memory_array = _shared_memory_array_force64

# POSEIDON imports
from POSEIDON.constants import R_Sun, R_E, M_E
from POSEIDON.core import (
    create_star,
    create_planet,
    load_data,
    wl_grid_constant_R,
    define_model,
    set_priors,
    read_opacities,
)
from POSEIDON.retrieval import run_retrieval
from POSEIDON.utility import read_retrieved_spectrum, plot_collection
from POSEIDON.visuals import plot_spectra_retrieved
from POSEIDON.corner import generate_cornerplot


def main():
    """Run the chemistry-only retrieval on the contaminated spectrum."""
    t_start = time.time()

    os.chdir(BASE_DIR)

    if _rank == 0:
        print(
            ">> Retrieval case: "
            f"N_TRANSITS={N_TRANSITS}, "
            f"F_SPOT_CASE={F_SPOT_CASE:.2f}, "
            f"F_FAC_CASE={F_FAC_CASE:.2f}, "
            f"SOURCE={INJECTION_LABEL}",
            flush=True,
        )

    # --- System Definition ---
    # TRAPPIST-1 stellar properties
    r_star = 0.1192 * R_Sun
    t_star = 2566.0
    metallicity_star = 0.00
    log_g_star = 5.2396

    # Photosphere-only stellar model for the no-contamination baseline
    star = create_star(
        r_star, t_star, log_g_star, metallicity_star,
        stellar_grid="phoenix",
    )

    # TRAPPIST-1e planetary properties
    planet_name = "Trappist-1e"
    r_planet = 0.917985 * R_E
    m_planet = 0.6356 * M_E
    t_eq_planet = 255.0

    planet = create_planet(planet_name, r_planet, mass=m_planet, T_eq=t_eq_planet)

    # --- Data and Instrument Setup ---
    # High-resolution wavelength grid used internally by POSEIDON
    wl_min, wl_max, res = 0.4, 6.0, 4000
    wl = wl_grid_constant_R(wl_min, wl_max, res)

    # Retrieval input: SPHINX-injected spectrum interpreted with a chemistry-only model
    data_dir = PARENT_DIR / "observations_sphinx"
    obs_file = build_observation_filename(
        N_TRANSITS,
        F_SPOT_CASE,
        F_FAC_CASE,
        source_label=INJECTION_LABEL,
    )
    ensure_dataset_exists(data_dir, obs_file)
    datasets = [obs_file]
    instruments = ["JWST_NIRSpec_PRISM"]

    data = load_data(str(data_dir), datasets, instruments, wl)

    # --- Atmospheric Model ---
    # Chemistry-only model used as the contamination-agnostic baseline
    model_name = build_model_name(
        MODE_LABEL,
        N_TRANSITS,
        F_SPOT_CASE,
        F_FAC_CASE,
        source_label=INJECTION_LABEL,
    )
    bulk_species = ["N2"]
    param_species = ["H2O", "CH4", "CO2", "O3"]

    model = define_model(
        model_name,
        bulk_species,
        param_species,
        PT_profile="isotherm",
        cloud_model="cloud-free",
    )

    if _rank == 0:
        print(f"Free parameters: {model['param_names']}", flush=True)

    # --- Priors ---
    prior_types = {
        "T": "uniform",
        "R_p_ref": "uniform",
        "log_H2O": "uniform",
        "log_CH4": "uniform",
        "log_CO2": "uniform",
        "log_O3": "uniform",
    }

    prior_ranges = {
        "T": [200.0, 400.0],
        "R_p_ref": [0.85 * r_planet, 1.15 * r_planet],
        "log_H2O": [-8.0, -1.0],
        "log_CH4": [-8.0, -1.0],
        "log_CO2": [-5.0, -1.0],
        "log_O3": [-8.0, -1.0],
    }

    priors = set_priors(planet, star, model, data, prior_types, prior_ranges)

    # --- Opacity Interpolation Grid ---
    opacity_treatment = "opacity_sampling"

    # Temperature and pressure support for opacity interpolation
    t_fine = np.arange(200.0, 400.0 + 10.0, 10.0)
    log_p_fine = np.arange(-2.0, 2.0 + 0.2, 0.2)

    opac = read_opacities(model, wl, opacity_treatment, t_fine, log_p_fine)

    # --- Vertical Pressure Grid ---
    p_min, p_max, n_layers = 1.0e-2, 2.0, 100
    p_grid = np.logspace(np.log10(p_max), np.log10(p_min), n_layers)
    p_ref = 1.0  # reference pressure in bar

    # --- Retrieval ---
    if _rank == 0:
        print(">> Initializing retrieval with MultiNest...", flush=True)
    
    t_ret_start = time.time()
    
    run_retrieval(
        planet, star, model, opac, data, priors, wl, p_grid, p_ref, R=res,
        spectrum_type="transmission",
        sampling_algorithm="MultiNest",
        N_live=500,
        verbose=(_rank == 0),
    )
    
    dt_ret = time.time() - t_ret_start
    
    if _rank == 0:
        print(f">> Retrieval completed. Elapsed time: {dt_ret/60:.2f} min", flush=True)

    # Synchronize ranks before starting root-only plotting
    comm.Barrier()

    # --- Post-processing ---
    if _rank == 0:
        # Load the stored posterior spectrum envelopes
        wl_ret, s_low2, s_low1, s_med, s_high1, s_high2 = \
            read_retrieved_spectrum(planet_name, model_name)

        # Convert posterior envelopes into plotting collections
        spectra_med = plot_collection(s_med, wl_ret, collection=[])
        spectra_low1 = plot_collection(s_low1, wl_ret, collection=[])
        spectra_low2 = plot_collection(s_low2, wl_ret, collection=[])
        spectra_high1 = plot_collection(s_high1, wl_ret, collection=[])
        spectra_high2 = plot_collection(s_high2, wl_ret, collection=[])

        out_dir = Path("figures")
        out_dir.mkdir(parents=True, exist_ok=True)

        # Save the retrieved-spectrum summary figure
        fig_spec = plot_spectra_retrieved(
            spectra_med, spectra_low2, spectra_low1, spectra_high1, spectra_high2,
            planet_name, data, R_to_bin=100,
            data_labels=instruments,
            data_colour_list=["lime"],
        )
        fig_spec.savefig(
            out_dir / f"{planet_name}_{model_name}_retrieved_spectrum.png",
            dpi=180, bbox_inches="tight"
        )

        # Save the posterior corner plot when POSEIDON returns a figure handle
        corner_out = generate_cornerplot(planet, model)
        fig_corner = None
        if hasattr(corner_out, "savefig"):
            fig_corner = corner_out
        elif (isinstance(corner_out, tuple) and len(corner_out) > 0 and 
              hasattr(corner_out[0], "savefig")):
            fig_corner = corner_out[0]

        if fig_corner is not None:
            fig_corner.savefig(
                out_dir / f"{planet_name}_{model_name}_corner.png",
                dpi=180, bbox_inches="tight"
            )
        else:
            print("[corner] Plot was not returned. It might have been saved internally.")

        total_minutes = (time.time() - t_start) / 60.0
        append_timing_entry(
            TIMES_PATH,
            N_TRANSITS,
            F_SPOT_CASE,
            F_FAC_CASE,
            total_minutes,
            MODE_LABEL,
        )
        print(f">> Figures saved in: {out_dir.resolve()}", flush=True)
        print(f">> Total execution time: {total_minutes:.2f} min", flush=True)
        print(f">> Timing entry appended to: {TIMES_PATH.resolve()}", flush=True)

    return True


if __name__ == "__main__":
    main()
