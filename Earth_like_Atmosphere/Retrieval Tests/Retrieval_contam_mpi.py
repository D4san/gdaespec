#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POSEIDON — Retrieval Pipeline for Contam strategy (stellar contamination).
This script is configured for MPI parallelization and models stellar spots.

Usage example:
    mpirun -n 8 python -u Retrieval_contam_mpi.py
"""

import os
import time
import warnings
from pathlib import Path

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
    """Run the joint chemistry and stellar-contamination retrieval."""
    t_start = time.time()

    # --- System Definition ---
    # TRAPPIST-1 stellar properties
    r_star = 0.1192 * R_Sun
    t_star = 2566.0
    metallicity_star = 0.00
    log_g_star = 5.2396

    # Stellar model including spots and faculae for the contamination retrieval
    star = create_star(
        r_star, t_star, log_g_star, metallicity_star,
        stellar_grid="phoenix",
        stellar_contam="two_spots",
        f_spot=0.01, T_spot=0.86 * t_star,
        f_fac=0.08,  T_fac=t_star + 100.0,
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

    # Retrieval input: contaminated spectrum modeled jointly with stellar heterogeneity
    data_dir = Path("observations")
    obs_file = "pandexo_output_10transits_fspot0.00_ffac0.00.dat"
    datasets = [obs_file]
    instruments = ["JWST_NIRSpec_PRISM"]

    data = load_data(str(data_dir), datasets, instruments, wl)

    # --- Atmospheric Model ---
    # Chemistry model augmented with the two-spots stellar contamination parameterization
    model_name = "contam_10T_0.00spot-0.00fac"
    bulk_species = ["N2"]
    param_species = ["H2O", "CH4", "CO2", "O3"]

    model = define_model(
        model_name,
        bulk_species,
        param_species,
        PT_profile="isotherm",
        cloud_model="cloud-free",
        stellar_contam="two_spots",
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
        "f_spot": "uniform",
        "f_fac": "uniform",
        "T_phot": "uniform",
        "T_fac": "uniform",
        "T_spot": "uniform",
    }

    prior_ranges = {
        "T": [200.0, 400.0],
        "R_p_ref": [0.85 * r_planet, 1.15 * r_planet],
        "log_H2O": [-8.0, -1.0],
        "log_CH4": [-8.0, -1.0],
        "log_CO2": [-5.0, -1.0],
        "log_O3": [-8.0, -1.0],
        "f_spot": [0.0, 0.26],
        "f_fac": [0.0, 0.7],
        "T_phot": [0.9 * t_star, 1.1 * t_star],
        "T_fac": [t_star, t_star + 150.0],
        "T_spot": [0.8 * t_star, 0.9 * t_star],
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

        print(f">> Figures saved in: {out_dir.resolve()}", flush=True)
        print(f">> Total execution time: {(time.time()-t_start)/60:.2f} min", flush=True)

    return True


if __name__ == "__main__":
    main()
