#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POSEIDON — Pipeline de retrieval (contaminación estelar 'two_spots') listo para MPI.

Cómo lanzar (ejemplo):
    mpirun -n 8 python -u Retrieval_analysis_mpi.py

Notas:
- Backend matplotlib no interactivo para evitar conflictos en workers.
- Se fuerza float64 en los buffers de memoria compartida de POSEIDON con un
  monkey-patch cuya firma es compatible con tu versión (sin keyword dtype=).
- Parche robusto para MPICH: usa Shared_query(MPI.PROC_NULL) y evita assert.
"""

# ------------------------- Ajustes de entorno seguros -------------------------
import os
os.environ.setdefault("MPLBACKEND", "Agg")  # backend no interactivo para figuras
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# (Opcional) silenciar warnings de pysynphot
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module=r"pysynphot")

# --------------------------------- Imports -----------------------------------
import time
import numpy as np
from pathlib import Path
from mpi4py import MPI

# Comunicator global y rank
comm = MPI.COMM_WORLD
_rank = comm.Get_rank()

# -------- Parche: shared_memory_array compatible que fuerza float64 -----------
import POSEIDON.utility as _U

def _shared_memory_array_force64(node_rank, node_comm, shape):
    """
    Sustituye a POSEIDON.utility.shared_memory_array, preservando la firma original.
    Siempre asigna memoria compartida como float64 y devuelve (array, window).
    Implementación robusta en MPICH: usa Shared_query(MPI.PROC_NULL).
    """
    dtype = np.float64
    itemsize = np.dtype(dtype).itemsize

    # Solo el rank 0 del *nodo* reserva bytes; los demás 0 y consultan
    nbytes = int(np.prod(shape)) * itemsize if node_rank == 0 else 0

    # Crear ventana de memoria compartida
    win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=node_comm)

    # Consultar el segmento con tamaño no-cero (PROC_NULL = "el que tenga datos")
    buf, _disp_unit = win.Shared_query(MPI.PROC_NULL)

    # Construir el ndarray con dtype= float64 sobre ese buffer
    arr = np.ndarray(buffer=buf, dtype=dtype, shape=shape)

    return arr, win

# Activar el parche
_U.shared_memory_array = _shared_memory_array_force64

# ------------------------ Resto de imports de POSEIDON ------------------------
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
    t0 = time.time()

    # =========================== Estrella y Planeta ============================
    R_s   = 0.1192 * R_Sun
    T_s   = 2566.0
    Met_s = 0.00
    log_g_s = 5.2396

    star = create_star(
        R_s, T_s, log_g_s, Met_s,
        stellar_grid="phoenix",
    )

    planet_name = "Trappist-1e"
    R_p = 0.917985 * R_E
    M_p = 0.6356 * M_E
    T_eq = 255.0

    planet = create_planet(planet_name, R_p, mass=M_p, T_eq=T_eq)

    # ============================== Datos e instrumentos =======================
    wl_min, wl_max, R = 0.4, 6.0, 4000
    wl = wl_grid_constant_R(wl_min, wl_max, R)

    data_dir = Path("observations")
    observation = "pandexo_output_100transits_fspot0.26_ffac0.70.dat"
    datasets = [observation]
    instruments = ["JWST_NIRSpec_PRISM"]

    data = load_data(str(data_dir), datasets, instruments, wl)

    # ================================ Modelo ==================================
    model_name = "uncontam_100T_0.26spot-0.70fac"

    bulk_species  = ["N2"]
    param_species = ["H2O", "CH4", "CO2", "O3"]

    model = define_model(
        model_name,
        bulk_species,
        param_species,
        PT_profile="isotherm",
        cloud_model="cloud-free",
    )

    if _rank == 0:
        print("Free parameters:", model["param_names"], flush=True)

    # ================================ Priors ==================================
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
        "R_p_ref": [0.85 * R_p, 1.15 * R_p],
        "log_H2O": [-8.0, -1.0],
        "log_CH4": [-8.0, -1.0],
        "log_CO2": [-5.0, -1.0],
        "log_O3": [-8.0, -1.0],
    }

    priors = set_priors(planet, star, model, data, prior_types, prior_ranges)

    # ============================== Opacidades =================================
    opacity_treatment = "opacity_sampling"

    T_fine = np.arange(200.0, 400.0 + 10.0, 10.0)         # K
    log_P_fine = np.arange(-2.0, 2.0 + 0.2, 0.2)          # log10(P/bar)

    opac = read_opacities(model, wl, opacity_treatment, T_fine, log_P_fine)

    # ============================= Config. atmósfera ===========================
    P_min, P_max, N_layers = 1.0e-2, 2.0, 100
    P = np.logspace(np.log10(P_max), np.log10(P_min), N_layers)
    P_ref = 1.0  # bar

    # =============================== Retrieval ================================
    if _rank == 0:
        print(">> Iniciando retrieval con MultiNest...", flush=True)
    t1 = time.time()
    run_retrieval(
        planet, star, model, opac, data, priors, wl, P, P_ref, R=R,
        spectrum_type="transmission",
        sampling_algorithm="MultiNest",
        N_live=500,
        verbose=(_rank == 0),
        # Ejemplos de tuning (descomentar si los necesitas):
        # const_efficiency_mode=True,
        # sampling_efficiency=0.05,  # si const_efficiency_mode=True
        # evidence_tolerance=0.5,
    )
    dt = time.time() - t1
    if _rank == 0:
        print(f">> Retrieval terminado. Δt = {dt/60:.2f} min", flush=True)

    # ============================= Post-procesado ==============================
    # Asegura que MultiNest terminó de escribir todo
    comm.Barrier()

    result_spec_median = None
    if _rank == 0:
        # Leer espectros recuperados
        wl_ret, spec_low2, spec_low1, spec_median, spec_high1, spec_high2 = \
            read_retrieved_spectrum(planet_name, model_name)

        # Colección para el plot de espectros
        spectra_median = plot_collection(spec_median, wl_ret, collection=[])
        spectra_low1   = plot_collection(spec_low1,   wl_ret, collection=[])
        spectra_low2   = plot_collection(spec_low2,   wl_ret, collection=[])
        spectra_high1  = plot_collection(spec_high1,  wl_ret, collection=[])
        spectra_high2  = plot_collection(spec_high2,  wl_ret, collection=[])

        out_dir = Path("figures")
        out_dir.mkdir(parents=True, exist_ok=True)

        # Espectro recuperado
        fig_spec = plot_spectra_retrieved(
            spectra_median, spectra_low2, spectra_low1, spectra_high1, spectra_high2,
            planet_name, data, R_to_bin=100,
            data_labels=instruments,
            data_colour_list=["lime"],
        )
        fig_spec.savefig(out_dir / f"{planet_name}_{model_name}_retrieved_spectrum.png",
                         dpi=180, bbox_inches="tight")

        # Corner plot: puede devolver None, o (fig, axes), o fig
        corner_out = generate_cornerplot(planet, model)
        fig_corner = None
        if hasattr(corner_out, "savefig"):
            fig_corner = corner_out
        elif isinstance(corner_out, tuple) and len(corner_out) > 0 and hasattr(corner_out[0], "savefig"):
            fig_corner = corner_out[0]

        if fig_corner is not None:
            fig_corner.savefig(out_dir / f"{planet_name}_{model_name}_corner.png",
                               dpi=180, bbox_inches="tight")
        else:
            print("[corner] La función no devolvió Figure; es posible que ya haya guardado el plot internamente.")

        print(f">> Figuras guardadas en: {out_dir.resolve()}", flush=True)
        print(f">> Tiempo total: {(time.time()-t0)/60:.2f} min", flush=True)

        result_spec_median = spec_median

    return result_spec_median


if __name__ == "__main__":
    main()
