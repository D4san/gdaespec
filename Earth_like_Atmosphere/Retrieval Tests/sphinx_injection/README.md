# SPHINX Injection Retrievals

This subdirectory contains the retrieval scripts for the SPHINX injection experiment.

## Purpose

These runs are designed to test a model-mismatch scenario:

- the synthetic observations are generated with SPHINX stellar contamination curves
- the retrieval itself still uses PHOENIX stellar models

This allows a controlled comparison against the standard PHOENIX-based workflow.

## Scripts

- `Retrieval_no-contam_mpi.py`
  Chemistry-only retrieval of the SPHINX-injected observation.
- `Retrieval_contam_mpi.py`
  Joint atmospheric and stellar-contamination retrieval of the SPHINX-injected observation.
- `Retrieval_G-DAE_mpi.py`
  Retrieval of the G-DAE reconstructed SPHINX-injected observation.

## Analysis Notebooks

- `Retrieval_analysis no-contam_sphinx.ipynb`
  Notebook for analyzing the chemistry-only retrieval on SPHINX-injected observations.
- `Retrieval_analysis_contam_sphinx.ipynb`
  Notebook for analyzing the contamination-aware retrieval on SPHINX-injected observations.
- `Retrieval_analysis G-DAE sphinx.ipynb`
  Notebook for the G-DAE branch. It includes the preprocessing section that reconstructs the SPHINX-injected observation files into their `*_recon.dat` counterparts before the retrieval analysis.

## Inputs

These scripts expect the observation files to already exist in:

- [../observations_sphinx](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/observations_sphinx)

Those files are generated from the final export section of:

- [../../02_G-DAE_Analysis.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/02_G-DAE_Analysis.ipynb)

The G-DAE analysis notebook in this folder can also regenerate the reconstructed `*_recon.dat` files from `../observations_sphinx/` if you want the preprocessing step documented inside the retrieval-analysis branch itself.

## Default Case

The scripts are currently configured for:

- `n_transits = 10`
- `f_spot_case = 0.26`
- `f_fac_case = 0.70`

You can change the case in two equivalent ways:

- edit the module defaults inside the selected script
- export `RETRIEVAL_N_TRANSITS`, `RETRIEVAL_F_SPOT_CASE`, and `RETRIEVAL_F_FAC_CASE` before launching the MPI job

This keeps the three SPHINX retrieval strategies synchronized on the same case
without requiring separate file edits for each run.

## Runtime Log

Completed SPHINX runs append their total wall-clock time to:

- [Times](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Times)

The format mirrors the baseline `Retrieval Tests/Times` file, but it is stored
inside `sphinx_injection/` so the mismatch experiment stays self-contained.

## Running

Example:

```bash
mpirun -n <cores> python -u "Earth_like_Atmosphere/Retrieval Tests/sphinx_injection/Retrieval_no-contam_mpi.py"
```

Each script changes its working directory to this folder before running, so the generated retrieval outputs stay isolated from the baseline retrieval branch.
