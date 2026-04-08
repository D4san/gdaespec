# Retrieval Tests

This directory contains the POSEIDON retrieval workflows used to test how stellar contamination affects the recovery of atmospheric parameters.

## Retrieval Strategies

The standard comparison uses three strategies:

1. `No-contam`
   Retrieval of the contaminated spectrum with no stellar contamination model.
2. `Contam`
   Retrieval of the contaminated spectrum while explicitly fitting stellar contamination.
3. `G-DAE`
   Retrieval of the G-DAE reconstructed spectrum with a chemistry-only atmospheric model.

## Dependencies

- `POSEIDON`
- `MPI` such as MPICH or OpenMPI

## Directory Layout

- `observations/`
  Retrieval-ready PHOENIX-based observations and their reconstructed `*_recon.dat` companions.
- `observations_sphinx/`
  Retrieval-ready observations injected with SPHINX `epsilon(lambda)`.
- `sphinx_injection/`
  MPI retrieval scripts and companion analysis notebooks for the SPHINX injection experiment.
- `Retrieval_*_mpi.py`
  The original retrieval scripts for the baseline PHOENIX workflow.
- `Retrieval_analysis*.ipynb`
  Notebooks used to inspect retrieval outputs and compare strategies.

## Baseline PHOENIX Workflow

1. Generate or update the PHOENIX-based observation files in `observations/`.
2. Run one of the scripts:
   - `Retrieval_no-contam_mpi.py`
   - `Retrieval_contam_mpi.py`
   - `Retrieval_G-DAE_mpi.py`
3. Inspect the output with the corresponding retrieval analysis notebook.

### Posterior-propagated diagnostics

If you want to evaluate `MSE` and reduced `chi^2` using posterior sampled
spectra instead of only the retrieved median spectrum, run:

```bash
python "Earth_like_Atmosphere/Retrieval Tests/posterior_diagnostics.py" --mode uncontam --branch phoenix
```

This helper will:

1. locate the raw MultiNest posterior for the requested retrieval case
2. regenerate a configurable number of spectra with `retrieved_samples(...)`
3. save those propagated spectra to `posterior_diagnostics/<model>_posterior_samples.npz`
4. update `chi2_log_posterior.csv` with:
   - the legacy median-spectrum metrics (`MSE`, `chi2`, `chi2_reduced`)
   - new `posterior_*` summary columns for the propagated metric distribution

Useful flags:

- `--mode {uncontam,contam,recon}`
- `--branch {phoenix,sphinx}`
- `--n-transits`
- `--f-spot`
- `--f-fac`
- `--posterior-samples`

The same script is importable from notebooks via:

```python
from posterior_diagnostics import build_case_config, run_case
```

Typical MPI usage:

```bash
mpirun -n <cores> python -u Retrieval_no-contam_mpi.py
```

The three baseline MPI scripts now share the same case parameters:

- `RETRIEVAL_N_TRANSITS`
- `RETRIEVAL_F_SPOT_CASE`
- `RETRIEVAL_F_FAC_CASE`

If those environment variables are not defined, the scripts fall back to their
module defaults.

## SPHINX Injection Workflow

This branch is designed for the mismatch experiment:

- the synthetic observation is generated with SPHINX `epsilon(lambda)`
- the retrieval still uses PHOENIX stellar models

### Step 1: Generate SPHINX observations

Open [02_G-DAE_Analysis.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/02_G-DAE_Analysis.ipynb) and run the final section:

- `Retrieval Observations for PHOENIX and SPHINX Injection`

That section will:

1. read the contamination curves from [stellar_contamination](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination)
2. build PandExo input spectra for the selected cases
3. write the SPHINX-injected observations to `observations_sphinx/`
4. generate the matching `*_recon.dat` files for the G-DAE strategy

### Step 2: Run the SPHINX retrievals

Use the scripts in [sphinx_injection](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection):

- [Retrieval_no-contam_mpi.py](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Retrieval_no-contam_mpi.py)
- [Retrieval_contam_mpi.py](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Retrieval_contam_mpi.py)
- [Retrieval_G-DAE_mpi.py](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Retrieval_G-DAE_mpi.py)

Run them from a cluster or MPI-enabled environment, for example:

```bash
mpirun -n <cores> python -u "Earth_like_Atmosphere/Retrieval Tests/sphinx_injection/Retrieval_no-contam_mpi.py"
```

Each script:

- reads its input from `observations_sphinx/`
- keeps the stellar model in PHOENIX
- writes outputs relative to `sphinx_injection/` so the SPHINX branch stays separated from the baseline retrieval outputs
- appends the total wall-clock time of each completed run to `sphinx_injection/Times`

The posterior diagnostics helper also works on this branch:

```bash
python "Earth_like_Atmosphere/Retrieval Tests/posterior_diagnostics.py" --mode recon --branch sphinx
```

### Step 3: Analyze the SPHINX retrievals

The companion notebooks live in [sphinx_injection](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection):

- [Retrieval_analysis no-contam_sphinx.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Retrieval_analysis%20no-contam_sphinx.ipynb)
- [Retrieval_analysis_contam_sphinx.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Retrieval_analysis_contam_sphinx.ipynb)
- [Retrieval_analysis G-DAE sphinx.ipynb](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection/Retrieval_analysis%20G-DAE%20sphinx.ipynb)

In particular, `Retrieval_analysis G-DAE sphinx.ipynb` preserves the preprocessing block needed to:

- read the SPHINX-injected observation files
- reconstruct the `*_recon.dat` products from `../observations_sphinx/`
- evaluate the retrieval outputs against the clean reference spectrum

## Retrieval Cases

The currently documented coverage pairs for retrieval are:

- `(0.00, 0.00)`
- `(0.01, 0.08)`
- `(0.08, 0.54)`
- `(0.26, 0.70)`

You can switch between cases either by editing the defaults inside the script or
by exporting:

- `RETRIEVAL_N_TRANSITS`
- `RETRIEVAL_F_SPOT_CASE`
- `RETRIEVAL_F_FAC_CASE`

before launching the MPI run.
