# Retrieval Tests

This directory contains the POSEIDON retrieval workflow used to compare two
ways of handling stellar contamination in synthetic TRAPPIST-1e transmission
observations.

The active workflow is the five-observation campaign in `campaign_5obs/`.

## Scientific Comparison

Each campaign test uses a synthetic 10-transit observation generated from the
same clean atmospheric spectrum and one contamination case.

The noisy spectra are generated with PandExo in `campaign_observations.py`.
That script starts from `pandexo_spec.txt`, applies the selected
`epsilon(lambda)` contamination curve from `../stellar_contamination/`, and
writes the noisy observation files under `campaign_5obs/test_XX/<branch>/observations/`.
It also writes the matching G-DAE reconstructions as `*_recon.dat`.
For the G-DAE branch, the third column is the physical-space mean of the
MC-dropout ensemble. The fourth column follows the analysis-notebook contract:
it combines the fixed-input epistemic spread with half the PandExo instrumental
uncertainty in quadrature.
The two retrieval strategies are:

1. `gdae`
   The observed spectrum is preprocessed with the trained G-DAE. POSEIDON then
   retrieves the atmospheric parameters from the reconstructed `*_recon.dat`
   spectrum.
2. `contam`
   POSEIDON retrieves atmospheric parameters and stellar-contamination
   parameters jointly from the raw contaminated observation.

The campaign includes two contamination branches:

- `phoenix`: contamination curves from the PHOENIX-based stellar model grid.
- `sphinx`: observations injected with SPHINX contamination curves while the
  retrieval still uses the PHOENIX stellar model prescription.

## Files

| Path | Purpose |
| --- | --- |
| `campaign_common.py` | Shared campaign configuration, paths, case definitions, and CSV helpers. |
| `campaign_setup.py` | Creates the campaign directory tree and empty CSV headers. |
| `campaign_observations.py` | Generates PandExo observations and G-DAE reconstructions. |
| `campaign_retrieval_mpi.py` | Runs one POSEIDON retrieval case with MPI. |
| `campaign_metrics.py` | Computes MSE and reduced chi-square for completed retrievals. |
| `campaign_plot_aggregates.py` | Plots aggregate timing and metric summaries from campaign CSV files. |
| `campaign_plot_parameters.py` | Plots retrieved atmospheric parameters from POSEIDON result files. |
| `campaign_plot_observations.py` | Plots raw observations and G-DAE reconstructions across tests. |
| `campaign_run_gdae_queue.py` | Runs the missing `gdae` campaign jobs. |
| `campaign_run_contam_queue.py` | Runs the missing `contam` campaign jobs. |
| `pandexo_spec.txt` | Clean spectrum used as the PandExo input baseline. |
| `campaign_5obs/` | Campaign CSV summaries and generated per-test products. |

## Campaign Layout

`campaign_5obs/` contains:

- `metrics.csv`: aggregate metric table with
  `id,branch,f_spot,f_fac,strategy,MSE,chi2_reduced`.
- `times.csv`: aggregate run-time table with
  `id,branch,f_spot,f_fac,strategy,delta_time`.
- `test_01/` to `test_05/`: generated observations, retrieval outputs, and
  figures for each synthetic observation.
- `plots/`: aggregate figures generated from the campaign outputs.

The observations, samples, MultiNest state, POSEIDON products, figures, and
logs are preserved together with the aggregate tables.

## Inputs

The campaign expects these project files:

- `pandexo_spec.txt`
- `../Models/G-DAE.keras`
- `../stellar_contamination/`

Regenerating noisy observations requires PandExo/Pandeia. Running retrievals
requires a working POSEIDON installation with MultiNest, MPI, `mpirun`, and the
corresponding opacity and stellar-model data available locally.

## Typical Workflow

The included retrieval record can be inspected without rerunning POSEIDON:

```bash
python campaign_plot_aggregates.py
python campaign_plot_parameters.py
```

Regenerating the full experiment is an advanced workflow. First create the
layout and observations, then launch the two queues:

```bash
python campaign_setup.py
python campaign_observations.py --test-id test_01 --branch all
python campaign_run_gdae_queue.py --nproc 12 --include-test01 --keep-going
python campaign_run_contam_queue.py --nproc 12 --include-test01 --keep-going
```

Repeat the observation command for `test_02` through `test_05`. The queue
commands invoke `mpirun` and skip cases whose result file already exists.
Use the direct `campaign_retrieval_mpi.py` command only to run one selected
case.
