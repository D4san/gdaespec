# campaign_5obs

This folder stores the five-observation retrieval campaign for the Earth-like
TRAPPIST-1e analog.

Each `test_XX` folder is one independent 10-transit synthetic observation. For
each test, the campaign compares:

- `gdae`: retrieval on the G-DAE reconstructed spectrum.
- `contam`: retrieval on the raw contaminated spectrum with stellar
  contamination fitted inside POSEIDON.

The aggregate campaign tables are:

- `times.csv`: `id,branch,f_spot,f_fac,strategy,delta_time`
- `metrics.csv`: `id,branch,f_spot,f_fac,strategy,MSE,chi2_reduced`

The observations, POSEIDON products, samples, state files, figures, logs, and
plots are preserved as the retrieval record.

Typical commands are documented in the parent
[`README.md`](../README.md).
