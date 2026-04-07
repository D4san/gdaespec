# Experiment 0: G-DAE Performance on an Earth-like Atmosphere (TRAPPIST-1e Analog)

This directory contains the workflow used to train and evaluate the G-DAE on Earth-like transmission spectra for a TRAPPIST-1e analog.

## Data Source & Reference

The base spectral dataset used in this experiment is derived from [Duque-Castaño et al. (2025)](https://academic.oup.com/mnras/article/539/2/1528/8109635). The original source repository is:

- [MultiREx-public/examples/papers/DZF-MLBiosignatureClassification](https://github.com/D4san/MultiREx-public/tree/main/examples/papers/DZF-MLBiosignatureClassification)

The SPHINX stellar grids used to generate alternative contamination curves for TRAPPIST-1 were downloaded separately by the user and are documented inside [stellar_contamination](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination).

## Core Workflow

The main experiment is organized in the following order:

1. `01_G-DAE.ipynb`
   Train the G-DAE on the Earth-like spectral dataset. This notebook already ingests the PHOENIX-based and SPHINX-based `epsilon` families when they are present in `stellar_contamination/`.
2. `02_G-DAE_Analysis.ipynb`
   Evaluate the trained G-DAE and compute reconstruction metrics. The main analysis cells remain tied to the original PHOENIX contamination workflow.
3. `Retrieval Tests/`
   Run atmospheric retrieval experiments on synthetic observations derived from the spectra.

## SPHINX Workflow

The SPHINX branch is intentionally kept separate from the main analysis flow:

1. Generate the SPHINX `epsilon(lambda)` files in [stellar_contamination](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/stellar_contamination).
2. Open `02_G-DAE_Analysis.ipynb` and run only the final export section named `Retrieval Observations for PHOENIX and SPHINX Injection`.
3. That section writes:
   - PHOENIX-based observations to [Retrieval Tests/observations](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/observations)
   - SPHINX-injected observations to [Retrieval Tests/observations_sphinx](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/observations_sphinx)
   - the corresponding `*_recon.dat` files for the G-DAE strategy
4. Run the retrieval scripts in [Retrieval Tests/sphinx_injection](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/sphinx_injection) to perform retrievals on SPHINX-injected observations while keeping the retrieval star model in PHOENIX.

## Practical Note

`02_G-DAE_Analysis.ipynb` now serves two different purposes:

- the main notebook body evaluates the original G-DAE experiment
- the final export section prepares retrieval-ready observations for both PHOENIX and SPHINX injection cases

This separation is deliberate so the analysis notebook does not mix SPHINX into the uncertainty tests unless that is explicitly desired later.
