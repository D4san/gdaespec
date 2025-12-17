
# G-DAESpec
<div align="center">
  <img src="G-DAE_logo.png" alt="G-DAESpec Logo" width="600"/>
</div>



[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Powered by MultiREx](https://img.shields.io/badge/Powered%20by-MultiREx-blueviolet)](https://github.com/D4san/MultiREx-public)
[![Powered by POSEIDON](https://img.shields.io/badge/Powered%20by-POSEIDON-007ec6)](https://github.com/MartianColonist/POSEIDON)
[![Powered by TauREx 3](https://img.shields.io/badge/Powered%20by-TauREx%203-red)](https://github.com/ucl-exoplanets/TauREx3_public)

## Overview

**G-DAESpec** implements a **G-DAE** (General Denoising AutoEncoder) designed to correct the effects of stellar contamination in exoplanet transmission spectra.

This work validates the G-DAE architecture using two distinct planetary cases:
1.  **Rocky Planets**: Using a **TRAPPIST-1e** analogue.
2.  **Sub-Neptunes**: Using **K2-18b** analogues.

The project leverages the **MultiREx** library for spectral handling and **POSEIDON** for atmospheric retrieval and validation.

## Features

- **Stellar Contamination Correction**: Removes noise introduced by stellar spots and faculae to recover the true planetary spectrum.
- **Deep Learning Architecture**: Utilizes a Denoising AutoEncoder (DAE) trained on simulated datasets.
- **Broad Applicability**: Tested on both terrestrial and sub-Neptune atmospheric regimes.
- **Validation**: Results are verified against standard Bayesian atmospheric retrieval methods (Nested Sampling via POSEIDON).

## Repository Structure

The repository is organized by planetary case studies:

### 1. Earth-like Atmospheres (`Earth_like_Atmosphere/`)
Focuses on the TRAPPIST-1e analogue.
- **Data Generation**: `02_spec_data.ipynb` (Spectra), `02_stellar_contamination_epsilon.ipynb` (Contamination factors).
- **Models**: `Models/` directory containing the trained AutoEncoder (`AE.keras`) and training notebooks (`03_AE.ipynb`).
- **Validation**: `Retreival Tests/` containing POSEIDON retrieval outputs and comparison plots.
- **Uncertainties**: `AE pandexo_incertidumbres_final.ipynb` for instrument noise simulation (JWST).

### 2. Sub-Neptune Atmospheres (`Sub_Neptune_Atmosphere/`)
Focuses on the K2-18b analogue.
- **Data**: `specs/` folder with generated spectra (`.joblib`, `.parquet`).
- **Models**: `AE.keras` and training scripts (`AE_l2.ipynb`).
- **Analysis**: `RetrievalAnal/` for detailed retrieval performance assessment.

## Requirements

The following key libraries are used in this project:

- **[MultiREx](https://github.com/D4san/MultiREx-public)**: For handling exoplanet spectra.
- **[POSEIDON](https://github.com/MartianColonist/POSEIDON)**: For atmospheric retrieval and spectra generation.
- **[TauREx 3](https://github.com/ucl-exoplanets/TauREx3_public)**: Bayesian retrieval framework.
- **TensorFlow / Keras**: For building and training the G-DAE models.
- **Pandeia / PandExo**: For JWST instrument noise simulation.
- **Scikit-learn, NumPy, Pandas, Matplotlib**: For data processing and visualization.

### Opacity Data & Chemical Species

This project utilizes opacity data and chemical species line lists from the following sources:
- **[ExoMol](https://www.exomol.com/)**: Molecular line lists for exoplanet and other hot atmospheres.
- **[ExoTransmit](https://github.com/elizakempton/Exo_Transmit)**: Transmission spectra calculation.

Please ensure to cite these works appropriately when using the data.

## Citation

If you use this code in your research, please cite the associated paper (link/DOI to be added).

## License

[MIT License](LICENSE)
