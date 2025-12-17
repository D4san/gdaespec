# MultiREx: Machine-Assisted Biosignature Classification in Exoplanet Atmospheres

## Abstract

This repository hosts the computational framework and datasets associated with the research on classifying potential biosignatures and atmospheric compositions in exoplanets. The study leverages deep learning techniques, specifically AutoEncoders (AE), to analyze low signal-to-noise ratio (SNR) transmission spectra of both Earth-like and Sub-Neptune planets. By reducing the dimensionality of spectral data, we aim to efficiently identify atmospheric components (e.g., $H_2O$, $CH_4$, $O_3$) amidst stellar contamination and instrumental noise.

## Description

The project focuses on:
1.  **Spectra Generation**: Simulating transmission spectra for Earth-like (e.g., TRAPPIST-1e) and Sub-Neptune (e.g., K2-18b) atmospheres.
2.  **Stellar Contamination**: Modeling the impact of stellar spots and faculae on the observed spectra.
3.  **Dimensionality Reduction**: Training AutoEncoders to compress spectral data while retaining essential features for classification.
4.  **Retrieval Analysis**: Validating the machine learning approach against traditional atmospheric retrieval methods (using POSEIDON).

## Repository Structure

The project is divided into two main domains based on the planet type:

### 1. Earth-like Atmospheres (`Earth_like_Atmosphere/`)
Focused on terrestrial planets with potential biosignatures.

- **Data Preparation**:
  - `02_spec_data.ipynb`: Generates transmission spectra dataset.
  - `02_stellar_contamination_epsilon.ipynb`: Computes contamination factors ($\epsilon$).
- **Modeling**:
  - `Models/`: Contains AutoEncoder models (`AE.keras`) and training notebooks.
- **Validation**:
  - `Retreival Tests/`: POSEIDON retrieval results (MultiNest) and corner plots.

### 2. Sub-Neptune Atmospheres (`Sub_Neptune_Atmosphere/`)
Focused on gas-rich sub-Neptunes like K2-18b.

- **Spectra & Data**:
  - `specs/`: Contains spectral datasets (`K2-18b_data.joblib`, parquet files) and generation notebooks (`Spectra_Generatio.ipynb`).
  - `stellar_contamination.ipynb`: Stellar contamination analysis for this case.
- **Modeling**:
  - `AE.keras`, `AE_l2.ipynb`: AutoEncoder models and training scripts.
- **Analysis**:
  - `RetrievalAnal/`: Comprehensive retrieval analysis including `Retrieval_analysis.ipynb` and `POSEIDON_output/`.
  - `AE pandexo_incertidumbres.ipynb`: Uncertainty analysis with Pandexo.

### Shared Resources
- **`AE pandexo_incertidumbres_final.ipynb`** (in Earth-like folder): General uncertainty analysis framework.

## Requirements

To replicate the experiments, the following Python packages are required:

- `tensorflow` / `keras` (for AutoEncoders)
- `multirex` (for spectra generation)
- `poseidon` (for atmospheric retrieval)
- `pandexo` (for instrument simulation)
- `scikit-learn`
- `plotly`
- `numpy`, `pandas`, `matplotlib`

## Citation

If you use this code or data in your research, please refer to the associated publication (insert link/DOI here).

## License

[MIT License](LICENSE)
