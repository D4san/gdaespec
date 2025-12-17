# MultiREx: Machine-Assisted Biosignature Classification in Earth-like Exoplanets

## Abstract

This repository hosts the computational framework and datasets associated with the research on classifying potential biosignatures in Earth-like exoplanets. The study leverages deep learning techniques, specifically AutoEncoders (AE), to analyze low signal-to-noise ratio (SNR) transmission spectra. By reducing the dimensionality of spectral data, we aim to efficiently identify atmospheric components such as $H_2O$, $CH_4$, and $O_3$ amidst stellar contamination and instrumental noise.

## Description

The project focuses on:
1.  **Spectra Generation**: Simulating transmission spectra for Earth-like atmospheres with varying concentrations of biosignatures.
2.  **Stellar Contamination**: Modeling the impact of stellar spots and faculae on the observed spectra.
3.  **Dimensionality Reduction**: Training AutoEncoders to compress spectral data while retaining essential features for classification.
4.  **Retrieval Analysis**: Validating the machine learning approach against traditional atmospheric retrieval methods (using POSEIDON).

## Repository Structure

The core resources are located in the `Earth_like_Atmosphere/` directory:

### Data Preparation
- **`02_spec_data.ipynb`**: Notebook for generating the transmission spectra dataset used for training and testing.
- **`02_stellar_contamination_epsilon.ipynb`**: Computes the contamination factors ($\epsilon$) arising from stellar heterogeneity (spots and faculae).

### Modeling (Deep Learning)
- **`Models/`**: Directory containing:
  - `03_AE.ipynb`: The training pipeline for the AutoEncoder architecture.
  - Saved Keras models (`AE.keras`, `AE_l2.keras`).

### Uncertainty & Simulation
- **`AE pandexo_incertidumbres_final.ipynb`**: Analysis incorporating Pandexo simulations to account for observational uncertainties and instrument noise (e.g., JWST NIRSpec/MIRI).

### Validation
- **`Retreival Tests/`**: Contains outputs from atmospheric retrieval experiments, including nested sampling results (MultiNest) and corner plots comparing retrieved parameters against ground truths.

## Requirements

To replicate the experiments, the following Python packages are required:

- `tensorflow` / `keras` (for AutoEncoders)
- `pandexo` (for instrument simulation)
- `poseidon` (for atmospheric retrieval)
- `scikit-learn`
- `plotly`
- `numpy`, `pandas`, `matplotlib`

## Citation

If you use this code or data in your research, please refer to the associated publication (insert link/DOI here).

## License

[MIT License](LICENSE)
