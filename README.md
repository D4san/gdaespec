# gdaeaspec

## Overview

This repository contains tools and models for the analysis of astronomical spectra, with a focus on Earth-like atmospheres and stellar contamination. The project utilizes Machine Learning techniques, including Autoencoders (AE) and Random Forests (RF), to process and analyze spectral data.

## Features

- **Spectral Analysis**: Tools for processing and analyzing spectral data (e.g., `01_pandexo_spec_analysis.ipynb`, `02_spec_data.ipynb`).
- **Stellar Contamination**: Analysis of stellar contamination (e.g., `02_stellar_contamination_epsilon.ipynb`).
- **Machine Learning Models**:
  - **Autoencoders (AE)**: Trained models and notebooks for dimensionality reduction and feature extraction (e.g., `03_AE.ipynb`, `AE.keras`).
  - **Random Forests (RF)**: Models for regression or classification tasks involving specific molecules like CH4, H2O, and O3 (e.g., `04_CH4_RF.ipynb`, `04_H2O_RF.ipynb`, `04_O3_RF.ipynb`).
- **Data**: Includes synthetic or observational spectral data (e.g., Phoenix models).

## Structure

The main work is located in the `Earth_like_Atmosphere/` directory, which contains:
- `Models/`: Saved Keras models.
- `Phoenix/`: Stellar spectra files (BT-Settl models).
- Notebooks for data analysis, model training, and evaluation.

## Author

**David Santiago Duque Castaño**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
