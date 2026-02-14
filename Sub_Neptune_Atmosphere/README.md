# G-DAE Training on Sub-Neptunes (K2-18b Analog)

This directory contains the experimental workflow for generating data, training, and evaluating the G-DAE model on a synthetic library of Sub-Neptune atmospheres, specifically modeled after K2-18b.

## Workflow Structure

The experimental pipeline is designed to be executed in the following order:

1.  **`01_Spectra_Generation.ipynb`**: Generation of the synthetic spectral library for the K2-18b analog.
2.  **`02_Stellar_Contamination.ipynb`**: Calculation of the stellar contamination factor $\epsilon(\lambda)$ for various stellar activity configurations (spots and faculae).
3.  **`03_AE_Training.ipynb`**: Training of the G-DAE model using the generated spectral data and the computed contamination profiles.
4.  **`04_G-DAE_Evaluation.ipynb`**: Evaluation of the G-DAE model's performance on the test set, analyzing its ability to reconstruct clean spectra.
