# Experiment 0: G-DAE Performance on an Earth-like Atmosphere (TRAPPIST-1e Analog)

This directory contains the experimental workflow exploring the capacity of the G-DAE to mitigate stellar contamination in transmission spectra. The study focuses on an Earth-like atmosphere using a TRAPPIST-1e planetary analog.

## Data Source & Reference

The data used in this experiment is derived from the work of [Duque-Castaño et al. (2025)](https://academic.oup.com/mnras/article/539/2/1528/8109635). The original datasets can be found in the following repository:  
[MultiREx-public/examples/papers/DZF-MLBiosignatureClassification](https://github.com/D4san/MultiREx-public/tree/main/examples/papers/DZF-MLBiosignatureClassification)

## Workflow Structure

The experimental pipeline is designed to be executed in the following order:

1.  **`01_G-DAE.ipynb`**: Training regarding the G-DAE model.
2.  **`02_G-DAE_Analysis.ipynb`**: Analysis of the G-DAE output and performance metrics.
3.  **`Retrieval Tests/`**: Atmospheric retrieval experiments performed on the resulting data to validate the decontamination strategy.
