# Stellar Contamination Assets

This directory stores the stellar contamination curves used in the Earth-like TRAPPIST-1 workflow.

## Contents

- `TRAPPIST-1_contam_fspot*_ffac*.txt`
  Original contamination curves used in the PHOENIX-based workflow.
- `sphinx_TRAPPIST-1_contam_fspot*_ffac*.txt`
  Alternative contamination curves generated with SPHINX stellar spectra for the same spot and facula coverage fractions.
- `generate_sphinx_trappist1_contamination.ipynb`
  Notebook used to generate the SPHINX curves.
- `bitacora_generacion_sphinx_trappist1.md`
  Step-by-step generation log for the SPHINX contamination files.
- `sphinx_data/`
  Local SPHINX stellar grids used as input to generate the SPHINX `epsilon(lambda)` curves.

## Role in the Workflow

This folder feeds two different stages:

1. `01_G-DAE.ipynb`
   The training dataset can include both the original PHOENIX contamination curves and the SPHINX ones.
2. `02_G-DAE_Analysis.ipynb`
   The final export section reads these files to generate retrieval-ready observations for:
   - the standard PHOENIX contamination branch
   - the SPHINX injection branch

## SPHINX Branch

The SPHINX contamination files are meant to support the following experiment:

- inject synthetic observations using SPHINX `epsilon(lambda)`
- run atmospheric retrievals with PHOENIX stellar models
- compare the behavior of the retrieval under model mismatch

The retrieval-ready observation files produced from these curves are written outside this folder, in:

- [Retrieval Tests/observations_sphinx](/C:/Proyectos/Astro/gdaespec/Earth_like_Atmosphere/Retrieval%20Tests/observations_sphinx)

## References

- Base Earth-like spectral dataset: [MultiREx-public - DZF-MLBiosignatureClassification](https://github.com/D4san/MultiREx-public/tree/main/examples/papers/DZF-MLBiosignatureClassification)
- Associated paper: [Monthly Notices of the Royal Astronomical Society (MNRAS), 2025](https://academic.oup.com/mnras/article/539/2/1528/8109635)
- SPHINX grid archive used as user-provided input: [Zenodo record 7416042](https://zenodo.org/records/7416042)
