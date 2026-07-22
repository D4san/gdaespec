# Examples

This folder contains short notebooks that show how to apply the trained models
without walking through the full experiment workflow.

| Notebook | Purpose |
| --- | --- |
| `G-DAE-Example.ipynb` | Earth-like G-DAE example using one campaign observation and its clean reference. |

The input files are `data/observation.dat` and `data/reference.dat`. The clean
reference is used only to evaluate the reconstruction and prepare the plots.

## Model contract

The published models use `log(depth / ref_flat)` as their input and output
representation. Predictions are returned to physical transit-depth units with
`exp(prediction) * ref_flat`. The scalar reference depths, channel count,
wavelength order, and file hashes are recorded in `data/model-contracts.json`.
