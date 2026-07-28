# Corrected G-DAE Retrieval Campaign Archive

This manifest records the preservation copy created before consolidating the
G-DAE retrieval branches.

## Campaign

- Campaign ID: `gdae_campaign_5obs_corrected_20260724_v1`
- Preservation date: `2026-07-28`
- Completed retrievals: `35/35`
- Archived source: `/home/dasan/gdae_campaign_5obs_corrected_20260724_v1`
- Source files verified: `707/707`
- Result files: `35`
- Reconstructed observations: `35`
- MultiNest raw files: `280`
- Log files: `38`
- PNG/PDF products: `119`

## Preservation copy

- Archive:
  `/home/dasan/gdae_retrieval_archives/gdae_campaign_5obs_corrected_20260724_v1_preserved_20260728.tar.gz`
- Archive SHA-256:
  `90da17d000dfbe7c739ca77e4bf514dbd74d3a5ccbb510ccf56b2337cea53e63`
- Per-file checksum manifest:
  `/home/dasan/gdae_retrieval_archives/gdae_campaign_5obs_corrected_20260724_v1_preserved_20260728_FILES.sha256`
- Archive checksum file:
  `/home/dasan/gdae_retrieval_archives/gdae_campaign_5obs_corrected_20260724_v1_preserved_20260728.tar.gz.sha256`

The compressed archive passed `gzip -t`, its SHA-256 checksum was verified, and
all 707 source files matched the per-file checksum manifest.

## Reconstruction provenance

- Earth-like G-DAE checkpoint SHA-256:
  `2eef4becfa00b7eafcb0cc0ba012b9654d5b613bc7c64ab918c3a170af256d0f`
- Input transform: `log(depth / ref_flat)`
- Physical-space inverse: `exp(pred_log) * ref_flat`
- Model wavelength order: strictly descending
- Central reconstruction: physical-space mean of the MC-dropout predictions
- Epistemic uncertainty: physical-space standard deviation of the MC-dropout
  predictions
- Total diagonal uncertainty:
  `sqrt(sigma_epistemic**2 + (0.5 * sigma_instrumental)**2)`

The archive is intentionally stored outside the Git working tree. It must not
be deleted until it has been copied to a second durable storage location.
