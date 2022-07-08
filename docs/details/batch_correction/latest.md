# Batch correction (latest) {#details-batch_correction}

## Overview

We expect a `batch_correction` HDF5 group at the root of the file, containing the parameters and results of batch correction across samples.
The group itself contains the `parameters` and `results` subgroups.

A separate batch correction step was not used prior to version 2.0 of the format, so the `batch_correction` group may be absent in pre-v2.0 files.
In such cases, we can check the `corrected` dataset in the `pca/results` group, see `pca::validate()` for more details.

## Parameters

`parameters` should contain:

- `num_neighbors`: a scalar integer specifying the number of neighbors to use for mutual nearest neighbors correction.
- `approximate`: a scalar integer to be treated as a boolean, indicating whether an approximate neighbor search was used.
- `method`: a scalar string specifying the correction method to use, either `"none"` or `"mnn"`.

## Results

If `method = "mnn"` and there are multiple samples (i.e., `num_samples > 1` in `batch_correction::validate()`), `results` should contain:

- `corrected`: a 2-dimensional float dataset containing the corrected PCs in a row-major layout.
  Each row corresponds to a cell and each column corresponds to a dimension.
  (The expected number of each is determined by `num_cells` and `num_dims`, respectively, from `batch_correction::validate()`.)

Otherwise, correction is assumed to be a no-op and `results` may be empty.
Downstream steps should instead fetch coordinates from the `combined` dataset in `combine_embeddings::validate()`.

## Changelog

- Version 2.0: added the `batch_correction` step.
