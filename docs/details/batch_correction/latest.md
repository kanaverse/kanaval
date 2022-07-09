# Batch correction (latest) {#details-batch_correction}

## Overview

We expect a `batch_correction` HDF5 group at the root of the file, containing the parameters and results of batch correction across samples.
The group itself contains the `parameters` and `results` subgroups.

A separate batch correction step was not used prior to version 2.0 of the format, so the `batch_correction` group may be absent in pre-v2.0 files.
In such cases, we can check for the presence of a `corrected` dataset in the [`pca`](../pca/v1_1.md) step.

**Definitions:**

- `num_samples`: the number of samples in the dataset.
  This is typically determined from the [`inputs`](../inputs/latest.md) step.
- `num_cells`: the number of cells remaining after quality control.
  This is typically determined from the [`cell_filtering`](../cell_filtering/latest.md) step.b
- `num_dims`: the expected number of dimensions in the original (uncorrected) embedding.
  This is typically determined from the [`combine_embeddings`](../combine_embeddings/latest.md) step.

## Parameters

`parameters` should contain:

- `num_neighbors`: a scalar integer specifying the number of neighbors to use for mutual nearest neighbors correction.
- `approximate`: a scalar integer to be treated as a boolean, indicating whether an approximate neighbor search was used.
- `method`: a scalar string specifying the correction method to use, either `"none"` or `"mnn"`.

## Results

If `method = "mnn"` and `num_samples > 1`, `results` should contain:

- `corrected`: a 2-dimensional float dataset containing the corrected PCs in a row-major layout.
  Each row corresponds to a cell and each column corresponds to a dimension.

Otherwise, correction is assumed to be a no-op and `results` may be empty.
Downstream steps should instead fetch coordinates from the `combined` dataset in the [`combine_embeddings`](../combine_embeddings/latest.md) step.

## Changelog

- Version 2.0: added the `batch_correction` step.
