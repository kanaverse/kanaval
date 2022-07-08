# Combine embeddings (latest) {#details-combine_embeddings}

## Overview

We expect a `combine_embedding` HDF5 group at the root of the file, describing how embeddings are combined across multiple modalities.
The group itself contains the `parameters` and `results` subgroups.

Only a single embedding was generated prior to version 2.0 of the format, so the `combine_embeddings` group may be absent in pre-v2.0 files.
No concept of multi-modality exists in earlier versions so downstream steps should use the PCs directly from `pca::validate()`.

## Parameters

`parameters` should contain:

- `approximate`: an integer scalar to be interpreted as a boolean,
  indicating whether an approximate neighbor search was used to compute the per-embedding scaling factors.
- `weights`: a group containing the (scaling) weights to apply to each modality.
  If empty, weights are implicitly assumed to be equal to unity for all modalities.
  Otherwise, the group should contain a float scalar dataset named after each modality (see `modalities` in `combine_embedding::validate()`), 
  containing the weight to be applied to each modality.

## Results

If multiple modalities have available embeddings (i.e., `modalities.size() > 1` in `combine_embedding::validate()`), `results` should contain:

- `combined`: a 2-dimensional float dataset containing the combined embeddings in a row-major layout.
  Each row corresponds to a cell and each column corresponds to a dimension.
  The expected number of dimensions should be the sum of the number of dimensions for each embedding-relevant modality, i.e., `total_dims` in `combine_embedding::validate()`.
  (The expected number of cells is determined from `num_cells` in `combine_embedding::validate()`).

Otherwise, `combined` may be missing, in which case it is assumed that only a single embedding exists in the analysis and no combining is necessary.
Downstream steps should instead use the coordinates from the PCA group of the available modality, see `pca::validate()` or `adt_pca::validate()`.

## Changelog

- Version 2.0: added the `combine_embedding` group.
