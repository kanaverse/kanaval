# UMAP (latest) {#details-umap}

## Overview

We expect an `umap` HDF5 group at the root of the file, containing the UMAP parameters and results.
The `umap` group itself contains the `parameters` and `results` subgroups.

**Definitions:**

- `num_cells`: number of cells remaining after QC filtering.
  This is usually determined from the [`cell_filtering`](../cell_filtering/latest.md) step.

## Parameters

`parameters` should contain:

- `num_epochs`: a scalar integer containing the number of epochs to perform.
- `num_neighbors`: a scalar integer containing the number of nearest neighbors to use when constructing the sets.
- `min_dist`: a scalar float specifying the minimum distance between points.
- `animate`: a scalar integer to be interpreted as a boolean, indicating whether an animation should be performed.

## Results

`results` should contain:

- `x`: a float dataset of length equal to `num_cells`, containing the x-coordinates for each cell.
- `y`: a float dataset of length equal to `num_cells`, containing the y-coordinates for each cell.

## Changelog

- Version 1.0: added the `umap` group.
