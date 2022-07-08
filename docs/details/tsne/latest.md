# t-SNE (latest) {#details-tsne}

## Overview

We expect an `tsne` HDF5 group at the root of the file, containing the t-SNE parameters and results.
The group itself contains the `parameters` and `results` subgroups.

## Parameters

`parameters` should contain:

- `perplexity`: a scalar float specifying the t-SNE perplexity.
- `iterations`: a scalar integer specifying the t-SNE iterations.
- `animate`: a scalar integer to be interpreted as a boolean, indicating whether an animation should be performed.

## Results

`results` should contain:

- `x`: a float dataset of length equal to the number of cells (after QC filtering), containing the x-coordinates for each cell.
- `y`: a float dataset of length equal to the number of cells (after QC filtering), containing the y-coordinates for each cell.

## Changelog

- Version 1.0: added the `tsne` group.
