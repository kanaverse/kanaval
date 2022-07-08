# Feature selection (latest) {#details-feature_selection}

## Overview

We expect a `feature_selection` HDF5 group at the root of the file, containing statistics for highly variable gene detection from the RNA data.
The group itself contains the `parameters` and `results` subgroups.
In multi-modal contexts, the results here are only relevant for the RNA modality.

## Parameters

`parameters` should contain:

- `span`: a scalar float specifying the span to use for the LOWESS smoother.

## Results

`results` should contain:

- `means`: a 1-dimensional float dataset of length equal to the number of genes,
  containing the mean log-expression of each gene.
- `vars`: a 1-dimensional float dataset of length equal to the number of genes,
  containing the variance in log-expression of each gene.
- `fitted`: a 1-dimensional float dataset of length equal to the number of genes,
  containing the fitted value of the trend for each gene.
- `resids`: a 1-dimensional float dataset of length equal to the number of genes,
  containing the residuals from the trend for each gene.

## Changelog

- Version 1.0: added the `feature_selection` group.
