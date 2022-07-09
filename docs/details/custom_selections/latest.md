# Custom selections (latest) {#details-custom_selections}

## Overview

We expect a `custom_selections` HDF5 group at the root of the file, containing information about the custom selections.
The group itself contains the `parameters` and `results` subgroups.

**Definitions:**

- `num_cells`: number of cells remaining after QC filtering.
  This is typically determined from the [`cell_filtering`](../cell_filtering/latest.md) step.
- `modalities`: an array of names of all available modalities.
  This is typically determined from the [`inputs`](../inputs/latest.md) step.
- `num_features`: an array of length equal to `modalities`, containing the number of features present in each modality.
  This is typically determined from the [`inputs`](../inputs/latest.md) step.

## Parameters

`parameters` should contain:

- `selections`: a group defining the custom selections.
  Each child is named after a user-created selection.
  Each child is an integer dataset of arbitrary length containing the indices of the selected cells.
  Note that indices refer to the dataset after QC filtering and should be less than `num_cells`.
  Indices in each dataset should also be unique and sorted.

## Results

`results` should contain `per_selection`, a group containing the marker results for each selection after a comparison to a group containing all other cells.
Each child of `per_selection` is named after its selection and is itself a group containing children named after `modalities`.
Each modality-specific child is yet another group containing the statistics for that modality:

- `means`: a float dataset of length equal to the number of features in that modality (as determined from `num_features`), containing the mean expression of each gene in the selection.
- `detected`: a float dataset of length equal to the number of features in that modality, containing the proportion of cells with detected expression of each gene in the selection.
- `lfc`: a float dataset of length equal to the number of features in that modality, containing the log-fold change in the selection compared to all other cells.
- `delta_detected`: same as `lfc`, but for the delta-detected, i.e., difference in the percentage of detected expression.
- `cohen`: same as `lfc`, but for Cohen's d.
- `auc`: same as `lfc`, but for the AUCs.

## Changelog

- Version 2.1: require that the indices be sorted.
- [Version 2.0](v2_0.md): support multiple modalities in the results.
- [Version 1.0](v1_0.md): added the `custom_selections` group.
