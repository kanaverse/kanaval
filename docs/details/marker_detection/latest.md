# Marker detection (latest) {#details-marker_detection}

## Overview

We expect a `marker_detection` HDF5 group at the root of the file, containing marker statistics for each cluster.
The group itself contains the `parameters` and `results` subgroups.

**Definitions:**

- `num_clusters`: number of clusters in the analysis.
  This is typically determined from [`choose_clustering`](../choose_clustering/latest.md).
- `modalities`: an array of names of all available modalities.
  This is typically determined from the [`inputs`](../inputs/latest.md) step.
- `num_features`: an array of length equal to `modalities`, containing the number of features present in each modality.
  This is typically determined from the [`inputs`](../inputs/latest.md) step.

## Parameters

`parameters` should be empty.

## Results

`results` should contain `per_cluster`, a group containing the marker results for each cluster.
Each child of `per_cluster` is named after a cluster index from 0 to `num_clusters - 1`, and is itself a group containing children named according to `modalities`.
Each modality-specific child is yet another group containing the statistics for that modality:

- `means`: a float dataset of length equal to the number of features for this modality (as determined from `num_features`), containing the mean expression of each feature in the current cluster.
- `detected`: a float dataset of length equal to the number of features, containing the proportion of cells with detected expression of each feature in the current cluster.
- `lfc`: an group containing statistics for the log-fold changes from all pairwise comparisons involving the current cluster.
  This contains:
  - `min`: a float dataset of length equal to the number of features, containing the minimum log-fold change across all pairwise comparisons for each feature.
  - `mean`: a float dataset of length equal to the number of features, containing the mean log-fold change across all pairwise comparisons for each feature.
  - `min_rank`: a float dataset of length equal to the number of features, containing the minimum rank of the log-fold changes across all pairwise comparisons for each feature.
- `delta_detected`: same as `lfc`, but for the delta-detected (i.e., difference in the percentage of detected expression).
- `cohen`: same as `lfc`, but for Cohen's d.
- `auc`: same as `lfc`, but for the AUCs.

## Changelog

- Version 2.0: support markers for multiple modalities.
- [Version 1.0](v1_0.md): added the `marker_detection` group.
