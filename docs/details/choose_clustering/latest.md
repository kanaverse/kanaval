# Choose clustering (latest) {#details-choose_clustering}

## Overview

We expect a `choose_clustering` HDF5 group at the root of the file, which defines the choice of clustering method.
The `choose_clustering` group itself contains the `parameters` and `results` subgroups.

## Parameters

`parameters` should contain:

- `method`: a scalar string specifying the clustering method to use.
  This is currently either `snn_graph` or `kmeans`.

## Results

`results` is empty.
 
Depending on the `method`, [`snn_graph_cluster`](../snn_graph_cluster/latest.md) or [`kmeans_cluster`](../kmeans_cluster/latest.md) must have non-empty `results`.
Both may also be non-empty, in which case the appropriate clustering is chosen based on `method`.

## Changelog

- Version 1.0: added the `choose_clustering` group.
