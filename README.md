# Validating `kana` files

## Overview

The `kana` file format contains an embedded HDF5 file that captures the analysis state of the [**kana**](https://github.com/jkanche/kana) application.
This embedded state file stores the parameters and results for each step in a simple single-cell RNA-seq analysis.
By storing this state, we can easily reload existing analyses into **kana** without recomputation.
It is also straightforward to extract results from this state in other data analysis frameworks (e.g., R/Bioconductor).

The **kanaval** repository contains a specification of the expected structure and content of the state file.
We use a specification-as-code approach that enforces the specification with a validator library, implemented with header-only C++ for portablity to any system that supports a foreign function interface. 
It is thus possible to create `kana` files from other languages, validate them, and upload them to **kana**.

## Background on the `kana` layout

This section describes version 1 of the `kana` export format.

The first 8 bytes define an unsigned 64-bit integer in little-endian, specifying the format type.
This is used to denote whether the input data files are embedded (0) or linked (1);
the former is used for export to a standalone file while the latter is used to save the state to the browser's cache.

The next 8 bytes define another unsigned 64-bit integer describing the format version.
We use semantic versioning where each version number is described by 3 digits, i.e., `XXXYYYZZZ`.
This document will only consider version 1, i.e., `001YYYZZZ`.

The next 8 bytes define another unsigned 64-bit integer specifying the size of the HDF5 file containing the analysis state.
Let's call this value `state_nbytes`.

The next `state_nbytes` bytes contain a HDF5 state file.
Each analysis step is represented by a HDF5 group that contains the parameters and results.
See the next section for details on the expected groups.

The remaining bytes contain the embedded input files when dealing with an embedded format type.
Each file can be excised by reading the offsets and sizes in the `inputs` group in the state file.

<img src="https://raw.githubusercontent.com/LTLA/kanaval/master/docs/layout.png" width="500">

## Structure of the HDF5 state file

Inside the HDF5 state file, each analysis step is represented by a HDF5 group.

**Version 3.0:**

- [Inputs](docs/details/inputs/v3_0.md)
- [RNA quality control](docs/details/rna_quality_control/v3_0.md)
- [ADT quality control](docs/details/adt_quality_control/v3_0.md)
- [CRISPR quality control](docs/details/crispr_quality_control/v3_0.md)
- [Cell filtering](docs/details/cell_filtering/v3_0.md)
- [RNA normalization](docs/details/rna_normalization/v3_0.md)
- [ADT normalization](docs/details/adt_normalization/v2_0.md)
- [CRISPR normalization](docs/details/crispr_normalization/v3_0.md)
- [Feature selection](docs/details/feature_selection/v1_0.md)
- [RNA PCA](docs/details/rna_pca/v3_0.md)
- [ADT PCA](docs/details/adt_pca/v2_0.md)
- [CRISPR PCA](docs/details/crispr_pca/v3_0.md)
- [Combine embeddings](docs/details/combine_embeddings/v3_0.md)
- [Batch correction](docs/details/batch_correction/v2_0.md)
- [Neighbor index](docs/details/neighbor_index/v2_0.md)
- [k-means clustering](docs/details/kmeans_cluster/v1_0.md)
- [SNN graph clustering](docs/details/snn_graph_cluster/v3_0.md)
- [Choose clustering](docs/details/choose_clustering/v1_0.md)
- [Marker detection](docs/details/marker_detection/v3_0.md)
- [Custom selections](docs/details/custom_selections/v3_0.md)
- [Cell labelling](docs/details/cell_labelling/v1_0.md)
- [t-SNE](docs/details/tsne/v1_0.md)
- [UMAP](docs/details/umap/v1_0.md)

**Version 2.1:**

- [Inputs](docs/details/inputs/v2_1.md)
  - [v2.0](docs/details/inputs/v2_0.md)
- [RNA quality control](docs/details/quality_control/v2_1.md)
  - [v2.0](docs/details/quality_control/v1_0.md)
- [ADT quality control](docs/details/adt_quality_control/v2_1.md)
  - [v2.0](docs/details/adt_quality_control/v2_0.md)
- [Cell filtering](docs/details/cell_filtering/v2_0.md)
- [RNA normalization](docs/details/normalization/v1_0.md)
- [ADT normalization](docs/details/adt_normalization/v2_0.md)
- [Feature selection](docs/details/feature_selection/v1_0.md)
- [RNA PCA](docs/details/pca/v2_0.md)
- [ADT PCA](docs/details/adt_pca/v2_0.md)
- [Combine embeddings](docs/details/combine_embeddings/v2_0.md)
- [Batch correction](docs/details/batch_correction/v2_0.md)
- [Neighbor index](docs/details/neighbor_index/v2_0.md)
- [k-means clustering](docs/details/kmeans_cluster/v1_0.md)
- [SNN graph clustering](docs/details/snn_graph_cluster/v1_0.md)
- [Choose clustering](docs/details/choose_clustering/v1_0.md)
- [Marker detection](docs/details/marker_detection/v2_0.md)
- [Custom selections](docs/details/custom_selections/v2_1.md)
  - [v2.0](docs/details/custom_selections/v2_0.md)
- [Cell labelling](docs/details/cell_labelling/v1_0.md)
- [t-SNE](docs/details/tsne/v1_0.md)
- [UMAP](docs/details/umap/v1_0.md)

**Version 1.2:**

- [Inputs](docs/details/inputs/v1_2.md)
  - [v1.1](docs/details/inputs/v1_1.md)
  - [v1.0](docs/details/inputs/v1_0.md)
- [Quality control](docs/details/quality_control/v1_0.md)
- [Normalization](docs/details/normalization/v1_0.md)
- [Feature selection](docs/details/feature_selection/v1_0.md)
- [PCA](docs/details/pca/v1_1.md)
  - [v1.0](docs/details/v1_0.md)
- [k-means clustering](docs/details/kmeans_cluster/v1_0.md)
- [SNN graph clustering](docs/details/snn_graph_cluster/v1_0.md)
- [Choose clustering](docs/details/choose_clustering/v1_0.md)
- [Marker detection](docs/details/marker_detection/v1_0.md)
- [Custom selections](docs/details/custom_selections/v1_0.md)
- [Cell labelling](docs/details/cell_labelling/v1_0.md)
- [t-SNE](docs/details/tsne/v1_0.md)
- [UMAP](docs/details/umap/v1_0.md)

## Running the validator

Calling the [`validate()`](https://ltla.github.io/kanaval/validate_8hpp.html) function will validate the state file,
which will throw a reasonably informative error if there are any problems.

```cpp
#include "H5Cpp.h"
#include "kanaval/validate.hpp"

H5::H5File handle(path, H5F_ACC_RDONLY);
kanaval::validate(handle, embedded, version);
```

