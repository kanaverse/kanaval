# Inputs (latest) {#details-inputs}

## Overview

An `inputs` group should be present at the root of the HDF5 file, containing all information relating to the data inputs.
This includes the definition of the input files as well as summary statistics like the number of cells and features.
The group itself contains the `parameters` and `results` subgroups.

In this section, a "matrix" refers to one or more files describing a single (count) matrix.
This should be exactly one file for HDF5-based formats, or multiple files for MatrixMarket formats, e.g., to include feature information - see below for details.
Multiple matrices may be supplied, in which case each matrix is assumed to contain data for one sample.
Alternatively, a single matrix may contain data for one or more samples as partitioned by a sample factor.

The loaded dataset refers to the in-memory representation of the matrix (for single matrix inputs) or the combined matrices (for multiple inputs).
The identities of the rows of the loaded dataset may be a permutation or subset of the rows in the input matrices.
This is especially true for multiple inputs where the loaded dataset only contains the intersection of features across inputs.
For multiple matrices, the loaded dataset is assumed to be a column-wise concatenation of the individual matrices after sorting them by `sample_names`.

**Definitions:**

- `embedded`: whether the input files were embedded in the `*.kana` file.
  If `false`, the input files are assumed to be linked instead.

## Parameters

`parameters` should contain:

- `format`: a scalar string specifying the file format for a single matrix.
  This is usually either `"MatrixMarket"`, for a MatrixMarket file with possible feature/barcode annotations;
  `"10X"`, for the 10X Genomics HDF5 matrix format;
  or `"H5AD"`, for the H5AD format.
  Other values are allowed but their interpretation is implementation-defined (e.g., for custom resources). 
  For multiple matrices, `format` should instead be a 1-dimensional string dataset of length equal to the number of uploads.
  Each element of the dataset is usually one of `"MatrixMarket"`, `"10X"` or `"H5AD"`; 
  different values can be present for mixed input formats.
- `files`: a group of groups representing an array of input file information.
  Each inner group is named by their positional index in the array and contains information about a file in an upload.
  Each inner group should contain:
  - `type`: a scalar string specifying the type of the file.
    For multiple matrices, the constraints below apply to all files corresponding to a single matrix.
    - If `format = "MatrixMarket"`, there should be exactly one `type = "mtx"` corresponding to the (possibly Gzipped) `*.mtx` file.
      There may be zero or one `type = "genes"`, containing a (possibly Gzipped) TSV file with the Ensembl and gene symbols for each row.
      There may be zero or one `type = "annotations"`, containing a (possibly Gzipped) TSV file with the annotations for each column.
    - If `format = "10X"` or `"H5AD"`, there should be exactly one `type = "h5"`.
    - For other `format`s, any `type` can be used, typically for custom resources.
  - `name`: a scalar string specifying the file name as it was provided to **kana**.

  If `embedded = true`, we additionally expect:
  - `offset`: a scalar integer specifying where the file starts as an offset from the start of the remaining bytes section.
    The offset for the first file should be zero, and entries in `files` should be ordered by increasing `offset`.
  - `size`: a non-negative scalar integer specifying the number of bytes in the file.
    The offset of each file should be equal to the sum of `size` and `offset` for the preceding file.

  If `embedded = false`, we expect:
  - `id`: a scalar string containing some unique identifier for this file.
    The interpretation of `id` is application-specific but usually refers to some cache or database.

Optionally, `parameters` may contain a `subset` group.
This may in turn contain a `cells` group, specifying the subset of cells to be used in downstream steps.
The `subset/cells` group should contain one of:

- `indices`, a 1-dimensional integer dataset specifying the indices of the cells to retain.
  Indices should be relative to the loaded dataset before subsetting.
  The length of this dataset should be equal to the number of cells in `results/num_cells`.
  Indices should be unique and sorted.
- `field`, a string scalar specifying the field of the column annotations to use for filtering.
  This field can either be a continuous or categorical variable.
  - If categorical, the group should also contain `values`, 
    a 1-dimensional string dataset specifying the allowable values for this field.
    The subset is defined as all cells with a `field` value in `values`.
  - If continuous, the group should also contain `ranges`,
    a 2-dimensional float dataset specifying the ranges of allowable values.
    Each row specifies the range `[a, b)` where `a` is the value in the first column and `b` is the value in the second column.
    Ranges should be non-overlapping and sorted in order of increasing `a`.
    The subset is defined as all cells that fall in any of the specified ranges.

For multiple matrices, `parameters` should also contain:

- `sample_groups`: an integer dataset of length equal to the number of samples.
  Each entry specifies the number of files in `files` that belong to a sample.
  (All files from the same sample are assumed to be contiguous in the array represented by `files`;
  so a `sample_groups` of `[3, 2, 1]` would mean that the first three files belong to the first sample, 
  the next 2 files belong to the second sample, and the last file belongs to the third sample.)
- `sample_names`: a string dataset of length equal to the number of samples.
  Each value contains the name for the sample defined by the corresponding entry of `sample_groups`.
  All names should be unique and sorted.

For single matrix inputs, `parameters` may also contain:

- `sample_factor`: a string scalar specifying the field in the per-cell annotation that contains the sample blocking factor. 
  If present, it is assumed that the matrix contains data for multiple samples.

## Results

`results` should contain:

- `num_cells`: an integer scalar specifying the number of cells in the loaded dataset.
- `num_features`: a group containing integer scalar datasets, each named after a modality (typically `"RNA"` and/or `"ADT"`).
  Each dataset contains the number of features for its named modality.
  When dealing with multiple matrix inputs, the number of features is defined as the intersection across all matrices.
- `identities`: a group containin 1-dimensional integer datasets, each named after a modality in `num_features`.
  Each dataset is of length equal to the number of features listed in its corresponding `num_features` dataset.
  Each dataset contains the identities of the rows in the loaded dataset corresponding to its named modality.
  If a single input was provided, `identities` identifies each row in terms of its index in the "original" input matrix (i.e., if it were loaded without any modification).
  If multiple inputs were provided, `identities` contains the intersection of features across inputs, and each value refers to the row index in the original matrix of the _first_ input.
  Row identities are parallel to the per-feature results in subsequent analysis steps.

All steps that generate per-gene results should use `identities` (or, for older formats, `permutation` or `indices`) to identify the genes corresponding to the statistics.
This includes [`feature_selection`](../feature_selection/latest.md), [`marker_detection`](../marker_detection/latest.md) and [`custom_selections`](../custom_selections/latest.md).

## Changelog

- Version 2.1: support arbitrary subsetting of the loaded matrix by cells.
- [Version 2.0](v2_0.md): store the number of features for multiple modalities.
- [Version 1.2](v1_2.md): streamlined results for multiple matrices and samples. 
- [Version 1.1](v1_1.md): added support for multiple matrices and samples. 
- [Version 1.0](v1_0.md): added the `inputs` group. 
