# PCA (latest) {#details-pca}

## Overview

We expect an `pca` HDF5 group at the root of the file, describing how the PCA was performed on the RNA expression matrix.
The group itself contains the `parameters` and `results` subgroups.

**Definitions:**

- `num_cells`: number of cells remaining after QC filtering.
  This is typically determined from the [`cell_filtering`](../cell_filtering/latest.md) step.

## Parameters 

`parameters` should contain:

- `num_hvgs`: a scalar integer containing the number of highly variable genes to use to compute the PCA.
- `num_pcs`: a scalar integer containing the number of PCs to compute.
- `block_method`: a scalar string specifying the method to use when dealing with multiple blocks in the dataset.
  This may be `"none"`, `"regress"` or `"weight"`.

## Results

`results` should contain:

- `pcs`: a 2-dimensional float dataset containing the PC coordinates in a row-major layout.
  Each row corresponds to a cell (after QC filtering) and each column corresponds to a PC.
  The total number of rows and columns is defined as `num_cells` and `num_pcs`, respectively.

  If `block_method = "weight"`, the PCs will be computed using a weighted method that adjusts for differences in the number of cells across blocks.
  If `block_method = "regress"`, the PCs will be computed on the residuals after regressing out the block-wise effects.
- `var_exp`: a float dataset of length equal to the number of PCs, containing the percentage of variance explained by each PC.

## Changelog

- Version 2.0: moved the MNN-corrected PCs to [`batch_correction`](../batch_correction/latest.md).
- [Version 1.1](v1_1.md): added support for various forms of batch correction.
- [Version 1.0](v1_0.md): added the `pca` group.
