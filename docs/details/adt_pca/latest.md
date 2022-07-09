# ADT PCA (latest) {#details-adt_pca}

## Overview

We expect an `adt_pca` group at the root of the file, containing information about the principal components analysis of the ADT count matrix.
The group itself contains the `parameters` and `results` subgroups.

No ADT data was available prior to version 2.0 of the format, so the `adt_pca` group may be absent in pre-v2.0 files.

**Definitions:**

- `adt_in_use`: whether ADTs are present in the dataset.
  This is typically determined by examining the [`inputs`](../inputs/latest.md). 
- `num_cells`: number of cells remaining after QC filtering.
  This is typically determined from the [`cell filtering`](../cell_filtering/latest.md) step.

## Parameters

`parameters` should contain:

- `num_pcs`: a scalar integer containing the maximum number of PCs to compute.
- `block_method`: a scalar string specifying the method to use when blocking on sample identity.
  This may be `"none"`, `"regress"` or `"weight"`.

## Results

If `adt_in_use = false`, `results` should be empty.

If `adt_in_use = true`, `results` should contain:

- `pcs`: a 2-dimensional float dataset containing the PC coordinates in a row-major layout.
  Each row corresponds to a cell after QC filtering, and the the total number of rows should be equal to `num_cells`.
  Each column corresponds to a PC.
  The number of PCs should be no greater than `num_pcs`, and may be less if not enough PCs are available in the original dataset.

  If `block_method = "weight"`, the PCs will be computed using a weighted method that adjusts for differences in the number of cells across blocks.
  If `block_method = "regress"`, the PCs will be computed on the residuals after regressing out the block-wise effects.
- `var_exp`: a float dataset of length equal to the number of PCs, containing the percentage of variance explained by each PC.

## Changelog

- Version 2.0: added the `adt_pca` group.
