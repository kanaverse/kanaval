# ADT normalization (latest) {#details-adt_normalization}

## Overview

We expect an `adt_normalization` group at the root of the file, containing information about the normalization of the ADT count matrix.
The group itself contains the `parameters` and `results` subgroups.

No ADT data was available prior to version 2.0 of the format, so the `adt_normalization` group may be absent in pre-v2.0 files.

**Definitions:**

- `adt_in_use`: whether ADTs are present in the dataset.
  This is typically determined by examining the [`inputs`](../inputs/latest.md) step. 
- `num_cells`: number of cells remaining after QC filtering.
  This is typically determined from the [`cell_filtering`](../cell_filtering/latest.md) step.

## Parameters

No contents are mandated for `parameters`.

## Results

If `adt_in_use = false`, `results` should be empty.

If `adt_in_use = true`, `results` should contain:

- `size_factors`, a 1-dimensional float dataset of length equal to `num_cells`, containing the size factor for each cell.

## Changelog

- Version 2.0: added the `adt_normalization` group.
