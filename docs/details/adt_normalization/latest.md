# ADT normalization (latest) {#details-adt_normalization}

## Overview

We expect an `adt_normalization` group at the root of the file, containing information about the normalization of the ADT count matrix.
The group itself contains the `parameters` and `results` subgroups.

No ADT data was available prior to version 2.0 of the format, so the `adt_normalization` group may be absent in pre-v2.0 files.

## Parameters

No contents are mandated for `parameters`.

## Results

If ADTs are not present (i.e., `adt_in_use = false` in `adt_normalization::validate()`), `results` should be empty.

If `adt_in_use = true`, `results` should contain:

- `size_factors`, a float dataset of length equal to the number of cells after QC filtering (see `num_cells` in `adt_normalization::validate()`).
  This contains the size factor for each cell.

## Changelog

- Version 2.0: added the `adt_normalization` group.
