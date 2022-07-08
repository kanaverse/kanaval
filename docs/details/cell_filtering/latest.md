# Cell filtering (latest) {#details-cell_filtering}

## Overview

We expect a `cell_filtering` HDF5 group at the root of the file, containing information about the filtering of low-quality cells.
The itself contains the `parameters` and `results` subgroups.

Cell filtering was part of the `quality_control` step prior to version 2.0 of the format, so the `cell_filtering` group may be absent in pre-v2.0 files.
For such files, the `discards` vector is implicitly defined as the one from the `quality_control` group, see `quality_control::validate()` for details.

## Parameters

`parameters` should be empty.

## Results

If there are multiple QC-relevant modalities (i.e., `num_modalities > 1` in `cell_filtering::validate()`), `results` should contain:

- `discards`: an integer dataset of length equal to the number of cells (i.e., `num_cells` in `cell_filtering::validate()`).
  Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the filter thresholds.
  This is typically some function of the per-modality `discards` from `quality_control::validate()` and `adt_quality_control::validate()`.

Otherwise, `discards` may be absent, in which case the discard dataset is implicitly defined as the `discards` from the single modality.
See `quality_control::validate()` or `adt_quality_control::validate()` for more detials.

## Changelog

- Version 2.0: added the `cell_filtering` group.
