# Cell filtering (latest) {#details-cell_filtering}

## Overview

We expect a `cell_filtering` HDF5 group at the root of the file, containing information about the filtering of low-quality cells.
The itself contains the `parameters` and `results` subgroups.

Cell filtering was part of the `quality_control` step prior to version 2.0 of the format, so the `cell_filtering` group may be absent in pre-v2.0 files.
For such files, the `discards` vector is implicitly defined as the one from the [`quality_control`](../quality_control/latest.md) step.

**Definitions:**

- `num_modalities`: number of modalities that were used for quality control (i.e., not skipped).
  This is typically determined by examining [`inputs`](../inputs/latest.md) along with the [`quality_control`](../quality_control/latest.md) and [`adt_quality_control`](../adt_quality_control/latest.md) steps.
- `num_cells`: number of cells in the unfiltered dataset.
  This is typically determined from the [`inputs`](../inputs/latest.md) step.

## Parameters

`parameters` should be empty.

## Results

If `num_modalities > 1`, `results` should contain:

- `discards`: an integer dataset of length equal to `num_cells`.
  Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the filter thresholds.
  This is typically some function of the per-modality `discards` from the [`quality_control`](../quality_control/latest.md) or [`adt_quality_control`](../adt_quality_control/latest.md) steps.

If `num_modalities = 1`, `discards` may be absent.
In this case, the discard dataset is implicitly defined as the `discards` from the single modality that was used for QC;
see the [`quality_control`](../quality_control/latest.md) or [`adt_quality_control`](../adt_quality_control/latest.md) steps for more details.

If `num_modalities = 0`, `discards` may be absent.
In this case, the discard dataset is implicitly defined as an all-zero vector of length `num_cells`, i.e., no removal of any cells.

## Changelog

- Version 2.0: added the `cell_filtering` group.
