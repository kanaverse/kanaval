# ADT quality control (latest) {#details-adt_quality_control}

## Overview

We expect an `adt_quality_control` HDF5 group at the root of the file, containing information about the quality control metrics and filters derived from the ADT counts.
The group itself contains the `parameters` and `results` subgroups.

No ADT data was available prior to version 2.0 of the format, so the `adt_quality_control` group may be absent in such files.

## Parameters

`parameters` should contain:

- `igg_prefix`: a scalar string containing the expected prefix for IgG features.
- `nmads`: a scalar float specifying the number of MADs to use to define the QC thresholds.
- `min_detected_drop`: a scalar float specifying the minimum relative drop in the number of detected features before a cell is considered to be low-quality.
- `skip`: a scalar integer to be interpreted as a boolean, specifying whether to skip the QC for the ADTs.

## Results

If ADTs are not available (i.e., `adt_in_use = false` in `adt_quality_control::validate()`), `results` should be empty.

If `adt_in_use = true`, `results` should contain:

- `metrics`, a group containing per-cell QC metrics derived from the RNA count data.
  This contains:
  - `sums`: a float dataset of length equal to the number of cells (i.e., `num_cells` in `adt_quality_control::validate()`), containing the total count for each cell.
  - `detected`:  an integer dataset of length equal to the number of cells, containing the total number of detected features for each cell.
  - `igg_total`: a float dataset of length equal to the number of cells, containing the total count in IgG features.
- `thresholds`, a group containing thresholds on the metrics for each sample.
  (The number of samples is defined from `num_samples` in `adt_quality_control::validate()`.)
  This group contains:
  - `detected`:  a float dataset of length equal to the number of samples, containing the threshold on the total number of detected features for each sample.
  - `igg_total`: a float dataset of length equal to the number of samples, containing the threshold on the total counts in IgG features for each sample.
- `discards`: an integer dataset of length equal to the number of cells.
  Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the ADT-based filter thresholds.

If `adt_in_use = true` and `skip = true`, `results` may be an empty group.
However, if any of `metrics`, `thresholds` or `discards` is present, they should follow the requirements listed above.

## Changelog

- Version 2.1: allow the QC step to be skipped.
- [Version 2.0](v2_0.md): added the `adt_quality_control` group.
