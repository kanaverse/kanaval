#ifndef KANAVAL_CELL_FILTERING_HPP
#define KANAVAL_CELL_FILTERING_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file cell_filtering.hpp
 *
 * @brief Validate cell filtering contents.
 */

namespace kanaval {

namespace cell_filtering {

/**
 * @cond
 */
inline int validate_results(const H5::Group& qhandle, int num_cells, int num_modalities) {
    auto rhandle = qhandle.openGroup("results");
    int remaining = -1;
    if (num_modalities > 1) {
        remaining = quality_control::check_discard_vector(rhandle, num_cells);
    }
    return remaining;
}
/**
 * @endcond
 */

/**
 * Check contents for the cell filtering step.
 * Contents are stored inside a `cell_filtering` HDF5 group at the root of the file.
 * The `cell_filtering` group itself contains the `parameters` and `results` subgroups.
 *
 * Cell filtering was rolled into the `quality_control` step prior to version 2.0 of the format, so the `cell_filtering` group may be absent in pre-v2.0 files.
 *
 * <HR>
 * `parameters` should be empty.
 *
 * <HR>
 * If `num_modalities > 1`, `results` should contain:
 *
 * - `discards`: an integer dataset of length equal to the number of cells.
 *   Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the filter thresholds.
 *   This is usually a union of the discarded sets from `quality_control::validate()` and `adt_quality_control::validate()`.
 *   
 * Otherwise, `discards` may be absent, in which case the discard dataset is implicitly defined as the `discards` from the QC group of the available modality,
 * i.e., `quality_control::validate()` or `adt_quality_control::validate()`.
 * 
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality filtering is applied.
 * @param num_modalities Number of (QC-relevant) modalities present in the dataset.
 * @param version Version of the state file.
 *
 * @return The number of cells remaining after filtering, or -1 if `discards` is absent.
 * If the format is invalid, an error is raised instead.
 */
inline int validate(const H5::H5File& handle, int num_cells, int num_modalities, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return -1;
    }

    auto qhandle = utils::check_and_open_group(handle, "cell_filtering");
    int remaining = 0;
    try {
        remaining = validate_results(qhandle, num_cells, num_modalities);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'quality_control'");
    }

    return remaining;
}

}

}

#endif
