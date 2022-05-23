#ifndef KANAVAL_QUALITY_CONTROL_HPP
#define KANAVAL_QUALITY_CONTROL_HPP

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
inline int validate_results(const H5::Group& qhandle, int num_cells, bool adt_in_use) {
    auto rhandle = qhandle.open("results");
    int remaining = -1;
    if (!adt_in_use) {
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
 * <HR>
 * `parameters` should be empty.
 *
 * <HR>
 * If `adt_in_use = true`, `results` should contain:
 *
 * - `discards`: an integer dataset of length equal to the number of cells.
 *   Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the filter thresholds.
 *   This is usually a union of the discarded sets from the upstream `quality_control` groups.
 *   
 * If `adt_in_use = false`, `discards` may be absent, in which case the discard dataset is implicitly defined as the `discards` from the upstream `quality_control` groups;
 * see `quality_control::validate()` and `adt_uality_control::validate()` for more details.
 * 
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality filtering is applied.
 * @param num_samples Number of batches in the dataset.
 * @param adt_in_use Whether ADTs are being used in this dataset.
 *
 * @return The number of cells remaining after filtering, or -1 if `discards` is absent.
 * If the format is invalid, an error is raised instead.
 */
inline int validate(const H5::H5File& handle, int num_cells, bool adt_in_use, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return -1;
    }

    auto qhandle = utils::check_and_open_group(handle, "cell_filtering");
    int remaining = 0;
    try {
        remaining = validate_results(qhandle, num_cells, adt_in_use);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'quality_control'");
    }

    return remaining;
}

}

}

#endif
