#ifndef KANAVAL_CELL_FILTERING_HPP
#define KANAVAL_CELL_FILTERING_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

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
    auto rhandle = utils::check_and_open_group(qhandle, "results");
    int remaining = -1;
    if (num_modalities > 1) {
        remaining = quality_control::check_discard_vector(rhandle, num_cells, false);
    }
    return remaining;
}
/**
 * @endcond
 */

/**
 * Check contents for the cell filtering step, see [here](@ref details-cell_filtering) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality filtering is applied.
 * @param num_modalities Number of modalities present in the dataset.
 * Note that only QC-relevant modalities need to be counted here;
 * if QC was skipped for a modality, it should not be included in this count.
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
    utils::check_and_open_group(qhandle, "parameters");

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
