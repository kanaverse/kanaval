#ifndef KANAVAL_NEIGHBOR_INDEX_HPP
#define KANAVAL_NEIGHBOR_INDEX_HPP

#include "H5Cpp.h"
#include "utils.hpp"

/**
 * @file neighbor_index.hpp
 *
 * @brief Validate neighbor index contents.
 */

namespace kanaval {

namespace neighbor_index {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    utils::check_and_open_scalar(phandle, "approximate", H5T_INTEGER);
    return;
}

inline void validate_results(const H5::Group& handle) {
    utils::check_and_open_group(handle, "results");
    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the neighbor index step.
 * Contents are stored inside an `neighbor_index` HDF5 group at the root of the file.
 * The `neighbor_index` group itself contains the `parameters` and `results` subgroups.
 *
 * <HR>
 * `parameters` should contain:
 * 
 * - `approximate`: an integer scalar to be interpreted as a boolean, specifying whether an approximate nearest neighbor search should be performed.
 * 
 * <HR>
 * No contents are mandated for `results`.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * 
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "neighbor_index");

    try {
        validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'neighbor_index'");
    }

    try {
        validate_results(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'neighbor_index'");
    }

    return;
}

}

}

#endif
