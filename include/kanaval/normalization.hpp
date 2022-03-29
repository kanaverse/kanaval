#ifndef KANAVAL_NORMALIZATION_HPP
#define KANAVAL_NORMALIZATION_HPP

#include "H5Cpp.h"
#include "utils.hpp"

/**
 * @file normalization.hpp
 *
 * @brief Validate normalization contents.
 */

namespace kanaval {

namespace normalization {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
    utils::check_and_open_group(handle, "parameters");
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
 * Check contents for the normalization step.
 * Contents are stored inside a `normalization` HDF5 group at the root of the file.
 * The `normalization` group itself contains the `parameters` and `results` subgroups.
 *
 * <HR>
 * No contents are mandated for `parameters`.
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
    auto nhandle = utils::check_and_open_group(handle, "normalization");

    try {
        validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'normalization'");
    }

    try {
        validate_results(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'normalization'");
    }

    return;
}

}

}

#endif
