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
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto npcs = utils::load_integer_scalar(phandle, "num_pcs");
    if (npcs < 0){ 
        throw std::runtime_error("number of PCs used for ADT normalization should be positive");
    }

    auto nclust = utils::load_integer_scalar(phandle, "num_clusters");
    if (nclust < 0){ 
        throw std::runtime_error("number of clusters used for ADT normalization should be positive");
    }

    return;
}

inline void validate_results(const H5::Group& handle, size_t num_cells) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    if (rhandle.exists("size_factors")) {
        utils::check_and_open_dataset(rhandle, "size_factors", H5T_FLOAT, { num_cells });
    }
    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the ADT normalization step.
 * Contents are stored inside a `normalization` HDF5 group at the root of the file.
 * The `normalization` group itself contains the `parameters` and `results` subgroups.
 *
 * <HR>
 * No contents are mandated for `parameters`.
 *
 * <HR>
 * `results` should contain:
 *
 * - `size_factors`, a float dataset of length equal to the number of cells, containing the size factor for each cell.
 *
 * Alternatively, if the dataset does not contain ADTs, `size_factors` may be absent.
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset.
 * 
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, size_t num_cells) {
    auto nhandle = utils::check_and_open_group(handle, "normalization");

    try {
        validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'normalization'");
    }

    try {
        validate_results(nhandle, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'normalization'");
    }

    return;
}

}

}

#endif
