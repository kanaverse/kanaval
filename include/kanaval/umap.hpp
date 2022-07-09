#ifndef KANAVAL_UMAP_HPP
#define KANAVAL_UMAP_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file umap.hpp
 *
 * @brief Validate UMAP contents.
 */

namespace kanaval {

/**
 * Validation for UMAP results
 */
namespace umap {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto nn = utils::load_integer_scalar<>(phandle, "num_neighbors");
    if (nn <= 0) {
        throw std::runtime_error("'num_neighbors' value should be positive");
    }

    auto iters = utils::load_integer_scalar<>(phandle, "num_epochs");
    if (iters <= 0) {
        throw std::runtime_error("'num_epochs' should be positive");
    }

    auto md = utils::load_float_scalar<>(phandle, "min_dist");
    if (md <= 0) {
        throw std::runtime_error("'min_dist' should be positive");
    }

    utils::load_integer_scalar<>(phandle, "animate");
    return; 
}

inline void validate_results(const H5::Group& handle, int num_cells) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    std::vector<size_t> dims { static_cast<size_t>(num_cells) };
    utils::check_and_open_dataset(rhandle, "x", H5T_FLOAT, dims);
    utils::check_and_open_dataset(rhandle, "y", H5T_FLOAT, dims);

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the UMAP step, see [here](@ref details-umap) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, int num_cells) {
    auto thandle = utils::check_and_open_group(handle, "umap");

    try {
        validate_parameters(thandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'umap'");
    }

    try {
        validate_results(thandle, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'umap'");
    }

    return;
}

}

}

#endif
