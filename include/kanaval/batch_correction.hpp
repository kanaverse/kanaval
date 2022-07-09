#ifndef KANAVAL_BATCH_CORRECTION_HPP
#define KANAVAL_BATCH_CORRECTION_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file batch_correction.hpp
 *
 * @brief Validate batch correction contents.
 */

namespace kanaval {

/** 
 * Validation for the batch correction
 */
namespace batch_correction {

/**
 * @cond
 */
inline std::string validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto nneighbors = utils::load_integer_scalar<>(phandle, "num_neighbors");
    if (nneighbors <= 0) {
        throw std::runtime_error("number of neighbors should be positive in 'num_neighbors'");
    }

    utils::check_and_open_scalar<>(phandle, "approximate", H5T_INTEGER);

    std::string method = utils::load_string(phandle, "method");
    if (method != "none" && method != "mnn") {
        throw std::runtime_error("unrecognized value '" + method + "' for the 'method'");
    }

    return method;
}

inline void validate_results(const H5::Group& handle, int num_dims, int num_cells, int num_samples, const std::string& method) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    if (method == "mnn" && num_samples > 1) {
        std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(num_dims) };
        utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, pdims);
    }

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the batch correction step, see [here](@ref details-batch_correction.md) for details.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_dims Number of dimensions of interest.
 * @param num_cells Number of high-quality cells in the dataset, i.e., after any quality-based filtering has been applied.
 * @param num_samples Number of batches available in the dataset.
 * @param version Version of the state file.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, int num_dims, int num_cells, int num_samples, int version) {
    if (version < 2000000) {
        return;
    }

    auto phandle = utils::check_and_open_group(handle, "batch_correction");

    std::string method;
    try {
        method = validate_parameters(phandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'batch_correction'");
    }

    try {
        validate_results(phandle, num_dims, num_cells, num_samples, method);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'batch_correction'");
    }

    return;
}

}

}

#endif
