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

inline void validate_results(const H5::Group& handle, int num_pcs, int num_cells, int num_samples, const std::string& method) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    if (method == "mnn" && num_samples > 1) {
        std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(num_pcs) };
        utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, pdims);
    }

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the batch correction step.
 * Contents are stored inside an `batch_correction` HDF5 group at the root of the file.
 * The `batch_correction` group itself contains the `parameters` and `results` subgroups.
 *
 * A separate batch correction step was not used prior to version 2.0 of the format, so the `batch_correction` group may be absent in pre-v2.0 files.
 *
 * <HR>
 * `parameters` should contain:
 *
 * - `num_neighbors`: a scalar integer specifying the number of neighbors to use for mutual nearest neighbors correction.
 * - `approximate`: a scalar integer to be treated as a boolean, indicating whether an approximate neighbor search was used.
 * - `method`: a scalar string specifying the correction method to use, either `"none"` or `"mnn"`.
 *
 * <HR>
 * If `method = "mnn"` and `num_samples > 1`, `results` should contain:
 *
 * - `corrected`: a 2-dimensional float dataset containing the corrected PCs in a row-major layout.
 *   Each row corresponds to a cell (after QC filtering) and each column corresponds to a PC.
 *
 * Otherwise, correction is assumed to be a no-op and `results` may be empty.
 * Downstream steps should instead fetch coordinates from the `combined` dataset in `combine_embeddings::validate()`.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_pcs Number of PCs used in upstream steps.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 * @param num_samples Number of batches available in the dataset.
 * @param version Version of the state file.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, int num_pcs, int num_cells, int num_samples, int version) {
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
        validate_results(phandle, num_pcs, num_cells, num_samples, method);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'batch_correction'");
    }

    return;
}

}

}

#endif
