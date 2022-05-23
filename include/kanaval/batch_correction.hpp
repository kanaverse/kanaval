#ifndef KANAVAL_BATCH_CORRECTION_HPP
#define KANAVAL_BATCH_CORRECTION_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file pca.hpp
 *
 * @brief Validate batch correction contents.
 */

namespace kanaval {

namespace pca {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
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

    return;
}

inline void validate_results(const H5::Group& handle, int num_pcs, int num_cells) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(num_pcs) };
    utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, pdims);

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the PCA step.
 * Contents are stored inside an `pca` HDF5 group at the root of the file.
 * The `pca` group itself contains the `parameters` and `results` subgroups.
 *
 * <HR>
 * `parameters` should contain:
 *
 * - `num_hvgs`: a scalar integer containing the number of highly variable genes to use to compute the PCA.
 * - `num_pcs`: a scalar integer containing the number of PCs to compute.
 * - \v1_1{\[**since version 1.1**\] `block_method`: a scalar string specifying the method to use when dealing with multiple blocks in the dataset.
 *   This may be `"none"`, `"regress"` or `"mnn"`.}
 *
 * <HR>
 * `results` should contain:
 *
 * - `pcs`: a 2-dimensional float dataset containing the PC coordinates in a row-major layout.
 *   Each row corresponds to a cell (after QC filtering) and each column corresponds to a PC.
 *   Note that this is deliberately transposed from the Javascript/Wasm representation for easier storage.
 *   \v1_1{\[**since version 1.1**\] If `block_type = "mnn"`, the PCs will be computed using a weighted method that adjusts for differences in the number of cells across blocks.
 *   If `block_type = "regress"`, the PCs will be computed on the residuals after regressing out the block-wise effects.}
 * - `var_exp`: a float dataset of length equal to the number of PCs, containing the percentage of variance explained by each PC.
 *
 * \v1_1{\[**since version 1.1**\] If `block_type = "mnn"`, the `results` group will also contain:}
 *
 * - \v1_1{`corrected`, a float dataset with the same dimensions as `pcs`, containing the MNN-corrected PCs for each cell.}
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 * @param version Version of the state file.
 *
 * @return The number of cells remaining after QC filtering.
 * If the format is invalid, an error is raised instead.
 */
inline void validate(const H5::H5File& handle, int num_pcs, int num_cells, int version) {
    if (version < 2000000) {
        return;
    }

    auto phandle = utils::check_and_open_group(handle, "batch_correction");

    try {
        validate_parameters(phandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'batch_correction'");
    }

    try {
        validate_results(phandle, num_pcs, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'batch_correction'");
    }

    return;
}

}

}

#endif
