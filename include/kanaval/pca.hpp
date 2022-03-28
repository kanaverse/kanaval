#ifndef KANAVAL_FEATURE_SELECTION_HPP
#define KANAVAL_FEATURE_SELECTION_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file pca.hpp
 *
 * @brief Validate PCA contents.
 */

namespace kanaval {

namespace pca {

/**
 * @cond
 */
inline std::pair<int, std::string> validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto nhvgs = utils::load_integer_scalar<>(phandle, "num_hvgs");
    if (nhvgs <= 0) {
        throw std::runtime_error("number of HVGs must be positive in 'num_hvgs'");
    }

    auto npcs = utils::load_integer_scalar<>(phandle, "num_pcs");
    if (npcs <= 0) {
        throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
    }

    auto method = utils::load_string(phandle, "block_method");
    if (method != "none" && method != "regress" && method != "mnn") {
        throw std::runtime_error("unrecognized value '" + method + "' for the 'block_method'");
    }

    return std::make_pair(npcs, method);
}

inline void validate_results(const H5::Group& handle, int num_pcs, std::string block_method, int num_cells) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(num_pcs) };
    utils::check_and_open_dataset(rhandle, "pcs", H5T_FLOAT, pdims);

    std::vector<size_t> vdims { static_cast<size_t>(num_pcs) };
    utils::check_and_open_dataset(rhandle, "var_exp", H5T_FLOAT, vdims);

    if (block_method == "mnn") {
        utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, pdims);
    }

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the PCA step.
 *
 * `handle` should contain a `pca` group, itself containing the `parameters` and `results` subgroups.
 *
 * `parameters` should contain:
 *
 * - `num_hvgs`: a scalar integer containing the number of highly variable genes to use to compute the PCA.
 * - `num_pcs`: a scalar integer containing the number of PCs to compute.
 * - `block_method`: a scalar string specifying the method to use when dealing with multiple blocks in the dataset.
 *   This may be `"none"`, `"regress"` or `"mnn"`.
 *
 * `results` should contain:
 *
 * - `pcs`: a 2-dimensional float dataset containing the PC coordinates in a row-major layout.
 *   Each row corresponds to a cell (after QC filtering) and each column corresponds to a PC.
 *   Note that this is deliberately transposed from the Javascript/Wasm representation for easier storage.
 *   - If `block_type = "mnn"`, the PCs will be computed using a weighted method that adjusts for differences in the number of cells across blocks.
 *   - If `block_type = "regress"`, the PCs will be computed on the residuals after regressing out the block-wise effects.
 * - `var_exp`: a float dataset of length equal to the number of PCs, containing the percentage of variance explained by each PC.
 *
 * If `block_type = "mnn"`, the `results` group will also contain:
 *
 * - `corrected`, a float dataset with the same dimensions as `pcs`, containing the MNN-corrected PCs for each cell.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 *
 * @return The number of cells remaining after QC filtering.
 * If the format is invalid, an error is raised instead.
 */
inline void validate(const H5::H5File& handle, int num_cells) {
    auto phandle = utils::check_and_open_group(handle, "pca");

    int npcs;
    std::string bmethod;
    try {
        auto pout = validate_parameters(phandle);
        npcs = pout.first;
        bmethod = pout.second;
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'pca'");
    }

    try {
        validate_results(phandle, npcs, bmethod, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'pca'");
    }

    return;
}

}

}

#endif
