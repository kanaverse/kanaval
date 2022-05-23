#ifndef KANAVAL_ADT_PCA_HPP
#define KANAVAL_ADT_PCA_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file adt_pca.hpp
 *
 * @brief Validate ADT PCA contents.
 */

namespace kanaval {

namespace adt_pca {

/**
 * @cond
 */
inline int validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto npcs = utils::load_integer_scalar<>(phandle, "num_pcs");
    if (npcs <= 0) {
        throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
    }

    std::string method = utils::load_string(phandle, "block_method");
    if (method != "none" && method != "regress" && method != "weight") {
        throw std::runtime_error("unrecognized value '" + method + "' for the 'block_method'");
    }

    return npcs;
}

inline int validate_results(const H5::Group& handle, int num_pcs, int num_cells, bool adt_in_use) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    int obs_pcs = -1;

    if (adt_in_use) {
        auto vhandle = utils::check_and_open_dataset(rhandle, "var_exp", H5T_FLOAT);
        auto dspace = dhandle.getSpace();
        if (dspace.getSimpleExtentNdims() != 0) {
            throw std::runtime_error("'var_exp' dataset does not have the expected dimensions");
        }

        hsize_t observed;
        dspace.getSimpleExtentDims(&observed);
        if (static_cast<int>(observed) <= num_pcs) {
            throw std::runtime_error("length of 'var_exp' dataset exceeds the requested number of PCs");
        }
        
        utils::check_and_open_dataset(rhandle, "pcs", H5T_FLOAT, { static_cast<size_t>(num_cells), static_cast<size_t>(observed) });
        obs_pcs = observed;
    }

    return obs_pcs;
}
/**
 * @endcond
 */

/**
 * Check contents for the PCA step on the ADT log-normalized matrix.
 * Contents are stored inside an `adt_pca` HDF5 group at the root of the file.
 * The `adt_pca` group itself contains the `parameters` and `results` subgroups.
 *
 * <HR>
 * `parameters` should contain:
 *
 * - `num_pcs`: a scalar integer containing the maximum number of PCs to compute.
 * - `block_method`: a scalar string specifying the method to use when dealing with multiple blocks in the dataset.
 *   This may be `"none"`, `"regress"` or `"mnn"`.
 *
 * <HR>
 * `results` should contain:
 *
 * - `pcs`: a 2-dimensional float dataset containing the PC coordinates in a row-major layout.
 *   Each row corresponds to a cell (after QC filtering) and each column corresponds to a PC.
 *   The number of PCs should be no greater than `num_pcs` (and may be less if not enough PCs are available in the original dataset).
 *   If `block_method = "weight"`, the PCs will be computed using a weighted method that adjusts for differences in the number of cells across blocks.
 *   If `block_method = "regress"`, the PCs will be computed on the residuals after regressing out the block-wise effects.
 * - `var_exp`: a float dataset of length equal to the number of PCs, containing the percentage of variance explained by each PC.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 * @param version Version of the state file.
 *
 * @return The available number of PCs, or -1 if no ADTs were available for this dataset.
 * If the format is invalid, an error is raised instead.
 */
inline int validate(const H5::H5File& handle, int num_cells, bool adt_in_use, int version) {
    if (version < 2000000) {
        return -1;        
    }

    auto phandle = utils::check_and_open_group(handle, "adt_pca");

    int num_pcs;
    try {
        num_pcs = validate_parameters(phandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_pca'");
    }

    int obs_pcs;
    try {
        obs_pcs = validate_results(phandle, num_pcs, num_cells, adt_in_use);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_pca'");
    }

    return obs_pcs;
}

}

}

#endif
