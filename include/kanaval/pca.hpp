#ifndef KANAVAL_PCA_HPP
#define KANAVAL_PCA_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

/**
 * @file pca.hpp
 *
 * @brief Validate PCA contents.
 */

namespace kanaval {

/**
 * Validation for the RNA-based PCA
 */
namespace pca {

/**
 * @cond
 */
inline std::pair<int, std::string> validate_parameters(const H5::Group& handle, int version) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto nhvgs = utils::load_integer_scalar<>(phandle, "num_hvgs");
    if (nhvgs <= 0) {
        throw std::runtime_error("number of HVGs must be positive in 'num_hvgs'");
    }

    auto npcs = utils::load_integer_scalar<>(phandle, "num_pcs");
    if (npcs <= 0) {
        throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
    }

    std::string method;
    if (version >= 1001000) {
        method = utils::load_string(phandle, "block_method");
        check_block_method(method, version);
    }

    return std::make_pair(npcs, method);
}

inline int validate_results(const H5::Group& handle, int max_pcs, std::string block_method, int num_cells, int version) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    int obs_pcs = check_pca_contents(rhandle, max_pcs, num_cells);

    if (version >= 1001000 && version < 2000000) {
        if (block_method == "mnn") {
            utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, { static_cast<size_t>(num_cells), static_cast<size_t>(obs_pcs) });
        }
    }

    return obs_pcs;
}
/**
 * @endcond
 */

/**
 * Check contents for the PCA step on the RNA log-expression matrix, see [here](@ref details-pca) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 * @param version Version of the state file.
 *
 * @return The number of computed PCs.
 * If the format is invalid, an error is raised instead.
 */
inline int validate(const H5::H5File& handle, int num_cells, int version = 1001000) {
    auto phandle = utils::check_and_open_group(handle, "pca");

    int npcs;
    std::string bmethod;
    try {
        auto pout = validate_parameters(phandle, version);
        npcs = pout.first;
        bmethod = pout.second;
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'pca'");
    }

    int obs_pcs;
    try {
        obs_pcs = validate_results(phandle, npcs, bmethod, num_cells, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'pca'");
    }

    return obs_pcs;
}

}

}

#endif
