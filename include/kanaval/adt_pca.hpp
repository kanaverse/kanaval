#ifndef KANAVAL_ADT_PCA_HPP
#define KANAVAL_ADT_PCA_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

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
inline int validate_parameters(const H5::Group& handle, int version) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto npcs = utils::load_integer_scalar<>(phandle, "num_pcs");
    if (npcs <= 0) {
        throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
    }

    std::string method = utils::load_string(phandle, "block_method");
    pca::check_block_method(method, version);

    return npcs;
}

inline int validate_results(const H5::Group& handle, int max_pcs, int num_cells, bool adt_in_use) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    int obs_pcs = -1;
    if (adt_in_use) {
        obs_pcs = pca::check_pca_contents(rhandle, max_pcs, num_cells);
    }

    return obs_pcs;
}
/**
 * @endcond
 */

/**
 * Check contents for the PCA step on the ADT log-normalized abundance matrix.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of high-quality cells in the dataset, i.e., after any quality-based filtering has been applied.
 * @param adt_in_use Whether ADTs are being used in this dataset.
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
        num_pcs = validate_parameters(phandle, version);
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
