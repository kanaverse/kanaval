#ifndef KANAVAL_ADT_NORMALIZATION_HPP
#define KANAVAL_ADT_NORMALIZATION_HPP

#include "H5Cpp.h"
#include "utils.hpp"

/**
 * @file adt_normalization.hpp
 *
 * @brief Validate ADT normalization contents.
 */

namespace kanaval {

/** 
 * Validation for ADT normalization
 */
namespace adt_normalization {

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

inline void validate_results(const H5::Group& handle, size_t num_cells, bool adt_in_use) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    if (adt_in_use) {
        utils::check_and_open_dataset(rhandle, "size_factors", H5T_FLOAT, { num_cells });
    }
    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the normalization step on the ADT count matrix, see [here](@ref details-adt_normalization) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of high-quality cells in the dataset, i.e., after any quality-based filtering has been applied.
 * @param adt_in_use Whether ADTs are being used in this dataset.
 * @param version Version of the format.
 * 
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, size_t num_cells, bool adt_in_use, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return;
    }

    auto nhandle = utils::check_and_open_group(handle, "adt_normalization");

    try {
        validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_normalization'");
    }

    try {
        validate_results(nhandle, num_cells, adt_in_use);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_normalization'");
    }

    return;
}

}

}

#endif
