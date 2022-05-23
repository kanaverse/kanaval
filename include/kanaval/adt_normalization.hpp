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

inline bool validate_results(const H5::Group& handle, size_t num_cells, bool adt_in_use) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    if (adt_in_use) {
        utils::check_and_open_dataset(rhandle, "size_factors", H5T_FLOAT, { num_cells });
    }
    return in_use;
}
/**
 * @endcond
 */

/**
 * Check contents for the ADT normalization step.
 * Contents are stored inside a `adt_normalization` HDF5 group at the root of the file.
 * The `adt_normalization` group itself contains the `parameters` and `results` subgroups.
 *
 * No ADT data was available prior to version 2.0 of the format, so the `adt_normalization` group may be absent in pre-v2.0 files.
 *
 * <HR>
 * No contents are mandated for `parameters`.
 *
 * <HR>
 * If `adt_in_use = false`, `results` should be empty.
 *
 * If `adt_in_use = true`, `results` should contain:
 *
 * - `size_factors`, a float dataset of length equal to the number of cells, containing the size factor for each cell.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset.
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
