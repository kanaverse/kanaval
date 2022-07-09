#ifndef KANAVAL_FEATURE_SELECTION_HPP
#define KANAVAL_FEATURE_SELECTION_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file feature_selection.hpp
 *
 * @brief Validate feature selection contents.
 */

namespace kanaval {

/**
 * Validation for feature selection
 */
namespace feature_selection {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    auto span = utils::load_float_scalar<>(phandle, "span");
    if (span < 0 || span > 1) {
        throw std::runtime_error("LOWESS span should lie in [0, 1]");
    }
}

inline void validate_results(const H5::Group& handle, int num_genes) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    std::vector<size_t> dims{ static_cast<size_t>(num_genes) };
    utils::check_and_open_dataset(rhandle, "means", H5T_FLOAT, dims);
    utils::check_and_open_dataset(rhandle, "vars", H5T_FLOAT, dims);
    utils::check_and_open_dataset(rhandle, "fitted", H5T_FLOAT, dims);
    utils::check_and_open_dataset(rhandle, "resids", H5T_FLOAT, dims);

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the feature selection step on the log-expression matrix, see [here](@ref details-feature_selection) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_genes Number of genes in the dataset.
 * For multi-modal contexts, this should refer to the RNA modality only.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::Group& handle, int num_genes) {
    auto fhandle = utils::check_and_open_group(handle, "feature_selection");

    try {
        validate_parameters(fhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'feature_selection'");
    }

    try {
        validate_results(fhandle, num_genes);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'feature_selection'");
    }

    return;
}

}

}

#endif
