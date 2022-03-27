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
 * Check contents for the feature selection step.
 * 
 * @param handle An open HDF5 file handle.
 * @param num_genes Number of genes in the dataset.
 *
 * @return If the format is invalid, an error is raised.
 *
 * @description
 * `handle` should contain a `feature_selection` group, itself containing the `parameters` and `results` subgroups.
 *
 * `parameters` should contain:
 * 
 * - `span`: a scalar float specifying the span to use for the LOWESS smoother.
 * 
 * `results` should contain:
 * 
 * - `means`: a 1-dimensional float dataset of length equal to the number of genes,
 *   containing the mean log-expression of each gene.
 * - `vars`: a 1-dimensional float dataset of length equal to the number of genes,
 *   containing the variance in log-expression of each gene.
 * - `fitted`: a 1-dimensional float dataset of length equal to the number of genes,
 *   containing the fitted value of the trend for each gene.
 * - `resids`: a 1-dimensional float dataset of length equal to the number of genes,
 *   containing the residuals from the trend for each gene.
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
