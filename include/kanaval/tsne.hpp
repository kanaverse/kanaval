#ifndef KANAVAL_TSNE_HPP
#define KANAVAL_TSNE_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file tsne.hpp
 *
 * @brief Validate t-SNE contents.
 */

namespace kanaval {

/**
 * Validation for t-SNE results
 */
namespace tsne {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    auto perp = utils::load_float_scalar<>(phandle, "perplexity");
    if (perp <= 0) {
        throw std::runtime_error("'perplexity' value should be positive");
    }

    auto iters = utils::load_integer_scalar<>(phandle, "iterations");
    if (iters <= 0) {
        throw std::runtime_error("'iterations' should be positive");
    }

    utils::load_integer_scalar<>(phandle, "animate");
    return; 
}

inline void validate_results(const H5::Group& handle, int num_cells) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    std::vector<size_t> dims { static_cast<size_t>(num_cells) };
    utils::check_and_open_dataset(rhandle, "x", H5T_FLOAT, dims);
    utils::check_and_open_dataset(rhandle, "y", H5T_FLOAT, dims);

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the t-SNE step, see [here](@ref details-tsne) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, int num_cells) {
    auto thandle = utils::check_and_open_group(handle, "tsne");

    try {
        validate_parameters(thandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'tsne'");
    }

    try {
        validate_results(thandle, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'tsne'");
    }

    return;
}

}

}

#endif
