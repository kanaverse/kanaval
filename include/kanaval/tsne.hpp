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
 * Check contents for the t-SNE step.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 *
 * @return If the format is invalid, an error is raised.
 *
 * @details
 * `handle` should contain a `tsne` group, itself containing the `parameters` and `results` subgroups.
 *
 * `parameters` should contain:
 *
 * - `perplexity`: a scalar float specifying the t-SNE perplexity.
 * - `iterations`: a scalar integer specifying the t-SNE iterations.
 * - `animate`: a scalar integer to be interpreted as a boolean, indicating whether an animation should be performed.
 *
 * `results` should contain:
 *
 * - `x`: a float dataset of length equal to the number of cells (after QC filtering), containing the x-coordinates for each cell.
 * - `y`: a float dataset of length equal to the number of cells (after QC filtering), containing the y-coordinates for each cell.
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
