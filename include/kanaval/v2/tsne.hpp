#ifndef KANAVAL_TSNE_V2_HPP
#define KANAVAL_TSNE_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_tsne(const H5::H5File& handle, int num_cells) {
    auto thandle = utils::check_and_open_group(handle, "tsne");

    try {
        auto phandle = utils::check_and_open_group(thandle, "parameters");

        auto perp = utils::load_float_scalar<>(phandle, "perplexity");
        if (perp <= 0) {
            throw std::runtime_error("'perplexity' value should be positive");
        }

        auto iters = utils::load_integer_scalar<>(phandle, "iterations");
        if (iters <= 0) {
            throw std::runtime_error("'iterations' should be positive");
        }

        utils::load_integer_scalar<>(phandle, "animate");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'tsne'");
    }

    try {
        auto rhandle = utils::check_and_open_group(thandle, "results");

        std::vector<size_t> dims { static_cast<size_t>(num_cells) };
        utils::check_and_open_dataset(rhandle, "x", H5T_FLOAT, dims);
        utils::check_and_open_dataset(rhandle, "y", H5T_FLOAT, dims);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'tsne'");
    }

    return;
}

}

}

#endif
