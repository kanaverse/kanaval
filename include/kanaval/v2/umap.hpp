#ifndef KANAVAL_UMAP_V2_HPP
#define KANAVAL_UMAP_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_umap(const H5::H5File& handle, int num_cells) {
    auto thandle = utils::check_and_open_group(handle, "umap");

    try {
        auto phandle = utils::check_and_open_group(thandle, "parameters");

        auto nn = utils::load_integer_scalar<>(phandle, "num_neighbors");
        if (nn <= 0) {
            throw std::runtime_error("'num_neighbors' value should be positive");
        }

        auto iters = utils::load_integer_scalar<>(phandle, "num_epochs");
        if (iters <= 0) {
            throw std::runtime_error("'num_epochs' should be positive");
        }

        auto md = utils::load_float_scalar<>(phandle, "min_dist");
        if (md <= 0) {
            throw std::runtime_error("'min_dist' should be positive");
        }

        utils::load_integer_scalar<>(phandle, "animate");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'umap'");
    }

    try {
        auto rhandle = utils::check_and_open_group(thandle, "results");

        std::vector<size_t> dims { static_cast<size_t>(num_cells) };
        utils::check_and_open_dataset(rhandle, "x", H5T_FLOAT, dims);
        utils::check_and_open_dataset(rhandle, "y", H5T_FLOAT, dims);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'umap'");
    }

    return;
}

}

}

#endif
