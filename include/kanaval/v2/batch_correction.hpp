#ifndef KANAVAL_BATCH_CORRECTION_V2_HPP
#define KANAVAL_BATCH_CORRECTION_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_batch_correction(const H5::H5File& handle, int num_dims, int num_cells, int num_samples, int version) {
    if (version < 2000000) {
        return;
    }

    auto xhandle = utils::check_and_open_group(handle, "batch_correction");

    std::string method;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");

        auto nneighbors = utils::load_integer_scalar<>(phandle, "num_neighbors");
        if (nneighbors <= 0) {
            throw std::runtime_error("number of neighbors should be positive in 'num_neighbors'");
        }

        utils::check_and_open_scalar<>(phandle, "approximate", H5T_INTEGER);

        method = utils::load_string(phandle, "method");
        if (method != "none" && method != "mnn") {
            throw std::runtime_error("unrecognized value '" + method + "' for the 'method'");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'batch_correction'");
    }

    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");

        if (method == "mnn" && num_samples > 1) {
            std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(num_dims) };
            utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, pdims);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'batch_correction'");
    }

    return;
}

}

}

#endif
