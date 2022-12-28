#ifndef KANAVAL_COMBINE_EMBEDDINGS_V2_HPP
#define KANAVAL_COMBINE_EMBEDDINGS_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include <numeric>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_combine_embeddings(const H5::H5File& handle, int num_cells, const std::vector<std::string>& modalities, int total_dims, int version) {
    if (version < 2000000) {
        return;
    }

    auto chandle = utils::check_and_open_group(handle, "combine_embeddings");

    try {
        auto phandle = utils::check_and_open_group(chandle, "parameters");
        utils::check_and_open_scalar(phandle, "approximate", H5T_INTEGER);

        auto whandle = utils::check_and_open_group(phandle, "weights");
        if (whandle.getNumObjs()) { // if non-zero weights, all modalities should be present.
            for (const auto& m : modalities) {
                utils::check_and_open_scalar(whandle, m, H5T_FLOAT);
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'combine_embeddings'");
    }

    try {
        auto rhandle = utils::check_and_open_group(chandle, "results");

        if (modalities.size() > 1) {
            std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(total_dims) };
            utils::check_and_open_dataset(rhandle, "combined", H5T_FLOAT, pdims);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'combine_embeddings'");
    }

    return;
}

}

}

#endif
