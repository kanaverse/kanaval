#ifndef KANAVAL_PCA_V2_HPP
#define KANAVAL_PCA_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "misc.hpp"

namespace kanaval {

namespace v2 {

inline int validate_pca(const H5::H5File& handle, int num_cells, int version = 1001000) {
    auto xhandle = utils::check_and_open_group(handle, "pca");

    int npcs;
    std::string block_method;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");

        auto nhvgs = utils::load_integer_scalar<>(phandle, "num_hvgs");
        if (nhvgs <= 0) {
            throw std::runtime_error("number of HVGs must be positive in 'num_hvgs'");
        }

        npcs = utils::load_integer_scalar<>(phandle, "num_pcs");
        if (npcs <= 0) {
            throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
        }

        if (version >= 1001000) {
            block_method = utils::load_string(phandle, "block_method");
            pca::check_block_method(block_method, version);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'pca'");
    }

    int obs_pcs;
    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");
        obs_pcs = pca::check_pca_contents(rhandle, npcs, num_cells);

        if (version >= 1001000 && version < 2000000) {
            if (block_method == "mnn") {
                utils::check_and_open_dataset(rhandle, "corrected", H5T_FLOAT, { static_cast<size_t>(num_cells), static_cast<size_t>(obs_pcs) });
            }
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'pca'");
    }

    return obs_pcs;
}

}

}

#endif
