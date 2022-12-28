#ifndef KANAVAL_ADT_PCA_V2_HPP
#define KANAVAL_ADT_PCA_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "misc.hpp"

namespace kanaval {

namespace v2 {

inline int validate_adt_pca(const H5::H5File& handle, int num_cells, bool adt_in_use, int version) {
    if (version < 2000000) {
        return -1;        
    }

    auto ahandle = utils::check_and_open_group(handle, "adt_pca");

    int num_pcs;
    try {
        auto phandle = utils::check_and_open_group(ahandle, "parameters");

        num_pcs = utils::load_integer_scalar<>(phandle, "num_pcs");
        if (num_pcs <= 0) {
            throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
        }

        std::string method = utils::load_string(phandle, "block_method");
        pca::check_block_method(method, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_pca'");
    }

    int obs_pcs = -1;
    try {
        auto rhandle = utils::check_and_open_group(ahandle, "results");
        if (adt_in_use) {
            obs_pcs = pca::check_pca_contents(rhandle, num_pcs, num_cells);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_pca'");
    }

    return obs_pcs;
}

}

}

#endif
