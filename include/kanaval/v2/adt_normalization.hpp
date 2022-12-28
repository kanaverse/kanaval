#ifndef KANAVAL_ADT_NORMALIZATION_V2_HPP
#define KANAVAL_ADT_NORMALIZATION_V2_HPP

#include "H5Cpp.h"
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_adt_normalization(const H5::H5File& handle, size_t num_cells, bool adt_in_use, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return;
    }

    auto nhandle = utils::check_and_open_group(handle, "adt_normalization");

    try {
        auto phandle = utils::check_and_open_group(nhandle, "parameters");

        auto npcs = utils::load_integer_scalar(phandle, "num_pcs");
        if (npcs < 0){ 
            throw std::runtime_error("number of PCs used for ADT normalization should be positive");
        }

        auto nclust = utils::load_integer_scalar(phandle, "num_clusters");
        if (nclust < 0){ 
            throw std::runtime_error("number of clusters used for ADT normalization should be positive");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_normalization'");
    }

    try {
        auto rhandle = utils::check_and_open_group(nhandle, "results");
        if (adt_in_use) {
            utils::check_and_open_dataset(rhandle, "size_factors", H5T_FLOAT, { num_cells });
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_normalization'");
    }

    return;
}

}

}

#endif
