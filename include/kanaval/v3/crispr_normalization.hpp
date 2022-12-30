#ifndef KANAVAL_CRISPR_NORMALIZATION_V3_HPP
#define KANAVAL_CRISPR_NORMALIZATION_V3_HPP

#include "H5Cpp.h"
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

inline void validate_crispr_normalization(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "crispr_normalization");

    try {
        utils::check_and_open_group(nhandle, "parameters");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'crispr_normalization'");
    }

    try {
        utils::check_and_open_group(nhandle, "results");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'crispr_normalization'");
    }

    return;
}

}

}

#endif
