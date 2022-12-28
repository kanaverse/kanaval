#ifndef KANAVAL_NORMALIZATION_V2_HPP
#define KANAVAL_NORMALIZATION_V2_HPP

#include "H5Cpp.h"
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_normalization(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "normalization");

    try {
        utils::check_and_open_group(nhandle, "parameters");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'normalization'");
    }

    try {
        utils::check_and_open_group(nhandle, "results");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'normalization'");
    }

    return;
}

}

}

#endif
