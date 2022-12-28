#ifndef KANAVAL_NEIGHBOR_INDEX_V2_HPP
#define KANAVAL_NEIGHBOR_INDEX_V2_HPP

#include "H5Cpp.h"
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_neighbor_index(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "neighbor_index");

    try {
        auto phandle = utils::check_and_open_group(nhandle, "parameters");
        utils::check_and_open_scalar(phandle, "approximate", H5T_INTEGER);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'neighbor_index'");
    }

    try {
        utils::check_and_open_group(nhandle, "results");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'neighbor_index'");
    }

    return;
}

}

}

#endif
