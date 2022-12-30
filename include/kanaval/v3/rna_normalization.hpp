#ifndef KANAVAL_RNA_NORMALIZATION_V3_HPP
#define KANAVAL_RNA_NORMALIZATION_V3_HPP

#include "H5Cpp.h"
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

inline void validate_rna_normalization(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "rna_normalization");

    try {
        utils::check_and_open_group(nhandle, "parameters");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'rna_normalization'");
    }

    try {
        utils::check_and_open_group(nhandle, "results");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'rna_normalization'");
    }

    return;
}

}

}

#endif
