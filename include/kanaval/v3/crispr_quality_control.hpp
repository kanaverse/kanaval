#ifndef KANAVAL_CRISPR_QUALITY_CONTROL_V3_HPP
#define KANAVAL_CRISPR_QUALITY_CONTROL_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "quality_control.hpp"

namespace kanaval {

namespace v3 {

inline int validate_crispr_quality_control(const H5::H5File& handle, int num_cells, int num_blocks, bool crispr_available, int version) {
    auto xhandle = utils::check_and_open_group(handle, "crispr_quality_control");

    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        auto nmads = utils::load_float_scalar<>(phandle, "nmads");
        if (nmads < 0) {
            throw std::runtime_error("number of MADs in 'nmads' should be non-negative");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'crispr_quality_control'");
    }

    int remaining;
    try {
        remaining = quality_control::validate_results(
            xhandle, 
            num_cells, 
            num_blocks, 
            { { "sums", H5T_FLOAT }, { "detected", H5T_INTEGER }, { "max_proportion", H5T_FLOAT }, { "max_index", H5T_INTEGER }},
            { "max_count" },
            crispr_available
        );
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'crispr_quality_control'");
    }

    return remaining;
}

}

}

#endif
