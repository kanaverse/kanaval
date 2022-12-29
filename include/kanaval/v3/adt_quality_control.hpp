#ifndef KANAVAL_ADT_QUALITY_CONTROL_V3_HPP
#define KANAVAL_ADT_QUALITY_CONTROL_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "quality_control.hpp"

namespace kanaval {
    
namespace v3 {

inline int validate_adt_quality_control(const H5::H5File& handle, int num_cells, int num_blocks, bool adt_in_use, int version) {
    auto xhandle = utils::check_and_open_group(handle, "adt_quality_control");

    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        utils::check_and_open_scalar(phandle, "igg_prefix", H5T_STRING);

        auto nmads = utils::load_float_scalar<>(phandle, "nmads");
        if (nmads < 0) {
            throw std::runtime_error("number of MADs in 'nmads' should be non-negative");
        }

        auto mindrop = utils::load_float_scalar<>(phandle, "min_detected_drop");
        if (mindrop < 0 || mindrop >= 1) {
            throw std::runtime_error("minimum detected drop should lie in [0, 1)");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_quality_control'");
    }

    int remaining;
    try {
        remaining = quality_control::validate_results(
            xhandle, 
            num_cells, 
            num_blocks, 
            { { "sums", H5T_FLOAT }, { "detected", H5T_INTEGER }, { "igg_total", H5T_FLOAT }},
            { "detected", "igg_total" },
            adt_in_use
        );
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_quality_control'");
    }

    return remaining;
}

}

}

#endif
