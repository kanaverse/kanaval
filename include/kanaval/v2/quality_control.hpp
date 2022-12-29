#ifndef KANAVAL_QUALITY_CONTROL_V2_HPP
#define KANAVAL_QUALITY_CONTROL_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "misc.hpp"

namespace kanaval {

namespace v2 {

inline std::pair<bool, int> validate_quality_control(const H5::H5File& handle, int num_cells, int num_samples, int version) {
    auto qhandle = utils::check_and_open_group(handle, "quality_control");

    bool skip = false;
    try {
        auto phandle = utils::check_and_open_group(qhandle, "parameters");
        utils::check_and_open_scalar(phandle, "use_mito_default", H5T_INTEGER);
        utils::check_and_open_scalar(phandle, "mito_prefix", H5T_STRING);

        auto nmads = utils::load_float_scalar<>(phandle, "nmads");
        if (nmads < 0) {
            throw std::runtime_error("number of MADs in 'nmads' should be non-negative");
        }

        if (version >= 2001000) {
            skip = utils::load_integer_scalar(phandle, "skip");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'quality_control'");
    }

    int remaining;
    try {
        remaining = quality_control::validate_results(
            qhandle,
            num_cells,
            num_samples,
            { { "sums", H5T_FLOAT }, { "detected", H5T_INTEGER }, { "proportion", H5T_FLOAT } },
            { "sums", "detected", "proportion" },
            true, /* RNA is always assumed to be in use */
            skip
        );
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'quality_control'");
    }

    return std::make_pair(skip, remaining);
}

}

}

#endif
