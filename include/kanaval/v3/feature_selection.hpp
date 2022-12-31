#ifndef KANAVAL_FEATURE_SELECTION_V3_HPP
#define KANAVAL_FEATURE_SELECTION_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

inline void validate_feature_selection(const H5::Group& handle, int num_genes, bool rna_available, int version) {
    auto fhandle = utils::check_and_open_group(handle, "feature_selection");

    try {
        auto phandle = utils::check_and_open_group(fhandle, "parameters");
        auto span = utils::load_float_scalar<>(phandle, "span");
        if (span < 0 || span > 1) {
            throw std::runtime_error("LOWESS span should lie in [0, 1]");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'feature_selection'");
    }

    try {
        auto rhandle = utils::check_and_open_group(fhandle, "results");

        if (rna_available) {
            std::vector<size_t> dims{ static_cast<size_t>(num_genes) };
            utils::check_and_open_dataset(rhandle, "means", H5T_FLOAT, dims);
            utils::check_and_open_dataset(rhandle, "vars", H5T_FLOAT, dims);
            utils::check_and_open_dataset(rhandle, "fitted", H5T_FLOAT, dims);
            utils::check_and_open_dataset(rhandle, "resids", H5T_FLOAT, dims);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'feature_selection'");
    }

    return;
}

}

}

#endif
