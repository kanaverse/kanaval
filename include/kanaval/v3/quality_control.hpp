#ifndef KANAVAL_QUALITY_CONTROL_V3_HPP
#define KANAVAL_QUALITY_CONTROL_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

namespace quality_control {

template<class Object>
int check_discard_vector(const Object& rhandle, size_t num_cells) {
    int remaining = 0;

    try {
        std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
        auto dihandle = utils::check_and_open_dataset(rhandle, "discards", H5T_INTEGER, dims);

        std::vector<int> discards(num_cells);
        dihandle.read(discards.data(), H5::PredType::NATIVE_INT);
        for (auto d : discards) {
            remaining += (d == 0);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve discard information from 'results'");
    }

    return remaining;
}

inline int validate_results(
    const H5::Group& handle, 
    int num_cells, 
    int num_blocks, 
    const std::vector<std::pair<std::string, H5T_class_t> >& metrics,
    const std::vector<std::string>& thresholds,
    bool in_use)
{
    auto rhandle = utils::check_and_open_group(handle, "results");
    int remaining = -1;

    if (in_use) {
        try {
            auto mhandle = utils::check_and_open_group(rhandle, "metrics");
            std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
            for (const auto& m : metrics) {
                utils::check_and_open_dataset(mhandle, m.first, m.second, dims);
            }
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve metrics from 'results'");
        }

        try {
            auto thandle = utils::check_and_open_group(rhandle, "thresholds");
            std::vector<size_t> dims{ static_cast<size_t>(num_blocks) };
            for (const auto& t : thresholds) {
                utils::check_and_open_dataset(thandle, t, H5T_FLOAT, dims);
            }
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve thresholds from 'results'");
        }

        remaining = check_discard_vector(rhandle, num_cells);
    }

    return remaining;
}

}

}

}

#endif
