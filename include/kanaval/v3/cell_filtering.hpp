#ifndef KANAVAL_CELL_FILTERING_V3_HPP
#define KANAVAL_CELL_FILTERING_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include <unordered_map>
#include "../utils.hpp"
#include "quality_control.hpp"

namespace kanaval {

namespace v3 {

inline int validate_cell_filtering(const H5::H5File& handle, int num_cells, const std::unordered_map<std::string, int>& modalities, int version) {
    auto qhandle = utils::check_and_open_group(handle, "cell_filtering");

    // Checking parameters.
    int modalities_for_filtering = 0;
    int last_found = 0;
    try {
        auto phandle = utils::check_and_open_group(qhandle, "parameters");

        auto find_modality = [&](const std::string& flag, const std::string& target) -> void {
            if (utils::load_integer_scalar(phandle, flag)) {
                auto it = modalities.find(target);
                if (it != modalities.end()) {
                    ++modalities_for_filtering;
                    last_found = it->second;
                }
            }
        };

        find_modality("use_rna", "RNA");
        find_modality("use_adt", "ADT");
        find_modality("use_crispr", "CRISPR");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'cell_filtering'");
    }

    // Checking if there's a discard vector.
    int remaining;
    try {
        auto rhandle = utils::check_and_open_group(qhandle, "results");
        if (modalities_for_filtering > 1) {
            remaining = quality_control::check_discard_vector(rhandle, num_cells);
        } else if (modalities_for_filtering == 1) {
            remaining = last_found;
        } else {
            remaining = num_cells;
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'cell_filtering'");
    }

    return remaining;
}

}

}

#endif
