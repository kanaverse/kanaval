#ifndef KANAVAL_COMBINE_EMBEDDINGS_V3_HPP
#define KANAVAL_COMBINE_EMBEDDINGS_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include <numeric>
#include <unordered_map>
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

inline int validate_combine_embeddings(const H5::H5File& handle, int num_cells, const std::unordered_map<std::string, int>& modalities, int version) {
    auto xhandle = utils::check_and_open_group(handle, "combine_embeddings");

    // Checking the parameters.
    int total_dims = 0;
    int num_modalities = 0;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        utils::check_and_open_scalar(phandle, "approximate", H5T_INTEGER);

        auto check_modality = [&](const std::string& weight, const std::string& target) -> void {
            if (utils::load_float_scalar(phandle, weight) > 0) {
                auto it = modalities.find(target);
                if (it != modalities.end()) {
                    ++num_modalities;
                    total_dims += it->second;
                }
            }
        };

        check_modality("rna_weight", "RNA");
        check_modality("adt_weight", "ADT");
        check_modality("crispr_weight", "CRISPR");

        if (num_modalities == 0) {
            throw std::runtime_error("could not find any available modality with non-zero weight");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'combine_embeddings'");
    }

    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");

        if (num_modalities > 1) {
            std::vector<size_t> pdims { static_cast<size_t>(num_cells), static_cast<size_t>(total_dims) };
            utils::check_and_open_dataset(rhandle, "combined", H5T_FLOAT, pdims);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'combine_embeddings'");
    }

    return total_dims;
}

}

}

#endif
