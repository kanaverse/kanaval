#ifndef KANAVAL_CUSTOM_SELECTIONS_V3_HPP
#define KANAVAL_CUSTOM_SELECTIONS_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include <unordered_map>
#include "../utils.hpp"
#include "markers.hpp"

namespace kanaval {

namespace v3 {

inline void validate_custom_selections(const H5::Group& handle, int num_cells, const std::unordered_map<std::string, int>& modalities, int version) {
    auto cshandle = utils::check_and_open_group(handle, "custom_selections");

    // Checking the parameters.
    std::vector<std::string> selections;
    bool has_auc;
    try {
        auto phandle = utils::check_and_open_group(cshandle, "parameters");
        has_auc = markers::check_common_parameters(phandle);

        auto shandle = utils::check_and_open_group(phandle, "selections");
        for (hsize_t i = 0; i < shandle.getNumObjs(); ++i) {
            auto name = shandle.getObjnameByIdx(i);
            selections.push_back(name);

            auto involved = utils::load_integer_vector(shandle, name);
            for (auto i : involved) {
                if (i < 0 || i >= num_cells) {
                    throw std::runtime_error("indices out of range for selection '" + selections.back() + "'");
                }
            }

            if (!utils::is_unique_and_sorted(involved)) {
                throw std::runtime_error("indices should be sorted and unique for selection '" + selections.back() + "'");
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'custom_selections'");
    }

    // Checking the results.
    try {
        auto rhandle = utils::check_and_open_group(cshandle, "results");
        auto mhandle = utils::check_and_open_group(rhandle, "per_selection");
        if (mhandle.getNumObjs() != selections.size()) {
            throw std::runtime_error("number of groups in 'per_selection' is not consistent with the expected number of selections");
        }

        for (const auto& s : selections) {
            try {
                auto shandle = utils::check_and_open_group(mhandle, s);
                for (const auto& mod : modalities) {
                    std::vector<size_t> dims{ static_cast<size_t>(mod.second) };

                    try {
                        auto ahandle = utils::check_and_open_group(shandle, mod.first);
                        utils::check_and_open_dataset(ahandle, "means", H5T_FLOAT, dims);
                        utils::check_and_open_dataset(ahandle, "detected", H5T_FLOAT, dims);

                        for (const auto& eff : markers::effects) {
                            if (!has_auc && eff == "auc") {
                                continue;
                            }
                            utils::check_and_open_dataset(ahandle, eff, H5T_FLOAT, dims);
                        }
                    } catch (std::exception& e) {
                        throw utils::combine_errors(e, "failed to retrieve statistics for modality '" + mod.first + "'");
                    }
                }
            } catch (std::exception& e) {
                throw utils::combine_errors(e, "failed to retrieve statistics for selection '" + s + "' in 'results/per_selection'");
            }
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'custom_selections'");
    }

    return;
}

}

}

#endif
