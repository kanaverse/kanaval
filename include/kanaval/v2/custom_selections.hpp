#ifndef KANAVAL_CUSTOM_SELECTIONS_V2_HPP
#define KANAVAL_CUSTOM_SELECTIONS_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "misc.hpp"

namespace kanaval {

namespace v2 {

namespace custom_selections {

inline void validate_custom_markers(const H5::Group& shandle, int num_features) {
    std::vector<size_t> dims{ static_cast<size_t>(num_features) };
    utils::check_and_open_dataset(shandle, "means", H5T_FLOAT, dims);
    utils::check_and_open_dataset(shandle, "detected", H5T_FLOAT, dims);
    for (const auto& eff : markers::effects) {
        utils::check_and_open_dataset(shandle, eff, H5T_FLOAT, dims);
    }
}

}

inline void validate_custom_selections(const H5::Group& handle, int num_cells, const std::vector<std::string>& modalities, const std::vector<int>& num_features, int version) {
    auto cshandle = utils::check_and_open_group(handle, "custom_selections");

    std::vector<std::string> selections;
    try {
        auto phandle = utils::check_and_open_group(cshandle, "parameters");
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

            if (version >= 2001000) {
                if (!utils::is_unique_and_sorted(involved)) {
                    throw std::runtime_error("indices should be sorted and unique for selection '" + selections.back() + "'");
                }
            } else {
                std::sort(involved.begin(), involved.end());
                if (!utils::is_unique_and_sorted(involved)) {
                    throw std::runtime_error("indices should be unique for selection '" + selections.back() + "'");
                }
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'custom_selections'");
    }

    try {
        auto rhandle = utils::check_and_open_group(cshandle, "results");

        if (version >= 2000000) {
            auto mhandle = utils::check_and_open_group(rhandle, "per_selection");
            if (mhandle.getNumObjs() != selections.size()) {
                throw std::runtime_error("number of groups in 'per_selection' is not consistent with the expected number of selections");
            }

            for (const auto& s : selections) {
                try {
                    auto shandle = utils::check_and_open_group(mhandle, s);
                    for (size_t a = 0; a < modalities.size(); ++a) {
                        try {
                            auto ahandle = utils::check_and_open_group(shandle, modalities[a]);
                            custom_selections::validate_custom_markers(ahandle, num_features[a]);
                        } catch (std::exception& e) {
                            throw utils::combine_errors(e, "failed to retrieve statistics for modality '" + modalities[a] + "'");
                        }
                    }
                } catch (std::exception& e) {
                    throw utils::combine_errors(e, "failed to retrieve statistics for selection '" + s + "' in 'results/per_selection'");
                }
            }

        } else {
            auto mhandle = utils::check_and_open_group(rhandle, "markers");
            if (mhandle.getNumObjs() != selections.size()) {
                throw std::runtime_error("number of groups in 'markers' is not consistent with the expected number of selections");
            }

            for (const auto& s : selections) {
                try {
                    auto shandle = utils::check_and_open_group(mhandle, s);
                    custom_selections::validate_custom_markers(shandle, num_features[0]);
                } catch (std::exception& e) {
                    throw utils::combine_errors(e, "failed to retrieve statistics for selection '" + s + "' in 'results/markers'");
                }
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
