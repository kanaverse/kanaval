#ifndef KANAVAL_MARKER_DETECTION_V3_HPP
#define KANAVAL_MARKER_DETECTION_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include <unordered_map>
#include "../utils.hpp"
#include "markers.hpp"

namespace kanaval {

namespace v3 {

inline void validate_marker_detection(const H5::Group& handle, int num_clusters, const std::unordered_map<std::string, int>& modalities, int version) {
    auto xhandle = utils::check_and_open_group(handle, "marker_detection");

    // Checking the parameters.
    bool has_auc;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        has_auc = markers::check_common_parameters(phandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'marker_detection'");
    }

    // Checking the results.
    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");
        auto chandle = utils::check_and_open_group(rhandle, "per_cluster");
        for (const auto& mod : modalities) {
            auto mohandle = utils::check_and_open_group(chandle, mod.first);
            if (mohandle.getNumObjs() != num_clusters) {
                throw std::runtime_error("number of groups in 'per_cluster/" + mod.first + "' is not consistent with the expected number of clusters");
            }

            std::vector<size_t> dims{ static_cast<size_t>(mod.second) };
            for (int i = 0; i < num_clusters; ++i) {
                try {
                    auto ihandle = utils::check_and_open_group(mohandle, std::to_string(i));
                    utils::check_and_open_dataset(ihandle, "means", H5T_FLOAT, dims);
                    utils::check_and_open_dataset(ihandle, "detected", H5T_FLOAT, dims);

                    for (const auto& eff : markers::effects) {
                        if (!has_auc && eff == "auc") {
                            continue;
                        }

                        try {
                            auto ehandle = utils::check_and_open_group(ihandle, eff);
                            utils::check_and_open_dataset(ehandle, "mean", H5T_FLOAT, dims);
                            utils::check_and_open_dataset(ehandle, "min", H5T_FLOAT, dims);
                            utils::check_and_open_dataset(ehandle, "min_rank", H5T_FLOAT, dims);
                        } catch (std::exception& e) {
                            throw utils::combine_errors(e, "failed to retrieve summary statistic for '" + eff + "'");
                        }
                    }
                } catch (std::exception& e) {
                    throw utils::combine_errors(e, "failed to retrieve statistics for cluster " + std::to_string(i) + " in 'per_cluster/" + mod.first + "'");
                }
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'marker_detection'");
    }

    return;
}

}

}

#endif
