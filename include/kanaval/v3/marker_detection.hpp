#ifndef KANAVAL_MARKER_DETECTION_V3_HPP
#define KANAVAL_MARKER_DETECTION_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "../v2/misc.hpp"

namespace kanaval {

namespace v3 {

inline void validate_marker_detection(const H5::Group& handle, int num_clusters, const std::vector<std::string>& modalities, const std::vector<int>& num_features, int version) {
    auto xhandle = utils::check_and_open_group(handle, "marker_detection");

    // Checking the parameters.
    bool has_auc;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        has_auc = utils::load_integer_scalar<>(phandle, "compute_auc");

        auto lfc_threshold = utils::load_float_scalar<>(phandle, "lfc_threshold");
        if (lfc_threshold < 0) {
            throw std::runtime_error("'lfc_threshold' must be non-negative");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'marker_detection'");
    }

    // Checking the results.
    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");
        auto chandle = utils::check_and_open_group(rhandle, "per_cluster");
        for (size_t m = 0; m < modalities.size(); ++m) {
            auto mohandle = utils::check_and_open_group(chandle, modalities[m]);
            if (mohandle.getNumObjs() != num_clusters) {
                throw std::runtime_error("number of groups in 'per_cluster/" + modalities[m] + "' is not consistent with the expected number of clusters");
            }

            std::vector<size_t> dims{ static_cast<size_t>(num_features[m]) };
            for (int i = 0; i < num_clusters; ++i) {
                try {
                    auto ihandle = utils::check_and_open_group(mohandle, std::to_string(i));
                    utils::check_and_open_dataset(ihandle, "means", H5T_FLOAT, dims);
                    utils::check_and_open_dataset(ihandle, "detected", H5T_FLOAT, dims);

                    for (const auto& eff : v2::markers::effects) {
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
                    throw utils::combine_errors(e, "failed to retrieve statistics for cluster " + std::to_string(i) + " in 'per_cluster/" + modalities[m] + "'");
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
