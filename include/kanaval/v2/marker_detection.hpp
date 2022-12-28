#ifndef KANAVAL_MARKER_DETECTION_V2_HPP
#define KANAVAL_MARKER_DETECTION_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "misc.hpp"

namespace kanaval {

namespace v2 {

namespace marker_detection {

inline void validate_markers(const H5::Group& chandle, int num_features, int num_clusters, std::string parent) {
    if (chandle.getNumObjs() != num_clusters) {
        throw std::runtime_error("number of groups in '" + parent + "' is not consistent with the expected number of clusters");
    }

    for (int i = 0; i < num_clusters; ++i) {
        try {
            auto ihandle = utils::check_and_open_group(chandle, std::to_string(i));
            std::vector<size_t> dims{ static_cast<size_t>(num_features) };
            utils::check_and_open_dataset(ihandle, "means", H5T_FLOAT, dims);
            utils::check_and_open_dataset(ihandle, "detected", H5T_FLOAT, dims);

            for (const auto& eff : markers::effects) {
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
            throw utils::combine_errors(e, "failed to retrieve statistics for cluster " + std::to_string(i) + " in '" + parent + "'");
        }
    }
    return;
}

}

inline void validate_marker_detection(const H5::Group& handle, int num_clusters, const std::vector<std::string>& modalities, const std::vector<int>& num_features, int version) {
    auto xhandle = utils::check_and_open_group(handle, "marker_detection");

    try {
        utils::check_and_open_group(xhandle, "parameters");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'marker_detection'");
    }

    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");
        if (version >= 2000000) {
            auto chandle = utils::check_and_open_group(rhandle, "per_cluster");
            for (size_t m = 0; m < modalities.size(); ++m) {
                auto mohandle = utils::check_and_open_group(chandle, modalities[m]);
                marker_detection::validate_markers(mohandle, num_features[m], num_clusters, "per_cluster/" + modalities[m]);
            }
        } else {
            auto chandle = utils::check_and_open_group(rhandle, "clusters");
            marker_detection::validate_markers(chandle, num_features[0], num_clusters, "clusters");
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'marker_detection'");
    }

    return;
}

}

}

#endif
