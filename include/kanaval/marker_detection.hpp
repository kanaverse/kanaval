#ifndef KANAVAL_MARKER_DETECTION_HPP
#define KANAVAL_MARKER_DETECTION_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

/**
 * @file marker_detection.hpp
 *
 * @brief Validate marker detection contents.
 */

namespace kanaval {

namespace marker_detection {

/**
 * @cond
 */
inline void validate_parameters(const H5::Group& handle) {
    utils::check_and_open_group(handle, "parameters");
}

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
/**
 * @endcond
 */

/**
 * Check contents for the marker detection step, see [here](@ref details-marker_detection) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_clusters Number of clusters produced by previous steps.
 * @param modalities Modalities available in the dataset, should be some combination of `"RNA"` or `"ADT"`.
 * If `version < 2000000`, this is ignored and an RNA modality is always assumed.
 * @param num_features Number of features for each modality in `modalities`.
 * If `version < 2000000`, only the first value is used and is assumed to refer to the number of genes for the RNA modality.
 * @param modalities Available modalities in the dataset.
 * Ignored for `version < 2000000`.
 * @param version Version of the format.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::Group& handle, int num_clusters, const std::vector<std::string>& modalities, const std::vector<int>& num_features, int version) {
    auto mhandle = utils::check_and_open_group(handle, "marker_detection");

    try {
        validate_parameters(mhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'marker_detection'");
    }

    try {
        auto rhandle = utils::check_and_open_group(mhandle, "results");
        if (version >= 2000000) {
            auto chandle = utils::check_and_open_group(rhandle, "per_cluster");
            for (size_t m = 0; m < modalities.size(); ++m) {
                auto mohandle = utils::check_and_open_group(chandle, modalities[m]);
                validate_markers(mohandle, num_features[m], num_clusters, "per_cluster/" + modalities[m]);
            }
        } else {
            auto chandle = utils::check_and_open_group(rhandle, "clusters");
            validate_markers(chandle, num_features[0], num_clusters, "clusters");
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'marker_detection'");
    }

    return;
}

}

}

#endif
