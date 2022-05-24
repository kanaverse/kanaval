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

inline void validate_markers(const H5::Group& chandle, int num_genes, int num_clusters, std::string parent) {
    if (chandle.getNumObjs() != num_clusters) {
        throw std::runtime_error("number of groups in '" + parent + "' is not consistent with the expected number of clusters");
    }

    for (int i = 0; i < num_clusters; ++i) {
        try {
            auto ihandle = utils::check_and_open_group(chandle, std::to_string(i));
            std::vector<size_t> dims{ static_cast<size_t>(num_genes) };
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
 * Check contents for the marker detection step.
 * Contents are stored inside an `marker_detection` HDF5 group at the root of the file.
 * The `marker_detection` group itself contains the `parameters` and `results` subgroups.
 *
 * <HR>
 * `parameters` should be empty.
 * 
 * <HR>
 * `results` should contain `per_cluster`, a group containing the marker results for each cluster.
 * Each child of `per_cluster` is named after a cluster index from 0 to `num_clusters - 1`, and is itself a group containing children named according to `modalities`.
 * Each modality-specific child is yet another group containing the statistics for that modality:
 *
 * - `means`: a float dataset of length equal to the number of genes, containing the mean expression of each gene in the current cluster.
 * - `detected`: a float dataset of length equal to the number of genes, containing the proportion of cells with detected expression of each gene in the current cluster.
 * - `lfc`: an group containing statistics for the log-fold changes from all pairwise comparisons involving the current cluster.
 *   This contains:
 *   - `min`: a float dataset of length equal to the number of genes, containing the minimum log-fold change across all pairwise comparisons for each gene.
 *   - `mean`: a float dataset of length equal to the number of genes, containing the mean log-fold change across all pairwise comparisons for each gene.
 *   - `min_rank`: a float dataset of length equal to the number of genes, containing the minimum rank of the log-fold changes across all pairwise comparisons for each gene.
 * - `delta_detected`: same as `lfc`, but for the delta-detected (i.e., difference in the percentage of detected expression).
 * - `cohen`: same as `lfc`, but for Cohen's d.
 * - `auc`: same as `lfc`, but for the AUCs.
 *
 * <DIV style="color:blue">
 * <details>
 * <summary>For versions 1.0-1.2</summary>
 * `results` should contain:
 *
 * - `clusters`: a group representing an array of length equal to the number of clusters.
 *   Each child is another group that is named by the cluster index from 0 to `num_clusters - 1`, containing the marker details for that cluster.
 *   Each child group contains:
 *   - `means`: a float dataset of length equal to the number of genes, containing the mean expression of each gene in the current cluster.
 *   - `detected`: a float dataset of length equal to the number of genes, containing the proportion of cells with detected expression of each gene in the current cluster.
 *   - `lfc`: an group containing statistics for the log-fold changes from all pairwise comparisons involving the current cluster.
 *     This contains:
 *     - `min`: a float dataset of length equal to the number of genes, containing the minimum log-fold change across all pairwise comparisons for each gene.
 *     - `mean`: a float dataset of length equal to the number of genes, containing the mean log-fold change across all pairwise comparisons for each gene.
 *     - `min_rank`: a float dataset of length equal to the number of genes, containing the minimum rank of the log-fold changes across all pairwise comparisons for each gene.
 *   - `delta_detected`: same as `lfc`, but for the delta-detected (i.e., difference in the percentage of detected expression).
 *   - `cohen`: same as `lfc`, but for Cohen's d.
 *   - `auc`: same as `lfc`, but for the AUCs.
 * </details>
 * </DIV>
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_clusters Number of clusters produced by previous steps.
 * @param modalities Available modalities in the dataset.
 * Ignored for `version < 2000000`.
 * @param num_genes Number of features for each modality.
 * Should only contain a single value for `version < 2000000`.
 * @param version Version of the format.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::Group& handle, int num_clusters, const std::vector<std::string>& modalities, const std::vector<int>& num_genes, int version) {
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
                validate_markers(mohandle, num_genes[m], num_clusters, "per_cluster/" + modalities[m]);
            }
        } else {
            auto chandle = utils::check_and_open_group(rhandle, "clusters");
            validate_markers(chandle, num_genes[0], num_clusters, "clusters");
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'marker_detection'");
    }

    return;
}

}

}

#endif
