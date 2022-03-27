#ifndef KANAVAL_KMEANS_CLUSTER_HPP
#define KANAVAL_KMEANS_CLUSTER_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file kmeans_cluster.hpp
 *
 * @brief Validate k-means clustering contents.
 */

namespace kanaval {

namespace kmeans_cluster { 

/**
 * @cond
 */
inline int validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    auto k = utils::load_integer_scalar<>(phandle, "k");
    if (k <= 0) {
        throw std::runtime_error("number of clusters 'k' must be positive");
    }
    return k;
}

inline int validate_results(const H5::Group& handle, int k, int num_cells, bool in_use) {
    auto phandle = utils::check_and_open_group(handle, "results");

    int nclusters = 0;
    if (phandle.exists("clusters") || in_use) {
        std::vector<size_t> dim { static_cast<size_t>(num_cells) };
        auto clushandle = utils::check_and_open_dataset(phandle, "clusters", H5T_INTEGER, dim);
        std::vector<int> clusters(num_cells);
        clushandle.read(clusters.data(), H5::PredType::NATIVE_INT);

        if (num_cells) {
            auto minned = *std::min_element(clusters.begin(), clusters.end());
            auto maxed = *std::max_element(clusters.begin(), clusters.end());
            if (minned < 0 || maxed >= k) {
                throw std::runtime_error("entries in 'clusters' are out of range for the given 'k'");
            }

            nclusters = maxed + 1;
            std::vector<int> counts(nclusters);
            for (auto c : clusters) {
                ++counts[c];
            }
            for (auto c : counts) {
                if (c == 0) {
                    throw std::runtime_error("each cluster must be represented at least once in 'clusters'");
                }
            }

        }
    }

    return nclusters;
}
/**
 * @endcond
 */

/**
 * Check contents for the k-means clustering step.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset after any quality filtering is applied.
 * @param in_use Was k-means clustering used by downstream steps?
 *
 * @return The total number of clusters.
 * If the format is invalid, an error is raised.
 *
 * @description
 * `handle` should contain a `kmeans_cluster` group, itself containing the `parameters` and `results` subgroups.
 * 
 * `parameters` should contain:
 * 
 * - `k`: a scalar integer specifying the number of clusters to create.
 * 
 * If `in_use = true`, `results` should contain:
 * 
 * - `clusters`: an integer dataset of length equal to `num_cells`.
 *   This contains the cluster assignment for each cell, which should lie in $[0, k)$.
 *   The total number of clusters may be less than $k$, e.g., when there are too few cells.
 *   For $N$ clusters, there should be at least one occurrence of each integer in $[0, N)$.
 * 
 * If `in_use = false`, `clusters` may be absent.
 * If it is present, it should follow the constraints listed above.
 */
inline int validate(const H5::H5File& handle, int num_cells, bool in_use = true) {
    auto nhandle = utils::check_and_open_group(handle, "kmeans_cluster");

    int k;
    try {
        k = validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'kmeans_cluster'");
    }

    int nclusters;
    try {
        nclusters = validate_results(nhandle, k, num_cells, in_use);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'kmeans_cluster'");
    }

    return nclusters;
}

}

}

#endif
