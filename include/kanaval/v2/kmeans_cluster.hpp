#ifndef KANAVAL_KMEANS_CLUSTER_V2_HPP
#define KANAVAL_KMEANS_CLUSTER_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline int validate_kmeans_cluster(const H5::H5File& handle, int num_cells, bool in_use = true) {
    auto nhandle = utils::check_and_open_group(handle, "kmeans_cluster");

    int k;
    try {
        auto phandle = utils::check_and_open_group(nhandle, "parameters");
        k = utils::load_integer_scalar<>(phandle, "k");
        if (k <= 0) {
            throw std::runtime_error("number of clusters 'k' must be positive");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'kmeans_cluster'");
    }

    int nclusters = 0;
    try {
        auto phandle = utils::check_and_open_group(nhandle, "results");

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
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'kmeans_cluster'");
    }

    return nclusters;
}

}

}

#endif
