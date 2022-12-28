#ifndef KANAVAL_SNN_GRAPH_CLUSTER_V2_HPP
#define KANAVAL_SNN_GRAPH_CLUSTER_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline int validate_snn_graph_cluster(const H5::H5File& handle, int num_cells, bool in_use = true) {
    auto xhandle = utils::check_and_open_group(handle, "snn_graph_cluster");

    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");

        auto k = utils::load_integer_scalar<>(phandle, "k");
        if (k <= 0) {
            throw std::runtime_error("number of neighbors 'k' must be positive");
        }

        auto scheme = utils::load_string(phandle, "scheme");
        if (scheme != "rank" && scheme != "jaccard" && scheme != "number") {
            throw std::runtime_error("'scheme' must be one of 'rank', 'jaccard' or 'number'");
        }

        auto res = utils::load_float_scalar<>(phandle, "resolution");
        if (res < 0) {
            throw std::runtime_error("'resolution' must be non-negative");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'snn_graph_cluster'");
    }

    int nclusters = 0;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "results");

        if (phandle.exists("clusters") || in_use) {
            std::vector<size_t> dim { static_cast<size_t>(num_cells) };
            auto clushandle = utils::check_and_open_dataset(phandle, "clusters", H5T_INTEGER, dim);
            std::vector<int> clusters(num_cells);
            clushandle.read(clusters.data(), H5::PredType::NATIVE_INT);

            if (num_cells) {
                auto minned = *std::min_element(clusters.begin(), clusters.end());
                if (minned < 0) {
                    throw std::runtime_error("entries in 'clusters' should be non-negative");
                }

                auto maxed = *std::max_element(clusters.begin(), clusters.end());
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
        throw utils::combine_errors(e, "failed to retrieve results from 'snn_graph_cluster'");
    }

    return nclusters;
}

}

}

#endif
