#ifndef KANAVAL_SNN_GRAPH_CLUSTER_V3_HPP
#define KANAVAL_SNN_GRAPH_CLUSTER_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

inline int validate_snn_graph_cluster(const H5::H5File& handle, int num_cells, bool in_use = true) {
    auto xhandle = utils::check_and_open_group(handle, "snn_graph_cluster");

    // Checking the parameters.
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

        auto algorithm = utils::load_string(phandle, "algorithm");
        if (algorithm != "multilevel" && algorithm != "walktrap" && algorithm != "leiden") {
            throw std::runtime_error("'algorithm' must be one of 'multilevel', 'walktrap' or 'leiden'");
        }

        auto mres = utils::load_float_scalar<>(phandle, "multilevel_resolution");
        if (mres < 0) {
            throw std::runtime_error("'multilevel_resolution' must be non-negative");
        }

        auto lres = utils::load_float_scalar<>(phandle, "leiden_resolution");
        if (lres < 0) {
            throw std::runtime_error("'leiden_resolution' must be non-negative");
        }

        auto wsteps = utils::load_integer_scalar<>(phandle, "walktrap_steps");
        if (wsteps < 0) {
            throw std::runtime_error("'walktrap_steps' must be non-negative");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'snn_graph_cluster'");
    }

    // Checking the results.
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
