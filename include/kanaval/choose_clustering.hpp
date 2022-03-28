#ifndef KANAVAL_CHOOSE_CLUSTERING_HPP
#define KANAVAL_CHOOSE_CLUSTERING_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file choose_clustering.hpp
 *
 * @brief Validate the clustering choice.
 */

namespace kanaval {

namespace choose_clustering { 

/**
 * @cond
 */
inline std::string validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    auto method = utils::load_string(phandle, "method");
    if (method != "kmeans" && method != "snn_graph") {
        throw std::runtime_error("'method' should be either 'kmeans' or 'snn_graph'");
    }
    return method;
}

inline void validate_results(const H5::Group& handle) {
    utils::check_and_open_group(handle, "results");
    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the clustering choice.
 *
 * @param handle An open HDF5 file handle.
 *
 * @return The clustering method of choice. 
 * If the format is invalid, an error is raised.
 *
 * @details
 * `handle` should contain a `choose_clustering` group, itself containing the `parameters` and `results` subgroups.
 * 
 * `parameters` should contain:
 * 
 * - `method`: a scalar string specifying the clustering method to use.
 *   This is currently either `snn_graph` or `kmeans`.
 * 
 * `results` is empty.
 * 
 * Depending on the `method`, `snn_graph_cluster` or `kmeans_cluster` must have non-empty `results`.
 * Both may also be non-empty, in which case the appropriate clustering is chosen based on `method`.
 */
inline std::string validate(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "choose_clustering");

    std::string output;
    try {
        output = validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'choose_clustering'");
    }

    try {
        validate_results(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'choose_clustering'");
    }

    return output;
}

}

}

#endif
