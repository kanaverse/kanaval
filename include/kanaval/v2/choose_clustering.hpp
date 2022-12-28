#ifndef KANAVAL_CHOOSE_CLUSTERING_V2_HPP
#define KANAVAL_CHOOSE_CLUSTERING_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline std::string validate_choose_clustering(const H5::H5File& handle) {
    auto nhandle = utils::check_and_open_group(handle, "choose_clustering");

    std::string method;
    try {
        auto phandle = utils::check_and_open_group(nhandle, "parameters");
        method = utils::load_string(phandle, "method");
        if (method != "kmeans" && method != "snn_graph") {
            throw std::runtime_error("'method' should be either 'kmeans' or 'snn_graph'");
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'choose_clustering'");
    }

    try {
        utils::check_and_open_group(nhandle, "results");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'choose_clustering'");
    }

    return method;
}

}

}

#endif
