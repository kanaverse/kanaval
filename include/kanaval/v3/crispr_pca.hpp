#ifndef KANAVAL_CRISPR_PCA_V3_HPP
#define KANAVAL_CRISPR_PCA_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include <stdexcept>
#include <string>
#include "../utils.hpp"
#include "pca.hpp"

namespace kanaval {

namespace v3 {

inline int validate_crispr_pca(const H5::H5File& handle, int num_cells, int version) {
    auto xhandle = utils::check_and_open_group(handle, "crispr_pca");

    int npcs;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        npcs = pca::check_parameters(phandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'crispr_pca'");
    }

    int obs_pcs;
    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");
        obs_pcs = pca::check_pca_contents(rhandle, npcs, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'crispr_pca'");
    }

    return obs_pcs;
}

}

}

#endif
