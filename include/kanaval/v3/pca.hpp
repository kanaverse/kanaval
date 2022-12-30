#ifndef KANAVAL_PCA_V3_HPP
#define KANAVAL_PCA_V3_HPP

#include "H5Cpp.h"
#include <vector>
#include <string>
#include <stdexcept>
#include "../utils.hpp"

namespace kanaval {

namespace v3 {

namespace pca {

inline int check_parameters(const H5::Group& phandle) {
    auto block_method = utils::load_string(phandle, "block_method");
    std::vector<std::string> options{ "none", "regress", "weight" };
    bool found = false;
    for (const auto& o : options) {
        if (o == block_method) {
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("unrecognized value '" + block_method + "' for the 'block_method'");
    }

    auto npcs = utils::load_integer_scalar<>(phandle, "num_pcs");
    if (npcs <= 0) {
        throw std::runtime_error("number of PCs must be positive in 'num_pcs'");
    }
    return npcs;
}

template<class Object>
size_t check_pca_contents(const Object& rhandle, int max_pcs, int num_cells) {
    auto vhandle = utils::check_and_open_dataset(rhandle, "var_exp", H5T_FLOAT);
    auto dspace = vhandle.getSpace();
    if (dspace.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("'var_exp' dataset does not have the expected dimensions");
    }

    hsize_t observed;
    dspace.getSimpleExtentDims(&observed);
    if (static_cast<int>(observed) > max_pcs) {
        throw std::runtime_error("length of 'var_exp' dataset exceeds the requested number of PCs");
    }
    
    utils::check_and_open_dataset(rhandle, "pcs", H5T_FLOAT, { static_cast<size_t>(num_cells), static_cast<size_t>(observed) });
    return observed;
}

}

}

}

#endif
