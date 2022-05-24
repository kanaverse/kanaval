#ifndef KANAVAL_MISC_HPP
#define KANAVAL_MISC_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

namespace kanaval {

namespace pca {

template<class Object>
size_t check_pca_contents(const Object& rhandle, int max_pcs, int num_cells) {
    auto vhandle = utils::check_and_open_dataset(rhandle, "var_exp", H5T_FLOAT);
    auto dspace = vhandle.getSpace();
    if (dspace.getSimpleExtentNdims() != 0) {
        throw std::runtime_error("'var_exp' dataset does not have the expected dimensions");
    }

    hsize_t observed;
    dspace.getSimpleExtentDims(&observed);
    if (static_cast<int>(observed) <= max_pcs) {
        throw std::runtime_error("length of 'var_exp' dataset exceeds the requested number of PCs");
    }
    
    utils::check_and_open_dataset(rhandle, "pcs", H5T_FLOAT, { static_cast<size_t>(num_cells), static_cast<size_t>(observed) });
    return observed;
}

}

namespace quality_control {

template<class Object>
int check_discard_vector(const Object& rhandle, size_t num_cells) {
    int remaining = 0;
    try {
        std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
        auto dihandle = utils::check_and_open_dataset(rhandle, "discards", H5T_INTEGER, dims);
        std::vector<int> discards(num_cells);
        dihandle.read(discards.data(), H5::PredType::NATIVE_INT);
        for (auto d : discards) {
            remaining += (d == 0);
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve discard information from 'results'");
    }
    return remaining;
}

}

namespace markers {
    
inline const std::vector<std::string> effects { "lfc", "delta_detected", "cohen", "auc" };

}

}

#endif
