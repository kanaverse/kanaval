#ifndef KANAVAL_MISC_V2_HPP
#define KANAVAL_MISC_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

namespace pca {

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

inline void check_block_method(const std::string& method, int version) {
    std::vector<std::string> options{ "none", "regress" };
    if (version < 2000000) {
        options.push_back("mnn");
    } else {
        options.push_back("weight");
    }

    bool found = false;
    for (const auto& o : options) {
        if (o == method) {
            found = true;
            break;
        }
    }

    if (!found) {
        throw std::runtime_error("unrecognized value '" + method + "' for the 'block_method'");
    }
}

}

namespace quality_control {

template<class Object>
int check_discard_vector(const Object& rhandle, size_t num_cells, bool skip) {
    int remaining = 0;

    if (!skip || rhandle.exists("discards")) {
        try {
            std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
            auto dihandle = utils::check_and_open_dataset(rhandle, "discards", H5T_INTEGER, dims);

            if (!skip) {
                std::vector<int> discards(num_cells);
                dihandle.read(discards.data(), H5::PredType::NATIVE_INT);
                for (auto d : discards) {
                    remaining += (d == 0);
                }
            } else {
                remaining = num_cells;
            }

        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve discard information from 'results'");
        }
    } else {
        remaining = num_cells;
    }

    return remaining;
}

inline int validate_results(
    const H5::Group& handle, 
    int num_cells, 
    int num_samples, 
    const std::vector<std::pair<std::string, H5T_class_t> >& metrics,
    const std::vector<std::string>& thresholds,
    bool in_use,
    bool skip)
{
    auto rhandle = utils::check_and_open_group(handle, "results");
    int remaining = -1;

    if (in_use) {
        if (!skip || rhandle.exists("metrics")) {
            try {
                auto mhandle = utils::check_and_open_group(rhandle, "metrics");
                std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
                for (const auto& m : metrics) {
                    utils::check_and_open_dataset(mhandle, m.first, m.second, dims);
                }
            } catch (std::exception& e) {
                throw utils::combine_errors(e, "failed to retrieve metrics from 'results'");
            }
        }

        if (!skip || rhandle.exists("thresholds")) {
            try {
                auto thandle = utils::check_and_open_group(rhandle, "thresholds");
                std::vector<size_t> dims{ static_cast<size_t>(num_samples) };
                for (const auto& t : thresholds) {
                    utils::check_and_open_dataset(thandle, t, H5T_FLOAT, dims);
                }
            } catch (std::exception& e) {
                throw utils::combine_errors(e, "failed to retrieve thresholds from 'results'");
            }
        }

        remaining = check_discard_vector(rhandle, num_cells, skip);
    }

    return remaining;
}

}

namespace markers {
    
inline const std::vector<std::string> effects { "lfc", "delta_detected", "cohen", "auc" };

}

}

}

#endif
