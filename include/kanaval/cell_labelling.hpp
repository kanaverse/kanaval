#ifndef KANAVAL_CELL_LABELLING_HPP
#define KANAVAL_CELL_LABELLING_HPP

#include "H5Cpp.h"
#include <vector>
#include <unordered_set>
#include "utils.hpp"

/**
 * @file cell_labelling.hpp
 *
 * @brief Validate cell labelling contents.
 */

namespace kanaval {

namespace cell_labelling {

/**
 * @cond
 */
inline std::unordered_set<std::string> validate_parameters(const H5::Group& handle) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    
    std::unordered_set<std::string> output;
    for (size_t i = 0; i < 2; ++i) {
        std::string species = (i ? "mouse_references" : "human_references");
        auto refs = utils::load_string_vector(phandle, species);
        for (const auto& r : refs) {
            if (output.find(r) != output.end()) {
                throw std::runtime_error("duplicated reference '" + r + "' in '" + species + "'");
            }
            output.insert(r);
        }
    }

    return output;
}

inline void validate_results(const H5::Group& handle, const std::unordered_set<std::string>& available, int num_clusters) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    auto perhandle = utils::check_and_open_group(rhandle, "per_reference");

    size_t nchilds = perhandle.getNumObjs();
    std::vector<size_t> dims { static_cast<size_t>(num_clusters) };
    for (size_t a = 0; a < nchilds; ++a) {
        std::string name = perhandle.getObjnameByIdx(a);
        if (available.find(name) == available.end()) {
            throw std::runtime_error("reference '" + name + "' in 'results/per_reference' not listed in the parameters");
        }
        utils::check_and_open_dataset(perhandle, name, H5T_STRING, dims);
    }

    if (nchilds > 1) {
        auto integrated = utils::load_string_vector(rhandle, "integrated");
        if (integrated.size() != num_clusters) {
            throw std::runtime_error("'integrated' should have length equal to the number of clusters"); 
        }
        for (const auto& i : integrated) {
            if (available.find(i) == available.end()) {
                throw std::runtime_error("reference '" + i + "' not listed in the parameters");
            }
        }
    }

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the cell labelling step, see [here](@ref details-cell_labelling) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_clusters Number of clusters to be labelled.
 * 
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::H5File& handle, int num_clusters) {
    auto nhandle = utils::check_and_open_group(handle, "cell_labelling");

    std::unordered_set<std::string> refs;
    try {
        refs = validate_parameters(nhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'cell_labelling'");
    }

    try {
        validate_results(nhandle, refs, num_clusters);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'cell_labelling'");
    }

    return;
}

}

}

#endif
