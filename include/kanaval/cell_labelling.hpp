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

    for (const auto& a : available) {
        if (!perhandle.exists(a)) {
            throw std::runtime_error("reference '" + a + "' in parameters not listed in 'results/per_reference'");
        }
    }

    size_t nchilds = perhandle.getNumObjs();
    std::vector<size_t> dims { static_cast<size_t>(num_clusters) };
    for (size_t a = 0; a < nchilds; ++a) {
        std::string name = perhandle.getObjnameByIdx(a);
        if (available.find(name) == available.end()) {
            throw std::runtime_error("reference '" + name + "' in 'results/per_reference' not listed in the parameters");
        }
        utils::check_and_open_dataset(perhandle, name, H5T_STRING, dims);
    }

    if (available.size() > 1) {
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
 * Check contents for the cell labelling step.
 *
 * `handle` should contain a `cell_labelling` group, itself containing the `parameters` and `results` subgroups.
 *
 * `parameters` should contain:
 * 
 * - `human_references`: a string dataset defining the human reference datasets used for labelling.
 *   Each entry contains the name of a reference dataset, e.g., `"BlueprintEncode"`.
 * - `mouse_references`: a string dataset defining the mouse reference datasets used for labelling.
 *   Each entry contains the name of a reference dataset, e.g., `"ImmGen"`.
 * 
 * `results` should contain:
 * 
 * - `per_reference`: a group containing the label assignments for each cluster in each reference.
 *   Each child is named after its corresponding reference, and is a string dataset of length equal to the number of clusters.
 *   Entries of the dataset contain the assigned label for the corresponding cluster.
 * 
 * For multiple references of the relevant species, `results` will also contain:
 * 
 * - `integrated`: a string dataset of length equal to the number of clusters.
 *   This specifies the reference with the top-scoring label for each cluster, after integrating the results of all per-reference classifications.
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
