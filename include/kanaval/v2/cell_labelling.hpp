#ifndef KANAVAL_CELL_LABELLING_V2_HPP
#define KANAVAL_CELL_LABELLING_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include <unordered_set>
#include "../utils.hpp"

namespace kanaval {

namespace v2 {

inline void validate_cell_labelling(const H5::H5File& handle, int num_clusters) {
    auto nhandle = utils::check_and_open_group(handle, "cell_labelling");

    std::unordered_set<std::string> refs;
    try {
        auto phandle = utils::check_and_open_group(nhandle, "parameters");
        
        for (size_t i = 0; i < 2; ++i) {
            std::string species = (i ? "mouse_references" : "human_references");
            auto species_refs = utils::load_string_vector(phandle, species);
            for (const auto& r : species_refs) {
                if (refs.find(r) != refs.end()) {
                    throw std::runtime_error("duplicated reference '" + r + "' in '" + species + "'");
                }
                refs.insert(r);
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'cell_labelling'");
    }

    try {
        auto rhandle = utils::check_and_open_group(nhandle, "results");
        auto perhandle = utils::check_and_open_group(rhandle, "per_reference");

        size_t nchilds = perhandle.getNumObjs();
        std::vector<size_t> dims { static_cast<size_t>(num_clusters) };
        for (size_t a = 0; a < nchilds; ++a) {
            std::string name = perhandle.getObjnameByIdx(a);
            if (refs.find(name) == refs.end()) {
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
                if (refs.find(i) == refs.end()) {
                    throw std::runtime_error("reference '" + i + "' not listed in the parameters");
                }
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'cell_labelling'");
    }

    return;
}

}

}

#endif
