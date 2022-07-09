#ifndef KANAVAL_CUSTOM_SELECTIONS_HPP
#define KANAVAL_CUSTOM_SELECTIONS_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

/**
 * @file custom_selections.hpp
 *
 * @brief Validate custom selection contents.
 */

namespace kanaval {

/**
 * Validation for custom selections
 */
namespace custom_selections {

/**
 * @cond
 */
inline std::vector<std::string> validate_parameters(const H5::Group& handle, int num_cells) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    auto shandle = utils::check_and_open_group(phandle, "selections");

    std::vector<std::string> output;
    for (hsize_t i = 0; i < shandle.getNumObjs(); ++i) {
        auto name = shandle.getObjnameByIdx(i);
        output.push_back(name);

        auto involved = utils::load_integer_vector(shandle, name);
        for (auto i : involved) {
            if (i < 0 || i >= num_cells) {
                throw std::runtime_error("indices out of range for selection '" + output.back() + "'");
            }
        }

        if (!utils::is_unique_and_sorted(involved)) {
            throw std::runtime_error("indices should be sorted and unique for selection '" + output.back() + "'");
        }
    }

    return output;
}

inline void validate_custom_markers(const H5::Group& shandle, int num_features) {
    std::vector<size_t> dims{ static_cast<size_t>(num_features) };
    utils::check_and_open_dataset(shandle, "means", H5T_FLOAT, dims);
    utils::check_and_open_dataset(shandle, "detected", H5T_FLOAT, dims);
    for (const auto& eff : markers::effects) {
        utils::check_and_open_dataset(shandle, eff, H5T_FLOAT, dims);
    }
}

inline void validate_results(const H5::Group& handle, const std::vector<std::string>& selections, int num_features) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    auto mhandle = utils::check_and_open_group(rhandle, "markers");
    if (mhandle.getNumObjs() != selections.size()) {
        throw std::runtime_error("number of groups in 'markers' is not consistent with the expected number of selections");
    }

    for (const auto& s : selections) {
        try {
            auto shandle = utils::check_and_open_group(mhandle, s);
            validate_custom_markers(shandle, num_features);
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve statistics for selection '" + s + "' in 'results/markers'");
        }
    }
    return;
}

inline void validate_results(const H5::Group& handle, const std::vector<std::string>& selections, const std::vector<std::string>& modalities, const std::vector<int>& num_features) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    auto mhandle = utils::check_and_open_group(rhandle, "per_selection");
    if (mhandle.getNumObjs() != selections.size()) {
        throw std::runtime_error("number of groups in 'per_selection' is not consistent with the expected number of selections");
    }

    for (const auto& s : selections) {
        try {
            auto shandle = utils::check_and_open_group(mhandle, s);
            for (size_t a = 0; a < modalities.size(); ++a) {
                try {
                    auto ahandle = utils::check_and_open_group(shandle, modalities[a]);
                    validate_custom_markers(ahandle, num_features[a]);
                } catch (std::exception& e) {
                    throw utils::combine_errors(e, "failed to retrieve statistics for modality '" + modalities[a] + "'");
                }
            }
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve statistics for selection '" + s + "' in 'results/per_selection'");
        }
    }
    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the custom selections step, see [here](@ref details-custom_selections) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of high-quality cells in the dataset, i.e., after any quality-based filtering has been applied.
 * @param modalities Modalities available in the dataset, should be some combination of `"RNA"` or `"ADT"`.
 * If `version < 2000000`, this is ignored and an RNA modality is always assumed.
 * @param num_features Number of features for each modality in `modalities`.
 * If `version < 2000000`, only the first value is used and is assumed to refer to the number of genes for the RNA modality.
 * @param version Version of the format.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::Group& handle, int num_cells, const std::vector<std::string>& modalities, const std::vector<int>& num_features, int version) {
    auto mhandle = utils::check_and_open_group(handle, "custom_selections");

    std::vector<std::string> collected;
    try {
        collected = validate_parameters(mhandle, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'custom_selections'");
    }

    try {
        if (version >= 2000000) {
            validate_results(mhandle, collected, modalities, num_features);
        } else {
            validate_results(mhandle, collected, num_features[0]);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'custom_selections'");
    }

    return;
}

}

}

#endif
