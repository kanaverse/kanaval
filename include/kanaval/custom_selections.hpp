#ifndef KANAVAL_CUSTOM_SELECTIONS_HPP
#define KANAVAL_CUSTOM_SELECTIONS_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"

/**
 * @file custom_selections.hpp
 *
 * @brief Validate custom selection contents.
 */

namespace kanaval {

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
    }

    return output;
}

inline void validate_results(const H5::Group& handle, const std::vector<std::string>& selections, int num_genes) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    
    auto mhandle = utils::check_and_open_group(rhandle, "markers");
    if (mhandle.getNumObjs() != selections.size()) {
        throw std::runtime_error("number of groups in 'markers' is not consistent with the expected number of selections");
    }

    for (const auto& s : selections) {
        try {
            auto shandle = utils::check_and_open_group(mhandle, s);
            std::vector<size_t> dims{ static_cast<size_t>(num_genes) };
            utils::check_and_open_dataset(shandle, "means", H5T_FLOAT, dims);
            utils::check_and_open_dataset(shandle, "detected", H5T_FLOAT, dims);
            for (const auto& eff : markers::effects) {
                utils::check_and_open_dataset(shandle, eff, H5T_FLOAT, dims);
            }
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve statistics for selection '" + s + "' in 'results/markers'");
        }
    }
    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the custom selections step.
 * 
 * `handle` should contain a `custom_selections` group, itself containing the `parameters` and `results` subgroups.
 *
 * `parameters` should contain:
 *
 * - `selections`: a group defining the custom selections.
 *   Each child is named after a user-created selection.
 *   Each child is an integer dataset of arbitrary length containing the indices of the selected cells.
 *   Note that indices refer to the dataset after QC filtering and should be less than `num_cells`.
 * 
 * `results` should contain:
 *
 * - `markers`: a group containing the marker results for each selection after a comparison to a group containing all other cells.
 *   Each child is named after its selection and is a group containing:
 *   - `means`: a float dataset of length equal to the number of genes, containing the mean expression of each gene in the selection.
 *   - `detected`: a float dataset of length equal to the number of genes, containing the proportion of cells with detected expression of each gene in the selection.
 *   - `lfc`: a float dataset of length equal to the number of genes, containing the log-fold change in the selection compared to all other cells.
 *   - `delta_detected`: same as `lfc`, but for the delta-detected (i.e., difference in the percentage of detected expression).
 *   - `cohen`: same as `lfc`, but for Cohen's d.
 *   - `auc`: same as `lfc`, but for the AUCs.
 *
 * @param handle An open HDF5 file handle.
 * @param num_genes Number of genes in the dataset.
 * @param num_cells Number of cells in the dataset after QC filtering.
 *
 * @return If the format is invalid, an error is raised.
 */
inline void validate(const H5::Group& handle, int num_genes, int num_cells) {
    auto mhandle = utils::check_and_open_group(handle, "custom_selections");

    std::vector<std::string> collected;
    try {
        collected = validate_parameters(mhandle, num_cells);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'custom_selections'");
    }

    try {
        validate_results(mhandle, collected, num_genes);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'custom_selections'");
    }

    return;
}

}

}

#endif
