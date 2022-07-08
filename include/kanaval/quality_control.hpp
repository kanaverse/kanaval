#ifndef KANAVAL_QUALITY_CONTROL_HPP
#define KANAVAL_QUALITY_CONTROL_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

/**
 * @file quality_control.hpp
 *
 * @brief Validate quality control contents.
 */

namespace kanaval {

namespace quality_control {

/**
 * @cond
 */
inline bool validate_parameters(const H5::Group& qhandle, int version) {
    auto phandle = utils::check_and_open_group(qhandle, "parameters");
    utils::check_and_open_scalar(phandle, "use_mito_default", H5T_INTEGER);
    utils::check_and_open_scalar(phandle, "mito_prefix", H5T_STRING);

    auto nmads = utils::load_float_scalar<>(phandle, "nmads");
    if (nmads < 0) {
        throw std::runtime_error("number of MADs in 'nmads' should be non-negative");
    }

    bool skip = false;
    if (version >= 2001000) {
        skip = utils::load_integer_scalar(phandle, "skip");
    }
    return skip;
}

inline int validate_results(const H5::Group& qhandle, int num_cells, int num_samples, bool skip) {
    auto rhandle = utils::check_and_open_group(qhandle, "results");

    if (!skip || rhandle.exists("metrics")) {
        try {
            auto mhandle = utils::check_and_open_group(rhandle, "metrics");
            std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
            utils::check_and_open_dataset(mhandle, "sums", H5T_FLOAT, dims);
            utils::check_and_open_dataset(mhandle, "detected", H5T_INTEGER, dims);
            utils::check_and_open_dataset(mhandle, "proportion", H5T_FLOAT, dims);
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve metrics from 'results'");
        }
    }

    if (!skip || rhandle.exists("thresholds")) {
        try {
            auto thandle = utils::check_and_open_group(rhandle, "thresholds");

            std::vector<size_t> dims{ static_cast<size_t>(num_samples) };
            utils::check_and_open_dataset(thandle, "sums", H5T_FLOAT, dims);
            utils::check_and_open_dataset(thandle, "detected", H5T_FLOAT, dims);
            utils::check_and_open_dataset(thandle, "proportion", H5T_FLOAT, dims);

        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve thresholds from 'results'");
        }
    }

    int remaining = check_discard_vector(rhandle, num_cells, skip);

    return remaining;
}
/**
 * @endcond
 */

/**
 * Check contents for the quality control step on the RNA count matrix, see [here](@ref details-quality_control) for more details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality filtering is applied.
 * @param num_samples Number of batches in the dataset.
 * @param version Version of the kana file.
 * 
 * @return A pair where the first element specifies whether this QC step was skipped,
 * and the second element contains the number of high-quality cells according to the RNA-based metrics.
 * If the format is invalid, an error is raised instead.
 */ 
inline std::pair<bool, int> validate(const H5::H5File& handle, int num_cells, int num_samples, int version) {
    auto qhandle = utils::check_and_open_group(handle, "quality_control");

    bool skip;
    try {
        skip = validate_parameters(qhandle, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'quality_control'");
    }

    int remaining;
    try {
        remaining = validate_results(qhandle, num_cells, num_samples, skip);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'quality_control'");
    }

    return std::make_pair(skip, remaining);
}

}

}

#endif
