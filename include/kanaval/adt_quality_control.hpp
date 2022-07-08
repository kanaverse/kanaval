#ifndef KANAVAL_ADT_QUALITY_CONTROL_HPP
#define KANAVAL_ADT_QUALITY_CONTROL_HPP

#include "H5Cpp.h"
#include <vector>
#include "utils.hpp"
#include "misc.hpp"

/**
 * @file adt_quality_control.hpp
 *
 * @brief Validate ADT quality control contents.
 */

namespace kanaval {

namespace adt_quality_control {

/**
 * @cond
 */
inline bool validate_parameters(const H5::Group& qhandle, int version) {
    auto phandle = utils::check_and_open_group(qhandle, "parameters");
    utils::check_and_open_scalar(phandle, "igg_prefix", H5T_STRING);

    auto nmads = utils::load_float_scalar<>(phandle, "nmads");
    if (nmads < 0) {
        throw std::runtime_error("number of MADs in 'nmads' should be non-negative");
    }

    auto mindrop = utils::load_float_scalar<>(phandle, "min_detected_drop");
    if (mindrop < 0 || mindrop >= 1) {
        throw std::runtime_error("minimum detected drop should lie in [0, 1)");
    }

    bool skip = false;
    if (version >= 2001000) {
        skip = utils::load_integer_scalar(phandle, "skip");
    }
    return skip;
}

inline int validate_results(const H5::Group& qhandle, int num_cells, int num_samples, bool skip, bool adt_in_use) {
    auto rhandle = utils::check_and_open_group(qhandle, "results");

    int remaining = -1;
    if (adt_in_use) {
        if (!skip || rhandle.exists("metrics")) {
            try {
                auto mhandle = utils::check_and_open_group(rhandle, "metrics");
                std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
                utils::check_and_open_dataset(mhandle, "sums", H5T_FLOAT, dims);
                utils::check_and_open_dataset(mhandle, "detected", H5T_INTEGER, dims);
                utils::check_and_open_dataset(mhandle, "igg_total", H5T_FLOAT, dims);
            } catch (std::exception& e) {
                throw utils::combine_errors(e, "failed to retrieve metrics from 'results'");
            }
        }

        if (!skip || rhandle.exists("thresholds")) {
            try {
                auto thandle = utils::check_and_open_group(rhandle, "thresholds");

                std::vector<size_t> dims{ static_cast<size_t>(num_samples) };
                utils::check_and_open_dataset(thandle, "detected", H5T_FLOAT, dims);
                utils::check_and_open_dataset(thandle, "igg_total", H5T_FLOAT, dims);

            } catch (std::exception& e) {
                throw utils::combine_errors(e, "failed to retrieve thresholds from 'results'");
            }
        }

        remaining = quality_control::check_discard_vector(rhandle, num_cells, skip);
    }

    return remaining;
}
/**
 * @endcond
 */

/**
 * Check contents for the quality control step on the ADT count matrix, see [here](@ref details-adt_quality_control) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality control is applied.
 * @param num_samples Number of samples in the dataset.
 * @param adt_in_use Whether ADTs are being used in this dataset.
 * @param version Version of the kana file.
 * 
 * @return A pair where the first element specifies whether this QC step was skipped,
 * and the second element contains the number of high-quality cells according to the ADT-based metrics.
 * If ADTs were not present in the dataset, -1 is returned in the second element instead.
 * An error is raised if the format is invalid.
 */ 
inline std::pair<bool, int> validate(const H5::H5File& handle, int num_cells, int num_samples, bool adt_in_use, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return std::make_pair(false, -1);
    }

    auto qhandle = utils::check_and_open_group(handle, "adt_quality_control");

    bool skip;
    try {
        skip = validate_parameters(qhandle, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_quality_control'");
    }

    int remaining;
    try {
        remaining = validate_results(qhandle, num_cells, num_samples, skip, adt_in_use);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_quality_control'");
    }

    return std::make_pair(skip, remaining);
}

}

}

#endif
