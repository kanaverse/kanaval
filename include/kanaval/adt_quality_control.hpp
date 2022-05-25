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
inline void validate_parameters(const H5::Group& qhandle) {
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
    return;
}

inline int validate_results(const H5::Group& qhandle, int num_cells, int num_samples, bool adt_in_use) {
    auto rhandle = utils::check_and_open_group(qhandle, "results");

    int remaining = -1;
    if (adt_in_use) {
        try {
            auto mhandle = utils::check_and_open_group(rhandle, "metrics");
            std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
            utils::check_and_open_dataset(mhandle, "sums", H5T_FLOAT, dims);
            utils::check_and_open_dataset(mhandle, "detected", H5T_INTEGER, dims);
            utils::check_and_open_dataset(mhandle, "igg_total", H5T_FLOAT, dims);
        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve metrics from 'results'");
        }

        try {
            auto thandle = utils::check_and_open_group(rhandle, "thresholds");

            std::vector<size_t> dims{ static_cast<size_t>(num_samples) };
            utils::check_and_open_dataset(thandle, "detected", H5T_FLOAT, dims);
            utils::check_and_open_dataset(thandle, "igg_total", H5T_FLOAT, dims);

        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to retrieve thresholds from 'results'");
        }

        remaining = quality_control::check_discard_vector(rhandle, num_cells);
    }

    return remaining;
}
/**
 * @endcond
 */

/**
 * Check contents for the quality control step on the ADT count matrix.
 * Contents are stored inside an `adt_quality_control` HDF5 group at the root of the file.
 * The `adt_quality_control` group itself contains the `parameters` and `results` subgroups.
 * 
 * No ADT data was available prior to version 2.0 of the format, so the `adt_quality_control` group may be absent in such files.
 *
 * <HR>
 * `parameters` should contain:
 * 
 * - `igg_prefix`: a scalar string containing the expected prefix for IgG features.
 * - `nmads`: a scalar float specifying the number of MADs to use to define the QC thresholds.
 * - `min_detected_drop`: a scalar float specifying the minimum relative drop in the number of detected features before a cell is considered to be low-quality.
 * 
 * <HR>
 * If `adt_in_use = false`, `results` should be empty.
 *
 * If `adt_in_use = true`, `results` should contain:
 * 
 * - `metrics`, a group containing per-cell QC metrics derived from the RNA count data.
 *   This contains:
 *   - `sums`: a float dataset of length equal to the number of cells, containing the total count for each cell.
 *   - `detected`:  an integer dataset of length equal to the number of cells, containing the total number of detected features for each cell.
 *   - `igg_total`: a float dataset of length equal to the number of cells, containing the total count in IgG features.
 * - `thresholds`, a group containing thresholds on the metrics for each sample.
 *   This contains:
 *   - `detected`:  a float dataset of length equal to the number of samples, containing the threshold on the total number of detected features for each sample.
 *   - `igg_total`: a float dataset of length equal to the number of samples, containing the threshold on the total counts in IgG features for each sample.
 * - `discards`: an integer dataset of length equal to the number of cells.
 *   Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the ADT-based filter thresholds.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality control is applied.
 * @param num_samples Number of samples in the dataset.
 * @param adt_in_use Whether ADTs are being used in this dataset.
 * @param version Version of the format.
 * 
 * @return The number of high-quality cells, according to the ADT-based metrics.
 * Alternatively, if ADTs were not present in the dataset, -1 is returned instead.
 * An error is raised if the format is invalid.
 */ 
inline int validate(const H5::H5File& handle, int num_cells, int num_samples, bool adt_in_use, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return -1;
    }

    auto qhandle = utils::check_and_open_group(handle, "adt_quality_control");

    try {
        validate_parameters(qhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'adt_quality_control'");
    }

    int remaining = 0;
    try {
        remaining = validate_results(qhandle, num_cells, num_samples, adt_in_use);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'adt_quality_control'");
    }

    return remaining;
}

}

}

#endif
