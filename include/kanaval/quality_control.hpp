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
inline void validate_parameters(const H5::Group& qhandle) {
    auto phandle = utils::check_and_open_group(qhandle, "parameters");
    utils::check_and_open_scalar(phandle, "use_mito_default", H5T_INTEGER);
    utils::check_and_open_scalar(phandle, "mito_prefix", H5T_STRING);
    auto nmads = utils::load_float_scalar<>(phandle, "nmads");
    if (nmads < 0) {
        throw std::runtime_error("number of MADs in 'nmads' should be non-negative");
    }
    return;
}

inline int validate_results(const H5::Group& qhandle, int num_cells, int num_samples) {
    auto rhandle = utils::check_and_open_group(qhandle, "results");

    try {
        auto mhandle = utils::check_and_open_group(rhandle, "metrics");
        std::vector<size_t> dims{ static_cast<size_t>(num_cells) };
        utils::check_and_open_dataset(mhandle, "sums", H5T_FLOAT, dims);
        utils::check_and_open_dataset(mhandle, "detected", H5T_INTEGER, dims);
        utils::check_and_open_dataset(mhandle, "proportion", H5T_FLOAT, dims);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve metrics from 'results'");
    }

    try {
        auto thandle = utils::check_and_open_group(rhandle, "thresholds");

        std::vector<size_t> dims{ static_cast<size_t>(num_samples) };
        utils::check_and_open_dataset(thandle, "sums", H5T_FLOAT, dims);
        utils::check_and_open_dataset(thandle, "detected", H5T_FLOAT, dims);
        utils::check_and_open_dataset(thandle, "proportion", H5T_FLOAT, dims);

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve thresholds from 'results'");
    }

    int remaining = check_discard_vector(rhandle, num_cells);

    return remaining;
}
/**
 * @endcond
 */

/**
 * Check contents for the RNA-based quality control step.
 * Contents are stored inside a `quality_control` HDF5 group at the root of the file.
 * The `quality_control` group itself contains the `parameters` and `results` subgroups.
 * 
 * <HR>
 * `parameters` should contain:
 * 
 * - `use_mito_default`: a scalar integer to be interpreted as a boolean.
 *   This specifies whether to use the default mitochondrial gene list.
 * - `mito_prefix`: a scalar string containing the expected prefix for mitochondrial gene symbols.
 * - `nmads`: a scalar float specifying the number of MADs to use to define the QC thresholds.
 * 
 * <HR>
 * `results` should contain:
 * 
 * - `metrics`, a group containing per-cell QC metrics derived from the RNA count data.
 *   This contains:
 *   - `sums`: a float dataset of length equal to the number of cells, containing the total count for each cell.
 *   - `detected`:  an integer dataset of length equal to the number of cells, containing the total number of detected genes for each cell.
 *   - `proportion`: a float dataset of length equal to the number of cells, containing the percentage of counts in (mitochondrial) genes.
 * - `thresholds`, a group containing thresholds on the metrics for each batch.
 *   This contains:
 *   - `sums`: a float dataset of length equal to the number of batches, containing the total count threshold for each batch.
 *   - `detected`:  an integer dataset of length equal to the number of batches, containing the threshold on the total number of detected genes for each batch.
 *   - `proportion`: a float dataset of length equal to the number of batches, containing the threshold on the percentage of counts in (mitochondrial) genes for each batch.
 * - `discards`: an integer dataset of length equal to the number of cells.
 *   Each value is interpreted as a boolean and specifies whether the corresponding cell would be discarded by the RNA-based filter thresholds.
 *
 * <HR>
 * @param handle An open HDF5 file handle.
 * @param num_cells Number of cells in the dataset before any quality filtering is applied.
 * @param num_samples Number of batches in the dataset.
 * 
 * @return The number of high-quality cells, according to the RNA-based metrics.
 * If the format is invalid, an error is raised instead.
 */ 
inline int validate(const H5::H5File& handle, int num_cells, int num_samples) {
    auto qhandle = utils::check_and_open_group(handle, "quality_control");

    try {
        validate_parameters(qhandle);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'quality_control'");
    }

    int remaining = 0;
    try {
        remaining = validate_results(qhandle, num_cells, num_samples);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'quality_control'");
    }

    return remaining;
}

}

}

#endif
