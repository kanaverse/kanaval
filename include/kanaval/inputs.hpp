#ifndef KANAVAL_INPUTS_HPP
#define KANAVAL_INPUTS_HPP

#include "H5Cpp.h"
#include "utils.hpp"
#include <stdexcept>
#include <vector>
#include <string>
#include <numeric>
#include <unordered_set>
#include <unordered_map>

/**
 * @file inputs.hpp
 *
 * @brief Validate input contents.
 */

namespace kanaval {

/** 
 * Validation of inputs
 */
namespace inputs {

/**
 * @brief Details about the dataset.
 */
struct Details {
    /**
     * Available modalities in the dataset.
     * Currently, this may contain `"RNA"` and/or `"ADT"`.
     */
    std::vector<std::string> modalities;

    /**
     * Number of features for each modality listed in `modalities`.
     * For datasets containing multiple samples, each entry contains the number of features in the intersection of feature spaces across all samples.
     */
    std::vector<int> num_features;

    /**
     * Number of cells.
     * For multi-sample datasets, this considers the total number of cells in all samples.
     */
    int num_cells;

    /**
     * Number of samples.
     * Note that a single matrix may contain multiple samples.
     */
    int num_samples;
};

/**
 * @cond
 */
struct ParamDump {
    int num_matrices;
    bool multi_matrix;
    bool multi_sample;
    int subset_num;
};

inline ParamDump validate_parameters(const H5::Group& handle, bool embedded, int version) {
    auto phandle = utils::check_and_open_group(handle, "parameters");
    ParamDump output;

    // Formats can either be a scalar... or not.
    std::vector<std::string> formats;
    output.multi_matrix = false;
    {
        auto fhandle = utils::check_and_open_dataset(phandle, "format", H5T_STRING);
        auto fspace = fhandle.getSpace();
        if (fspace.getSimpleExtentNdims() == 0) {
            formats.push_back(utils::load_string(fhandle));
        } else {
            if (version < 1001000) {
                throw std::runtime_error("'format' should be a scalar string in version 1.0");
            }
            output.multi_matrix = true;
            formats = utils::load_string_vector(fhandle);
        }
    }
    output.num_matrices = formats.size();

    auto fihandle = utils::check_and_open_group(phandle, "files");
    auto nfiles = fihandle.getNumObjs();

    // Checking the runs.
    std::vector<int> runs;
    if (output.multi_matrix) {
        runs = utils::load_integer_vector(phandle, "sample_groups");
        if (runs.size() != formats.size()) {
            throw std::runtime_error("'sample_groups' and 'format' should have the same length");
        }

        int total_files = std::accumulate(runs.begin(), runs.end(), 0);
        if (total_files != static_cast<int>(nfiles)) {
            throw std::runtime_error("sum of 'sample_groups' is not equal to the length of 'files'");
        }

        // Checking that everyone has unique and sorted names.
        auto names = utils::load_string_vector(phandle, "sample_names");
        if (names.size() != formats.size()) {
            throw std::runtime_error("'sample_names' and 'format' should have the same length");
        }

        if (version >= 2001000) {
            if (!utils::is_unique_and_sorted(names)) {
                throw std::runtime_error("duplicated or unsorted values in 'sample_names'");
            }
        } else {
            std::sort(names.begin(), names.end());
            if (!utils::is_unique_and_sorted(names)) {
                throw std::runtime_error("duplicated values in 'sample_names'");
            }
        }
    } else {
        runs.push_back(nfiles);
    }

    // Checking the files.
    int sofar = 0;
    std::vector<std::pair<hsize_t, hsize_t> > bytes;
    for (size_t r = 0; r < runs.size(); ++r) {
        auto curf = formats[r];
        std::vector<std::string> types;

        for (int s = 0; s < runs[r]; ++s, ++sofar) {
            std::string current = std::to_string(sofar);
            try {
                auto curfihandle = utils::check_and_open_group(fihandle, current);

                utils::check_and_open_dataset(curfihandle, "name", H5T_STRING, {});
                types.push_back(utils::load_string(curfihandle, "type"));

                if (embedded) {
                    bytes.emplace_back(
                        utils::load_integer_scalar<hsize_t>(curfihandle, "offset"),
                        utils::load_integer_scalar<hsize_t>(curfihandle, "size")
                    );
                } else {
                    utils::check_and_open_dataset(curfihandle, "id", H5T_STRING, {});
                }
            } catch (std::exception& e) {
                throw utils::combine_errors(e, "failed to retrieve information for file " + current);
            }
        }

        if (curf == "MatrixMarket") {
            std::unordered_map<std::string, int> expected;
            expected["mtx"] = 0;
            expected["genes"] = 0;
            expected["annotations"] = 0;

            for (auto t : types) {
                auto it = expected.find(t);
                if (it == expected.end()) {
                    throw std::runtime_error("unknown file type '" + t + "' when format is 'MatrixMarket'");
                }
                ++(it->second);
            }

            if (expected["mtx"] != 1) {
                throw std::runtime_error("expected exactly one 'mtx' file when format is 'MatrixMarket'");
            }
            if (expected["genes"] > 1) {
                throw std::runtime_error("expected no more than one 'genes' file when format is 'MatrixMarket'");
            }
            if (expected["annotations"] > 1) {
                throw std::runtime_error("expected no more than one 'annotation' file when format is 'MatrixMarket'");
            }
            
        } else if (curf == "10X") {
            if (types.size() != 1 || types.front() != "h5") {
                throw std::runtime_error("expected exactly one 'h5' file when format is '10X'");
            }
        } else if (curf == "H5AD") {
            if (types.size() != 1 || types.front() != "h5") {
                throw std::runtime_error("expected exactly one 'h5' file when format is 'H5AD'");
            }
        }
    }

    // Checking the files make sense.
    if (embedded) {
        hsize_t sofar = 0;
        for (const auto& b : bytes) {
            if (b.first != sofar) {
                throw std::runtime_error("offsets and sizes of 'files' are not sorted and contiguous");
            }
            sofar += b.second;
        }
    }

    // Checking if there's a batch variable.
    if (!output.multi_matrix && phandle.exists("sample_factor")) {
        utils::check_and_open_dataset(phandle, "sample_factor", H5T_STRING, {});
        output.multi_sample = true;
    }  else {
        output.multi_sample = output.multi_matrix;
    }

    // Checking if there's a subset group.
    output.subset_num = -1;

    if (phandle.exists("subset")) {
        auto subhandle = utils::check_and_open_group(phandle, "subset");

        if (subhandle.exists("indices")) {
            auto subidx = utils::load_integer_vector(subhandle, "indices");

            for (auto i : subidx) {
                if (i < 0) {
                    throw std::runtime_error("indices in 'subset/indices' should be non-negative");
                }
            }

            if (!utils::is_unique_and_sorted(subidx)) {
                throw std::runtime_error("indices in 'subset/indices' should be unique and sorted");
            }

            output.subset_num = subidx.size();
        } else {
            utils::check_and_open_dataset(subhandle, "field", H5T_STRING, {});
            auto vhandle = utils::check_and_open_dataset(subhandle, "values", H5T_STRING);
            auto vdims = utils::load_dataset_dimensions(vhandle);
            if (vdims.size() != 1) {
                throw std::runtime_error("'subset/values' should be a 1-dimensional string dataset");
            }
        }
    }

    return output;
}

inline Details validate_results(const H5::Group& handle, const ParamDump& param_info, int version) {
    auto rhandle = utils::check_and_open_group(handle, "results");
    Details output;

    // Pulling out the dimensions and modalities.
    if (version < 2000000) {
        output.modalities.push_back("RNA");

        auto dims = utils::load_integer_vector<int>(rhandle, "dimensions");
        if (dims.size() != 2) {
            throw std::runtime_error("'dimensions' should be a dataset of length 2");
        }
        if (dims[0] < 0 || dims[1] < 0) {
            throw std::runtime_error("'dimensions' should contain non-negative integers");
        }
        output.num_features.push_back(dims[0]);
        output.num_cells = dims[1];
    } else {
        output.num_cells = utils::load_integer_scalar<>(rhandle, "num_cells");

        auto fhandle = utils::check_and_open_group(rhandle, "num_features");
        size_t nmodals = fhandle.getNumObjs();
        if (nmodals == 0) {
            throw std::runtime_error("number of modalities should be positive");
        }

        for (hsize_t idx = 0; idx < nmodals; ++idx) {
            std::string modality = fhandle.getObjnameByIdx(idx);
            output.modalities.push_back(modality);
            output.num_features.push_back(utils::load_integer_scalar<>(fhandle, modality));
        }
    }

    // Checking the number of samples.
    output.num_samples = 1;
    if (rhandle.exists("num_samples")) {
        output.num_samples = utils::load_integer_scalar<int>(rhandle, "num_samples");
    }
    if (param_info.multi_matrix) {
        if (output.num_samples != param_info.num_matrices) {
            throw std::runtime_error("'num_samples' should be equal to the number of matrices");
        }
    } else {
        if (!param_info.multi_sample && output.num_samples != 1) {
            throw std::runtime_error("'num_samples' should be 1 for single matrix inputs without 'sample_factor'");
        }
    }

    auto check_unique = [](const std::vector<int>& idx, const std::string& msg) -> void {
        std::unordered_set<int> used;
        used.reserve(idx.size());
        for (auto i : idx) {
            if (i < 0) {
                throw std::runtime_error(msg + " contains negative values");
            } else if (used.find(i) != used.end()) {
                throw std::runtime_error(msg + " contains duplicate values");
            }
            used.insert(i);
        }
    };

    if (version >= 2000000) {
        auto ihandle = utils::check_and_open_group<>(rhandle, "identities");
        for (size_t m = 0; m < output.modalities.size(); ++m) {
            auto idx = utils::load_integer_vector<int>(ihandle, output.modalities[m]);
            if (idx.size() != static_cast<size_t>(output.num_features[m])) {
                throw std::runtime_error("'identities' for modality '" + output.modalities[m] + "' should have length equal to its number of features");
            }
            check_unique(idx, "'identities' for modality '" + output.modalities[m] + "'");
        }

    } else if (version >= 1002000) {
        auto idx = utils::load_integer_vector<int>(rhandle, "identities");
        if (idx.size() != static_cast<size_t>(output.num_features[0])) {
            throw std::runtime_error("'identities' should have length equal to the number of genes");
        }
        check_unique(idx, "'identities'");
        
    } else {
        if (param_info.multi_matrix) {
            auto idx = utils::load_integer_vector<int>(rhandle, "indices");
            if (idx.size() != static_cast<size_t>(output.num_features[0])) {
                throw std::runtime_error("'indices' should have length equal to the number of genes");
            }
            check_unique(idx, "'indices'");

        } else {
            auto perms = utils::load_integer_vector<int>(rhandle, "permutation");
            if (perms.size() != static_cast<size_t>(output.num_features[0])) {
                throw std::runtime_error("'permutation' should have length equal to the number of genes");
            }

            // Note that the code below implies that all consecutive entries are present,
            // otherwise we would see duplicates.
            std::vector<unsigned char> used(perms.size());
            for (auto p : perms) {
                if (p < 0 || static_cast<size_t>(p) >= perms.size()) {
                    throw std::runtime_error("'permutation' contains out-of-range values");
                } else if (used[p]) {
                    throw std::runtime_error("duplicated index in 'permutation'");
                }
                used[p] = 1;
            }
        }
    }

    if (param_info.subset_num > -1) {
        if (param_info.subset_num != output.num_cells) {
            throw std::runtime_error("inconsistent number of cells in 'parameters/subset/indices' and 'results/num_cells'");
        }
    }

    return output;
}
/**
 * @endcond
 */

/**
 * Check contents for the input step, see [here](@ref details-inputs) for details.
 *
 * @param handle An open HDF5 file handle.
 * @param embedded Whether the data files are embedded or linked.
 * @param version Version of the state file.
 *
 * @return Details about the dataset.
 * If the format is invalid, an error is raised instead.
 */
inline Details validate(const H5::Group& handle, bool embedded, int version) {
    auto ihandle = utils::check_and_open_group(handle, "inputs");

    ParamDump dump; 
    try {
        dump = validate_parameters(ihandle, embedded, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'inputs'");
    }

    Details output;
    try {
        output = validate_results(ihandle, dump, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'inputs'");
    }

    return output;
}

}

}

#endif
