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

namespace inputs {

/**
 * @cond
 */
inline std::pair<bool, bool> validate_parameters(const H5::Group& handle, bool embedded, int version = 1001000) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    // Formats can either be a scalar... or not.
    std::vector<std::string> formats;
    bool multifile = false;
    {
        auto fhandle = utils::check_and_open_dataset(phandle, "format", H5T_STRING);
        auto fspace = fhandle.getSpace();
        if (fspace.getSimpleExtentNdims() == 0) {
            formats.push_back(utils::load_string(fhandle));
        } else {
            if (version < 1001000) {
                throw std::runtime_error("'format' should be a scalar string in version 1.0");
            }
            multifile = true;
            formats = utils::load_string_vector(fhandle);
        }
    }

    auto fihandle = utils::check_and_open_group(phandle, "files");
    auto nfiles = fihandle.getNumObjs();

    // Checking the runs.
    std::vector<int> runs;
    if (multifile) {
        runs = utils::load_integer_vector(phandle, "sample_groups");
        if (runs.size() != formats.size()) {
            throw std::runtime_error("'sample_groups' and 'format' should have the same length");
        }

        int total_files = std::accumulate(runs.begin(), runs.end(), 0);
        if (total_files != static_cast<int>(nfiles)) {
            throw std::runtime_error("sum of 'sample_groups' is not equal to the length of 'files'");
        }

        // Checking that everyone has unique groups.
        auto names = utils::load_string_vector(phandle, "sample_names");
        if (names.size() != formats.size()) {
            throw std::runtime_error("'sample_names' and 'format' should have the same length");
        }

        std::unordered_set<std::string> stuff;
        for (auto s : names) {
            if (stuff.find(s) != stuff.end()) {
                throw std::runtime_error("duplicated sample name '" + s + "' in 'sample_names'");
            }
            stuff.insert(s);
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
            expected["annotation"] = 0;

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
            if (expected["annotation"] > 1) {
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
        std::sort(bytes.begin(), bytes.end());
        hsize_t sofar = 0;
        for (const auto& b : bytes) {
            if (b.first != sofar) {
                throw std::runtime_error("offsets and sizes of 'files' are not contiguous");
            }
            sofar += b.second;
        }
    }

    // Checking if there's a batch variable.
    bool multisample = multifile;
    if (phandle.exists("sample_factor")) {
        utils::check_and_open_dataset(phandle, "sample_factor", H5T_STRING, {});
        multisample = true;
    } 
    
    return std::make_pair(multifile, multisample);
}

inline void validate_results(const H5::Group& handle, bool blocked, int version = 1001000) {
    auto rhandle = utils::check_and_open_group(handle, "results");

    auto dims = utils::load_integer_vector<int>(rhandle, "dimensions");
    if (dims.size() != 2) {
        throw std::runtime_error("'dimensions' should be a dataset of length 2");
    }
    if (dims[0] < 0 || dims[1] < 0) {
        throw std::runtime_error("'dimensions' should contain non-negative integers");
    }

    if (blocked) {
        auto idx = utils::load_integer_vector<int>(rhandle, "indices");
        if (idx.size() != static_cast<size_t>(dims[0])) {
            throw std::runtime_error("'indices' should have length equal to the number of genes");
        }
        for (auto i : idx) {
            if (i < 0) {
                throw std::runtime_error("'indices' contains negative values");
            }
        }

    } else {
        auto perms = utils::load_integer_vector<int>(rhandle, "permutation");
        if (perms.size() != static_cast<size_t>(dims[0])) {
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

    return;
}
/**
 * @endcond
 */

/**
 * Check contents for the input step.
 * Contents are stored inside an `inputs` HDF5 group at the root of the file.
 * The `inputs` group itself contains the `parameters` and `results` subgroups.
 *
 * In this section, a "matrix" refers to one or more files describing a single (count) matrix.
 * This should be exactly one file for HDF5-based formats, or multiple files for MatrixMarket formats, e.g., to include feature information - see below for details.
 * A matrix may contain data for one or more samples.
 *
 * @section Parameters
 * `parameters` should contain:
 * 
 * - `format`: a scalar string specifying the file format for a single matrix.
 *   This is usually either `"MatrixMarket"`, for a MatrixMarket file with possible feature/barcode annotations;
 *   `"10X"`, for the 10X Genomics HDF5 matrix format;
 *   or `"H5AD"`, for the H5AD format.
 *   Other values are allowed but their interpretation is implementation-defined (e.g., for in-house resources). 
 *   @v1_1{\[**since version 1.1**\] For multiple matrices, `format` should instead be a 1-dimensional string dataset of length equal to the number of uploads.
 *   Each element of the dataset is usually one of `"MatrixMarket"`, `"10X"` or `"H5AD"`; 
 *   different values can be present for mixed input formats.}
 * - `files`: a group of groups representing an array of input file information.
 *   Each inner group is named by their positional index in the array and contains information about a file in an upload.
 *   Each inner group should contain:
 *   - `type`: a scalar string specifying the type of the file.
 *     @v1_1{\[**since version 1.1**\] For multiple matrices, the constraints below apply to all files corresponding to a single matrix.}
 *     - If `format = "MatrixMarket"`, there should be exactly one `type = "mtx"` corresponding to the (possibly Gzipped) `*.mtx` file.
 *       There may be zero or one `type = "genes"`, containing a (possibly Gzipped) TSV file with the Ensembl and gene symbols for each row.
 *       There may be zero or one `type = "annotation"`, containing a (possibly Gzipped) TSV file with the annotations for each column.
 *     - If `format = "10X"` or `"H5AD"`, there should be exactly one `type = "h5"`.
 *   - `name`: a scalar string specifying the file name as it was provided to **kana**.
 *   - `offset`: a scalar integer specifying where the file starts as an offset from the start of the remaining bytes section.
 *   - `size`: a scalar integer specifying the number of bytes in the file.
 *
 * @v1_1{\[**since version 1.1**\] For multiple matrices, `parameters` should also contain:}
 *
 * - @v1_1{`sample_groups`: an integer dataset of length equal to the number of samples.
 *   Each entry specifies the number of files in `files` that belong to a sample.
 *   (All files from the same sample are assumed to be contiguous in the array represented by `files`;
 *   so a `sample_groups` of `[3, 2, 1]` would mean that the first three files belong to the first sample, 
 *   the next 2 files belong to the second sample, and the last file belongs to the third sample.)}
 * - @v1_1{`sample_names`: a string dataset of length equal to the number of samples, containing the sample name.}
 *
 * @v1_1{\[**since version 1.1**\] For single matrix inputs, `parameters` may also contain:}
 *
 * - @v1_1{`sample_factor`: a string scalar specifying the field in the per-cell annotation that contains the sample blocking factor. 
 *   If present, it is assumed that the matrix contains data for multiple samples.}
 *
 * @section Results
 * `results` should contain:
 * 
 * - `dimensions`: an integer dataset of length 2,
 *   containing the number of features and the number of cells in the dataset.
 *   @v1_1{\[**since version 1.1**\] When dealing with multi-sample inputs, the first entry is instead defined as the size of the intersection of features across all samples.}
 *
 * If there is only a single matrix, `results` should also contain:
 *
 * - `permutation`: an integer dataset of length equal to the number of cells,
 *   describing the permutation to be applied to the per-gene results to recover the original row order.
 *
 * @v1_1{\[**since version 1.1**\] If there are multiple matrices, `results` should instead contain:}
 *
 * - @v1_1{`indices`: an integer dataset containing the row index of each feature in the intersection.
 *   For each entry, the gene is defined as the indexed row in the first sample _without permutation_.
 *   The `indices` are parallel to the per-gene results.}
 * 
 * All steps that generate per-gene results should use `permutation` or `indices` to identify the genes corresponding to the statistics.
 * See the documentation for the functions listed below:
 *
 * - `feature_selection::validate()`
 * - `marker_detection::validate()`
 * - `custom_selections::validate()`
 * 
 * @param handle An open HDF5 file handle.
 * @param version Version of the state file.
 *
 * @return Boolean indicating whether downstream analyses need to block on the sample.
 * If the format is invalid, an error is raised instead.
 */
inline bool validate(const H5::Group& handle, bool embedded = true, int version = 1001000) {
    auto ihandle = utils::check_and_open_group(handle, "inputs");

    std::pair<bool, bool> blocked;
    try {
        blocked = validate_parameters(ihandle, embedded, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'inputs'");
    }

    try {
        validate_results(ihandle, blocked.first, version);
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'inputs'");
    }

    return blocked.second;
}

}

}

#endif
