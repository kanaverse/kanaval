#ifndef KANAVAL_INPUTS_V3_HPP
#define KANAVAL_INPUTS_V3_HPP

#include "H5Cpp.h"
#include "../utils.hpp"
#include <stdexcept>
#include <vector>
#include <string>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

namespace kanaval {

namespace v3 {

namespace inputs {

inline size_t check_datasets(const H5::Group& phandle, bool embedded) {
    auto dhandle = utils::check_and_open_group(phandle, "datasets");
    size_t ndatasets = dhandle.getNumObjs();
    int last = 0;
    std::unordered_set<std::string> used_names;

    for (size_t d = 0; d < ndatasets; ++d) {
        try {
            auto curdhandle = utils::check_and_open_group(dhandle, std::to_string(d));
            utils::load_string(curdhandle, "format");
            
            auto dname = utils::load_string(curdhandle, "name");
            if (used_names.find(dname) != used_names.end()) {
                throw std::runtime_error("detected duplicate dataset name '" + dname + "'");
            }
            used_names.insert(dname);

            // Checking the files.
            auto fhandle = utils::check_and_open_group(curdhandle, "files");
            size_t nfiles = fhandle.getNumObjs();
            for (size_t f = 0; f < nfiles; ++f) {
                try {
                    auto curfhandle = utils::check_and_open_group(fhandle, std::to_string(f));
                    utils::load_string(curfhandle, "type");
                    utils::load_string(curfhandle, "name");

                    if (embedded) {
                        int offset = utils::load_integer_scalar(curfhandle, "offset");
                        if (offset < 0) {
                            throw std::runtime_error("offset should be non-negative");
                        }
                        if (offset != last) {
                            throw std::runtime_error("byte range is not contiguous with previous file");
                        }

                        int size = utils::load_integer_scalar(curfhandle, "size");
                        if (size < 0) {
                            throw std::runtime_error("size should be non-negative");
                        }
                        last += size;
                    } else {
                        utils::load_string(curfhandle, "id");
                    }

                } catch (std::exception& e) {
                    throw utils::combine_errors(e, "failed to inspect details for file " + std::to_string(f)); 
                }
            }

            // Checking the options.
            if (curdhandle.exists("options")) {
                auto ohandle = utils::check_and_open_group(curdhandle, "options");
                try {
                    size_t noptions = ohandle.getNumObjs();
                    for (size_t o = 0; o < noptions; ++o) {
                        if (ohandle.childObjType(o) != H5O_TYPE_DATASET) {
                            std::string optname = ohandle.getObjnameByIdx(o);
                            throw std::runtime_error("option '" + optname + "' should be a dataset");
                        }
                    }
                } catch (std::exception& e) {
                    throw utils::combine_errors(e, "failed to load options");
                }
            }

        } catch (std::exception& e) {
            throw utils::combine_errors(e, "failed to check parameters for dataset " + std::to_string(d)); 
        }
    }

    return ndatasets;
}

inline int check_subset_cells(const H5::Group& subhandle) {
    int subset_limit = -1;

    if (subhandle.exists("cells")) {
        auto subcellhandle = utils::check_and_open_group(subhandle, "cells");
        if (subcellhandle.exists("indices")) {
            auto subidx = utils::load_integer_vector(subcellhandle, "indices");

            for (auto i : subidx) {
                if (i < 0) {
                    throw std::runtime_error("indices in 'subset/indices' should be non-negative");
                }
            }

            if (!utils::is_unique_and_sorted(subidx)) {
                throw std::runtime_error("indices in 'subset/indices' should be unique and sorted");
            }

            subset_limit = subidx.size();
        } else {
            utils::check_and_open_dataset(subcellhandle, "field", H5T_STRING, {});

            if (subcellhandle.exists("values")) {
                auto vhandle = utils::check_and_open_dataset(subcellhandle, "values", H5T_STRING);
                auto vdims = utils::load_dataset_dimensions(vhandle);
                if (vdims.size() != 1) {
                    throw std::runtime_error("'subset/values' should be a 1-dimensional string dataset");
                }

            } else {
                auto rahandle = utils::check_and_open_dataset(subcellhandle, "ranges", H5T_FLOAT);

                auto radims = utils::load_dataset_dimensions(rahandle);
                if (radims.size() != 2) {
                    throw std::runtime_error("'subset/ranges' should be a 2-dimensional float dataset");
                }
                if (radims[1] != 2) {
                    throw std::runtime_error("'subset/ranges' should have two columns");
                }

                std::vector<double> loaded(radims[0] * radims[1]);
                rahandle.read(loaded.data(), H5::PredType::NATIVE_DOUBLE);

                // Should be sorted. Each row is a [start, end) pair,
                // defining an interval to retain.  Intervals should be
                // non-overlapping and sorted, and HDF5 stores its matrices
                // as row-order, so start1 <= end1 <= start2 <= end2 <= ...
                if (!std::is_sorted(loaded.begin(), loaded.end())) {
                    throw std::runtime_error("'subset/ranges' should specify sorted, non-overlapping intervals");
                }
            }
        }
    }

    return subset_limit;
}

inline std::unordered_map<std::string, int> check_feature_identities(const H5::Group& rhandle) {
    std::unordered_map<std::string, int> num_features;
    auto ihandle = utils::check_and_open_group(rhandle, "feature_identities");

    auto fill_identities = [&](const std::string& key) -> void {
        if (ihandle.exists(key)) {
            auto khandle = utils::check_and_open_dataset(ihandle, key);
            auto indices = utils::load_integer_vector(khandle);

            std::sort(indices.begin(), indices.end());
            int last = -1;
            for (auto i : indices) {
                if (i <= last || i < 0) {
                    throw std::runtime_error("identities for modality '" + key + "' should be unique and non-negative");
                }
                last = i;
            }

            num_features[key] =  indices.size();
        }
    };

    try {
        fill_identities("RNA");
        fill_identities("ADT");
        fill_identities("CRISPR");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to load 'feature_identities'");
    }

    if (num_features.size() == 0) {
        throw std::runtime_error("'feature_identities' should contain at least one recognized modality");
    }

    return num_features;
}

inline void check_feature_names(const H5::Group& rhandle, const std::unordered_map<std::string, int>& num_features) {
    if (!rhandle.exists("feature_names")) {
        return;
    }

    auto nhandle = utils::check_and_open_group(rhandle, "feature_names");
    try {
        for (const auto& mod : num_features) {
            const auto& key = mod.first;
            if (!nhandle.exists(key)) {
                continue;
            }

            auto khandle = utils::check_and_open_dataset(nhandle, key);
            if (khandle.getDataType().getClass() != H5T_STRING) {
                throw std::runtime_error("names for modality '" + key + "' should be a string dataset");
            }

            auto kdims = utils::load_dataset_dimensions(khandle);
            if (kdims.size() != 1) {
                throw std::runtime_error("names for modality '" + key + "' should be a 1-dimensional string dataset");
            }
            if (kdims[0] != mod.second) {
                throw std::runtime_error("names for modality '" + key + "' should be of length equal to 'feature_identities/" + key + "'");
            }
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to load 'feature_names'");
    }
}

struct Details {
    int num_cells;
    int num_blocks;
    std::unordered_map<std::string, int> num_features;
};

}

inline inputs::Details validate_inputs(const H5::Group& handle, bool embedded, int version) {
    auto xhandle = utils::check_and_open_group(handle, "inputs");

    // Checking parameters.
    int subsetted = -1;
    int ndatasets;
    bool has_block = false;
    try {
        auto phandle = utils::check_and_open_group(xhandle, "parameters");
        ndatasets = inputs::check_datasets(phandle, embedded);

        if (phandle.exists("subset")) {
            auto subhandle = utils::check_and_open_group(phandle, "subset");
            subsetted = inputs::check_subset_cells(subhandle);
        }

        if (ndatasets == 1 && phandle.exists("block_factor")) {
            utils::load_string(phandle, "block_factor");
            has_block = true;
        }

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve parameters from 'inputs'");
    }

    // Checking results.
    inputs::Details output;
    try {
        auto rhandle = utils::check_and_open_group(xhandle, "results");

        // Checking the number of cells.
        output.num_cells = utils::load_integer_scalar(rhandle, "num_cells");
        if (output.num_cells <= 0) {
            throw std::runtime_error("number of cells should be a positive integer");
        }
        if (subsetted >= 0 && output.num_cells != subsetted) { 
            throw std::runtime_error("number of cells should be equal to the length of 'parameters/subset/indices'");
        }

        // Checking the number of blocks.
        output.num_blocks = utils::load_integer_scalar(rhandle, "num_blocks");
        if (output.num_blocks <= 0) {
            throw std::runtime_error("number of blocks should be a positive integer");
        }
        if (ndatasets > 1 || !has_block) {
            if (ndatasets != output.num_blocks) {
                throw std::runtime_error("number of blocks should be equal to the number of datasets");
            }
        }

        // Checking the identities.
        output.num_features = inputs::check_feature_identities(rhandle);

        // Checking the names.
        inputs::check_feature_names(rhandle, output.num_features);

    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'inputs'");
    }

    return output;
}

}

}

#endif
