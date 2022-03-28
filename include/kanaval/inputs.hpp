#ifndef KANAVAL_INPUTS_HPP
#define KANAVAL_INPUTS_HPP

namespace kanaval {

namespace inputs {

inline bool validate_parameters(const H5::Group& handle, bool embedded) {
    auto phandle = utils::check_and_open_group(handle, "parameters");

    // Formats can either be a scalar... or not.
    std::vector<std::string> formats;
    bool multifile = false;
    {
        auto fhandle = utils::check_and_open_dataset(handle, "format", H5T_STRING);
        auto fspace = fhandle.getSpace();
        if (fspace.getSimpleExtentNdims() == 0) {
            formats.push_back(load_string(fhandle));
        } else {
            multifile = true;
            formats = load_string_vector(fhandle);
        }
    }

    auto fihandle = utils::check_and_open_group(handle, "files");
    auto nfiles = fhandle.getNumObjs();

    // Checking the runs.
    std::vector<int> runs;
    if (multifile) {
        runs = utils::load_integer_vector(handle, "sample_groups");
        if (runs.size() != formats.size()) {
            throw std::runtime_error("'sample_groups' and 'format' should have the same length");
        }

        int total_files = std::accumulate(runs.begin(), runs.end(), 0);
        if (total_files != static_cast<int>(nfiles)) {
            throw std::runtime_error("sum of 'sample_groups' is not equal to the length of 'files'");
        }

        // Checking that everyone has unique groups.
        auto names = utils::load_string_vector(handle, "sample_names");
        if (runs.size() != formats.size()) {
            throw std::runtime_error("'names' and 'format' should have the same length");
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
    std::vector<std::pair<size_t, size_t> > bytes;
    for (size_t r = 0; r < runs.size(); ++r) {
        auto curf = formats[r];
        std::vector<std::string> types;

        for (int s = 0; s < runs[r]; ++s) {
            std::string current = std::to_string(s + sofar);
            auto curfihandle = utils::check_and_open_group(fihandle, current);

            utils::check_and_open_dataset(current, "name", H5T_STRING, {});
            types.push_back(utils::load_string(current, "type"));

            if (embedded) {
                bytes.emplace_back(
                    utils::load_integer<hsize_t>(current, "start"),
                    utils::load_integer<hsize_t>(current, "offset")
                );
            } else {
                utils::check_and_open_dataset(current, "id", H5T_STRING, {});
            }
        }

        if (curf == "MatrixMarket") {
            std::unordered_map<std::string> expected;
            expected["mtx"] = 0;
            expected["genes"] = 0;
            expected["annotation"] = 0;

            for (auto t : types) {
                auto it = expected.find(t);
                if (it == expected.end()) {
                    throw std::runtime_error("unknown file type '" + t + "' when format is 'MatrixMarket'");
                }
                ++(it->second)
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
            if (types.size() != 1 && types.front() != "h5") {
                throw std::runtime_error("expected exactly one 'h5' file when format is '10X'");
            }
        } else if (curf == "H5AD") {
            if (types.size() != 1 && types.front() != "h5") {
                throw std::runtime_error("expected exactly one 'h5' file when format is 'H5AD'");
            }
        }
    }

    // Checking if there's a batch variable.
    if (phandle.exists("sample_factor")) {
        utils::check_and_open_dataset(phandle, "sample_factor", H5T_STRING, {});
        return true;
    } else {
        return multifile;
    }
}

}

}
