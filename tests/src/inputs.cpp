#include <gtest/gtest.h>
#include "kanaval/inputs.hpp"
#include "utils.h"
#include <iostream>

void add_single_matrix(H5::H5File& handle, std::string mode = "MatrixMarket") {
    auto ihandle = handle.createGroup("inputs");

    auto phandle = ihandle.createGroup("parameters");
    quick_write_dataset(phandle, "format", mode);
    auto fihandle = phandle.createGroup("files");

    if (mode == "MatrixMarket") {
        auto fhandle0 = fihandle.createGroup("0");
        quick_write_dataset(fhandle0, "type", "mtx");
        quick_write_dataset(fhandle0, "name", "foo.mtx");
        quick_write_dataset(fhandle0, "size", 1);
        quick_write_dataset(fhandle0, "offset", 0);

        auto fhandle1 = fihandle.createGroup("1");
        quick_write_dataset(fhandle1, "type", "genes");
        quick_write_dataset(fhandle1, "name", "genes.tsv");
        quick_write_dataset(fhandle1, "size", 2);
        quick_write_dataset(fhandle1, "offset", 1);
    } else {
        auto fhandle = fihandle.createGroup("0");
        quick_write_dataset(fhandle, "type", "h5");
        quick_write_dataset(fhandle, "name", "foo.h5");
        quick_write_dataset(fhandle, "size", 1);
        quick_write_dataset(fhandle, "offset", 0);
    }

    auto rhandle = ihandle.createGroup("results");
    quick_write_dataset(rhandle, "dimensions", std::vector<int>{1000, 100});

    std::vector<int> permutation(1000);
    std::iota(permutation.rbegin(), permutation.rend(), 0); // reversed... outta control, bruh.
    quick_write_dataset(rhandle, "permutation", permutation);

    return;
}

TEST(SingleInputs, AllOK) {
    const std::string path = "TEST_inputs.h5";

    // Trying with MatrixMarket
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_FALSE(kanaval::inputs::validate(handle));
    }

    // Trying with HDF5.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle, "10X");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_FALSE(kanaval::inputs::validate(handle));
    }

    // Trying with a sample factor.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle, "H5AD");
        auto phandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(phandle, "sample_factor", "FOO");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_TRUE(kanaval::inputs::validate(handle));
    }

    // Trying with links.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle, "H5AD");

        auto phandle = handle.openGroup("inputs/parameters");
        auto zerohandle = phandle.openGroup("files/0");
        zerohandle.unlink("size");
        zerohandle.unlink("offset");
        quick_write_dataset(zerohandle, "id", "FOO");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_FALSE(kanaval::inputs::validate(handle, false));
    }
}

void quick_input_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::inputs::validate(handle);
    }, msg);
}

TEST(SingleInputs, ParametersFail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/parameters/format");
    }
    quick_input_throw(path, "'format' dataset");

    // How does it handle an older version?
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        auto phandle = handle.openGroup("inputs/parameters");
        phandle.unlink("format");
        quick_write_dataset(phandle, "format", std::vector<std::string>{ "FOO", "BAR" });
    }
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::inputs::validate(handle, true, 1000000);
    }, "version 1.0");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/parameters/files/0");
    }
    quick_input_throw(path, "file 0");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/parameters/files/0/offset");
    }
    quick_input_throw(path, "'offset' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        auto fihandle = handle.openGroup("inputs/parameters/files");

        auto fhandle2 = fihandle.createGroup("2");
        quick_write_dataset(fhandle2, "type", "annotation");
        quick_write_dataset(fhandle2, "name", "anno.tsv");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "not contiguous");
}

TEST(SingleInputs, MatrixMarketFail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        auto fihandle = handle.openGroup("inputs/parameters/files");

        auto fhandle2 = fihandle.createGroup("2");
        quick_write_dataset(fhandle2, "type", "mtx");
        quick_write_dataset(fhandle2, "name", "blah.mtx");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "expected exactly one 'mtx' file");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        auto fihandle = handle.openGroup("inputs/parameters/files");

        auto fhandle2 = fihandle.createGroup("2");
        quick_write_dataset(fhandle2, "type", "whee");
        quick_write_dataset(fhandle2, "name", "blah.mtx");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "unknown file type 'whee'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        auto fihandle = handle.openGroup("inputs/parameters/files");

        auto fhandle2 = fihandle.createGroup("2");
        quick_write_dataset(fhandle2, "type", "genes");
        quick_write_dataset(fhandle2, "name", "genes.tsv");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "expected no more than one 'genes' file");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        auto fihandle = handle.openGroup("inputs/parameters/files");

        for (int i = 2; i < 4; ++i) {
            auto fhandle = fihandle.createGroup(std::to_string(i));
            quick_write_dataset(fhandle, "type", "annotation");
            quick_write_dataset(fhandle, "name", "annotation.tsv");
            quick_write_dataset(fhandle, "size", 2);
            quick_write_dataset(fhandle, "offset", 1);
        }
    }
    quick_input_throw(path, "expected no more than one 'annotation' file");
}

TEST(SingleInputs, HDF5Fail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle, "10X");
        auto fihandle = handle.openGroup("inputs/parameters/files");

        auto fhandle2 = fihandle.createGroup("1");
        quick_write_dataset(fhandle2, "type", "h5");
        quick_write_dataset(fhandle2, "name", "blah.h5");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "expected exactly one 'h5' file");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle, "H5AD");
        auto fihandle = handle.openGroup("inputs/parameters/files");

        auto fhandle2 = fihandle.createGroup("1");
        quick_write_dataset(fhandle2, "type", "h5");
        quick_write_dataset(fhandle2, "name", "blah.h5");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "expected exactly one 'h5' file");
}

TEST(SingleInputs, ResultsFail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/results");
    }
    quick_input_throw(path, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/results/dimensions");
        auto rhandle = handle.openGroup("inputs/results");
        quick_write_dataset(rhandle, "dimensions", std::vector<int>{1,2,3});
    }
    quick_input_throw(path, "dataset of length 2");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/results/permutation");
        auto rhandle = handle.openGroup("inputs/results");
        quick_write_dataset(rhandle, "permutation", std::vector<int>{1,2,3});
    }
    quick_input_throw(path, "length equal to the number of genes");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/results/permutation");
        auto rhandle = handle.openGroup("inputs/results");
        quick_write_dataset(rhandle, "permutation", std::vector<int>(1000, -1));
    }
    quick_input_throw(path, "out-of-range");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle);
        handle.unlink("inputs/results/permutation");
        auto rhandle = handle.openGroup("inputs/results");
        quick_write_dataset(rhandle, "permutation", std::vector<int>(1000));
    }
    quick_input_throw(path, "duplicated");
}

void add_multiple_matrices(H5::H5File& handle) {
    auto ihandle = handle.createGroup("inputs");

    auto phandle = ihandle.createGroup("parameters");
    quick_write_dataset(phandle, "format", std::vector<std::string>{ "10X", "MatrixMarket" });
    quick_write_dataset(phandle, "sample_groups", std::vector<int>{ 1, 2 });
    quick_write_dataset(phandle, "sample_names", std::vector<std::string>{ "A", "B" });
    auto fihandle = phandle.createGroup("files");
    
    auto fhandle = fihandle.createGroup("0");
    quick_write_dataset(fhandle, "type", "h5");
    quick_write_dataset(fhandle, "name", "foo.h5");
    quick_write_dataset(fhandle, "size", 1);
    quick_write_dataset(fhandle, "offset", 3);

    auto fhandle0 = fihandle.createGroup("1");
    quick_write_dataset(fhandle0, "type", "mtx");
    quick_write_dataset(fhandle0, "name", "foo.mtx");
    quick_write_dataset(fhandle0, "size", 1);
    quick_write_dataset(fhandle0, "offset", 0);

    auto fhandle1 = fihandle.createGroup("2");
    quick_write_dataset(fhandle1, "type", "genes");
    quick_write_dataset(fhandle1, "name", "genes.tsv");
    quick_write_dataset(fhandle1, "size", 2);
    quick_write_dataset(fhandle1, "offset", 1);

    auto rhandle = ihandle.createGroup("results");
    quick_write_dataset(rhandle, "dimensions", std::vector<int>{500, 100});

    std::vector<int> indices(500);
    std::iota(indices.begin(), indices.end(), 0); 
    quick_write_dataset(rhandle, "indices", indices);

    return;
}

TEST(MultipleInputs, AllOK) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_TRUE(kanaval::inputs::validate(handle));
    }
}

TEST(MultipleInputs, ParametersFailed) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
        handle.unlink("inputs/parameters/sample_names");
    }
    quick_input_throw(path, "sample_names");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
        handle.unlink("inputs/parameters/sample_names");
        auto phandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(phandle, "sample_names", std::vector<std::string>{ "A" });
    }
    quick_input_throw(path, "'sample_names' and 'format'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
        handle.unlink("inputs/parameters/sample_names");
        auto phandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(phandle, "sample_names", std::vector<std::string>{ "A", "A" });
    }
    quick_input_throw(path, "duplicated sample name");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
        handle.unlink("inputs/parameters/sample_groups");
        auto phandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(phandle, "sample_groups", std::vector<int>{ 1 });
    }
    quick_input_throw(path, "'sample_groups' and 'format'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
        handle.unlink("inputs/parameters/sample_groups");
        auto phandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(phandle, "sample_groups", std::vector<int>{ 1, 1 });
    }
    quick_input_throw(path, "sum of 'sample_groups'");
}

TEST(MultipleInputs, ResultsFail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_multiple_matrices(handle);
        handle.unlink("inputs/results");
    }
    quick_input_throw(path, "'results' group");
}
