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
        EXPECT_NO_THROW(kanaval::inputs::validate(handle));
    }

    // Trying with HDF5.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_single_matrix(handle, "10X");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::inputs::validate(handle));
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
