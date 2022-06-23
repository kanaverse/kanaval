#include <gtest/gtest.h>
#include "kanaval/pca.hpp"
#include "utils.h"
#include <iostream>

void add_pca(H5::H5File& handle, int num_pcs, int num_cells) {
    auto qhandle = handle.createGroup("pca");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "num_pcs", num_pcs);
    quick_write_dataset(phandle, "num_hvgs", 4000);
    quick_write_dataset(phandle, "block_method", "none");

    auto rhandle = qhandle.createGroup("results");

    H5::DataSpace space;
    std::vector<hsize_t> dims(2);
    dims[0] = num_cells;
    dims[1] = num_pcs;
    space.setExtentSimple(2, dims.data());
    rhandle.createDataSet("pcs", H5::PredType::NATIVE_DOUBLE, space);

    quick_write_dataset(rhandle, "var_exp", std::vector<double>(num_pcs));

    return;
}

TEST(PCA, AllOK) {
    const std::string path = "TEST_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::pca::validate(handle, 1000));
    }

    // Checking with MNN.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "mnn");

        auto rhandle = handle.openGroup("pca/results");
        auto dhandle = rhandle.openDataSet("pcs");
        H5::DataSpace space = dhandle.getSpace();
        rhandle.createDataSet("corrected", H5::PredType::NATIVE_DOUBLE, space);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::pca::validate(handle, 1000));
    }
}

void quick_pca_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::pca::validate(handle, num_cells);
    }, msg);
}

TEST(PCA, ParametersFailed) {
    const std::string path = "TEST_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        handle.unlink("pca/parameters");
    }
    quick_pca_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        handle.unlink("pca/parameters/num_hvgs");
    }
    quick_pca_throw(path, 1000, "'num_hvgs' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        handle.unlink("pca/parameters/num_hvgs");
    
        auto phandle = handle.openGroup("pca/parameters");
        quick_write_dataset(phandle, "num_hvgs", -1);
    }
    quick_pca_throw(path, 1000, "must be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        handle.unlink("pca/parameters/num_pcs");
    
        auto phandle = handle.openGroup("pca/parameters");
        quick_write_dataset(phandle, "num_pcs", -1);
    }
    quick_pca_throw(path, 1000, "must be positive");
}

TEST(PCA, ResultsFailed) {
    const std::string path = "TEST_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        handle.unlink("pca/results");
    }
    quick_pca_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
    }
    quick_pca_throw(path, 500, "'pcs' dataset");

    // Checking with MNN.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "mnn");
    }
    quick_pca_throw(path, 1000, "'corrected' dataset");

    quick_throw([&]() -> void { // Different error at a later version.
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::pca::validate(handle, 1000, 2000000);
    }, "unrecognized value 'mnn'");
}
