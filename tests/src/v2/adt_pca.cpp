#include <gtest/gtest.h>
#include "kanaval/v2/adt_pca.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_adt_pca(H5::H5File& handle, int num_pcs, int num_cells) {
    auto qhandle = handle.createGroup("adt_pca");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "num_pcs", num_pcs);
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

}

TEST(AdtPcaV2, AllOK) {
    const std::string path = "TEST_adt_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_pca(handle, 1000, true, latest));
    }

    // Checking with regression.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("adt_pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "regress");

        auto rhandle = handle.openGroup("adt_pca/results");
        auto dhandle = rhandle.openDataSet("pcs");
        H5::DataSpace space = dhandle.getSpace();
        rhandle.createDataSet("corrected", H5::PredType::NATIVE_DOUBLE, space);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_pca(handle, 1000, true, latest));
    }
}

TEST(AdtPcaV2, NoOp) {
    const std::string path = "TEST_adt_pca.h5";

    // No-ops if the version is too early.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_pca(handle, 100, true, 1000000));
    }

    // Skipping results if ADTs are not involved.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
        handle.unlink("adt_pca/results");
        handle.createGroup("adt_pca/results");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_pca(handle, 100, false, latest));
    }
}

static void quick_adt_pca_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_adt_pca(handle, num_cells, true, latest);
    }, msg);
}

TEST(AdtPcaV2, ParametersFailed) {
    const std::string path = "TEST_adt_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
        handle.unlink("adt_pca/parameters");
    }
    quick_adt_pca_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
        handle.unlink("adt_pca/parameters/num_pcs");
    
        auto phandle = handle.openGroup("adt_pca/parameters");
        quick_write_dataset(phandle, "num_pcs", -1);
    }
    quick_adt_pca_throw(path, 1000, "must be positive");
}

TEST(AdtPcaV2, ResultsFailed) {
    const std::string path = "TEST_adt_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
        handle.unlink("adt_pca/results");
    }
    quick_adt_pca_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
    }
    quick_adt_pca_throw(path, 500, "'pcs' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("adt_pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "RANDOM");
    }
    quick_adt_pca_throw(path, 1000, "unrecognized value 'RANDOM'");
}
