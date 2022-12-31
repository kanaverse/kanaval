#include <gtest/gtest.h>
#include "kanaval/v3/crispr_pca.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_crispr_pca(H5::H5File& handle, int num_pcs, int num_cells) {
    auto qhandle = handle.createGroup("crispr_pca");

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

TEST(CrisprPcaV3, AllOK) {
    const std::string path = "TEST_crispr_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_crispr_pca(handle, 1000, true, latest));
    }

    // Checking with other block_method.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("crispr_pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "regress");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_crispr_pca(handle, 1000, true, latest));
    }
}

static void quick_crispr_pca_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_crispr_pca(handle, num_cells, true, latest);
    }, msg);
}

TEST(CrisprPcaV3, ParametersFailed) {
    const std::string path = "TEST_crispr_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
        handle.unlink("crispr_pca/parameters");
    }
    quick_crispr_pca_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
        handle.unlink("crispr_pca/parameters/num_pcs");
    
        auto phandle = handle.openGroup("crispr_pca/parameters");
        quick_write_dataset(phandle, "num_pcs", -1);
    }
    quick_crispr_pca_throw(path, 1000, "must be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("crispr_pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "foo");
    }
    quick_crispr_pca_throw(path, 1000, "'block_method'");
}

TEST(CrisprPcaV3, ResultsFailed) {
    const std::string path = "TEST_crispr_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
        handle.unlink("crispr_pca/results");
    }
    quick_crispr_pca_throw(path, 1000, "'results' group");

    // Wrong dimensions.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
    }
    quick_crispr_pca_throw(path, 500, "'pcs' dataset");
}

TEST(CrisprPcaV3, NotInUse) {
    const std::string path = "TEST_crispr_pca.h5";

    // Not okay if CRISPR is available...
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_pca(handle, 10, 1000);
        auto rhandle = handle.openGroup("crispr_pca/results");
        rhandle.unlink("pcs");
    }
    quick_crispr_pca_throw(path, 500, "'pcs' dataset");

    // Okay if it's missing.
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v3::validate_crispr_pca(handle, 1000, false, latest), -1);
    }
}
