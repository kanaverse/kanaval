#include <gtest/gtest.h>
#include "kanaval/v3/rna_pca.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_rna_pca(H5::H5File& handle, int num_pcs, int num_cells) {
    auto qhandle = handle.createGroup("rna_pca");

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

}

TEST(RnaPcaV3, AllOK) {
    const std::string path = "TEST_rna_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_rna_pca(handle, 1000, true, latest));
    }

    // Checking with other block_method.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("rna_pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "regress");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_rna_pca(handle, 1000, true, latest));
    }
}

static void quick_rna_pca_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_rna_pca(handle, num_cells, true, latest);
    }, msg);
}

TEST(RnaPcaV3, ParametersFailed) {
    const std::string path = "TEST_rna_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        handle.unlink("rna_pca/parameters");
    }
    quick_rna_pca_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        handle.unlink("rna_pca/parameters/num_hvgs");
    }
    quick_rna_pca_throw(path, 1000, "'num_hvgs' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        handle.unlink("rna_pca/parameters/num_hvgs");
    
        auto phandle = handle.openGroup("rna_pca/parameters");
        quick_write_dataset(phandle, "num_hvgs", -1);
    }
    quick_rna_pca_throw(path, 1000, "must be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        handle.unlink("rna_pca/parameters/num_pcs");
    
        auto phandle = handle.openGroup("rna_pca/parameters");
        quick_write_dataset(phandle, "num_pcs", -1);
    }
    quick_rna_pca_throw(path, 1000, "must be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        auto phandle = handle.openGroup("rna_pca/parameters");
        phandle.unlink("block_method");
        quick_write_dataset(phandle, "block_method", "foo");
    }
    quick_rna_pca_throw(path, 1000, "'block_method'");
}

TEST(RnaPcaV3, ResultsFailed) {
    const std::string path = "TEST_rna_pca.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        handle.unlink("rna_pca/results");
    }
    quick_rna_pca_throw(path, 1000, "'results' group");

    // Wrong dimensions.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
    }
    quick_rna_pca_throw(path, 500, "'pcs' dataset");
}

TEST(RnaPcaV3, NotInUse) {
    const std::string path = "TEST_rna_pca.h5";

    // Not okay if RNA is available...
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_pca(handle, 10, 1000);
        auto rhandle = handle.openGroup("rna_pca/results");
        rhandle.unlink("pcs");
    }
    quick_rna_pca_throw(path, 500, "'pcs' dataset");

    // Okay if it's missing.
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v3::validate_rna_pca(handle, 1000, false, latest), -1);
    }
}
