#include <gtest/gtest.h>
#include "kanaval/umap.hpp"
#include "utils.h"
#include <iostream>

void add_umap(H5::H5File& handle, int num_cells) {
    auto qhandle = handle.createGroup("umap");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "num_neighbors", 15);
    quick_write_dataset(phandle, "num_epochs", 1000);
    quick_write_dataset(phandle, "min_dist", 0.1);
    quick_write_dataset(phandle, "animate", 1);

    auto rhandle = qhandle.createGroup("results");
    quick_write_dataset(rhandle, "x", std::vector<double>(num_cells));
    quick_write_dataset(rhandle, "y", std::vector<double>(num_cells));

    return;
}

TEST(UMAP, AllOK) {
    const std::string path = "TEST_umap.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::umap::validate(handle, 1000));
    }
}

void quick_umap_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::umap::validate(handle, num_cells);
    }, msg);
}

TEST(UMAP, ParametersFailed) {
    const std::string path = "TEST_umap.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
        handle.unlink("umap/parameters");
    }
    quick_umap_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
        handle.unlink("umap/parameters/num_neighbors");
    }
    quick_umap_throw(path, 1000, "'num_neighbors' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
        handle.unlink("umap/parameters/num_neighbors");
    
        auto phandle = handle.openGroup("umap/parameters");
        quick_write_dataset(phandle, "num_neighbors", -1);
    }
    quick_umap_throw(path, 1000, "should be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
        handle.unlink("umap/parameters/num_epochs");
    
        auto phandle = handle.openGroup("umap/parameters");
        quick_write_dataset(phandle, "num_epochs", -1);
    }
    quick_umap_throw(path, 1000, "should be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
        handle.unlink("umap/parameters/min_dist");
    
        auto phandle = handle.openGroup("umap/parameters");
        quick_write_dataset(phandle, "min_dist", -1.0);
    }
    quick_umap_throw(path, 1000, "should be positive");
}

TEST(UMAP, ResultsFailed) {
    const std::string path = "TEST_umap.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
        handle.unlink("umap/results");
    }
    quick_umap_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_umap(handle, 1000);
    }
    quick_umap_throw(path, 500, "'x' dataset");
}
