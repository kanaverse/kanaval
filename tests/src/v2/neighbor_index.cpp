#include <gtest/gtest.h>
#include "kanaval/v2/neighbor_index.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_neighbor_index(H5::H5File& handle) {
    auto qhandle = handle.createGroup("neighbor_index");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "approximate", 1);
    qhandle.createGroup("results");
    return;
}

}

TEST(NeighborIndex, AllOK) {
    const std::string path = "TEST_neighbor_index.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_neighbor_index(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_neighbor_index(handle));
    }
}

void quick_index_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_neighbor_index(handle);
    }, msg);
}

TEST(NeighborIndex, ParametersFailed) {
    const std::string path = "TEST_neighbor_index.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_neighbor_index(handle);
        handle.unlink("neighbor_index/parameters");
    }
    quick_index_throw(path, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_neighbor_index(handle);
        handle.unlink("neighbor_index/parameters/approximate");
    }
    quick_index_throw(path, "'approximate' dataset");
}

TEST(NeighborIndex, ResultsFailed) {
    const std::string path = "TEST_neighbor_index.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_neighbor_index(handle);
        handle.unlink("neighbor_index/results");
    }
    quick_index_throw(path, "'results' group");
}
