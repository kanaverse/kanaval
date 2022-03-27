#include <gtest/gtest.h>
#include "kanaval/normalization.hpp"
#include "utils.h"
#include <iostream>

void add_normalization(H5::H5File& handle) {
    auto qhandle = handle.createGroup("normalization");
    qhandle.createGroup("parameters");
    qhandle.createGroup("results");
    return;
}

TEST(Normalization, AllOK) {
    const std::string path = "TEST_normalization.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_normalization(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::normalization::validate(handle));
    }
}

void quick_norm_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::normalization::validate(handle);
    }, msg);
}

TEST(Normalization, ParametersFailed) {
    const std::string path = "TEST_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_normalization(handle);
        handle.unlink("normalization/parameters");
    }
    quick_norm_throw(path, "'parameters' group");
}

TEST(Normalization, ResultsFailed) {
    const std::string path = "TEST_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_normalization(handle);
        handle.unlink("normalization/results");
    }
    quick_norm_throw(path, "'results' group");
}
