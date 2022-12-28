#include <gtest/gtest.h>
#include "kanaval/v2/normalization.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_normalization(H5::H5File& handle) {
    auto qhandle = handle.createGroup("normalization");
    qhandle.createGroup("parameters");
    qhandle.createGroup("results");
    return;
}

}

TEST(NormalizationV2, AllOK) {
    const std::string path = "TEST_normalization.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_normalization(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_normalization(handle));
    }
}

static void quick_norm_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_normalization(handle);
    }, msg);
}

TEST(NormalizationV2, ParametersFailed) {
    const std::string path = "TEST_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_normalization(handle);
        handle.unlink("normalization/parameters");
    }
    quick_norm_throw(path, "'parameters' group");
}

TEST(NormalizationV2, ResultsFailed) {
    const std::string path = "TEST_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_normalization(handle);
        handle.unlink("normalization/results");
    }
    quick_norm_throw(path, "'results' group");
}
