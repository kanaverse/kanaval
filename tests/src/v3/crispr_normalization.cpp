#include <gtest/gtest.h>
#include "kanaval/v3/crispr_normalization.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_crispr_normalization(H5::H5File& handle) {
    auto qhandle = handle.createGroup("crispr_normalization");
    qhandle.createGroup("parameters");
    qhandle.createGroup("results");
    return;
}

}

TEST(CrisprNormalizationV3, AllOK) {
    const std::string path = "TEST_crispr_normalization.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_normalization(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_crispr_normalization(handle));
    }
}

static void quick_norm_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_crispr_normalization(handle);
    }, msg);
}

TEST(CrisprNormalizationV3, ParametersFailed) {
    const std::string path = "TEST_crispr_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_normalization(handle);
        handle.unlink("crispr_normalization/parameters");
    }
    quick_norm_throw(path, "'parameters' group");
}

TEST(CrisprNormalizationV3, ResultsFailed) {
    const std::string path = "TEST_crispr_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_crispr_normalization(handle);
        handle.unlink("crispr_normalization/results");
    }
    quick_norm_throw(path, "'results' group");
}
