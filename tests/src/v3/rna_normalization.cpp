#include <gtest/gtest.h>
#include "kanaval/v3/rna_normalization.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_rna_normalization(H5::H5File& handle) {
    auto qhandle = handle.createGroup("rna_normalization");
    qhandle.createGroup("parameters");
    qhandle.createGroup("results");
    return;
}

}

TEST(RnaNormalizationV3, AllOK) {
    const std::string path = "TEST_rna_normalization.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_normalization(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_rna_normalization(handle));
    }
}

static void quick_norm_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_rna_normalization(handle);
    }, msg);
}

TEST(RnaNormalizationV3, ParametersFailed) {
    const std::string path = "TEST_rna_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_normalization(handle);
        handle.unlink("rna_normalization/parameters");
    }
    quick_norm_throw(path, "'parameters' group");
}

TEST(RnaNormalizationV3, ResultsFailed) {
    const std::string path = "TEST_rna_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_normalization(handle);
        handle.unlink("rna_normalization/results");
    }
    quick_norm_throw(path, "'results' group");
}
