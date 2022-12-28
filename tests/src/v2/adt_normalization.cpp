#include <gtest/gtest.h>
#include "kanaval/v2/adt_normalization.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_adt_normalization(H5::H5File& handle, int ncells) {
    auto qhandle = handle.createGroup("adt_normalization");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "num_clusters", 20);
    quick_write_dataset(phandle, "num_pcs", 24);

    auto rhandle = qhandle.createGroup("results");
    quick_write_dataset(rhandle, "size_factors", std::vector<double>(ncells));
    return;
}

}

TEST(AdtNormalizationV2, AllOK) {
    const std::string path = "TEST_adt_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_normalization(handle, 100, true, latest));
    }
}

TEST(AdtNormalizationV2, NoOp) {
    const std::string path = "TEST_adt_normalization.h5";

    // No-ops if the version is too early.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_normalization(handle, 100, true, 1000000));
    }

    // Skipping results if ADTs are not involved.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
        handle.unlink("adt_normalization/results");
        handle.createGroup("adt_normalization/results");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_adt_normalization(handle, 100, false, latest));
    }
}

void quick_adt_norm_throw(const std::string& path, size_t ncells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_adt_normalization(handle, ncells, true, latest);
    }, msg);
}

TEST(AdtNormalizationV2, ParametersFailed) {
    const std::string path = "TEST_adt_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
        handle.unlink("adt_normalization/parameters");
    }
    quick_adt_norm_throw(path, 100, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
        handle.unlink("adt_normalization/parameters/num_pcs");
    }
    quick_adt_norm_throw(path, 100, "num_pcs");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
        handle.unlink("adt_normalization/parameters/num_clusters");
        quick_write_dataset(handle, "adt_normalization/parameters/num_clusters", -1);
    }
    quick_adt_norm_throw(path, 100, "positive");
}

TEST(AdtNormalizationV2, ResultsFailed) {
    const std::string path = "TEST_adt_normalization.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
        handle.unlink("adt_normalization/results");
    }
    quick_adt_norm_throw(path, 100, "'results' group");
    
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_normalization(handle, 100);
    }
    quick_adt_norm_throw(path, 99, "size_factors");
}
