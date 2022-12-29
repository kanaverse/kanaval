#include <gtest/gtest.h>
#include "kanaval/v2/adt_quality_control.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_adt_quality_control(H5::H5File& handle, int num_cells, int num_samples, int lost = 10) {
    auto qhandle = handle.createGroup("adt_quality_control");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "igg_prefix", "foobar");
    quick_write_dataset(phandle, "nmads", 3.0);
    quick_write_dataset(phandle, "min_detected_drop", 0.1);
    quick_write_dataset(phandle, "skip", 0);

    auto rhandle = qhandle.createGroup("results");

    auto mhandle = rhandle.createGroup("metrics");
    quick_write_dataset(mhandle, "sums", std::vector<double>(num_cells));
    quick_write_dataset(mhandle, "detected", std::vector<int>(num_cells));
    quick_write_dataset(mhandle, "igg_total", std::vector<double>(num_cells));

    auto thandle = rhandle.createGroup("thresholds");
    quick_write_dataset(thandle, "detected", std::vector<double>(num_samples));
    quick_write_dataset(thandle, "igg_total", std::vector<double>(num_samples));

    std::vector<int> discard(num_cells);
    std::fill(discard.begin(), discard.begin() + lost, 1);
    quick_write_dataset(rhandle, "discards", discard);
    return;
}

}

TEST(AdtQualityControlV2, AllOK) {
    const std::string path = "TEST_adt_quality_control.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 100, 1, true, latest);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, 90);
    }

    // Works with more batches.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 200, 2);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 200, 2, true, latest);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, 190);
    }
}

TEST(AdtQualityControlV2, NoOp) {
    const std::string path = "TEST_adt_quality_control.h5";

    // No-ops if the version is too early.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 100, 1, true, 1000000);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, -1);
    }

    // Skipping results if ADTs are not involved.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/results");
        handle.createGroup("adt_quality_control/results");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 200, 2, false, latest);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, -1);
    }
}

static void quick_adt_qc_throw(const std::string& path, int num_cells, int num_samples, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_adt_quality_control(handle, num_cells, num_samples, true, latest);
    }, msg);
}

TEST(AdtQualityControlV2, ParametersFailed) {
    const std::string path = "TEST_adt_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters");
    }
    quick_adt_qc_throw(path, 100, 1, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters/igg_prefix");
    }
    quick_adt_qc_throw(path, 100, 1, "'igg_prefix'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters/nmads");
        quick_write_dataset(handle, "adt_quality_control/parameters/nmads", -1.0);
    }
    quick_adt_qc_throw(path, 100, 1, "non-negative");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters/min_detected_drop");
        quick_write_dataset(handle, "adt_quality_control/parameters/min_detected_drop", -1.0);
    }
    quick_adt_qc_throw(path, 100, 1, "should lie in [0, 1)");
}

TEST(AdtQualityControlV2, ResultsFailed) {
    const std::string path = "TEST_adt_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/results");
    }
    quick_adt_qc_throw(path, 100, 1, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 200, 1);
    }
    quick_adt_qc_throw(path, 100, 1, "failed to retrieve metrics");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 2);
    }
    quick_adt_qc_throw(path, 100, 1, "failed to retrieve thresholds");
}

TEST(AdtQualityControlV2, SkipOK) {
    const std::string path = "TEST_adt_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 100, 1);
        auto phandle = handle.openGroup("adt_quality_control/parameters");
        phandle.unlink("skip");
        quick_write_dataset(phandle, "skip", 1);
    }

    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 100, 1, true, latest);
        EXPECT_TRUE(out.first);
        EXPECT_EQ(out.second, 100);
    }

    // Tolerates loss of all the values.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        handle.unlink("adt_quality_control/results/metrics");
        handle.unlink("adt_quality_control/results/thresholds");
        handle.unlink("adt_quality_control/results/discards");
    }

    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 100, 1, true, latest);
        EXPECT_TRUE(out.first);
        EXPECT_EQ(out.second, 100);
    }

    // Throws an error correctly if values are present but weird.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto rhandle = handle.openGroup("adt_quality_control/results");
        quick_write_dataset(rhandle, "discards", 987654321);
    }

    quick_adt_qc_throw(path, 100, 1, "discards");
}

TEST(AdtQualityControlV2, LegacySkipOK) {
    const std::string path = "TEST_adt_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_adt_quality_control(handle, 200, 2);
        handle.unlink("adt_quality_control/parameters/skip");
    }

    quick_adt_qc_throw(path, 200, 2, "skip");

    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_adt_quality_control(handle, 200, 2, true, 2000000);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, 190);
    }
}
