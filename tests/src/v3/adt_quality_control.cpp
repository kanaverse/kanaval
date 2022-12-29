#include <gtest/gtest.h>
#include "kanaval/v3/adt_quality_control.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_adt_quality_control(H5::H5File& handle, int num_cells, int num_samples, int lost = 10) {
    auto qhandle = handle.createGroup("adt_quality_control");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "igg_prefix", "foobar");
    quick_write_dataset(phandle, "nmads", 3.0);
    quick_write_dataset(phandle, "min_detected_drop", 0.1);

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

TEST(AdtQualityControlV3, AllOK) {
    const std::string path = "TEST_adt_quality_control.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v3::validate_adt_quality_control(handle, 100, 1, true, latest);
        EXPECT_EQ(out, 90);
    }

    // Works with more batches.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 200, 2);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v3::validate_adt_quality_control(handle, 200, 2, true, latest);
        EXPECT_EQ(out, 190);
    }
}

TEST(AdtQualityControlV3, NotInUse) {
    const std::string path = "TEST_adt_quality_control.h5";

    // Skipping results if ADTs are not involved.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/results");
        handle.createGroup("adt_quality_control/results");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v3::validate_adt_quality_control(handle, 200, 2, false, latest);
        EXPECT_EQ(out, -1);
    }
}

static void quick_adt_qc_throw(const std::string& path, int num_cells, int num_samples, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_adt_quality_control(handle, num_cells, num_samples, true, latest);
    }, msg);
}

TEST(AdtQualityControlV3, ParametersFailed) {
    const std::string path = "TEST_adt_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters");
    }
    quick_adt_qc_throw(path, 100, 1, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters/igg_prefix");
    }
    quick_adt_qc_throw(path, 100, 1, "'igg_prefix'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters/nmads");
        quick_write_dataset(handle, "adt_quality_control/parameters/nmads", -1.0);
    }
    quick_adt_qc_throw(path, 100, 1, "non-negative");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/parameters/min_detected_drop");
        quick_write_dataset(handle, "adt_quality_control/parameters/min_detected_drop", -1.0);
    }
    quick_adt_qc_throw(path, 100, 1, "should lie in [0, 1)");
}

TEST(AdtQualityControlV3, ResultsFailed) {
    const std::string path = "TEST_adt_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 1);
        handle.unlink("adt_quality_control/results");
    }
    quick_adt_qc_throw(path, 100, 1, "'results' group");

    // Mismatch in the number of cells.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 200, 1);
    }
    quick_adt_qc_throw(path, 100, 1, "failed to retrieve metrics");

    // Mismatch in the number of blocks.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_adt_quality_control(handle, 100, 2);
    }
    quick_adt_qc_throw(path, 100, 1, "failed to retrieve thresholds");
}
