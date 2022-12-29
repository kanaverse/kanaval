#include <gtest/gtest.h>
#include "kanaval/v2/quality_control.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_quality_control(H5::H5File& handle, int num_cells, int num_samples, int lost = 10) {
    auto qhandle = handle.createGroup("quality_control");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "use_mito_default", int(0));
    quick_write_dataset(phandle, "mito_prefix", "foobar");
    quick_write_dataset(phandle, "nmads", 3.0);
    quick_write_dataset(phandle, "skip", 0);

    auto rhandle = qhandle.createGroup("results");

    auto mhandle = rhandle.createGroup("metrics");
    quick_write_dataset(mhandle, "sums", std::vector<double>(num_cells));
    quick_write_dataset(mhandle, "detected", std::vector<int>(num_cells));
    quick_write_dataset(mhandle, "proportion", std::vector<double>(num_cells));

    auto thandle = rhandle.createGroup("thresholds");
    quick_write_dataset(thandle, "sums", std::vector<double>(num_samples));
    quick_write_dataset(thandle, "detected", std::vector<double>(num_samples));
    quick_write_dataset(thandle, "proportion", std::vector<double>(num_samples));

    std::vector<int> discard(num_cells);
    std::fill(discard.begin(), discard.begin() + lost, 1);
    quick_write_dataset(rhandle, "discards", discard);
    return;
}

}

TEST(QualityControlV2, AllOK) {
    const std::string path = "TEST_quality_control.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_quality_control(handle, 100, 1, latest);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, 90);
    }

    // Works with more batches.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 200, 2);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_quality_control(handle, 200, 2, latest);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, 190);
    }
}

static void quick_qc_throw(const std::string& path, int num_cells, int num_samples, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_quality_control(handle, num_cells, num_samples, latest);
    }, msg);
}

TEST(QualityControlV2, ParametersFailed) {
    const std::string path = "TEST_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters");
    }
    quick_qc_throw(path, 100, 1, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/mito_prefix");
    }
    quick_qc_throw(path, 100, 1, "'mito_prefix'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/use_mito_default");
    }
    quick_qc_throw(path, 100, 1, "'use_mito_default'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/nmads");
    }
    quick_qc_throw(path, 100, 1, "'nmads'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/nmads");
        quick_write_dataset(handle, "quality_control/parameters/nmads", -1.0);
    }
    quick_qc_throw(path, 100, 1, "non-negative");
}

TEST(QualityControlV2, ResultsFailed) {
    const std::string path = "TEST_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/results");
    }
    quick_qc_throw(path, 100, 1, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 200, 1);
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve metrics");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 2);
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve thresholds");
}

TEST(QualityControlV2, SkipOK) {
    const std::string path = "TEST_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 100, 1);
        auto phandle = handle.openGroup("quality_control/parameters");
        phandle.unlink("skip");
        quick_write_dataset(phandle, "skip", 1);
    }

    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_quality_control(handle, 100, 1, latest);
        EXPECT_TRUE(out.first);
        EXPECT_EQ(out.second, 100);
    }

    // Tolerates loss of all the values.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        handle.unlink("quality_control/results/metrics");
        handle.unlink("quality_control/results/thresholds");
        handle.unlink("quality_control/results/discards");
    }

    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_quality_control(handle, 100, 1, latest);
        EXPECT_TRUE(out.first);
        EXPECT_EQ(out.second, 100);
    }

    // Throws an error correctly if values are present but weird.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto rhandle = handle.openGroup("quality_control/results");
        quick_write_dataset(rhandle, "discards", 987654321);
    }

    quick_qc_throw(path, 100, 1, "discards");
}

TEST(QualityControlV2, LegacySkipOK) {
    const std::string path = "TEST_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_quality_control(handle, 200, 2);
        handle.unlink("quality_control/parameters/skip");
    }

    quick_qc_throw(path, 200, 2, "skip");

    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v2::validate_quality_control(handle, 200, 2, 2000000);
        EXPECT_FALSE(out.first);
        EXPECT_EQ(out.second, 190);
    }
}


