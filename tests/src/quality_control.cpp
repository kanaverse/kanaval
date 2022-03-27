#include <gtest/gtest.h>
#include "kanaval/quality_control.hpp"
#include "utils.h"
#include <iostream>

void add_quality_control(H5::H5File& handle, int num_cells, int num_batches) {
    auto qhandle = handle.createGroup("quality_control");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "use_mito_default", int(0));
    quick_write_dataset(phandle, "mito_prefix", "foobar");
    quick_write_dataset(phandle, "nmads", 3.0);

    auto rhandle = qhandle.createGroup("results");

    auto mhandle = rhandle.createGroup("metrics");
    quick_write_dataset(mhandle, "sums", std::vector<double>(num_cells));
    quick_write_dataset(mhandle, "detected", std::vector<int>(num_cells));
    quick_write_dataset(mhandle, "proportion", std::vector<double>(num_cells));

    auto thandle = rhandle.createGroup("thresholds");
    quick_write_dataset(thandle, "sums", std::vector<double>(num_batches));
    quick_write_dataset(thandle, "detected", std::vector<double>(num_batches));
    quick_write_dataset(thandle, "proportion", std::vector<double>(num_batches));

    quick_write_dataset(rhandle, "discards", std::vector<int>(num_cells));
    return;
}

TEST(QualityControl, AllOK) {
    const std::string path = "TEST_quality_control.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::quality_control::validate(handle, 100, 1));
    }

    // Works with more batches.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 200, 2);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::quality_control::validate(handle, 200, 2));
    }
}

void quick_qc_throw(const std::string& path, int num_cells, int num_batches, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::quality_control::validate(handle, num_cells, num_batches);
    }, msg);
}

TEST(QualityControl, ParametersFailed) {
    const std::string path = "TEST_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters");
    }
    quick_qc_throw(path, 100, 1, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/mito_prefix");
    }
    quick_qc_throw(path, 100, 1, "'mito_prefix'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/use_mito_default");
    }
    quick_qc_throw(path, 100, 1, "'use_mito_default'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/nmads");
    }
    quick_qc_throw(path, 100, 1, "'nmads'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/parameters/nmads");
        quick_write_dataset(handle, "quality_control/parameters/nmads", -1.0);
    }
    quick_qc_throw(path, 100, 1, "non-negative");
}

TEST(QualityControl, ResultsFailed) {
    const std::string path = "TEST_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 1);
        handle.unlink("quality_control/results");
    }
    quick_qc_throw(path, 100, 1, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 200, 1);
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve metrics");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_quality_control(handle, 100, 2);
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve thresholds");
}
