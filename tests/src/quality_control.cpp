#include <gtest/gtest.h>
#include "kanaval/quality_control.hpp"
#include "utils.h"

void add_quality_control(H5::H5File handle, int num_cells, int num_batches) {
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

TEST(QualityControl, ParametersOk) {
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



