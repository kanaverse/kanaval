#include <gtest/gtest.h>
#include "kanaval/v3/rna_quality_control.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_rna_quality_control(H5::H5File& handle, int num_cells, int num_blocks, int lost = 10) {
    auto qhandle = handle.createGroup("rna_quality_control");

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
    quick_write_dataset(thandle, "sums", std::vector<double>(num_blocks));
    quick_write_dataset(thandle, "detected", std::vector<double>(num_blocks));
    quick_write_dataset(thandle, "proportion", std::vector<double>(num_blocks));

    std::vector<int> discard(num_cells);
    std::fill(discard.begin(), discard.begin() + lost, 1);
    quick_write_dataset(rhandle, "discards", discard);
    return;
}

}

TEST(RnaQualityControlV3, AllOK) {
    const std::string path = "TEST_rna_quality_control.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v3::validate_rna_quality_control(handle, 100, 1, true, latest);
        EXPECT_EQ(out, 90);
    }

    // Works with more batches.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 200, 2);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v3::validate_rna_quality_control(handle, 200, 2, true, latest);
        EXPECT_EQ(out, 190);
    }
}

static void quick_qc_throw(const std::string& path, int num_cells, int num_blocks, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_rna_quality_control(handle, num_cells, num_blocks, true, latest);
    }, msg);
}

TEST(RnaQualityControlV3, ParametersFailed) {
    const std::string path = "TEST_rna_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        handle.unlink("rna_quality_control/parameters");
    }
    quick_qc_throw(path, 100, 1, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        handle.unlink("rna_quality_control/parameters/mito_prefix");
    }
    quick_qc_throw(path, 100, 1, "'mito_prefix'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        handle.unlink("rna_quality_control/parameters/use_mito_default");
    }
    quick_qc_throw(path, 100, 1, "'use_mito_default'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        handle.unlink("rna_quality_control/parameters/nmads");
    }
    quick_qc_throw(path, 100, 1, "'nmads'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        handle.unlink("rna_quality_control/parameters/nmads");
        quick_write_dataset(handle, "rna_quality_control/parameters/nmads", -1.0);
    }
    quick_qc_throw(path, 100, 1, "non-negative");
}

TEST(RnaQualityControlV3, ResultsFailed) {
    const std::string path = "TEST_rna_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        handle.unlink("rna_quality_control/results");
    }
    quick_qc_throw(path, 100, 1, "'results' group");

    // Mismatch in number of cells.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 200, 1);
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve metrics");

    // Mismatch in number of blocks.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 2);
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve thresholds");
}

TEST(RnaQualityControlV3, NotInUse) {
    const std::string path = "TEST_rna_quality_control.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_rna_quality_control(handle, 100, 1);
        auto rhandle = handle.openGroup("rna_quality_control/results");
        rhandle.unlink("thresholds");
        rhandle.unlink("discards");
        rhandle.unlink("metrics");
    }
    quick_qc_throw(path, 100, 1, "failed to retrieve metrics");

    // Tolerates absence of all results.
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto out = kanaval::v3::validate_rna_quality_control(handle, 100, 1, false, latest);
        EXPECT_EQ(out, -1);
    }
}
