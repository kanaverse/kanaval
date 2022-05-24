#include <gtest/gtest.h>
#include "kanaval/batch_correction.hpp"
#include "utils.h"
#include <iostream>

void add_batch_correction(H5::H5File& handle, int num_cells, int num_pcs) {
    auto qhandle = handle.createGroup("batch_correction");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "approximate", 1);
    quick_write_dataset(phandle, "num_neighbors", 20);
    quick_write_dataset(phandle, "method", "mnn");

    auto rhandle = qhandle.createGroup("results");
    H5::DataSpace space;
    std::vector<hsize_t> dims(2);
    dims[0] = num_cells;
    dims[1] = num_pcs;
    space.setExtentSimple(2, dims.data());
    rhandle.createDataSet("corrected", H5::PredType::NATIVE_DOUBLE, space);
}

TEST(BatchCorrection, AllOK) {
    const std::string path = "TEST_batch_correction.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::batch_correction::validate(handle, 10, 100, 2, latest));
    }

    // Works if there's only one sample and no PCs.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/results/corrected");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::batch_correction::validate(handle, 10, 100, 1, latest));
    }

    // Works with multiple samples and no correction.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/results/corrected");
        handle.unlink("batch_correction/parameters/method");
        quick_write_dataset(handle, "batch_correction/parameters/method", "none");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::batch_correction::validate(handle, 10, 100, 2, latest));
    }
}

void quick_correct_throw(const std::string& path, int num_pcs, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::batch_correction::validate(handle, num_pcs, num_cells, 2, latest);
    }, msg);
}

TEST(BatchCorrection, ParametersFailed) {
    const std::string path = "TEST_batch_correction.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/parameters");
    }
    quick_correct_throw(path, 10, 100, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/parameters/approximate");
    }
    quick_correct_throw(path, 10, 100, "approximate");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/parameters/method");
    }
    quick_correct_throw(path, 10, 100, "method");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/parameters/num_neighbors");
    }
    quick_correct_throw(path, 10, 100, "num_neighbors");
}

TEST(BatchCorrection, ResultsFailed) {
    const std::string path = "TEST_batch_correction.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/results");
    }
    quick_correct_throw(path, 10, 100, "'results' group");

    // Throws if multiple samples are requested but 'corrected' is absent.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/results/corrected");
    }
    quick_correct_throw(path, 10, 100, "'corrected'");

    // Wrong dimensions.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_batch_correction(handle, 100, 10);
        handle.unlink("batch_correction/results/corrected");
    }
    quick_correct_throw(path, 5, 100, "'corrected'");
}
