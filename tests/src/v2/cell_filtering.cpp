#include <gtest/gtest.h>
#include "kanaval/v2/cell_filtering.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_cell_filtering(H5::H5File& handle, int num_cells, int lost = 10) {
    auto qhandle = handle.createGroup("cell_filtering");
    qhandle.createGroup("parameters");

    auto rhandle = qhandle.createGroup("results");
    std::vector<int> discard(num_cells);
    std::fill(discard.begin(), discard.begin() + lost, 1);
    quick_write_dataset(rhandle, "discards", discard);
}

}

TEST(CellFilteringV2, AllOK) {
    const std::string path = "TEST_cell_filtering.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_cell_filtering(handle, 100);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v2::validate_cell_filtering(handle, 100, 2, latest), 90);
    }

    // Works if there's only one modality and no discards.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results/discards");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v2::validate_cell_filtering(handle, 100, 1, latest), -1);
    }
}

static void quick_filter_throw(const std::string& path, int num_cells, int num_modalities, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_cell_filtering(handle, num_cells, num_modalities, latest);
    }, msg);
}

TEST(CellFilteringV2, ParametersFailed) {
    const std::string path = "TEST_cell_filtering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/parameters");
    }
    quick_filter_throw(path, 100, 2, "'parameters' group");
}

TEST(CellFilteringV2, ResultsFailed) {
    const std::string path = "TEST_cell_filtering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results");
    }
    quick_filter_throw(path, 100, 2, "'results' group");

    // Throws if multiple modalities are requested but 'discards' is absent.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results/discards");
    }
    quick_filter_throw(path, 100, 2, "'discards' dataset");

    // Wrong number of discards.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_cell_filtering(handle, 100);
    }
    quick_filter_throw(path, 90, 2, "dimensions");
}
