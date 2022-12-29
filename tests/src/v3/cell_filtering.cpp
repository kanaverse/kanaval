#include <gtest/gtest.h>
#include "kanaval/v3/cell_filtering.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_cell_filtering(H5::H5File& handle, int num_cells, int lost = 10) {
    auto qhandle = handle.createGroup("cell_filtering");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "use_rna", 1);
    quick_write_dataset(phandle, "use_adt", 1);
    quick_write_dataset(phandle, "use_crispr", 0);

    auto rhandle = qhandle.createGroup("results");
    std::vector<int> discard(num_cells);
    std::fill(discard.begin(), discard.begin() + lost, 1);
    quick_write_dataset(rhandle, "discards", discard);
}

}

TEST(CellFilteringV3, AllOK) {
    const std::string path = "TEST_cell_filtering.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100, /* lost */ 10);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v3::validate_cell_filtering(handle, 100, { { "RNA", 77 }, { "ADT", 61 } }, latest), 90);
    }

    // Works if there's only one modality and no discards,
    // in which case that modality is used directly.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results/discards");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v3::validate_cell_filtering(handle, 100, { { "RNA", 77 } }, latest), 77);
    }

    // Works if there's no modalities that are also being used;
    // in which case, all cells are to be used.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results/discards");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_EQ(kanaval::v3::validate_cell_filtering(handle, 100, { { "CRISPR", 99 } }, latest), 100);
    }
}

static void quick_filter_throw(const std::string& path, int num_cells, const std::unordered_map<std::string, int>& modalities, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_cell_filtering(handle, num_cells, modalities, latest);
    }, msg);
}

TEST(CellFilteringV3, ParametersFailed) {
    const std::string path = "TEST_cell_filtering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/parameters/use_rna");
    }
    quick_filter_throw(path, 100, {}, "'use_rna'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100);
        auto phandle = handle.openGroup("cell_filtering/parameters");
        phandle.unlink("use_adt");
        quick_write_dataset(phandle, "use_adt", 1234.0);
    }
    quick_filter_throw(path, 100, {}, "'use_adt'");
}

TEST(CellFilteringV3, ResultsFailed) {
    const std::string path = "TEST_cell_filtering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results");
    }
    quick_filter_throw(path, 100, { { "RNA", 100 } }, "'results' group");

    // Throws if multiple modalities are requested but 'discards' is absent.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 100);
        handle.unlink("cell_filtering/results/discards");
    }
    quick_filter_throw(path, 100, { { "RNA", 81 }, { "ADT", 77 } }, "'discards' dataset");

    // Wrong number of discards.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_filtering(handle, 99);
    }
    quick_filter_throw(path, 100, { { "RNA", 81 }, { "ADT", 77 } }, "'discards' dataset");
}
