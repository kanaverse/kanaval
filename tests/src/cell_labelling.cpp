#include <gtest/gtest.h>
#include "kanaval/cell_labelling.hpp"
#include "utils.h"
#include <iostream>

void add_cell_labelling(H5::H5File& handle) {
    auto qhandle = handle.createGroup("cell_labelling");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "mouse_references", std::vector<std::string>{ "ImmGen", "MouseRNAseq" });
    quick_write_dataset(phandle, "human_references", std::vector<std::string>{ "BlueprintEncode", "DatabaseImmuneCellExpression" });

    auto rhandle = qhandle.createGroup("results");
    auto perhandle = rhandle.createGroup("per_reference");
    quick_write_dataset(perhandle, "ImmGen", std::vector<std::string>{ "something", "something", "something" });
    quick_write_dataset(perhandle, "MouseRNAseq", std::vector<std::string>{ "something", "something", "something" });
    quick_write_dataset(perhandle, "DatabaseImmuneCellExpression", std::vector<std::string>{ "something", "something", "something" });
    quick_write_dataset(perhandle, "BlueprintEncode", std::vector<std::string>{ "something", "something", "something" });

    quick_write_dataset(rhandle, "integrated", std::vector<std::string>{ "ImmGen", "ImmGen", "ImmGen" });
    return;
}

TEST(CellLabelling, AllOK) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::cell_labelling::validate(handle, 3));
    }
}

void quick_label_throw(const std::string& path, int num_clusters, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::cell_labelling::validate(handle, num_clusters);
    }, msg);
}

TEST(CellLabelling, ParametersFailed) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/parameters");
    }
    quick_label_throw(path, 3, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/parameters/mouse_references");

        auto phandle = handle.openGroup("cell_labelling/parameters");
        quick_write_dataset(phandle, "mouse_references", std::vector<std::string>{ "ImmGen", "ImmGen" });
    }
    quick_label_throw(path, 3, "duplicated");
}

TEST(CellLabelling, ResultsFailed) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/results");
    }
    quick_label_throw(path, 3, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/per_reference/ImmGen");
    }
    quick_label_throw(path, 3, "ImmGen");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        auto per_handle = handle.openGroup("cell_labelling/results/per_reference");
        quick_write_dataset(per_handle, "FooBar", std::vector<std::string>{ "something", "something", "something" });
    }
    quick_label_throw(path, 3, "FooBar");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
    }
    quick_label_throw(path, 4, "expected dimensions");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/integrated");
    }
    quick_label_throw(path, 3, "integrated");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/integrated");
        auto rhandle = handle.openGroup("cell_labelling/results");
        quick_write_dataset(rhandle, "integrated", std::vector<std::string>{ "something", "something" });
    }
    quick_label_throw(path, 3, "number of clusters");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/integrated");
        auto rhandle = handle.openGroup("cell_labelling/results");
        quick_write_dataset(rhandle, "integrated", std::vector<std::string>{ "something", "something", "something" });
    }
    quick_label_throw(path, 3, "not listed");
}
