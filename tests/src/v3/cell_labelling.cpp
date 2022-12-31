#include <gtest/gtest.h>
#include "kanaval/v3/cell_labelling.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_cell_labelling(H5::H5File& handle, int num_clusters = 3) {
    auto qhandle = handle.createGroup("cell_labelling");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "mouse_references", std::vector<std::string>{ "ImmGen", "MouseRNAseq" });
    quick_write_dataset(phandle, "human_references", std::vector<std::string>{ "BlueprintEncode", "DatabaseImmuneCellExpression" });

    auto rhandle = qhandle.createGroup("results");
    auto perhandle = rhandle.createGroup("per_reference");
    std::vector<std::string> dummy(num_clusters, "something");
    quick_write_dataset(perhandle, "ImmGen", dummy);
    quick_write_dataset(perhandle, "MouseRNAseq", dummy);
    quick_write_dataset(perhandle, "DatabaseImmuneCellExpression", dummy);
    quick_write_dataset(perhandle, "BlueprintEncode", dummy);

    std::vector<std::string> refs(num_clusters, "ImmGen");
    quick_write_dataset(rhandle, "integrated", refs);
    return;
}

}

TEST(CellLabellingV3, AllOK) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_cell_labelling(handle, 3, true, latest));
    }

    // Still works if one of the references isn't used; that's okay.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/per_reference/ImmGen");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_cell_labelling(handle, 3, true, latest));
    }
}

static void quick_label_throw(const std::string& path, int num_clusters, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_cell_labelling(handle, num_clusters, true, latest);
    }, msg);
}

TEST(CellLabellingV3, ParametersFailed) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/parameters");
    }
    quick_label_throw(path, 3, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/parameters/mouse_references");

        auto phandle = handle.openGroup("cell_labelling/parameters");
        quick_write_dataset(phandle, "mouse_references", std::vector<std::string>{ "ImmGen", "ImmGen" });
    }
    quick_label_throw(path, 3, "duplicated");
}

TEST(CellLabellingV3, ResultsFailed) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/results");
    }
    quick_label_throw(path, 3, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        auto per_handle = handle.openGroup("cell_labelling/results/per_reference");
        quick_write_dataset(per_handle, "FooBar", std::vector<std::string>{ "something", "something", "something" });
    }
    quick_label_throw(path, 3, "FooBar");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
    }
    quick_label_throw(path, 4, "expected dimensions");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/integrated");
    }
    quick_label_throw(path, 3, "integrated");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/integrated");
        auto rhandle = handle.openGroup("cell_labelling/results");
        quick_write_dataset(rhandle, "integrated", std::vector<std::string>{ "something", "something" });
    }
    quick_label_throw(path, 3, "number of clusters");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/integrated");
        auto rhandle = handle.openGroup("cell_labelling/results");
        quick_write_dataset(rhandle, "integrated", std::vector<std::string>{ "something", "something", "something" });
    }
    quick_label_throw(path, 3, "not listed");
}

TEST(CellLabellingV3, NotInUse) {
    const std::string path = "TEST_cell_labelling.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_cell_labelling(handle);
        handle.unlink("cell_labelling/results/per_reference");
    }
    quick_label_throw(path, 3, "'per_reference' group");

    // Okay if it's missing.
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_cell_labelling(handle, 3, false, latest));
    }
}
