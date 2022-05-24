#include <gtest/gtest.h>
#include "kanaval/custom_selections.hpp"
#include "utils.h"
#include <numeric>

std::vector<std::string> add_custom_selection_parameters(H5::Group& handle, int ncells = 10) {
    auto phandle = handle.createGroup("parameters");
    auto shandle = phandle.createGroup("selections");

    std::vector<std::string> available { "foo", "bar", "whee" };
    for (const auto& a : available) {
        std::vector<int> chosen(ncells);
        std::iota(chosen.begin(), chosen.end(), 0);
        quick_write_dataset(shandle, a, chosen); 
    }

    return available;
}

void add_custom_selection_result_base(H5::Group& handle, int ngenes) {
    quick_write_dataset(handle, "means", std::vector<double>(ngenes));
    quick_write_dataset(handle, "detected", std::vector<double>(ngenes));
    for (const auto& e : kanaval::markers::effects) {
        quick_write_dataset(handle, e, std::vector<double>(ngenes));
    }
}

void add_custom_selections(H5::H5File& handle, const std::vector<std::string>& modalities, const std::vector<int>& ngenes, int ncells = 10) {
    auto qhandle = handle.createGroup("custom_selections");
    auto selections = add_custom_selection_parameters(qhandle, ncells);

    auto rhandle = qhandle.createGroup("results");
    auto pshandle = rhandle.createGroup("per_selection");
    for (const auto& s : selections) {
        auto shandle = pshandle.createGroup(s);
        for (size_t m = 0; m < modalities.size(); ++m) {
            auto mhandle = shandle.createGroup(modalities[m]);
            add_custom_selection_result_base(mhandle, ngenes[m]);
        }
    }
}

void add_custom_selections(H5::H5File& handle, int ngenes, int ncells = 10) {
    add_custom_selections(handle, { "RNA" }, { ngenes }, ncells);
}

void add_custom_selections_legacy(H5::H5File& handle, int ngenes, int ncells = 10) {
    auto qhandle = handle.createGroup("custom_selections");
    auto selections = add_custom_selection_parameters(qhandle, ncells);

    auto rhandle = qhandle.createGroup("results");
    auto mhandle = rhandle.createGroup("markers");
    for (const auto& s : selections) {
        auto shandle = mhandle.createGroup(s);
        add_custom_selection_result_base(shandle, ngenes);
    }
}

TEST(CustomSelections, AllOK) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::custom_selections::validate(handle, 10, { "RNA" }, { 1000 }, latest));
    }
}

TEST(CustomSelections, MultiModalOK) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, { "RNA", "ADT" }, { 1000, 5 });
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::custom_selections::validate(handle, 10, { "RNA", "ADT" }, { 1000, 5 }, latest));
    }

    // Fails if multiple modalities were expected but not found.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 100, 5);
    }
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::custom_selections::validate(handle, 10, { "RNA", "ADT" }, { 100, 5 }, latest);
    }, "ADT");
}

void quick_custom_throw(const std::string& path, int ngenes, int ncells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::custom_selections::validate(handle, ncells, { "RNA" }, { ngenes }, latest);
    }, msg);
}

TEST(CustomSelections, ParametersFailed) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/parameters");
    }
    quick_custom_throw(path, 500, 10, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);

        auto shandle = handle.openGroup("custom_selections/parameters/selections");
        quick_write_dataset(shandle, "stuff", std::vector<int>(10, 10000)); 
    }
    quick_custom_throw(path, 500, 10, "out of range");
}

TEST(CustomSelections, ResultsFailed) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results");
    }
    quick_custom_throw(path, 500, 10, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo");
    }
    quick_custom_throw(path, 500, 10, "number of groups");

    // Same number of groups but different names.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo");
        handle.createGroup("custom_selections/results/per_selection/yay");
    }
    quick_custom_throw(path, 500, 10, "'foo' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo/RNA/detected");
    }
    quick_custom_throw(path, 500, 10, "detected");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo/RNA/lfc");
    }
    quick_custom_throw(path, 500, 10, "lfc");
}

TEST(CustomSelections, LegacyOK) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections_legacy(handle, 100);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::custom_selections::validate(handle, 10, {}, {100}, 1001000));
    }
}

