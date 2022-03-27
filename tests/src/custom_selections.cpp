#include <gtest/gtest.h>
#include "kanaval/custom_selections.hpp"
#include "utils.h"
#include <numeric>

void add_custom_selections(H5::H5File& handle, int ngenes) {
    auto qhandle = handle.createGroup("custom_selections");
    auto phandle = qhandle.createGroup("parameters");

    auto shandle = phandle.createGroup("selections");
    std::vector<std::string> available { "foo", "bar", "whee" };
    for (const auto& a : available) {
        std::vector<int> chosen(10);
        std::iota(chosen.begin(), chosen.end(), 0);
        quick_write_dataset(shandle, a, chosen); 
    }

    auto rhandle = qhandle.createGroup("results");
    auto mhandle = rhandle.createGroup("markers");
    for (const auto& a :available) {
        auto ihandle = mhandle.createGroup(a);
        quick_write_dataset(ihandle, "means", std::vector<double>(ngenes));
        quick_write_dataset(ihandle, "detected", std::vector<double>(ngenes));
        for (const auto& e : kanaval::markers::effects) {
            quick_write_dataset(ihandle, e, std::vector<double>(ngenes));
        }
    }

    return;
}

TEST(CustomSelections, AllOK) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::custom_selections::validate(handle, 1000, 10));
    }
}

void quick_custom_throw(const std::string& path, int ngenes, int ncells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::custom_selections::validate(handle, ngenes, ncells);
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
        handle.unlink("custom_selections/results/markers/foo");
    }
    quick_custom_throw(path, 500, 10, "number of groups");

    // Same number of groups but different names.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/markers/foo");
        handle.createGroup("custom_selections/results/markers/yay");
    }
    quick_custom_throw(path, 500, 10, "'foo' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/markers/foo/detected");
    }
    quick_custom_throw(path, 500, 10, "detected");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/markers/foo/lfc");
    }
    quick_custom_throw(path, 500, 10, "lfc");
}
