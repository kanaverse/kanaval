#include <gtest/gtest.h>
#include "kanaval/v3/custom_selections.hpp"
#include "../utils.h"
#include <numeric>

namespace v3 {

void add_custom_selections(H5::H5File& handle, const std::unordered_map<std::string, int>& modalities, int ncells = 10, bool has_auc = true, double lfc_threshold = 0) {
    auto qhandle = handle.createGroup("custom_selections");

    std::vector<std::string> available { "foo", "bar", "whee" };
    {
        auto phandle = qhandle.createGroup("parameters");
        quick_write_dataset(phandle, "compute_auc", static_cast<int>(has_auc));
        quick_write_dataset(phandle, "lfc_threshold", lfc_threshold);

        auto shandle = phandle.createGroup("selections");
        for (const auto& a : available) {
            std::vector<int> chosen(ncells);
            std::iota(chosen.begin(), chosen.end(), 0);
            quick_write_dataset(shandle, a, chosen); 
        }
    }

    auto rhandle = qhandle.createGroup("results");
    auto pshandle = rhandle.createGroup("per_selection");
    for (const auto& s : available) {
        auto shandle = pshandle.createGroup(s);

        for (const auto& mod : modalities) {
            auto mhandle = shandle.createGroup(mod.first);
            quick_write_dataset(mhandle, "means", std::vector<double>(mod.second));
            quick_write_dataset(mhandle, "detected", std::vector<double>(mod.second));

            for (const auto& e : kanaval::v3::markers::effects) {
                if (!has_auc && e == "auc") {
                    continue;
                }
                quick_write_dataset(mhandle, e, std::vector<double>(mod.second));
            }
        }
    }
}

void add_custom_selections(H5::H5File& handle, int ngenes, int ncells = 10, bool has_auc = true) {
    add_custom_selections(handle, { { "RNA", ngenes } }, ncells, has_auc);
}

}

TEST(CustomSelectionsV3, AllOK) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_custom_selections(handle, 10, { { "RNA", 1000 } }, latest));
    }
}

TEST(CustomSelectionsV3, MultiModalOK) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, { { "RNA", 1000 }, { "ADT", 5 } });
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_custom_selections(handle, 10, { { "RNA", 1000 }, { "ADT", 5 } }, latest));
    }

    // Fails if multiple modalities were expected but not found.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 100, 5);
    }
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_custom_selections(handle, 10, { { "RNA", 100 }, { "ADT", 5 } }, latest);
    }, "ADT");
}

static void quick_custom_throw(const std::string& path, int ngenes, int ncells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_custom_selections(handle, ncells, { { "RNA", ngenes } }, latest);
    }, msg);
}

TEST(CustomSelectionsV3, ParametersFailed) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/parameters");
    }
    quick_custom_throw(path, 500, 10, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/parameters/compute_auc");
    }
    quick_custom_throw(path, 500, 10, "'compute_auc'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/parameters/lfc_threshold");
        auto phandle = handle.openGroup("custom_selections/parameters");
        quick_write_dataset(phandle, "lfc_threshold", -1.0);
    }
    quick_custom_throw(path, 500, 10, "'lfc_threshold'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        auto shandle = handle.openGroup("custom_selections/parameters/selections");
        quick_write_dataset(shandle, "stuff", std::vector<int>(10, 10000)); 
    }
    quick_custom_throw(path, 500, 10, "out of range");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        auto shandle = handle.openGroup("custom_selections/parameters/selections");
        quick_write_dataset(shandle, "stuff", std::vector<int>(5, 5)); 
    }
    quick_custom_throw(path, 500, 10, "sorted and unique");
}

TEST(CustomSelectionsV3, ResultsFailed) {
    const std::string path = "TEST_custom_selections.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results");
    }
    quick_custom_throw(path, 500, 10, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo");
    }
    quick_custom_throw(path, 500, 10, "number of groups");

    // Same number of groups but different names.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo");
        handle.createGroup("custom_selections/results/per_selection/yay");
    }
    quick_custom_throw(path, 500, 10, "'foo' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo/RNA/detected");
    }
    quick_custom_throw(path, 500, 10, "detected");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 500);
        handle.unlink("custom_selections/results/per_selection/foo/RNA/lfc");
    }
    quick_custom_throw(path, 500, 10, "lfc");
}

TEST(CustomSelectionsV3, NoAucOK) {
    const std::string path = "TEST_custom_selection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_custom_selections(handle, 100, 5, false);
    }

    // Missing an AUC is okay...
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_custom_selections(handle, 20, { { "RNA", 100 } }, latest));
    }

    // Unless we force it to use AUCs.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto phandle = handle.openGroup("custom_selections/parameters");
        phandle.unlink("compute_auc");
        quick_write_dataset(phandle, "compute_auc", 1);
    }
    quick_custom_throw(path, 100, 20, "auc");
}
