#include <gtest/gtest.h>
#include "kanaval/v3/marker_detection.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_marker_detection(H5::H5File& handle, const std::unordered_map<std::string, int>& modalities, int nclusters, bool has_auc = true, double lfc_threshold = 0) {
    auto qhandle = handle.createGroup("marker_detection");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "compute_auc", static_cast<int>(has_auc));
    quick_write_dataset(phandle, "lfc_threshold", lfc_threshold);

    auto rhandle = qhandle.createGroup("results");
    auto chandle = rhandle.createGroup("per_cluster");
    for (const auto& p : modalities) {
        auto mohandle = chandle.createGroup(p.first);
        auto ngenes = p.second;

        for (int i = 0; i < nclusters; ++i) {
            auto istr = std::to_string(i);
            auto xhandle = mohandle.createGroup(istr);
            quick_write_dataset(xhandle, "means", std::vector<double>(ngenes));
            quick_write_dataset(xhandle, "detected", std::vector<double>(ngenes));

            for (const auto& e : kanaval::v3::markers::effects) {
                if (!has_auc && e == "auc") {
                    continue;
                }

                auto ehandle = xhandle.createGroup(e);
                quick_write_dataset(ehandle, "mean", std::vector<double>(ngenes));
                quick_write_dataset(ehandle, "min", std::vector<double>(ngenes));
                quick_write_dataset(ehandle, "min_rank", std::vector<double>(ngenes));
            }
        }
    }
}

void add_marker_detection(H5::H5File& handle, int ngenes, int nclusters, bool has_auc = true, double lfc_threshold = 0) {
    add_marker_detection(handle, { { "RNA", ngenes } }, nclusters, has_auc, lfc_threshold);
}

}

TEST(MarkerDetectionV3, AllOK) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 5);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_marker_detection(handle, 5, { { "RNA", 100 } }, latest));
    }
}

TEST(MarkerDetectionV3, MultiModalOK) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, { { "RNA", 100 }, { "ADT", 5 } }, 5);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_marker_detection(handle, 5, { { "RNA", 100 }, { "ADT", 5 } }, latest));
    }

    // Fails if multiple modalities were expected but not found.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 5);
    }
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_marker_detection(handle, 5, { { "RNA", 100 }, { "ADT", 5 } }, latest);
    }, "ADT");
}

static void quick_marker_throw(const std::string& path, int ngenes, int nclusters, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_marker_detection(handle, nclusters, { { "RNA", ngenes } }, latest);
    }, msg);
}

TEST(MarkerDetectionV3, ParametersFailed) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/parameters/compute_auc");
    }
    quick_marker_throw(path, 100, 7, "'compute_auc'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7, false, -1);
    }
    quick_marker_throw(path, 100, 7, "'lfc_threshold' must be non-negative");
}

TEST(MarkerDetectionV3, ResultsFailed) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results");
    }
    quick_marker_throw(path, 100, 7, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA");
    }
    quick_marker_throw(path, 100, 7, "RNA");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2");
    }
    quick_marker_throw(path, 100, 7, "number of clusters");

    // Same number of clusters as expected, but now the numbers aren't not consecutive.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2");
    }
    quick_marker_throw(path, 100, 6, "cluster 2");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2/detected");
    }
    quick_marker_throw(path, 100, 7, "detected");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2/lfc");
    }
    quick_marker_throw(path, 100, 7, "summary statistic");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2/auc/min");
    }
    quick_marker_throw(path, 100, 7, "auc");
}

TEST(MarkerDetectionV3, NoAucOK) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_marker_detection(handle, 100, 5, false);
    }

    // Missing an AUC is okay...
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_marker_detection(handle, 5, { { "RNA", 100 } }, latest));
    }

    // Unless we force it to use AUCs.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto phandle = handle.openGroup("marker_detection/parameters");
        phandle.unlink("compute_auc");
        quick_write_dataset(phandle, "compute_auc", 1);
    }
    quick_marker_throw(path, 100, 5, "auc");
}


