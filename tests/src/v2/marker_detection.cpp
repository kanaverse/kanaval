#include <gtest/gtest.h>
#include "kanaval/v2/marker_detection.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_marker_detection_base(H5::Group& handle, int ngenes, int nclusters) {
    for (int i = 0; i < nclusters; ++i) {
        auto istr = std::to_string(i);
        auto xhandle = handle.createGroup(istr);
        quick_write_dataset(xhandle, "means", std::vector<double>(ngenes));
        quick_write_dataset(xhandle, "detected", std::vector<double>(ngenes));
        for (const auto& e : kanaval::v2::markers::effects) {
            auto ehandle = xhandle.createGroup(e);
            quick_write_dataset(ehandle, "mean", std::vector<double>(ngenes));
            quick_write_dataset(ehandle, "min", std::vector<double>(ngenes));
            quick_write_dataset(ehandle, "min_rank", std::vector<double>(ngenes));
        }
    }
}

void add_marker_detection(H5::H5File& handle, const std::vector<int>& ngenes, int nclusters, const std::vector<std::string>& modality) {
    auto qhandle = handle.createGroup("marker_detection");
    qhandle.createGroup("parameters");
    auto rhandle = qhandle.createGroup("results");
    auto chandle = rhandle.createGroup("per_cluster");
    for (size_t m = 0; m < modality.size(); ++m) {
        auto mohandle = chandle.createGroup(modality[m]);
        v2::add_marker_detection_base(mohandle, ngenes[m], nclusters);
    }
}

void add_marker_detection(H5::H5File& handle, int ngenes, int nclusters) {
    v2::add_marker_detection(handle, { ngenes }, nclusters, { "RNA" });
}

void add_marker_detection_legacy(H5::H5File& handle, int ngenes, int nclusters) {
    auto qhandle = handle.createGroup("marker_detection");
    qhandle.createGroup("parameters");
    auto rhandle = qhandle.createGroup("results");
    auto chandle = rhandle.createGroup("clusters");
    v2::add_marker_detection_base(chandle, ngenes, nclusters);
}

}

TEST(MarkerDetectionV2, AllOK) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 5);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_marker_detection(handle, 5, { "RNA" }, { 100 }, latest));
    }
}

TEST(MarkerDetectionV2, MultiModalOK) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, { 100, 5 }, 5, { "RNA", "ADT" });
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_marker_detection(handle, 5, { "RNA", "ADT" }, { 100, 5 }, latest));
    }

    // Fails if multiple modalities were expected but not found.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 5);
    }
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_marker_detection(handle, 5, { "RNA", "ADT" }, { 100, 5 }, latest);
    }, "ADT");
}

void quick_marker_throw(const std::string& path, int ngenes, int nclusters, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_marker_detection(handle, nclusters, { "RNA" }, { ngenes }, latest);
    }, msg);
}

TEST(MarkerDetectionV2, ParametersFailed) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/parameters");
    }
    quick_marker_throw(path, 100, 7, "'parameters' group");
}

TEST(MarkerDetectionV2, ResultsFailed) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results");
    }
    quick_marker_throw(path, 100, 7, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA");
    }
    quick_marker_throw(path, 100, 7, "RNA");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2");
    }
    quick_marker_throw(path, 100, 7, "number of clusters");

    // Same number of clusters as expected, but now the numbers aren't not consecutive.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2");
    }
    quick_marker_throw(path, 100, 6, "cluster 2");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2/detected");
    }
    quick_marker_throw(path, 100, 7, "detected");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2/lfc");
    }
    quick_marker_throw(path, 100, 7, "summary statistic");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/per_cluster/RNA/2/auc/min");
    }
    quick_marker_throw(path, 100, 7, "auc");
}

TEST(MarkerDetectionV2, LegacyOk) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_marker_detection_legacy(handle, 100, 5);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_marker_detection(handle, 5, {}, {100}, 1001000));
    }
}

