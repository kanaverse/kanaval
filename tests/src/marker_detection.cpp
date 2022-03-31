#include <gtest/gtest.h>
#include "kanaval/marker_detection.hpp"
#include "utils.h"
#include <iostream>

void add_marker_detection(H5::H5File& handle, int ngenes, int nclusters) {
    auto qhandle = handle.createGroup("marker_detection");
    qhandle.createGroup("parameters");

    auto rhandle = qhandle.createGroup("results");
    auto chandle = rhandle.createGroup("clusters");

    for (int i = 0; i < nclusters; ++i) {
        auto ihandle = chandle.createGroup(std::to_string(i));
        quick_write_dataset(ihandle, "means", std::vector<double>(ngenes));
        quick_write_dataset(ihandle, "detected", std::vector<double>(ngenes));

        for (const auto& e : kanaval::markers::effects) {
            auto ehandle = ihandle.createGroup(e);
            quick_write_dataset(ehandle, "mean", std::vector<double>(ngenes));
            quick_write_dataset(ehandle, "min", std::vector<double>(ngenes));
            quick_write_dataset(ehandle, "min_rank", std::vector<double>(ngenes));
        }
    }

    return;
}

TEST(MarkerDetection, AllOK) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 5);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::marker_detection::validate(handle, 100, 5));
    }
}

void quick_marker_throw(const std::string& path, int nclusters, int ngenes, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::marker_detection::validate(handle, nclusters, ngenes);
    }, msg);
}

TEST(MarkerDetection, ParametersFailed) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/parameters");
    }
    quick_marker_throw(path, 100, 7, "'parameters' group");
}

TEST(MarkerDetection, ResultsFailed) {
    const std::string path = "TEST_marker_detection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results");
    }
    quick_marker_throw(path, 100, 7, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/clusters/2");
    }
    quick_marker_throw(path, 100, 7, "number of clusters");

    // Same number of clusters as expected, but now the numbers aren't not consecutive.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/clusters/2");
    }
    quick_marker_throw(path, 100, 6, "cluster 2");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/clusters/2/detected");
    }
    quick_marker_throw(path, 100, 7, "detected");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/clusters/2/lfc");
    }
    quick_marker_throw(path, 100, 7, "summary statistic");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_marker_detection(handle, 100, 7);
        handle.unlink("marker_detection/results/clusters/2/auc/min");
    }
    quick_marker_throw(path, 100, 7, "auc");
}
