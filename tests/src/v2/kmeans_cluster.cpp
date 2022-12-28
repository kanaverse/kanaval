#include <gtest/gtest.h>
#include "kanaval/v2/kmeans_cluster.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_kmeans_cluster(H5::H5File& handle, int num_cells, int num_clusters = 10) {
    auto qhandle = handle.createGroup("kmeans_cluster");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "k", 10);

    auto rhandle = qhandle.createGroup("results");
    std::vector<int> clusters(num_cells);
    for (int i = 0; i < num_cells; ++i) {
        clusters[i] = i % num_clusters;
    }
    quick_write_dataset(rhandle, "clusters", clusters);

    return;
}

}

TEST(KmeansClusterV2, AllOK) {
    const std::string path = "TEST_kmeans_cluster.h5";

    // Basic case.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_kmeans_cluster(handle, 1000));
    }

    // Allowed to be missing.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/results/clusters");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_kmeans_cluster(handle, 1000, false));
    }
}

static void quick_kmeans_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_kmeans_cluster(handle, num_cells);
    }, msg);
}

TEST(KmeansClusterV2, ParametersFailed) {
    const std::string path = "TEST_kmeans_cluster.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/parameters");
    }
    quick_kmeans_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/parameters/k");
    }
    quick_kmeans_throw(path, 1000, "'k' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/parameters/k");
        auto phandle = handle.openGroup("kmeans_cluster/parameters");
        quick_write_dataset(phandle, "k", 0);
    }
    quick_kmeans_throw(path, 1000, "positive");
}

TEST(KmeansClusterV2, ResultsFailed) {
    const std::string path = "TEST_kmeans_cluster.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/results");
    }
    quick_kmeans_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/results/clusters");
    }
    quick_kmeans_throw(path, 1000, "'clusters' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
    }
    quick_kmeans_throw(path, 500, "expected dimensions");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/results/clusters");
        auto rhandle = handle.openGroup("kmeans_cluster/results");
        quick_write_dataset(rhandle, "clusters", std::vector<int>(1000, 10));
    }
    quick_kmeans_throw(path, 1000, "out of range");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_kmeans_cluster(handle, 1000);
        handle.unlink("kmeans_cluster/results/clusters");
        auto rhandle = handle.openGroup("kmeans_cluster/results");
        quick_write_dataset(rhandle, "clusters", std::vector<int>(1000, 1));
    }
    quick_kmeans_throw(path, 1000, "represented at least once");
}
