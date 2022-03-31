#include <gtest/gtest.h>
#include "kanaval/snn_graph_cluster.hpp"
#include "utils.h"
#include <iostream>

void add_snn_graph_cluster(H5::H5File& handle, int num_cells, int num_clusters = 10) {
    auto qhandle = handle.createGroup("snn_graph_cluster");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "k", 10);
    quick_write_dataset(phandle, "scheme", "rank");
    quick_write_dataset(phandle, "resolution", 0.5);

    auto rhandle = qhandle.createGroup("results");
    std::vector<int> clusters(num_cells);
    for (int i = 0; i < num_cells; ++i) {
        clusters[i] = i % num_clusters;
    }
    quick_write_dataset(rhandle, "clusters", clusters);

    return;
}

TEST(SnnGraphCluster, AllOK) {
    const std::string path = "TEST_snn_graph_cluster.h5";

    // Basic case.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::snn_graph_cluster::validate(handle, 1000));
    }

    // Allowed to be missing.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::snn_graph_cluster::validate(handle, 1000, false));
    }
}

void quick_snn_graph_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::snn_graph_cluster::validate(handle, num_cells);
    }, msg);
}

TEST(SnnGraphCluster, ParametersFailed) {
    const std::string path = "TEST_snn_graph_cluster.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters");
    }
    quick_snn_graph_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/k");
    }
    quick_snn_graph_throw(path, 1000, "'k' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/k");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "k", 0);
    }
    quick_snn_graph_throw(path, 1000, "positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/scheme");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "scheme", "foo");
    }
    quick_snn_graph_throw(path, 1000, "scheme");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/scheme");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "scheme", "foo");
    }
    quick_snn_graph_throw(path, 1000, "scheme");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/resolution");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "resolution", -1);
    }
    quick_snn_graph_throw(path, 1000, "resolution");
}

TEST(SnnGraphCluster, ResultsFailed) {
    const std::string path = "TEST_snn_graph_cluster.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results");
    }
    quick_snn_graph_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
    }
    quick_snn_graph_throw(path, 1000, "'clusters' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
    }
    quick_snn_graph_throw(path, 500, "expected dimensions");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
        auto rhandle = handle.openGroup("snn_graph_cluster/results");
        quick_write_dataset(rhandle, "clusters", std::vector<int>(1000, -1));
    }
    quick_snn_graph_throw(path, 1000, "non-negative");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
        auto rhandle = handle.openGroup("snn_graph_cluster/results");
        quick_write_dataset(rhandle, "clusters", std::vector<int>(1000, 10));
    }
    quick_snn_graph_throw(path, 1000, "represented at least once");
}
