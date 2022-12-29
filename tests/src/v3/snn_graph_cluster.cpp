#include <gtest/gtest.h>
#include "kanaval/v3/snn_graph_cluster.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_snn_graph_cluster(H5::H5File& handle, int num_cells, int num_clusters = 10) {
    auto qhandle = handle.createGroup("snn_graph_cluster");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "k", 10);
    quick_write_dataset(phandle, "scheme", "rank");
    quick_write_dataset(phandle, "algorithm", "multilevel");
    quick_write_dataset(phandle, "multilevel_resolution", 0.5);
    quick_write_dataset(phandle, "leiden_resolution", 1.0);
    quick_write_dataset(phandle, "walktrap_steps", 4);

    auto rhandle = qhandle.createGroup("results");
    std::vector<int> clusters(num_cells);
    for (int i = 0; i < num_cells; ++i) {
        clusters[i] = i % num_clusters;
    }
    quick_write_dataset(rhandle, "clusters", clusters);

    return;
}

}

TEST(SnnGraphClusterV3, AllOK) {
    const std::string path = "TEST_snn_graph_cluster.h5";

    // Basic case.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_snn_graph_cluster(handle, 1000));
    }

    // Allowed to be missing.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_snn_graph_cluster(handle, 1000, false));
    }
}

static void quick_snn_graph_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_snn_graph_cluster(handle, num_cells);
    }, msg);
}

TEST(SnnGraphClusterV3, ParametersFailed) {
    const std::string path = "TEST_snn_graph_cluster.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters");
    }
    quick_snn_graph_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/k");
    }
    quick_snn_graph_throw(path, 1000, "'k' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/k");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "k", 0);
    }
    quick_snn_graph_throw(path, 1000, "positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/scheme");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "scheme", "foo");
    }
    quick_snn_graph_throw(path, 1000, "scheme");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/scheme");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "scheme", "foo");
    }
    quick_snn_graph_throw(path, 1000, "scheme");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/algorithm");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "algorithm", "foo");
    }
    quick_snn_graph_throw(path, 1000, "algorithm");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/multilevel_resolution");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "multilevel_resolution", -1);
    }
    quick_snn_graph_throw(path, 1000, "multilevel_resolution");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/leiden_resolution");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "leiden_resolution", -1);
    }
    quick_snn_graph_throw(path, 1000, "leiden_resolution");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/parameters/walktrap_steps");
        auto phandle = handle.openGroup("snn_graph_cluster/parameters");
        quick_write_dataset(phandle, "walktrap_steps", -1);
    }
    quick_snn_graph_throw(path, 1000, "walktrap_steps");
}

TEST(SnnGraphClusterV3, ResultsFailed) {
    const std::string path = "TEST_snn_graph_cluster.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results");
    }
    quick_snn_graph_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
    }
    quick_snn_graph_throw(path, 1000, "'clusters' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
    }
    quick_snn_graph_throw(path, 500, "expected dimensions");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
        auto rhandle = handle.openGroup("snn_graph_cluster/results");
        quick_write_dataset(rhandle, "clusters", std::vector<int>(1000, -1));
    }
    quick_snn_graph_throw(path, 1000, "non-negative");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_snn_graph_cluster(handle, 1000);
        handle.unlink("snn_graph_cluster/results/clusters");
        auto rhandle = handle.openGroup("snn_graph_cluster/results");
        quick_write_dataset(rhandle, "clusters", std::vector<int>(1000, 10));
    }
    quick_snn_graph_throw(path, 1000, "represented at least once");
}
