#include <gtest/gtest.h>
#include "kanaval/tsne.hpp"
#include "utils.h"
#include <iostream>

void add_tsne(H5::H5File& handle, int num_cells) {
    auto qhandle = handle.createGroup("tsne");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "perplexity", 30.0);
    quick_write_dataset(phandle, "iterations", 1000);
    quick_write_dataset(phandle, "animate", 1);

    auto rhandle = qhandle.createGroup("results");
    quick_write_dataset(rhandle, "x", std::vector<double>(num_cells));
    quick_write_dataset(rhandle, "y", std::vector<double>(num_cells));

    return;
}

TEST(TSNE, AllOK) {
    const std::string path = "TEST_tsne.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::tsne::validate(handle, 1000));
    }
}

void quick_tsne_throw(const std::string& path, int num_cells, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::tsne::validate(handle, num_cells);
    }, msg);
}

TEST(TSNE, ParametersFailed) {
    const std::string path = "TEST_tsne.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
        handle.unlink("tsne/parameters");
    }
    quick_tsne_throw(path, 1000, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
        handle.unlink("tsne/parameters/perplexity");
    }
    quick_tsne_throw(path, 1000, "'perplexity' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
        handle.unlink("tsne/parameters/perplexity");
    
        auto phandle = handle.openGroup("tsne/parameters");
        quick_write_dataset(phandle, "perplexity", -1.0);
    }
    quick_tsne_throw(path, 1000, "should be positive");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
        handle.unlink("tsne/parameters/iterations");
    
        auto phandle = handle.openGroup("tsne/parameters");
        quick_write_dataset(phandle, "iterations", -1);
    }
    quick_tsne_throw(path, 1000, "should be positive");
}

TEST(TSNE, ResultsFailed) {
    const std::string path = "TEST_tsne.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
        handle.unlink("tsne/results");
    }
    quick_tsne_throw(path, 1000, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_tsne(handle, 1000);
    }
    quick_tsne_throw(path, 500, "'x' dataset");
}
