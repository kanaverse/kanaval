#include <gtest/gtest.h>
#include "kanaval/choose_clustering.hpp"
#include "utils.h"
#include <iostream>

void add_choose_clustering(H5::H5File& handle) {
    auto qhandle = handle.createGroup("choose_clustering");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "method", "kmeans");
    qhandle.createGroup("results");
    return;
}

TEST(ChooseClustering, AllOK) {
    const std::string path = "TEST_choose_clustering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_choose_clustering(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::choose_clustering::validate(handle));
    }
}

void quick_choice_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::choose_clustering::validate(handle);
    }, msg);
}

TEST(ChooseClustering, ParametersFailed) {
    const std::string path = "TEST_choose_clustering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_choose_clustering(handle);
        handle.unlink("choose_clustering/parameters");
    }
    quick_choice_throw(path, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_choose_clustering(handle);
        handle.unlink("choose_clustering/parameters/method");

        auto phandle = handle.openGroup("choose_clustering/parameters");
        quick_write_dataset(phandle, "method", "asdasd");
    }
    quick_choice_throw(path, "'method' should be");
}

TEST(ChooseClustering, ResultsFailed) {
    const std::string path = "TEST_choose_clustering.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_choose_clustering(handle);
        handle.unlink("choose_clustering/results");
    }
    quick_choice_throw(path, "'results' group");
}
