#include <gtest/gtest.h>
#include "kanaval/v2/feature_selection.hpp"
#include "../utils.h"
#include <iostream>

namespace v2 {

void add_feature_selection(H5::H5File& handle, int num_genes) {
    auto qhandle = handle.createGroup("feature_selection");

    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "span", 0.5);

    auto rhandle = qhandle.createGroup("results");
    quick_write_dataset(rhandle, "means", std::vector<double>(num_genes));
    quick_write_dataset(rhandle, "vars", std::vector<double>(num_genes));
    quick_write_dataset(rhandle, "fitted", std::vector<double>(num_genes));
    quick_write_dataset(rhandle, "resids", std::vector<double>(num_genes));

    return;
}

}

TEST(FeatureSelectionV2, AllOK) {
    const std::string path = "TEST_feature_selection.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 123);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate_feature_selection(handle, 123));
    }
}

static void quick_feat_throw(const std::string& path, int num_genes, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v2::validate_feature_selection(handle, num_genes);
    }, msg);
}

TEST(FeatureSelectionV2, ParametersFailed) {
    const std::string path = "TEST_feature_selection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 123);
        handle.unlink("feature_selection/parameters");
    }
    quick_feat_throw(path, 123, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 123);
        handle.unlink("feature_selection/parameters/span");
    }
    quick_feat_throw(path, 123, "'span' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 123);
        handle.unlink("feature_selection/parameters/span");
    
        auto phandle = handle.openGroup("feature_selection/parameters");
        quick_write_dataset(phandle, "span", -1.0);
    }
    quick_feat_throw(path, 123, "[0, 1]");
}

TEST(FeatureSelectionV2, ResultsFailed) {
    const std::string path = "TEST_feature_selection.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 123);
        handle.unlink("feature_selection/results");
    }
    quick_feat_throw(path, 123, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 123);
        handle.unlink("feature_selection/results/means");
    }
    quick_feat_throw(path, 123, "'means' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v2::add_feature_selection(handle, 100);
    }
    quick_feat_throw(path, 200, "'means' dataset");
}
