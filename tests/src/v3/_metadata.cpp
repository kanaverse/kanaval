#include <gtest/gtest.h>
#include "kanaval/v3/_metadata.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

H5::Group add_metadata(H5::H5File& handle) {
    auto mhandle = handle.createGroup("_metadata");
    quick_write_dataset(mhandle, "format_version", latest);
    quick_write_dataset(mhandle, "application_name", "bakana");
    quick_write_dataset(mhandle, "application_version", "1.1.1");
    return mhandle;
}

}

TEST(MetadataV3, AllOk) {
    const std::string path = "TEST_metadata.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_metadata(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_metadata(handle, latest);
    }
}

static void quick_metadata_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_metadata(handle, latest);
    }, msg);
}

TEST(MetadataV3, Failed) {
    const std::string path = "TEST_metadata.h5";

    // Inconsistent format version.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        auto mhandle = v3::add_metadata(handle);
        mhandle.unlink("format_version");
        quick_write_dataset(mhandle, "format_version", 12345);
    }
    quick_metadata_throw(path, "format_version");

    // Non-string application name.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        auto mhandle = v3::add_metadata(handle);
        mhandle.unlink("application_name");
        quick_write_dataset(mhandle, "application_name", 1);
    }
    quick_metadata_throw(path, "application_name");

    // Missing application version.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        auto mhandle = v3::add_metadata(handle);
        mhandle.unlink("application_version");
    }
    quick_metadata_throw(path, "application_version");
}
