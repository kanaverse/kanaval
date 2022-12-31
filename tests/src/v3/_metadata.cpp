#include <gtest/gtest.h>
#include "kanaval/v3/_metadata.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

H5::Group add__metadata(H5::H5File& handle) {
    auto mhandle = handle.createGroup("_metadata");
    quick_write_dataset(mhandle, "format_version", latest);
    quick_write_dataset(mhandle, "application_name", "bakana");
    quick_write_dataset(mhandle, "application_version", "1.1.1");
    return mhandle;
}

}

TEST(MetadataV3, AllOk) {
    const std::string path = "TEST__metadata.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add__metadata(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate__metadata(handle, latest);
    }
}

static void quick__metadata_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate__metadata(handle, latest);
    }, msg);
}

TEST(MetadataV3, Failed) {
    const std::string path = "TEST__metadata.h5";

    // Inconsistent format version.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        auto mhandle = v3::add__metadata(handle);
        mhandle.unlink("format_version");
        quick_write_dataset(mhandle, "format_version", 12345);
    }
    quick__metadata_throw(path, "format_version");

    // Non-string application name.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        auto mhandle = v3::add__metadata(handle);
        mhandle.unlink("application_name");
        quick_write_dataset(mhandle, "application_name", 1);
    }
    quick__metadata_throw(path, "application_name");

    // Missing application version.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        auto mhandle = v3::add__metadata(handle);
        mhandle.unlink("application_version");
    }
    quick__metadata_throw(path, "application_version");
}
