#include <gtest/gtest.h>
#include "kanaval/combine_embeddings.hpp"
#include "utils.h"
#include <iostream>

void add_combine_embeddings(H5::H5File& handle, int num_cells, int total_pcs) {
    auto qhandle = handle.createGroup("combine_embeddings");
    auto phandle = qhandle.createGroup("parameters");
    phandle.createGroup("weights");
    quick_write_dataset(phandle, "approximate", 1);

    auto rhandle = qhandle.createGroup("results");
    H5::DataSpace space;
    std::vector<hsize_t> dims(2);
    dims[0] = num_cells;
    dims[1] = total_pcs;
    space.setExtentSimple(2, dims.data());
    rhandle.createDataSet("combined", H5::PredType::NATIVE_DOUBLE, space);
}

TEST(CombineEmbeddings, AllOK) {
    const std::string path = "TEST_combine_embeddings.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::combine_embeddings::validate(handle, 100, { "RNA", "ADT" }, 10, latest));
    }

    // Works if there's only one modality and no PCs.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 5);
        handle.unlink("combine_embeddings/results/combined");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::combine_embeddings::validate(handle, 100, { "RNA" }, 5, latest));
    }

    // Works if the modalities are weighted.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 11);
        quick_write_dataset(handle, "combine_embeddings/parameters/weights/RNA", 2.0);
        quick_write_dataset(handle, "combine_embeddings/parameters/weights/ADT", 0.5);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::combine_embeddings::validate(handle, 100, { "RNA", "ADT" }, 11, latest));
    }
}

void quick_combine_throw(const std::string& path, int num_cells, const std::vector<std::string>& modalities, int num_pcs, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::combine_embeddings::validate(handle, num_cells, modalities, num_pcs, latest);
    }, msg);
}

TEST(CombineEmbeddings, ParametersFailed) {
    const std::string path = "TEST_combine_embeddings.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters");
    }
    quick_combine_throw(path, 100, { "RNA" }, 10, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters/approximate");
    }
    quick_combine_throw(path, 100, { "RNA" }, 10, "approximate");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters/weights");
    }
    quick_combine_throw(path, 100, { "RNA" }, 10, "weights");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
        quick_write_dataset(handle, "combine_embeddings/parameters/weights/RNA", 2.0);
    }
    quick_combine_throw(path, 100, { "RNA", "ADT" }, 10, "ADT");
}

TEST(CombineEmbeddings, ResultsFailed) {
    const std::string path = "TEST_combine_embeddings.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/results");
    }
    quick_combine_throw(path, 100, { "RNA", "ADT" }, 10, "'results' group");
    
    // Throws if multiple modalities are requested but 'combined' is absent.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/results/combined");
    }
    quick_combine_throw(path, 100, { "RNA", "ADT" }, 10, "combined");

    // Wrong dimensions.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        add_combine_embeddings(handle, 100, 5);
        handle.unlink("combine_embeddings/results/combined");
    }
    quick_combine_throw(path, 100, { "RNA", "ADT" }, 10, "combined");
}
