#include <gtest/gtest.h>
#include "kanaval/v3/combine_embeddings.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_combine_embeddings(H5::H5File& handle, int num_cells, int total_pcs) {
    auto qhandle = handle.createGroup("combine_embeddings");
    auto phandle = qhandle.createGroup("parameters");
    quick_write_dataset(phandle, "rna_weight", 1.0);
    quick_write_dataset(phandle, "adt_weight", 1.2);
    quick_write_dataset(phandle, "crispr_weight", 0.0);
    quick_write_dataset(phandle, "approximate", 1);

    auto rhandle = qhandle.createGroup("results");
    H5::DataSpace space;
    std::vector<hsize_t> dims(2);
    dims[0] = num_cells;
    dims[1] = total_pcs;
    space.setExtentSimple(2, dims.data());
    rhandle.createDataSet("combined", H5::PredType::NATIVE_DOUBLE, space);
}

}

TEST(CombineEmbeddingsV3, AllOK) {
    const std::string path = "TEST_combine_embeddings.h5";

    // Works in simple mode.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 15);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_combine_embeddings(handle, 100, { { "RNA", 10 }, { "ADT", 5 } }, latest));
    }

    // Works if there's only one modality and no PCs.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 5);
        handle.unlink("combine_embeddings/results/combined");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate_combine_embeddings(handle, 100, { { "RNA", 5 } }, latest));

        // CRISPR is not used by default, so still just one modality.
        EXPECT_NO_THROW(kanaval::v3::validate_combine_embeddings(handle, 100, { { "RNA", 5 }, { "CRISPR", 100 } }, latest));
    }
}

static void quick_combine_throw(const std::string& path, int num_cells, const std::unordered_map<std::string, int>& modalities, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_combine_embeddings(handle, num_cells, modalities, latest);
    }, msg);
}

TEST(CombineEmbeddingsV3, ParametersFailed) {
    const std::string path = "TEST_combine_embeddings.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters");
    }
    quick_combine_throw(path, 100, {}, "'parameters' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters/approximate");
    }
    quick_combine_throw(path, 100, {}, "approximate");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters/rna_weight");
    }
    quick_combine_throw(path, 100, {}, "rna_weight");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/parameters/adt_weight");
        quick_write_dataset(handle, "combine_embeddings/parameters/adt_weight", "foo");
    }
    quick_combine_throw(path, 100, {}, "'adt_weight'");

    // Fails if no modality has non-zero weight.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
    }
    quick_combine_throw(path, 100, {}, "could not find");

}

TEST(CombineEmbeddingsV3, ResultsFailed) {
    const std::string path = "TEST_combine_embeddings.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/results");
    }
    quick_combine_throw(path, 100, { { "RNA", 10 } }, "'results' group");
    
    // Throws if multiple modalities are requested but 'combined' is absent.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
        handle.unlink("combine_embeddings/results/combined");
    }
    quick_combine_throw(path, 100, { { "RNA", 10 }, { "ADT", 5 } }, "combined");

    // Wrong dimensions.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_combine_embeddings(handle, 100, 10);
    }
    quick_combine_throw(path, 100, { { "RNA", 10 }, { "ADT", 5 } }, "dimensions");
}
