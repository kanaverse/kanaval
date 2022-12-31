#include <gtest/gtest.h>
#include "kanaval/v3/inputs.hpp"
#include "../utils.h"
#include <iostream>

namespace v3 {

void add_single_matrix(H5::H5File& handle, std::string mode = "MatrixMarket", int ngenes = 1000, int ncells = 100, int nblocks = 1) {
    auto ihandle = handle.createGroup("inputs");

    auto phandle = ihandle.createGroup("parameters");
    auto dhandle = phandle.createGroup("datasets");
    auto firsthandle = dhandle.createGroup("0");
    quick_write_dataset(firsthandle, "format", mode);
    quick_write_dataset(firsthandle, "name", "IMMA_FIRST");
    auto fihandle = firsthandle.createGroup("files");

    if (mode == "MatrixMarket") {
        auto fhandle0 = fihandle.createGroup("0");
        quick_write_dataset(fhandle0, "type", "mtx");
        quick_write_dataset(fhandle0, "name", "foo.mtx");
        quick_write_dataset(fhandle0, "size", 1);
        quick_write_dataset(fhandle0, "offset", 0);

        auto fhandle1 = fihandle.createGroup("1");
        quick_write_dataset(fhandle1, "type", "genes");
        quick_write_dataset(fhandle1, "name", "genes.tsv");
        quick_write_dataset(fhandle1, "size", 2);
        quick_write_dataset(fhandle1, "offset", 1);
    } else {
        auto fhandle = fihandle.createGroup("0");
        quick_write_dataset(fhandle, "type", "h5");
        quick_write_dataset(fhandle, "name", "foo.h5");
        quick_write_dataset(fhandle, "size", 1);
        quick_write_dataset(fhandle, "offset", 0);
    }

    auto rhandle = ihandle.createGroup("results");
    quick_write_dataset(rhandle, "num_cells", ncells);
    quick_write_dataset(rhandle, "num_blocks", nblocks);

    auto idhandle = rhandle.createGroup("feature_identities");
    std::vector<int> identities(1000);
    std::iota(identities.rbegin(), identities.rend(), 0); // reversed... outta control, bruh.
    quick_write_dataset(idhandle, "RNA", identities);

    return;
}

}

static void quick_input_throw(const std::string& path, std::string msg) {
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_inputs(handle, true, latest);
    }, msg);
}

TEST(InputsV3, SingleDatasetOK) {
    const std::string path = "TEST_inputs.h5";

    // Trying with MatrixMarket
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_features["RNA"], 1000);
        EXPECT_EQ(output.num_cells, 100);
        EXPECT_EQ(output.num_blocks, 1);
    }

    // Trying with HDF5.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle, "10X");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_features["RNA"], 1000);
        EXPECT_EQ(output.num_cells, 100);
        EXPECT_EQ(output.num_blocks, 1);
    }
}

TEST(InputsV3, Blocking) {
    const std::string path = "TEST_inputs.h5";

    // Fails at first, if we inject multiple blocks without a blocking factor.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle, "H5AD", 1000, 100, 5);
    }
    quick_input_throw(path, "should be equal to the number of datasets");

    // Adding back the blocking factor!
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto phandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(phandle, "block_factor", "FOO");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_features["RNA"], 1000);
        EXPECT_EQ(output.num_cells, 100);
        EXPECT_EQ(output.num_blocks, 5);
    }
}

TEST(InputsV3, Links) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle, "H5AD");
        auto zerohandle = handle.openGroup("inputs/parameters/datasets/0/files/0");
        zerohandle.unlink("size");
        zerohandle.unlink("offset");
    }
    quick_throw([&]() -> void {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_inputs(handle, false, latest);
    }, "id");

    // Adding back in the ID.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto zerohandle = handle.openGroup("inputs/parameters/datasets/0/files/0");
        quick_write_dataset(zerohandle, "id", "FOO");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, false, latest);
        EXPECT_EQ(output.num_features["RNA"], 1000);
        EXPECT_EQ(output.num_cells, 100);
        EXPECT_EQ(output.num_blocks, 1);
    }
}

TEST(InputsV3, MultiModal) {
    const std::string path = "TEST_inputs.h5";

    // Adding ADT.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results/feature_identities");
        quick_write_dataset(rhandle, "ADT", std::vector<int>{2,4,6,8});
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);

        EXPECT_EQ(output.num_features["RNA"], 1000);
        EXPECT_EQ(output.num_features["ADT"], 4);
        EXPECT_EQ(output.num_cells, 100);
        EXPECT_EQ(output.num_blocks, 1);
    }

    // Adding CRISPR.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto rhandle = handle.openGroup("inputs/results/feature_identities");
        quick_write_dataset(rhandle, "CRISPR", std::vector<int>{7,5,1,3,9});
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_features["RNA"], 1000);
        EXPECT_EQ(output.num_features["ADT"], 4);
        EXPECT_EQ(output.num_features["CRISPR"], 5);
    }
}

TEST(InputsV3, Options) {
    const std::string path = "TEST_inputs.h5";

    // Works okay for dataset inputs.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto zerohandle = handle.openGroup("inputs/parameters/datasets/0");
        quick_write_dataset(zerohandle, "options", "{ \"featureTypeRnaName\": \"Gene Expression\", \"assayIndex\": 0 }");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_inputs(handle, true, latest);
    }

    // Fails for invalid JSON.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto zerohandle = handle.openGroup("inputs/parameters/datasets/0");
        quick_write_dataset(zerohandle, "options", "asdasdasd");
    }
    quick_input_throw(path, "options should be a valid JSON string");

    // Fails for non-object JSON.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto zerohandle = handle.openGroup("inputs/parameters/datasets/0");
        quick_write_dataset(zerohandle, "options", "[1, 2, 3]");
    }
    quick_input_throw(path, "options should encode a JSON object");
}

TEST(InputsV3, ParametersFail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/parameters/datasets/0/format");
    }
    quick_input_throw(path, "'format' dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/parameters/datasets/0/name");
    }
    quick_input_throw(path, "'name'");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/parameters/datasets/0/files/0/offset");
    }
    quick_input_throw(path, "'offset' dataset");

    // Checking what happens if we throw in an extra file with overlapping offsets.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto fhandle2 = handle.createGroup("inputs/parameters/datasets/0/files/2");
        quick_write_dataset(fhandle2, "type", "annotations");
        quick_write_dataset(fhandle2, "name", "anno.tsv");
        quick_write_dataset(fhandle2, "size", 2);
        quick_write_dataset(fhandle2, "offset", 1);
    }
    quick_input_throw(path, "not contiguous");
}

TEST(InputsV3, ResultsFail) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/results");
    }
    quick_input_throw(path, "'results' group");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/results/num_cells");
        auto rhandle = handle.openGroup("inputs/results");
        quick_write_dataset(rhandle, "num_cells", -1);
    }
    quick_input_throw(path, "should be a positive integer");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/results/num_blocks");
        auto rhandle = handle.openGroup("inputs/results");
        quick_write_dataset(rhandle, "num_blocks", -1);
    }
    quick_input_throw(path, "should be a positive integer");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        handle.unlink("inputs/results/feature_identities/RNA");
    }
    quick_input_throw(path, "at least one recognized modality");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results/feature_identities");
        rhandle.unlink("RNA");
        quick_write_dataset(rhandle, "RNA", std::vector<int>{ 0, 1, -1 });
    }
    quick_input_throw(path, "non-negative");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results/feature_identities");
        rhandle.unlink("RNA");
        quick_write_dataset(rhandle, "RNA", std::vector<int>(1000));
    }
    quick_input_throw(path, "unique");
}

TEST(InputsV3, FeatureNames) {
    const std::string path = "TEST_inputs.h5";

    // Works as expected.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results");
        auto nhandle = rhandle.createGroup("feature_names");
        std::vector<std::string> identities(1000, "A");
        quick_write_dataset(nhandle, "RNA", identities);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_inputs(handle, true, latest);
    }

    // Still works for multiple modalities.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results");
        auto idhandle = rhandle.openGroup("feature_identities");
        quick_write_dataset(rhandle, "ADT", std::vector<int>{2,4,6,8});

        auto nhandle = rhandle.createGroup("feature_names");
        quick_write_dataset(nhandle, "RNA", std::vector<std::string>(1000, "A"));
        quick_write_dataset(nhandle, "ADT", std::vector<std::string>(4, "B"));
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        kanaval::v3::validate_inputs(handle, true, latest);
    }

    // Fails if the lengths don't match up.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results");
        auto nhandle = rhandle.createGroup("feature_names");
        std::vector<std::string> identities(999, "A");
        quick_write_dataset(nhandle, "RNA", identities);
    }
    quick_input_throw(path, "should be of length equal to 'feature_identities/RNA'");

    // Fails on invalid feature_names.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results");
        auto nhandle = rhandle.createGroup("feature_names");
        quick_write_dataset(nhandle, "RNA", std::vector<int>(10));
    }
    quick_input_throw(path, "should be a string dataset");

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto rhandle = handle.openGroup("inputs/results");
        auto nhandle = rhandle.createGroup("feature_names");
        quick_write_dataset(nhandle, "RNA", "FOOBAR");
    }
    quick_input_throw(path, "should be a 1-dimensional string dataset");
}

namespace v3 {

int add_multiple_matrices(H5::H5File& handle, int ngenes = 500, int ncells = 100) {
    auto ihandle = handle.createGroup("inputs");

    auto phandle = ihandle.createGroup("parameters");
    auto dhandle = phandle.createGroup("datasets");

    // Adding two input datasets.
    {
        auto dhandle0 = dhandle.createGroup("0");
        quick_write_dataset(dhandle0, "format", "10X");
        quick_write_dataset(dhandle0, "name", "A");
        auto fihandle = dhandle0.createGroup("files");

        auto fhandle = fihandle.createGroup("0");
        quick_write_dataset(fhandle, "type", "h5");
        quick_write_dataset(fhandle, "name", "foo.h5");
        quick_write_dataset(fhandle, "size", 3);
        quick_write_dataset(fhandle, "offset", 0);
    }

    {
        auto dhandle1 = dhandle.createGroup("1");
        quick_write_dataset(dhandle1, "format", "MatrixMarket");
        quick_write_dataset(dhandle1, "name", "B");
        auto fihandle = dhandle1.createGroup("files");

        auto fhandle0 = fihandle.createGroup("0");
        quick_write_dataset(fhandle0, "type", "mtx");
        quick_write_dataset(fhandle0, "name", "foo.mtx");
        quick_write_dataset(fhandle0, "size", 1);
        quick_write_dataset(fhandle0, "offset", 3);

        auto fhandle1 = fihandle.createGroup("1");
        quick_write_dataset(fhandle1, "type", "genes");
        quick_write_dataset(fhandle1, "name", "genes.tsv");
        quick_write_dataset(fhandle1, "size", 2);
        quick_write_dataset(fhandle1, "offset", 4);
    }

    auto rhandle = ihandle.createGroup("results");
    quick_write_dataset(rhandle, "num_cells", ncells);
    quick_write_dataset(rhandle, "num_blocks", 2);

    auto idhandle = rhandle.createGroup("feature_identities");
    std::vector<int> identities(ngenes);
    std::iota(identities.begin(), identities.end(), 0); 
    quick_write_dataset(idhandle, "RNA", identities);

    return 2;
}

}

TEST(InputsV3, MultiDatasetOK) {
    const std::string path = "TEST_inputs.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_multiple_matrices(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_features["RNA"], 500);
        EXPECT_EQ(output.num_cells, 100);
        EXPECT_EQ(output.num_blocks, 2);
    }

    // Fails if there's the wrong number of blocks.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_multiple_matrices(handle);
        auto rhandle = handle.openGroup("inputs/results");
        rhandle.unlink("num_blocks");
        quick_write_dataset(rhandle, "num_blocks", 3);
    }
    quick_input_throw(path, "should be equal to the number of datasets");

    // Fails if the names are not unique.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_multiple_matrices(handle);
        auto dhandle = handle.openGroup("inputs/parameters/datasets/1");
        dhandle.unlink("name");
        quick_write_dataset(dhandle, "name", "A");
    }
    quick_input_throw(path, "duplicate dataset name");
}

TEST(InputsV3, SubsettingIndices) {
    const std::string path = "TEST_inputs.h5";

    // Works correctly.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto ihandle = handle.openGroup("inputs");
        auto subhandle = ihandle.createGroup("parameters/subset");
        auto subcellhandle = subhandle.createGroup("cells");
        quick_write_dataset(subcellhandle, "indices", std::vector<int>{1,2,3});

        auto rhandle = ihandle.openGroup("results");
        rhandle.unlink("num_cells");
        quick_write_dataset(rhandle, "num_cells", 3);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_cells, 3);
    }

    // Fails if the number of cells is not right.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto rhandle = handle.openGroup("inputs/results");
        rhandle.unlink("num_cells");
        quick_write_dataset(rhandle, "num_cells", 100);
    }
    quick_input_throw(path, "should be equal to the length of 'parameters/subset/indices'");

    // Fails if the indices are negative.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto subhandle = handle.openGroup("inputs/parameters/subset/cells");
        subhandle.unlink("indices");
        quick_write_dataset(subhandle, "indices", std::vector<int>{-1,0,1});
    }
    quick_input_throw(path, "negative");
}

TEST(InputsV3, SubsettingFieldValues) {
    const std::string path = "TEST_inputs.h5";

    // No check on the number of cells if we're using field/value.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto ihandle = handle.openGroup("inputs");
        auto shandle = ihandle.createGroup("parameters/subset");
        auto subhandle = shandle.createGroup("cells");
        quick_write_dataset(subhandle, "field", "FOO");
        quick_write_dataset(subhandle, "values", std::vector<std::string>{"BAR"});
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_cells, 100);
    }

    // Fails if field or value are missing.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto subhandle = handle.openGroup("inputs/parameters/subset/cells");
        subhandle.unlink("field");
    }
    quick_input_throw(path, "field");
}

TEST(InputsV3, SubsettingRanges) {
    const std::string path = "TEST_inputs.h5";

    H5::DataSpace space;
    std::vector<hsize_t> nranges{3, 2};
    space.setExtentSimple(2, nranges.data());

    // No check on the number of cells if we're using field/ranges.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        v3::add_single_matrix(handle);
        auto ihandle = handle.openGroup("inputs");
        auto shandle = ihandle.createGroup("parameters/subset");
        auto subhandle = shandle.createGroup("cells");
        quick_write_dataset(subhandle, "field", "FOO");

        auto dhandle = subhandle.createDataSet("ranges", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> ranges { 0.2, 0.5, 1, 1.2, 1.25, 3 };
        dhandle.write(ranges.data(), H5::PredType::NATIVE_DOUBLE);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        auto output = kanaval::v3::validate_inputs(handle, true, latest);
        EXPECT_EQ(output.num_cells, 100);
    }

    // Fails if ranges is not the right shape.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto subhandle = handle.openGroup("inputs/parameters/subset/cells");
        subhandle.unlink("ranges");
        quick_write_dataset(subhandle, "ranges", std::vector<double>{ 1, 2, 3 });
    }
    quick_input_throw(path, "2-dimensional");

    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto subhandle = handle.openGroup("inputs/parameters/subset/cells");
        subhandle.unlink("ranges");

        H5::DataSpace space;
        std::vector<hsize_t> nranges{2, 3};
        space.setExtentSimple(2, nranges.data());

        auto dhandle = subhandle.createDataSet("ranges", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> ranges { 0.2, 0.5, 1, 1.2, 1.25, 3 };
        dhandle.write(ranges.data(), H5::PredType::NATIVE_DOUBLE);
    }
    quick_input_throw(path, "two columns");

    // Fails if ranges is not sorted.
    {
        H5::H5File handle(path, H5F_ACC_RDWR);
        auto subhandle = handle.openGroup("inputs/parameters/subset/cells");
        subhandle.unlink("ranges");

        auto dhandle = subhandle.createDataSet("ranges", H5::PredType::NATIVE_DOUBLE, space);
        std::vector<double> ranges { 0.2, 1.5, 1, 1.2, 1.25, 3 };
        dhandle.write(ranges.data(), H5::PredType::NATIVE_DOUBLE);
    }
    quick_input_throw(path, "sorted");
}
