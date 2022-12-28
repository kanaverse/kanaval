#include <gtest/gtest.h>
#include "kanaval/v2/validate.hpp"
#include "H5Cpp.h"
#include "../utils.h"
#include "helpers.h"

static void spawn(H5::H5File& handle, bool multi_matrix = false, bool include_adts = false) {
    int num_cells = 10;
    int num_genes = 1000;
    int filtered_cells = 8;
    int num_clusters = 5;
    int num_samples = 1;

    if (multi_matrix) {
        num_samples = v2::add_multiple_matrices(handle, num_genes, num_cells);
    } else {
        v2::add_single_matrix(handle, "MatrixMarket", num_genes, num_cells);
    }
    if (include_adts) {
        // Filling in ADTs.
        quick_write_dataset(handle, "inputs/results/num_features/ADT", 4);
        quick_write_dataset(handle, "inputs/results/identities/ADT", std::vector<int>{2,4,6,8});
    }

    v2::add_quality_control(handle, num_cells, num_samples, num_cells - filtered_cells);
    v2::add_adt_quality_control(handle, num_cells, num_samples, num_cells - filtered_cells - 1); // small difference to check for differences in handling between RNA/ADT.
    v2::add_cell_filtering(handle, num_cells, num_cells - filtered_cells);

    v2::add_normalization(handle);
    v2::add_adt_normalization(handle, filtered_cells);

    v2::add_feature_selection(handle, num_genes);

    int num_pcs = 20, num_adt_pcs = 10, total_pcs = (include_adts ? num_pcs + num_adt_pcs : num_pcs);
    v2::add_pca(handle, num_pcs, filtered_cells);
    v2::add_adt_pca(handle, num_adt_pcs, filtered_cells); 
    v2::add_combine_embeddings(handle, filtered_cells, total_pcs);
    v2::add_batch_correction(handle, filtered_cells, total_pcs);

    v2::add_neighbor_index(handle);
    v2::add_tsne(handle, filtered_cells);
    v2::add_umap(handle, filtered_cells);

    v2::add_choose_clustering(handle);
    v2::add_kmeans_cluster(handle, filtered_cells, num_clusters);
    v2::add_snn_graph_cluster(handle, filtered_cells, num_clusters);

    if (include_adts) {
        v2::add_marker_detection(handle, { num_genes, 4 }, num_clusters, { "RNA", "ADT" });
        v2::add_custom_selections(handle, { "RNA", "ADT" }, { num_genes, 4 }, filtered_cells);
    } else {
        v2::add_marker_detection(handle, num_genes, num_clusters);
        v2::add_custom_selections(handle, num_genes, filtered_cells);
    }

    v2::add_cell_labelling(handle, num_clusters);
}

TEST(Overall, Single) {
    const std::string path = "TEST_overall.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        spawn(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate(handle, true, latest));
    }
}

TEST(Overall, SingleAdts) {
    const std::string path = "TEST_overall.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        spawn(handle, false, true);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate(handle, true, latest));
    }
}

TEST(Overall, Multiple) {
    const std::string path = "TEST_overall.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        spawn(handle, true);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate(handle, true, latest));
    }
}

TEST(Overall, MultipleAdts) {
    const std::string path = "TEST_overall.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        spawn(handle, true, true);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v2::validate(handle, true, latest));
    }
}
