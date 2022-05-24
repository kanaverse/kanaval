#include <gtest/gtest.h>
#include "kanaval/validate.hpp"
#include "H5Cpp.h"
#include "utils.h"
#include "helpers.h"

TEST(Overall, Basic) {
    const std::string path = "TEST_overall.h5";
    H5::H5File handle(path, H5F_ACC_TRUNC);
    int num_cells = 10;
    int num_genes = 1000;
    int filtered_cells = 8;
    int num_samples = 1;
    int num_clusters = 5;

    add_single_matrix(handle, "MatrixMarket", num_genes, num_cells);

    add_quality_control(handle, num_cells, num_samples, num_cells - filtered_cells);
    add_adt_quality_control(handle, num_cells, num_samples, num_cells - filtered_cells);

    add_normalization(handle);
    add_feature_selection(handle, num_genes);
    add_pca(handle, 50, filtered_cells);

    add_neighbor_index(handle);
    add_tsne(handle, filtered_cells);
    add_umap(handle, filtered_cells);

    add_choose_clustering(handle);
    add_kmeans_cluster(handle, filtered_cells, num_clusters);
    add_snn_graph_cluster(handle, filtered_cells, num_clusters);

    add_marker_detection(handle, num_genes, num_clusters);
    add_cell_labelling(handle, num_clusters);
    add_custom_selections(handle, num_genes, filtered_cells);

    EXPECT_NO_THROW(kanaval::validate(handle, true, latest));
}

TEST(Overall, Multiple) {
    const std::string path = "TEST_overall.h5";
    H5::H5File handle(path, H5F_ACC_TRUNC);
    int num_cells = 10;
    int num_genes = 1000;
    int filtered_cells = 8;
    int num_clusters = 5;

    int num_samples = add_multiple_matrices(handle, num_genes, num_cells);
    add_quality_control(handle, num_cells, num_samples, num_cells - filtered_cells);
    add_normalization(handle);
    add_feature_selection(handle, num_genes);
    add_pca(handle, 50, filtered_cells);

    add_neighbor_index(handle);
    add_tsne(handle, filtered_cells);
    add_umap(handle, filtered_cells);

    add_choose_clustering(handle);
    add_kmeans_cluster(handle, filtered_cells, num_clusters);
    add_snn_graph_cluster(handle, filtered_cells, num_clusters);

    add_marker_detection(handle, num_genes, num_clusters);
    add_cell_labelling(handle, num_clusters);
    add_custom_selections(handle, num_genes, filtered_cells);

    EXPECT_NO_THROW(kanaval::validate(handle, true, latest));
}
