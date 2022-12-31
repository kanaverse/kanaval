#include <gtest/gtest.h>
#include "kanaval/v3/_validate.hpp"
#include "H5Cpp.h"
#include "../utils.h"
#include "helpers.h"
#include "../v2/helpers.h"

static void spawn_full(H5::H5File& handle, int num_blocks = 1) {
    int num_cells = 20;
    int num_genes = 1000;
    int filtered_cells = 15;
    int num_clusters = 5;

    v3::add_single_matrix(handle, "MatrixMarket", num_genes, num_cells, num_blocks);
    auto rhandle = handle.openGroup("inputs/results/feature_identities");
    quick_write_dataset(rhandle, "ADT", std::vector<int>{2,4,6,8});
    quick_write_dataset(rhandle, "CRISPR", std::vector<int>{1,3,5,7,9,11});

    // Small differences to check for differences in handling between modalities.
    v3::add_rna_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells - 2);
    v3::add_adt_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells - 1); 
    v3::add_crispr_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells - 3);
    v3::add_cell_filtering(handle, num_cells, num_cells - filtered_cells);

    v3::add_rna_normalization(handle);
    v2::add_adt_normalization(handle, filtered_cells);
    v3::add_crispr_normalization(handle);

    v2::add_feature_selection(handle, num_genes);

    // Again, some small differences to check for differences in handling between modalities.
    int num_rna_pcs = 20, num_adt_pcs = 10, num_crispr_pcs = 5;
    v3::add_rna_pca(handle, num_rna_pcs, filtered_cells);
    v2::add_adt_pca(handle, num_adt_pcs, filtered_cells); 
    v3::add_crispr_pca(handle, num_crispr_pcs, filtered_cells); 

    int total_pcs = num_rna_pcs + num_adt_pcs;
    v3::add_combine_embeddings(handle, filtered_cells, total_pcs);
    v2::add_batch_correction(handle, filtered_cells, total_pcs);

    v2::add_neighbor_index(handle);
    v2::add_tsne(handle, filtered_cells);
    v2::add_umap(handle, filtered_cells);

    v2::add_choose_clustering(handle);
    v2::add_kmeans_cluster(handle, filtered_cells, num_clusters);
    v3::add_snn_graph_cluster(handle, filtered_cells, num_clusters);

    std::unordered_map<std::string, int> num_features { { "RNA", num_genes }, { "ADT", 4 }, { "CRISPR", 6 }};
    v3::add_marker_detection(handle, num_features, num_clusters);
    v3::add_custom_selections(handle, num_features, filtered_cells);

    v2::add_cell_labelling(handle, num_clusters);

    v3::add_metadata(handle);
}

TEST(OverallV3, MultiModalOk) {
    const std::string path = "TEST_overall.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        spawn_full(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate(handle, true, latest));
    }
}

TEST(OverallV3, MissingFail) {
    const std::string path = "TEST_overall.h5";

    std::vector<std::string> group {
        "inputs",
        "rna_quality_control",
        "adt_quality_control",
        "crispr_quality_control",
        "cell_filtering",
        "rna_normalization",
        "adt_normalization",
        "crispr_normalization",
        "feature_selection",
        "rna_pca",
        "adt_pca",
        "crispr_pca",
        "combine_embeddings",
        "batch_correction",
        "neighbor_index",
        "tsne",
        "umap",
        "kmeans_cluster",
        "snn_graph_cluster",
        "marker_detection",
        "custom_selections",
        "cell_labelling",
        "_metadata"
    };

    // Checking that we correctly throw upon deleting each group.
    for (const auto& g : group) {
        {
            H5::H5File handle(path, H5F_ACC_TRUNC);
            spawn_full(handle);
            handle.unlink(g);
        }
        quick_throw([&]() -> void {
            H5::H5File handle(path, H5F_ACC_RDONLY);
            kanaval::v3::validate(handle, true, latest);
        }, g);
    }
}

TEST(OverallV3, UniModalOk) {
    const std::string path = "TEST_overall.h5";

    // Pure RNA experiment.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);

        int num_cells = 20;
        int num_genes = 1000;
        int filtered_cells = 15;
        int num_clusters = 5;
        int num_blocks = 1;

        v3::add_single_matrix(handle, "MatrixMarket", num_genes, num_cells);

        v3::add_rna_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells);
        v3::add_adt_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells); 
        handle.unlink("adt_quality_control/results/metrics");
        v3::add_crispr_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells);
        handle.unlink("crispr_quality_control/results/metrics");
        v3::add_cell_filtering(handle, num_cells, num_cells - filtered_cells);
        handle.unlink("cell_filtering/results/discards");

        v3::add_rna_normalization(handle);
        v2::add_adt_normalization(handle, filtered_cells);
        handle.unlink("adt_normalization/results/size_factors");
        v3::add_crispr_normalization(handle);

        v2::add_feature_selection(handle, num_genes);

        // Again, some small differences to check for differences in handling between modalities.
        int num_rna_pcs = 20, num_adt_pcs = 10, num_crispr_pcs = 5;
        v3::add_rna_pca(handle, num_rna_pcs, filtered_cells);
        v2::add_adt_pca(handle, num_adt_pcs, filtered_cells); 
        handle.unlink("adt_pca/results/pcs");
        v3::add_crispr_pca(handle, num_crispr_pcs, filtered_cells); 
        handle.unlink("crispr_pca/results/pcs");

        v3::add_combine_embeddings(handle, filtered_cells, num_rna_pcs);
        handle.unlink("combine_embeddings/results/combined");
        v2::add_batch_correction(handle, filtered_cells, num_rna_pcs);
        handle.unlink("batch_correction/results/corrected");

        v2::add_neighbor_index(handle);
        v2::add_tsne(handle, filtered_cells);
        v2::add_umap(handle, filtered_cells);

        v2::add_choose_clustering(handle);
        v2::add_kmeans_cluster(handle, filtered_cells, num_clusters);
        v3::add_snn_graph_cluster(handle, filtered_cells, num_clusters);

        std::unordered_map<std::string, int> num_features { { "RNA", num_genes } };
        v3::add_marker_detection(handle, num_features, num_clusters);
        v3::add_custom_selections(handle, num_features, filtered_cells);

        v2::add_cell_labelling(handle, num_clusters);

        v3::add_metadata(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate(handle, true, latest));
    }
}

TEST(OverallV3, BlockedOk) {
    const std::string path = "TEST_overall.h5";

    {
        H5::H5File handle(path, H5F_ACC_TRUNC);
        spawn_full(handle, /* num_blocks = */ 2);
        auto pihandle = handle.openGroup("inputs/parameters");
        quick_write_dataset(pihandle, "block_factor", "FOO");
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate(handle, true, latest));
    }
}

TEST(OverallV3, MultiDatasetOk) {
    const std::string path = "TEST_overall.h5";

    // RNA and ADT only.
    {
        H5::H5File handle(path, H5F_ACC_TRUNC);

        int num_cells = 20;
        int num_genes = 1000;
        int filtered_cells = 15;
        int num_clusters = 5;

        int num_blocks = v3::add_multiple_matrices(handle, num_genes, num_cells);
        auto rhandle = handle.openGroup("inputs/results/feature_identities");
        quick_write_dataset(rhandle, "ADT", std::vector<int>{2,4,6,8});

        v3::add_rna_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells);
        v3::add_adt_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells); 
        v3::add_crispr_quality_control(handle, num_cells, num_blocks, num_cells - filtered_cells);
        handle.unlink("crispr_quality_control/results/metrics");
        v3::add_cell_filtering(handle, num_cells, num_cells - filtered_cells);

        v3::add_rna_normalization(handle);
        v2::add_adt_normalization(handle, filtered_cells);
        v3::add_crispr_normalization(handle);

        v2::add_feature_selection(handle, num_genes);

        // Again, some small differences to check for differences in handling between modalities.
        int num_rna_pcs = 20, num_adt_pcs = 10, num_crispr_pcs = 5;
        v3::add_rna_pca(handle, num_rna_pcs, filtered_cells);
        v2::add_adt_pca(handle, num_adt_pcs, filtered_cells); 
        v3::add_crispr_pca(handle, num_crispr_pcs, filtered_cells); 
        handle.unlink("crispr_pca/results/pcs");

        int total_pcs = num_rna_pcs + num_adt_pcs;
        v3::add_combine_embeddings(handle, filtered_cells, total_pcs);
        v2::add_batch_correction(handle, filtered_cells, total_pcs);

        v2::add_neighbor_index(handle);
        v2::add_tsne(handle, filtered_cells);
        v2::add_umap(handle, filtered_cells);

        v2::add_choose_clustering(handle);
        v2::add_kmeans_cluster(handle, filtered_cells, num_clusters);
        v3::add_snn_graph_cluster(handle, filtered_cells, num_clusters);

        std::unordered_map<std::string, int> num_features { { "RNA", num_genes }, { "ADT", 4 } };
        v3::add_marker_detection(handle, num_features, num_clusters);
        v3::add_custom_selections(handle, num_features, filtered_cells);

        v2::add_cell_labelling(handle, num_clusters);

        v3::add_metadata(handle);
    }
    {
        H5::H5File handle(path, H5F_ACC_RDONLY);
        EXPECT_NO_THROW(kanaval::v3::validate(handle, true, latest));
    }
}

