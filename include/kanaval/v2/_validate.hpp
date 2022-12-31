#ifndef KANAVAL__VALIDATE_V2_HPP
#define KANAVAL__VALIDATE_V2_HPP

#include "inputs.hpp"

#include "quality_control.hpp"
#include "adt_quality_control.hpp"
#include "cell_filtering.hpp"

#include "normalization.hpp"
#include "adt_normalization.hpp"

#include "feature_selection.hpp"

#include "pca.hpp"
#include "adt_pca.hpp"
#include "combine_embeddings.hpp"
#include "batch_correction.hpp"

#include "neighbor_index.hpp"

#include "choose_clustering.hpp"
#include "kmeans_cluster.hpp"
#include "snn_graph_cluster.hpp"

#include "tsne.hpp"
#include "umap.hpp"

#include "marker_detection.hpp"
#include "custom_selections.hpp"
#include "cell_labelling.hpp"

#include <algorithm>

namespace kanaval {

namespace v2 {

void validate(const H5::H5File& handle, bool embedded, int version) {
    auto i_out = validate_inputs(handle, embedded, version);

    size_t rna_idx = std::find(i_out.modalities.begin(), i_out.modalities.end(), std::string("RNA")) - i_out.modalities.begin();
    size_t adt_idx = std::find(i_out.modalities.begin(), i_out.modalities.end(), std::string("ADT")) - i_out.modalities.begin();
    bool rna_in_use = rna_idx != i_out.modalities.size();
    bool adt_in_use = adt_idx != i_out.modalities.size();

    // Quality control.
    auto rna_filtered = validate_quality_control(handle, i_out.num_cells, i_out.num_samples, version);
    auto adt_filtered = validate_adt_quality_control(handle, i_out.num_cells, i_out.num_samples, adt_in_use, version);

    int num_qc_modalities = (rna_in_use && !rna_filtered.first) + (adt_in_use && !adt_filtered.first); // in use, not skipped.
    auto filtered_cells = validate_cell_filtering(handle, i_out.num_cells, num_qc_modalities, version);
    if (filtered_cells < 0) {
        filtered_cells = std::max(rna_filtered.second, adt_filtered.second);
    }

    // Normalization.
    validate_normalization(handle);
    validate_adt_normalization(handle, filtered_cells, adt_in_use, version);

    validate_feature_selection(handle, i_out.num_features[rna_idx]);

    // Dimensionality reduction.
    auto rna_pcs = validate_pca(handle, filtered_cells, version);
    auto adt_pcs = validate_adt_pca(handle, filtered_cells, adt_in_use, version);

    int total_pcs = (rna_in_use ? rna_pcs : 0) + (adt_in_use ? adt_pcs : 0);
    validate_combine_embeddings(handle, filtered_cells, i_out.modalities, total_pcs, version);
    validate_batch_correction(handle, total_pcs, filtered_cells, i_out.num_samples, version);

    validate_neighbor_index(handle);

    // Clustering.
    auto cluster_method = validate_choose_clustering(handle);
    int nclusters = 0;

    {
        bool is_snn = (cluster_method == "snn_graph");
        int snn_found = validate_snn_graph_cluster(handle, filtered_cells, is_snn);
        if (is_snn) {
            nclusters = snn_found;
        }
    }

    {
        bool is_kmeans = (cluster_method == "kmeans");
        int kmeans_found = validate_kmeans_cluster(handle, filtered_cells, is_kmeans);
        if (is_kmeans) {
            nclusters = kmeans_found;
        }
    }

    validate_tsne(handle, filtered_cells);
    validate_umap(handle, filtered_cells);

    validate_marker_detection(handle, nclusters, i_out.modalities, i_out.num_features, version);
    validate_custom_selections(handle, filtered_cells, i_out.modalities, i_out.num_features, version);
    validate_cell_labelling(handle, nclusters);
}

}

}

#endif
