#ifndef KANAVAL_VALIDATE_HPP
#define KANAVAL_VALIDATE_HPP

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

/**
 * @file validate.hpp
 *
 * @brief Validate the embedded state file of a kana file.
 */

/** 
 * @namespace kanaval
 * @brief Utilities for validating kana files.
 */
namespace kanaval {

/**
 * Validate the analysis state HDF5 file embedded inside a `*.kana` file.
 * This calls the validation functions for all of the individual analysis steps, namely:
 *
 * - `inputs::validate()`
 * - `quality_control::validate()`
 * - `adt_quality_control::validate()`
 * - `normalization::validate()`
 * - `adt_normalization::validate()`
 * - `feature_selection::validate()`
 * - `pca::validate()`
 * - `adt_pca::validate()`
 * - `combine_embeddings::validate()`
 * - `batch_correction::validate()`
 * - `neighbor_index::validate()`
 * - `choose_clustering::validate()`
 * - `kmeans_cluster::validate()`
 * - `snn_graph_cluster::validate()`
 * - `tsne::validate()`
 * - `umap::validate()`
 * - `marker_detection::validate()`
 * - `custom_selections::validate()`
 * - `cell_labelling::validate()`
 *
 * See the documentation for each individual function for more details on the expected structure of the state file.
 * An error is raised if an invalid structure is detected in any step.
 *
 * @param handle Open handle to a HDF5 file.
 * @param embedded Whether the data files are embedded.
 * @param version Version of the kana file.
 */
void validate(const H5::H5File& handle, bool embedded, int version) {
    auto i_out = inputs::validate(handle, embedded, version);

    size_t rna_idx = std::find(i_out.modalities.begin(), i_out.modalities.end(), std::string("RNA")) - i_out.modalities.begin();
    size_t adt_idx = std::find(i_out.modalities.begin(), i_out.modalities.end(), std::string("ADT")) - i_out.modalities.begin();
    bool rna_in_use = rna_idx != i_out.modalities.size();
    bool adt_in_use = adt_idx != i_out.modalities.size();

    // Quality control.
    auto rna_filtered = quality_control::validate(handle, i_out.num_cells, i_out.num_samples, version);
    auto adt_filtered = adt_quality_control::validate(handle, i_out.num_cells, i_out.num_samples, adt_in_use, version);

    int num_qc_modalities = (rna_in_use && !rna_filtered.first) + (adt_in_use && !adt_filtered.first); // in use, not skipped.
    auto filtered_cells = cell_filtering::validate(handle, i_out.num_cells, num_qc_modalities, version);
    if (filtered_cells < 0) {
        filtered_cells = std::max(rna_filtered.second, adt_filtered.second);
    }

    // Normalization.
    normalization::validate(handle);
    adt_normalization::validate(handle, filtered_cells, adt_in_use, version);

    feature_selection::validate(handle, i_out.num_features[rna_idx]);

    // Dimensionality reduction.
    auto rna_pcs = pca::validate(handle, filtered_cells, version);
    auto adt_pcs = adt_pca::validate(handle, filtered_cells, adt_in_use, version);

    int total_pcs = (rna_in_use ? rna_pcs : 0) + (adt_in_use ? adt_pcs : 0);
    kanaval::combine_embeddings::validate(handle, filtered_cells, i_out.modalities, total_pcs, version);
    batch_correction::validate(handle, total_pcs, filtered_cells, i_out.num_samples, version);

    neighbor_index::validate(handle);

    // Clustering.
    auto cluster_method = choose_clustering::validate(handle);
    int nclusters = 0;

    {
        bool is_snn = (cluster_method == "snn_graph");
        int snn_found = snn_graph_cluster::validate(handle, filtered_cells, is_snn);
        if (is_snn) {
            nclusters = snn_found;
        }
    }

    {
        bool is_kmeans = (cluster_method == "kmeans");
        int kmeans_found = kmeans_cluster::validate(handle, filtered_cells, is_kmeans);
        if (is_kmeans) {
            nclusters = kmeans_found;
        }
    }

    tsne::validate(handle, filtered_cells);
    umap::validate(handle, filtered_cells);

    marker_detection::validate(handle, nclusters, i_out.modalities, i_out.num_features, version);
    custom_selections::validate(handle, filtered_cells, i_out.modalities, i_out.num_features, version);
    cell_labelling::validate(handle, nclusters);
}

}

#endif
