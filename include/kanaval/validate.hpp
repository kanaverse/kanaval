#ifndef KANAVAL_VALIDATE_HPP
#define KANAVAL_VALIDATE_HPP

#include "inputs.hpp"
#include "quality_control.hpp"
#include "normalization.hpp"
#include "feature_selection.hpp"
#include "pca.hpp"
#include "neighbor_index.hpp"
#include "choose_clustering.hpp"
#include "kmeans_cluster.hpp"
#include "snn_graph_cluster.hpp"
#include "tsne.hpp"
#include "umap.hpp"
#include "marker_detection.hpp"
#include "custom_selections.hpp"
#include "cell_labelling.hpp"

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
 * - `normalization::validate()`
 * - `feature_selection::validate()`
 * - `pca::validate()`
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
 *
 * @param handle Open handle to a HDF5 file.
 * @param embedded Whether the data files are embedded.
 * @param version Version of the kana file.
 *
 * @return An error is raised if an invalid structure is detected.
 */
void validate(const H5::H5File& handle, bool embedded, int version = 1001000) {
    auto i_out = inputs::validate(handle, embedded, version);
    auto filtered_cells = quality_control::validate(handle, i_out.num_cells, i_out.num_samples);

    normalization::validate(handle);
    feature_selection::validate(handle, i_out.num_genes);
    pca::validate(handle, filtered_cells, version);
    neighbor_index::validate(handle);

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

    marker_detection::validate(handle, i_out.num_genes, nclusters);
    custom_selections::validate(handle, i_out.num_genes, filtered_cells);
    cell_labelling::validate(handle, nclusters);
}

}

#endif
