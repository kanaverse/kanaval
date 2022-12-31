#ifndef KANAVAL_VALIDATE_V3_HPP
#define KANAVAL_VALIDATE_V3_HPP

#include "inputs.hpp"

#include "rna_quality_control.hpp"
#include "adt_quality_control.hpp"
#include "crispr_quality_control.hpp"
#include "cell_filtering.hpp"

#include "rna_normalization.hpp"
#include "../v2/adt_normalization.hpp"
#include "crispr_normalization.hpp"

#include "feature_selection.hpp"

#include "rna_pca.hpp"
#include "../v2/adt_pca.hpp"
#include "crispr_pca.hpp"
#include "combine_embeddings.hpp"
#include "../v2/batch_correction.hpp"

#include "../v2/neighbor_index.hpp"

#include "../v2/choose_clustering.hpp"
#include "../v2/kmeans_cluster.hpp"
#include "snn_graph_cluster.hpp"

#include "../v2/tsne.hpp"
#include "../v2/umap.hpp"

#include "marker_detection.hpp"
#include "custom_selections.hpp"
#include "../v2/cell_labelling.hpp"

#include "_metadata.hpp"
#include <algorithm>
#include <unordered_map>

namespace kanaval {

namespace v3 {

void validate(const H5::H5File& handle, bool embedded, int version) {
    auto i_out = validate_inputs(handle, embedded, version);

    auto rnaIt = i_out.num_features.find("RNA");
    bool rna_available = rnaIt != i_out.num_features.end();
    bool adt_available = i_out.num_features.find("ADT") != i_out.num_features.end();
    bool crispr_available = i_out.num_features.find("CRISPR") != i_out.num_features.end();

    auto add_modalities = [](auto& host, auto rna_val, auto adt_val, auto crispr_val) -> void {
        if (rna_val >= 0) { host["RNA"] = rna_val; }
        if (adt_val >= 0) { host["ADT"] = adt_val; }
        if (crispr_val >= 0) { host["CRISPR"] = crispr_val; }
    };

    // Quality control.
    int filtered_cells;
    {
        std::unordered_map<std::string, int> survivors;
        add_modalities(survivors, 
            validate_rna_quality_control(handle, i_out.num_cells, i_out.num_blocks, rna_available, version),
            validate_adt_quality_control(handle, i_out.num_cells, i_out.num_blocks, adt_available, version),
            validate_crispr_quality_control(handle, i_out.num_cells, i_out.num_blocks, crispr_available, version)
        );
        filtered_cells = validate_cell_filtering(handle, i_out.num_cells, survivors, version);
    }

    // Normalization.
    validate_rna_normalization(handle);
    v2::validate_adt_normalization(handle, filtered_cells, adt_available, version);
    validate_crispr_normalization(handle);

    // Feature selection.
    validate_feature_selection(handle, (rna_available ? rnaIt->second : -1), rna_available, version);

    // Dimensionality reduction.
    {
        std::unordered_map<std::string, int> num_pcs;
        add_modalities(num_pcs,
            validate_rna_pca(handle, filtered_cells, rna_available, version),
            v2::validate_adt_pca(handle, filtered_cells, adt_available, version),
            validate_crispr_pca(handle, filtered_cells, crispr_available, version)
        );
        int total_pcs = validate_combine_embeddings(handle, filtered_cells, num_pcs, version);
        v2::validate_batch_correction(handle, total_pcs, filtered_cells, i_out.num_blocks, version);
    }

    v2::validate_neighbor_index(handle);

    // Clustering.
    auto cluster_method = v2::validate_choose_clustering(handle);
    int nclusters = 0;
    {
        bool is_snn = (cluster_method == "snn_graph");
        int snn_found = validate_snn_graph_cluster(handle, filtered_cells, is_snn);

        bool is_kmeans = (cluster_method == "kmeans");
        int kmeans_found = v2::validate_kmeans_cluster(handle, filtered_cells, is_kmeans);

        if (is_snn) {
            nclusters = snn_found;
        } else if (is_kmeans) {
            nclusters = kmeans_found;
        }
    }

    v2::validate_tsne(handle, filtered_cells);
    v2::validate_umap(handle, filtered_cells);

    validate_marker_detection(handle, nclusters, i_out.num_features, version);
    validate_custom_selections(handle, filtered_cells, i_out.num_features, version);
    v2::validate_cell_labelling(handle, nclusters);

    // Checking metadata.
    validate_metadata(handle, version);
}

}

}

#endif
