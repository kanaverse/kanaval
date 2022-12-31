#ifndef HELPERS_V3_H
#define HELPERS_V3_H

#include "H5Cpp.h"
#include <string>

namespace v3 {

void add_adt_quality_control(H5::H5File& handle, int num_cells, int num_samples, int lost = 10);

void add_cell_filtering(H5::H5File& handle, int num_cells, int lost = 10);

void add_combine_embeddings(H5::H5File& handle, int num_cells, int total_pcs);

void add_crispr_normalization(H5::H5File& handle);

void add_crispr_pca(H5::H5File& handle, int num_pcs, int num_cells);

void add_crispr_quality_control(H5::H5File& handle, int num_cells, int num_blocks, int lost = 10);

void add_custom_selections(H5::H5File& handle, const std::unordered_map<std::string, int>& modalities, int ncells = 10, bool has_auc = true, double lfc_threshold = 0);

void add_single_matrix(H5::H5File& handle, std::string mode = "MatrixMarket", int ngenes = 1000, int ncells = 100, int nblocks = 1);

int add_multiple_matrices(H5::H5File& handle, int ngenes = 500, int ncells = 100);

void add_marker_detection(H5::H5File& handle, const std::unordered_map<std::string, int>& modalities, int nclusters, bool has_auc = true, double lfc_threshold = 0);

void add_rna_normalization(H5::H5File& handle);

void add_rna_pca(H5::H5File& handle, int num_pcs, int num_cells);

void add_rna_quality_control(H5::H5File& handle, int num_cells, int num_blocks, int lost = 10);

void add_snn_graph_cluster(H5::H5File& handle, int num_cells, int num_clusters = 10);

H5::Group add__metadata(H5::H5File& handle);

}

#endif
