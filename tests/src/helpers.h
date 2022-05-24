#ifndef HELPERS_H
#define HELPERS_H

#include "H5Cpp.h"
#include <string>

void add_single_matrix(H5::H5File&, std::string = "MatrixMarket", int = 1000, int = 100);

int add_multiple_matrices(H5::H5File& handle, int = 1000, int = 100);

void add_cell_labelling(H5::H5File&, int = 3);

void add_choose_clustering(H5::H5File&);

void add_custom_selections(H5::H5File&, int, int = 10);
void add_custom_selections(H5::H5File&, const std::vector<std::string>&, const std::vector<int>&, int = 10);

void add_feature_selection(H5::H5File&, int);

void add_kmeans_cluster(H5::H5File&, int, int = 10);

void add_marker_detection(H5::H5File&, int, int);
void add_marker_detection(H5::H5File&, const std::vector<int>&, int, const std::vector<std::string>&);

void add_neighbor_index(H5::H5File&);

void add_normalization(H5::H5File&);

void add_adt_normalization(H5::H5File&, int);

void add_pca(H5::H5File&, int, int);

void add_adt_pca(H5::H5File&, int, int);

void add_combine_embeddings(H5::H5File&, int, int);

void add_batch_correction(H5::H5File&, int, int);

void add_quality_control(H5::H5File&, int, int, int = 10);

void add_adt_quality_control(H5::H5File&, int, int, int = 10);

void add_cell_filtering(H5::H5File&, int, int = 10);

void add_snn_graph_cluster(H5::H5File&, int, int = 10);

void add_tsne(H5::H5File&, int);

void add_umap(H5::H5File&, int);

#endif
