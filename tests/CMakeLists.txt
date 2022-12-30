include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    libtest 

    src/v2/inputs.cpp
    src/v2/quality_control.cpp
    src/v2/adt_quality_control.cpp
    src/v2/cell_filtering.cpp
    src/v2/normalization.cpp
    src/v2/adt_normalization.cpp
    src/v2/feature_selection.cpp
    src/v2/pca.cpp
    src/v2/adt_pca.cpp
    src/v2/combine_embeddings.cpp
    src/v2/batch_correction.cpp
    src/v2/neighbor_index.cpp
    src/v2/tsne.cpp
    src/v2/umap.cpp
    src/v2/kmeans_cluster.cpp
    src/v2/snn_graph_cluster.cpp
    src/v2/choose_clustering.cpp
    src/v2/marker_detection.cpp
    src/v2/custom_selections.cpp
    src/v2/cell_labelling.cpp
    src/v2/validate.cpp

    src/v3/rna_quality_control.cpp
    src/v3/adt_quality_control.cpp
    src/v3/cell_filtering.cpp
    src/v3/rna_normalization.cpp
    src/v3/rna_pca.cpp
    src/v3/combine_embeddings.cpp
    src/v3/marker_detection.cpp
    src/v3/custom_selections.cpp
    src/v3/snn_graph_cluster.cpp
)

target_link_libraries(
    libtest
    gtest_main
    kanaval
)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
target_link_libraries(libtest hdf5::hdf5 hdf5::hdf5_cpp)

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

include(GoogleTest)
gtest_discover_tests(libtest)
