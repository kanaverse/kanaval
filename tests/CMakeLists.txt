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
    src/inputs.cpp

    src/quality_control.cpp
    src/adt_quality_control.cpp

    src/normalization.cpp
    src/adt_normalization.cpp

    src/feature_selection.cpp

    src/pca.cpp
    src/adt_pca.cpp

    src/neighbor_index.cpp
    src/tsne.cpp
    src/umap.cpp
    src/kmeans_cluster.cpp
    src/snn_graph_cluster.cpp
    src/choose_clustering.cpp
    src/marker_detection.cpp
    src/custom_selections.cpp
    src/cell_labelling.cpp
    src/validate.cpp
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
