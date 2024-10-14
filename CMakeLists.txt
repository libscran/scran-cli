cmake_minimum_required(VERSION 3.14)

project(scrancli
    VERSION 2.0.0
    DESCRIPTION "A command-line interface for scran"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
add_subdirectory(extern)

# Building the application.
add_executable(scrancli src/main.cpp)

target_link_libraries(
    scrancli

    tatami
    tatami_mtx

    scran_qc
    scran_norm
    scran_variances
    scran_pca
    scran_markers

    kmeans

    knncolle
    knncolle_annoy

    qdtsne
    umappp

    igraph::igraph
    scran_graph_cluster

    CLI11::CLI11
)

set_target_properties(scrancli PROPERTIES OUTPUT_NAME "scran")

include(CheckIncludeFiles)
check_include_files(filesystem HAVE_CXX_FS)
if (NOT HAVE_CXX_FS) 
    target_link_libraries(scrancli stdc++fs)
endif()

find_package(ZLIB REQUIRED)
target_link_libraries(scrancli ZLIB::ZLIB)

find_package(igraph REQUIRED)
target_link_libraries(scrancli igraph::igraph) 
