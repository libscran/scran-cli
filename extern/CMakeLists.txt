cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG v2.4.2
)

FetchContent_Declare(
  aarand
  GIT_REPOSITORY https://github.com/LTLA/aarand
  GIT_TAG v1.0.2
)

FetchContent_Declare(
  subpar
  GIT_REPOSITORY https://github.com/LTLA/subpar
  GIT_TAG v0.3.1
)

FetchContent_Declare(
  tatami 
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG v3.1.1
)

FetchContent_Declare(
  tatami_stats
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_stats
  GIT_TAG v1.1.2
)

FetchContent_Declare(
  eminem
  GIT_REPOSITORY https://github.com/tatami-inc/eminem
  GIT_TAG v1.0.0
)

FetchContent_Declare(
  tatami_mtx
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_mtx
  GIT_TAG v2.0.0
)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen
  GIT_TAG 3.4.0
)

FetchContent_Declare(
  irlba
  GIT_REPOSITORY https://github.com/LTLA/CppIrlba
  GIT_TAG v2.0.0
)

FetchContent_Declare(
  kmeans 
  GIT_REPOSITORY https://github.com/LTLA/CppKmeans
  GIT_TAG v3.1.0
)

FetchContent_Declare(
  annoy
  GIT_REPOSITORY https://github.com/spotify/annoy
  GIT_TAG v1.17.2
)

FetchContent_Declare(
  knncolle
  GIT_REPOSITORY https://github.com/knncolle/knncolle
  GIT_TAG v2.3.2
)

FetchContent_Declare(
  knncolle_annoy
  GIT_REPOSITORY https://github.com/knncolle/knncolle_annoy
  GIT_TAG v0.1.1
)

FetchContent_Declare(
  scran_qc
  GIT_REPOSITORY https://github.com/libscran/scran_qc
  GIT_TAG master
)

FetchContent_Declare(
  scran_norm
  GIT_REPOSITORY https://github.com/libscran/scran_norm
  GIT_TAG v0.1.0
)

FetchContent_Declare(
  scran_blocks
  GIT_REPOSITORY https://github.com/libscran/scran_blocks
  GIT_TAG v0.1.0
)

FetchContent_Declare(
  WeightedLowess 
  GIT_REPOSITORY https://github.com/LTLA/CppWeightedLowess
  GIT_TAG v2.1.2
)

FetchContent_Declare(
  scran_variances
  GIT_REPOSITORY https://github.com/libscran/scran_variances
  GIT_TAG v0.1.3
)

FetchContent_Declare(
  scran_pca
  GIT_REPOSITORY https://github.com/libscran/scran_pca
  GIT_TAG v0.1.0
)

FetchContent_Declare(
  raiigraph
  GIT_REPOSITORY https://github.com/LTLA/raiigraph
  GIT_TAG v1.1.0
)

FetchContent_Declare(
  scran_graph_cluster
  GIT_REPOSITORY https://github.com/libscran/scran_graph_cluster
  GIT_TAG v0.1.2
)

FetchContent_Declare(
  scran_markers
  GIT_REPOSITORY https://github.com/libscran/scran_markers
  GIT_TAG v0.1.1
)

FetchContent_Declare(
  qdtsne 
  GIT_REPOSITORY https://github.com/libscran/qdtsne
  GIT_TAG v2.0.1
)

FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/libscran/umappp
  GIT_TAG v2.0.1
)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG v1.2.2
)

FetchContent_MakeAvailable(cli11)

FetchContent_MakeAvailable(aarand)
FetchContent_MakeAvailable(subpar)

FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(tatami_stats)

FetchContent_MakeAvailable(eminem)
FetchContent_MakeAvailable(tatami_mtx)

FetchContent_MakeAvailable(kmeans)

FetchContent_MakeAvailable(knncolle)
FetchContent_MakeAvailable(annoy)
FetchContent_MakeAvailable(knncolle_annoy)

FetchContent_MakeAvailable(scran_blocks)
FetchContent_MakeAvailable(scran_qc)
FetchContent_MakeAvailable(scran_norm)

FetchContent_MakeAvailable(WeightedLowess)
FetchContent_MakeAvailable(scran_variances)

FetchContent_MakeAvailable(Eigen)
FetchContent_MakeAvailable(irlba)
FetchContent_MakeAvailable(scran_pca)

FetchContent_MakeAvailable(raiigraph)
FetchContent_MakeAvailable(scran_graph_cluster)

FetchContent_MakeAvailable(scran_markers)

FetchContent_MakeAvailable(umappp)

FetchContent_MakeAvailable(qdtsne)

FetchContent_MakeAvailable(byteme)
