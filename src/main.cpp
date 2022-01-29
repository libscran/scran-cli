#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#include <map>
#include <memory>
#include <thread>
#include <iostream>
#include <chrono>
#include <fstream>

int nthreads_global = 1;

template<class Function>
void run_parallel(int total, Function fun) {
    int nworkers = nthreads_global;
    int jobs_per_worker = std::ceil(static_cast<double>(total)/nworkers);
    std::vector<std::thread> workers;
    workers.reserve(nworkers);
    int first = 0;

    for (int w = 0; w < nworkers && first < total; ++w, first += jobs_per_worker) {
        int last = std::min(first + jobs_per_worker, total);
        workers.emplace_back(fun, first, last);
    }

    for (auto& wrk : workers) {
        wrk.join();
    }
}

#define TATAMI_CUSTOM_PARALLEL run_parallel
#define SCRAN_CUSTOM_PARALLEL run_parallel
#define KNNCOLLE_CUSTOM_PARALLEL run_parallel

#include "scran/scran.hpp"
#include "tatami/ext/MatrixMarket_layered.hpp"

#include "qdtsne/qdtsne.hpp"
#include "umappp/Umap.hpp"

#ifdef __has_include
  #if __has_include(<filesystem>)
    #include <filesystem>
    namespace fs = std::filesystem;
  #else
    #include <experimental/filesystem>
    namespace fs = std::experimental::filesystem;
  #endif
#else
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#endif

int main(int argc, char* argv []) {
    /******************************
     *** Parsing the arguments. ***
     ******************************/

    CLI::App app{"Single-cell RNA-seq analyses on the command-line"};
    
    std::string path;
    app.add_option("path", path, "Path to the Matrix Market file")->required();

    double nthreads;
    app.add_option("-t,--nthreads", nthreads, "Number of threads to use (+2 for UMAP and t-SNE, which use their own threads)")
        ->default_val(1);
    
    std::string output;
    app.add_option("-o,--output", output, "Path to the output directory")
        ->default_val("output");

    bool skip_output;
    app.add_flag("--skip-output", skip_output, "Run the analysis but do not save results")
        ->default_val(false);

    double nmads;
    app.add_option("--qc-nmads", nmads, "Number of MADs to use for filtering")
        ->default_val(3);

    double span;
    app.add_option("--hvg-span", span, "LOWESS span for variance modelling")
        ->default_val(0.4);

    int nhvgs;
    app.add_option("--hvg-num", nhvgs, "Number of HVGs to use for PCA")
        ->default_val(2500);

    int npcs;
    app.add_option("--pca-num", npcs, "Number of PCs to keep")
        ->default_val(25);

    bool approx;
    app.add_option("--nn-approx", approx, "Whether to use an approximate neighbor search")
        ->default_val(true);

    int snn_neighbors;
    app.add_option("--snn-neighbors", snn_neighbors, "Number of neighbors to use for the SNN graph")
        ->default_val(10);

    scran::BuildSNNGraph::Scheme snn_scheme;
    std::map<std::string, scran::BuildSNNGraph::Scheme> scheme_map{
        {"ranked", scran::BuildSNNGraph::RANKED}, 
        {"number", scran::BuildSNNGraph::NUMBER}, 
        {"jaccard", scran::BuildSNNGraph::JACCARD}
    };
    app.add_option("--snn-scheme", snn_scheme, "Edge weighting scheme: ranked, number or jaccard")
        ->transform(CLI::CheckedTransformer(scheme_map, CLI::ignore_case))
        ->default_val(scran::BuildSNNGraph::RANKED);

    double snn_res;
    app.add_option("--snn-res", snn_res, "Resolution to use in multi-level community detection")
        ->default_val(1);

    double tsne_perplexity;
    app.add_option("--tsne-perplexity", tsne_perplexity, "Perplexity to use in t-SNE")
        ->default_val(30);

    int tsne_iterations;
    app.add_option("--tsne-iter", tsne_iterations, "Number of iterations to use in t-SNE")
        ->default_val(500);

    int umap_neighbors;
    app.add_option("--umap-neighbors", umap_neighbors, "Number of neighbors to use in the UMAP")
        ->default_val(15);

    double umap_mindist;
    app.add_option("--umap-mindist", umap_mindist, "Minimum distance to use in the UMAP")
        ->default_val(0.01);
    
    int umap_epochs;
    app.add_option("--umap-epochs", umap_epochs, "Number of epochs to use in the UMAP")
        ->default_val(500);

    CLI11_PARSE(app, argc, argv);

    /********************************
     *** Performing the analysis. ***
     ********************************/

    // Setting up some bits and pieces.
    nthreads_global = nthreads;

    auto declare = [&](const auto& x) -> void {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - x).count())/1000 << "s" << std::endl;
    };

    // Loading the data from an MatrixMarket file.
    std::cout << "Initializing matrix... " << std::flush;
    auto start = std::chrono::high_resolution_clock::now();
    size_t n = path.size();
    tatami::LayeredMatrixData<double, int> mat;
    if (n > 3 && path[n-3] == '.' && path[n-2] == 'g' && path[n-1] == 'z') {
        mat = tatami::MatrixMarket::load_layered_sparse_matrix_gzip(path.c_str());
    } else {
        mat = tatami::MatrixMarket::load_layered_sparse_matrix(path.c_str());
    }
    declare(start);

    // Filtering out low-quality cells. 
    std::cout << "Computing QC metrics... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto qc_res = scran::PerCellQCMetrics().run(mat.matrix.get(), { /* mito subset definitions go here */ });
    declare(start);

    std::cout << "Computing QC thresholds... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto qc_filters = scran::PerCellQCFilters().set_nmads(nmads).run(qc_res);
    declare(start);

    std::cout << "Filtering cells... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto filtered = scran::FilterCells().run(mat.matrix, qc_filters.overall_filter.data());
    declare(start);

    // Computing log-normalized expression values, re-using the total count from the QC step.
    std::cout << "Log-normalizing the counts... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto size_factors = scran::subset_vector<false>(qc_res.sums, qc_filters.overall_filter.data());
    auto normalized = scran::LogNormCounts().run(filtered, std::move(size_factors));
    declare(start);

    // Identifying highly variable genes.
    std::cout << "Mean-variance modelling... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto var_res = scran::ModelGeneVar().set_span(span).run(normalized.get());
    auto keep = scran::ChooseHVGs().set_top(nhvgs).run(var_res.residuals[0].size(), var_res.residuals[0].data());
    declare(start);

    // Performing a PCA on the HVGs. We transpose the output so cells are columns again.
    std::cout << "Principal components analysis... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto pca_res = scran::RunPCA().set_rank(npcs).run(normalized.get(), keep.data());
    pca_res.pcs.adjointInPlace();
    declare(start);

    // Building the nearest neighbor index.
    std::cout << "Building the neighbor index... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    std::unique_ptr<knncolle::Base<> > ptr;
    if (approx) {
        ptr.reset(new knncolle::AnnoyEuclidean<>(npcs, pca_res.pcs.cols(), pca_res.pcs.data()));
    } else {
        ptr.reset(new knncolle::VpTreeEuclidean<>(npcs, pca_res.pcs.cols(), pca_res.pcs.data()));
    }
    declare(start);

    // Finding all the neighbors.
    std::cout << "Finding neighbors for clustering... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto snn_nns = knncolle::find_nearest_neighbors_index_only<int>(ptr.get(), snn_neighbors);
    declare(start);

    std::cout << "Finding neighbors for t-SNE... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto tsne_nns = knncolle::find_nearest_neighbors<int, double>(ptr.get(), qdtsne::perplexity_to_k(tsne_perplexity));
    declare(start);

    std::cout << "Finding neighbors for UMAP... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto umap_nns = knncolle::find_nearest_neighbors<int, float>(ptr.get(), umap_neighbors);
    declare(start);

    // Running all remaining steps in parallel threads.
    std::vector<double> tsne_output = qdtsne::initialize_random<>(ptr->nobs());
    std::vector<float> umap_output(ptr->nobs() * 2);
    std::vector<int> best_clustering;
    scran::ScoreMarkers::Results<double> marker_res;

    std::mutex cout_lock; 

    std::thread first([&]() -> void {
        auto start = std::chrono::high_resolution_clock::now();
        auto snn_graph = scran::BuildSNNGraph().set_weighting_scheme(snn_scheme).run(snn_nns);
        cout_lock.lock();
        std::cout << "SNN graph construction... " << std::flush;
        declare(start);
        cout_lock.unlock();

        start = std::chrono::high_resolution_clock::now();
        auto clust_res = scran::ClusterSNNGraphMultiLevel().set_resolution(snn_res).run(snn_nns.size(), snn_graph);
        cout_lock.lock();
        std::cout << "Multi-level clustering... " << std::flush;
        declare(start);
        cout_lock.unlock();

        start = std::chrono::high_resolution_clock::now();
        best_clustering = clust_res.membership[clust_res.max];
        marker_res = scran::ScoreMarkers().run(normalized.get(), best_clustering.data());
        cout_lock.lock();
        std::cout << "Marker detection... " << std::flush;
        declare(start);
        cout_lock.unlock();
    });

    std::thread second([&]() -> void {
        auto start = std::chrono::high_resolution_clock::now();
        qdtsne::Tsne<>().set_perplexity(tsne_perplexity).set_max_iter(tsne_iterations).run(tsne_nns, tsne_output.data());
        cout_lock.lock();
        std::cout << "t-SNE calculation... " << std::flush;
        declare(start);
        cout_lock.unlock();
    });

    std::thread third([&]() -> void {
        auto start = std::chrono::high_resolution_clock::now();
        umappp::Umap<float>().set_min_dist(umap_mindist).set_num_epochs(umap_epochs).run(umap_nns, 2, umap_output.data());
        cout_lock.lock();
        std::cout << "UMAP calculation... " << std::flush;
        declare(start);
        cout_lock.unlock();
    });

    first.join();
    second.join();
    third.join();

    /***************************
     *** Saving the results. ***
     ***************************/

    if (skip_output) {
        return 0;
    }

    fs::create_directories(output);
    
    // QC results.
    {
        std::ofstream handle(output + "/qc_metrics.tsv");
        handle << "sums\tdetected\n";
        auto def = std::cout.precision();
        for (size_t i = 0; i < qc_res.sums.size(); ++i) {
            handle << std::setprecision(def) << qc_res.sums[i] << "\t" << qc_res.detected[i] << "\n"; // std::setprecision(4) << "\n";
        }
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/qc_thresholds.tsv");
        handle << std::setprecision(6);
        handle << "sums\t" << qc_filters.thresholds.sums[0] << "\n";
        handle << "detected\t" << qc_filters.thresholds.detected[0] << "\n";
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/qc_discard.tsv");
        handle << std::setprecision(6);
        handle << "sums\tdetected\toverall\n";
        for (size_t i = 0; i < qc_filters.filter_by_sums.size(); ++i) {
            handle << qc_filters.filter_by_sums[i] << "\t" 
                << qc_filters.filter_by_detected[i] << "\t" 
                << qc_filters.overall_filter[i] << "\n";
        }
        handle << std::flush;
    }

    // Mean-variance results.
    auto perm = mat.permutation;
    {
        std::ofstream handle(output + "/variances.tsv");
        handle << std::setprecision(6);
        handle << "mean\tvariance\tfitted\tresid\n";
        for (auto i : perm) {
            handle << var_res.means[0][i] << "\t" 
                << var_res.variances[0][i] << "\t"
                << var_res.fitted[0][i] << "\t"
                << var_res.residuals[0][i] << "\n";
        }
        handle << std::flush;
    }

    // PCA-related results.
    {
        std::ofstream handle(output + "/pca.tsv");
        handle << std::setprecision(6);
        for (size_t i = 0; i < pca_res.pcs.cols(); ++i) {
            handle << pca_res.pcs(0, i); 
            for (size_t j = 1; j < pca_res.pcs.rows(); ++j) {
                handle << "\t" << pca_res.pcs(j, i);
            }
            handle << "\n";
        }
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/pca_varexp.tsv");
        handle << std::setprecision(6);
        for (auto v : pca_res.variance_explained) {
            handle << v / pca_res.total_variance << "\n";
        }
        handle << std::flush;
    }

    // Clustering results.
    {
        std::ofstream handle(output + "/snn_cluster_multilevel.tsv");
        for (auto i : best_clustering) {
            handle << i << "\n";
        }
        handle << std::flush;
    }

    // Marker results.
    {
        fs::create_directory(output + "/markers");
        for (size_t i = 0; i < marker_res.means.size(); ++i) {
            std::ofstream handle(output + "/markers/" + std::to_string(i) + ".tsv");
            handle << "mean\tdetected\tlfc\tdelta_detected\tcohen\tauc\n";
            handle << std::setprecision(6);
            for (auto j : perm) {
                handle << marker_res.means[i][0][j] << "\t" 
                    << marker_res.detected[i][0][j] << "\t"
                    << marker_res.lfc[scran::differential_analysis::MEAN][i][j] << "\t"
                    << marker_res.delta_detected[scran::differential_analysis::MEAN][i][j] << "\t"
                    << marker_res.cohen[scran::differential_analysis::MEAN][i][j] << "\t"
                    << marker_res.auc[scran::differential_analysis::MEAN][i][j] << "\n";
            }
            handle << std::flush;
        }
    }

    // Dimensionality reduction results.
    {
        std::ofstream handle(output + "/tsne.tsv");
        handle << std::setprecision(6);
        for (size_t i = 0; i < ptr->nobs(); ++i) {
            handle << tsne_output[2*i] << "\t" << tsne_output[2*i + 1] << "\n";
        }
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/umap.tsv");
        handle << std::setprecision(6);
        for (size_t i = 0; i < ptr->nobs(); ++i) {
            handle << umap_output[2*i] << "\t" << umap_output[2*i + 1] << "\n";
        }
        handle << std::flush;
    }

    return 0;
}
