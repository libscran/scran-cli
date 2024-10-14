#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#include <map>
#include <memory>
#include <thread>
#include <iostream>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <algorithm>

#include "tatami/tatami.hpp"
#include "tatami_mtx/tatami_mtx.hpp"

#include "scran_qc/scran_qc.hpp"
#include "scran_norm/scran_norm.hpp"
#include "scran_variances/scran_variances.hpp"
#include "scran_pca/scran_pca.hpp"
#include "scran_graph_cluster/scran_graph_cluster.hpp"
#include "scran_markers/scran_markers.hpp"

#include "knncolle/knncolle.hpp"
#include "knncolle_annoy/knncolle_annoy.hpp"

#include "qdtsne/qdtsne.hpp"
#include "umappp/umappp.hpp"

std::vector<int> split_mito_list(const std::string& mito_list) {
    std::vector<int> output;
    int current = 0;
    bool empty = true;

    for (auto x : mito_list) {
        switch (x) {
            case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9': 
                current *= 10;
                current += (x - '0');
                empty = false;
                break;
            case ',': 
                if (empty) {
                    throw std::runtime_error("empty row index in the mitochondrial list");
                }
                output.push_back(current);
                current = 0;
                empty = true;
                break;
            default:
                throw std::runtime_error(std::string("unknown character '") + x + "' in the mitochondrial list");
        }
    }

    if (empty) {
        throw std::runtime_error("empty row index in the mitochondrial list");
    } else {
        output.push_back(current);
    }
    return output;
}

int main(int argc, char* argv []) {
    /******************************
     *** Parsing the arguments. ***
     ******************************/

    struct {
        std::string file_path;
        std::string output;

        std::string mito_list;
        double num_mads;

        double fit_span;
        int num_hvgs;
        int num_pcs;

        int nn_approx;
        std::string snn_scheme;
        int snn_neighbors;
        double snn_resolution;

        double tsne_perplexity;
        double tsne_iterations;

        int umap_neighbors;
        double umap_mindist;
        int umap_epochs;

        int num_threads;
    } all_opt;

    {
        CLI::App app{"Single-cell RNA-seq analyses on the command-line"};
        
        app.add_option("counts", all_opt.file_path, "Path to the MatrixMarket file containing the counts.")->required();
        app.add_option("-o,--output", all_opt.output, "Path to the output directory. If empty, results are not saved.")->default_val("output");

        app.add_option("--mito-list", all_opt.mito_list, "Comma-separated list of the 0-based row indices of the mitochondrial genes")->default_val("");
        app.add_option("--num-mads", all_opt.num_mads, "Number of MADs to use for defining QC thresholds.")->default_val(3);

        app.add_option("--fit-span", all_opt.fit_span, "LOWESS span for fitting the mean-variance trend.")->default_val(0.3);
        app.add_option("--num-hvgs", all_opt.num_hvgs, "Number of HVGs to use for PCA.")->default_val(2500);
        app.add_option("--num-pcs", all_opt.num_pcs, "Number of PCs to keep.")->default_val(25);

        app.add_option("--nn-approx", all_opt.nn_approx, "Whether to use an approximate neighbor search.")->default_val(true);

        app.add_option("--snn-neighbors", all_opt.snn_neighbors, "Number of neighbors to use for the SNN graph.")->default_val(10);
        app.add_option("--snn-scheme", all_opt.snn_scheme, "Edge weighting scheme for SNN graph construction.")
            ->check(CLI::IsMember({ "ranked", "number", "jaccard" }))
            ->default_val("ranked");
        app.add_option("--snn-res", all_opt.snn_resolution, "Resolution to use in multi-level community detection.")->default_val(0.5);

        app.add_option("--tsne-perplexity", all_opt.tsne_perplexity, "Perplexity to use in t-SNE.")->default_val(30);
        app.add_option("--tsne-iter", all_opt.tsne_iterations, "Number of iterations to use in t-SNE.")->default_val(500);
        app.add_option("--umap-neighbors", all_opt.umap_neighbors, "Number of neighbors to use in the UMAP.")->default_val(15);
        app.add_option("--umap-mindist", all_opt.umap_mindist, "Minimum distance to use in the UMAP.")->default_val(0.1);
        app.add_option("--umap-epochs", all_opt.umap_epochs, "Number of epochs to use in the UMAP.")->default_val(500);

        app.add_option("-t,--nthreads", all_opt.num_threads, "Number of threads to use (+2 for UMAP and t-SNE, which use their own threads).")->default_val(1);

        CLI11_PARSE(app, argc, argv);
    }

    // Setting up an easy timekeeping function.
    auto declare = [&](const auto& x) -> void {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - x).count())/1000 << "s" << std::endl;
    };

    /********************************
     *** Performing the analysis. ***
     ********************************/

    std::cout << "Reading matrix from file... " << std::flush;
    auto start = std::chrono::high_resolution_clock::now();
    auto mat = tatami_mtx::load_matrix_from_some_file<double, int>(all_opt.file_path.c_str(), [&]{
        tatami_mtx::Options opt;
        opt.parallel = (all_opt.num_threads > 1);
        return opt;
    }());
    declare(start);

    std::vector<const uint8_t*> mito_ptrs; // Deriving mitochondrial subsets.
    std::vector<uint8_t> mito;
    bool has_mito = (all_opt.mito_list != "");
    if (has_mito) {
        auto mito_list = split_mito_list(all_opt.mito_list);
        mito.resize(mat->nrow());
        for (auto x : mito_list) {
            if (x >= mat->nrow()) {
                throw std::runtime_error("mitochondrial row index " + std::to_string(x) + " is out of range");
            }
            mito[x] = 1;
        }
        mito_ptrs.push_back(mito.data());
    }

    std::cout << "Computing QC metrics... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto qc_res = scran_qc::compute_rna_qc_metrics(*mat, mito_ptrs, [&]{
        scran_qc::ComputeRnaQcMetricsOptions opt;
        opt.num_threads = all_opt.num_threads;
        return opt;
    }());
    declare(start);

    std::cout << "Computing QC thresholds... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto qc_filters = scran_qc::compute_rna_qc_filters(qc_res, [&]{
        scran_qc::ComputeRnaQcFiltersOptions opt;
        opt.sum_num_mads = all_opt.num_mads;
        opt.detected_num_mads = all_opt.num_mads;
        opt.subset_proportion_num_mads = all_opt.num_mads;
        return opt;
    }());
    declare(start);

    std::cout << "Filtering out low-quality cells... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto keep_cells = qc_filters.filter(qc_res);
    auto keep_cells_index = scran_qc::filter_index<int>(keep_cells.size(), keep_cells.data());
    auto filtered = tatami::make_DelayedSubset(mat, keep_cells_index, /* row = */ false);
    declare(start);

    std::cout << "Log-normalizing the counts... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    std::vector<double> size_factors; // Reusing the total counts for the high-quality cells as their size factors.
    size_factors.reserve(keep_cells_index.size());
    for (auto x : keep_cells_index) {
        size_factors.push_back(qc_res.sum[x]);
    }
    scran_norm::center_size_factors(size_factors.size(), size_factors.data(), NULL, scran_norm::CenterSizeFactorsOptions());
    auto normalized = scran_norm::normalize_counts(filtered, std::move(size_factors), scran_norm::NormalizeCountsOptions());
    declare(start);

    std::cout << "Mean-variance modelling... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto var_res = scran_variances::model_gene_variances(*normalized, [&]{
        scran_variances::ModelGeneVariancesOptions opt;
        opt.fit_variance_trend_options.span = all_opt.fit_span;
        opt.num_threads = all_opt.num_threads;
        return opt;
    }());
    auto hvgs = scran_variances::choose_highly_variable_genes_index(var_res.residuals.size(), var_res.residuals.data(), [&]{
        scran_variances::ChooseHighlyVariableGenesOptions opt;
        opt.top = all_opt.num_hvgs;
        return opt;
    }());
    declare(start);

    std::cout << "Principal components analysis... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto hvg_mat = tatami::make_DelayedSubset(normalized, std::move(hvgs), /* row = */ true);
    auto pca_res = scran_pca::simple_pca(*hvg_mat, [&]{
        scran_pca::SimplePcaOptions opt;
        opt.number = all_opt.num_pcs;
        opt.num_threads = all_opt.num_threads;
        return opt;
    }());
    declare(start);

    // Building the nearest neighbor index.
    std::cout << "Building the neighbor index... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    knncolle::SimpleMatrix<int, int, double> matview(pca_res.components.rows(), pca_res.components.cols(), pca_res.components.data());
    std::unique_ptr<knncolle::Builder<decltype(matview), double> > builder;
    if (all_opt.nn_approx) {
        builder.reset(new knncolle_annoy::AnnoyBuilder<Annoy::Euclidean, decltype(matview), double>);
    } else {
        builder.reset(new knncolle::VptreeBuilder<knncolle::EuclideanDistance, decltype(matview), double>);
    }
    auto prebuilt_index = builder->build_unique(matview);
    declare(start);

    // Finding all the neighbors. 
    std::cout << "Finding nearest neighbors... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    int tsne_neighbors = qdtsne::perplexity_to_k(all_opt.tsne_perplexity);
    auto most_nns = std::max({ all_opt.snn_neighbors, all_opt.umap_neighbors, tsne_neighbors });
    auto common_nn_res = knncolle::find_nearest_neighbors(*prebuilt_index, most_nns, all_opt.num_threads);
    size_t nobs = common_nn_res.size();
    declare(start);

    // Running the visualization steps in parallel threads.
    std::mutex cout_lock; 

    std::vector<double> tsne_output = qdtsne::initialize_random<2>(prebuilt_index->num_observations());
    std::thread second([&]() -> void {
        auto start = std::chrono::high_resolution_clock::now();

        knncolle::NeighborList<int, double> tsne_nns;
        tsne_nns.reserve(common_nn_res.size());
        for (const auto& nn : common_nn_res) {
            tsne_nns.emplace_back(nn.begin(), nn.begin() + tsne_neighbors);
        }

        auto status = qdtsne::initialize<2>(std::move(tsne_nns), [&]{
            qdtsne::Options opt;
            opt.perplexity = all_opt.tsne_perplexity;
            opt.num_threads = all_opt.num_threads;
            return opt;
        }());
        status.run(tsne_output.data());

        cout_lock.lock();
        std::cout << "t-SNE calculation... " << std::flush;
        declare(start);
        cout_lock.unlock();
    });

    std::vector<float> umap_output(nobs * 2);
    std::thread third([&]() -> void {
        auto start = std::chrono::high_resolution_clock::now();

        knncolle::NeighborList<int, float> umap_nns(nobs); // single-precision speeds up the UMAP.
        for (size_t i = 0; i < nobs; ++i) {
            const auto& nn = common_nn_res[i];
            auto& output = umap_nns[i];
            for (int j = 0; j < all_opt.umap_neighbors; ++j) {
                output.emplace_back(nn[j].first, nn[j].second);
            }
        }

        auto status = umappp::initialize(std::move(umap_nns), 2, umap_output.data(), [&]{
            umappp::Options opt;
            opt.min_dist = all_opt.umap_mindist;
            opt.num_epochs = all_opt.umap_epochs;
            opt.num_threads = all_opt.num_threads;
            return opt;
        }());
        status.run();

        cout_lock.lock();
        std::cout << "UMAP calculation... " << std::flush;
        declare(start);
        cout_lock.unlock();
    });

    // Running everything else on the main thread.
    cout_lock.lock();
    std::cout << "SNN graph construction... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int> > snn_nns(nobs);
    for (size_t i = 0; i < nobs; ++i) {
        const auto& nn = common_nn_res[i];
        auto& output = snn_nns[i];
        for (int j = 0; j < all_opt.snn_neighbors; ++j) {
            output.emplace_back(nn[j].first);
        }
    }
    auto snn_graph = scran_graph_cluster::build_snn_graph(std::move(snn_nns), [&]{
        scran_graph_cluster::BuildSnnGraphOptions opt;
        if (all_opt.snn_scheme == "ranked") {
            opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::RANKED;
        } else if (all_opt.snn_scheme == "jaccard") {
            opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::JACCARD;
        } else {
            opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::NUMBER;
        }
        opt.num_threads = all_opt.num_threads;
        return opt;
    }());
    auto graph = scran_graph_cluster::convert_to_graph(snn_graph);
    const auto& weights = snn_graph.weights; // edge weights
    declare(start);
    cout_lock.unlock();

    cout_lock.lock();
    std::cout << "Multi-level graph clustering... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    auto clust_res = scran_graph_cluster::cluster_multilevel(graph, weights, [&]{
        scran_graph_cluster::ClusterMultilevelOptions opt;
        opt.resolution = all_opt.snn_resolution;
        return opt;
    }());
    declare(start);
    cout_lock.unlock();

    cout_lock.lock();
    std::cout << "Marker detection... " << std::flush;
    start = std::chrono::high_resolution_clock::now();
    const auto& best_clustering = clust_res.membership;
    auto marker_res = scran_markers::score_markers_summary(*normalized, best_clustering.data(), [&]{
        scran_markers::ScoreMarkersSummaryOptions opt;
        opt.num_threads = all_opt.num_threads;
        return opt;
    }());
    declare(start);
    cout_lock.unlock();

    second.join();
    third.join();

    /***************************
     *** Saving the results. ***
     ***************************/

    const auto& output = all_opt.output;
    if (output == "") {
        return 0;
    }
    std::filesystem::create_directories(output);
    
    // QC results.
    {
        std::ofstream handle(output + "/qc_metrics.tsv");
        handle << "sums\tdetected";
        if (has_mito) {
            handle << "\tmito_proportions";
        }
        handle << "\n";
        auto def = std::cout.precision();
        for (size_t i = 0, ncells = qc_res.sum.size(); i < ncells; ++i) {
            handle << std::setprecision(def) << qc_res.sum[i] << "\t" << qc_res.detected[i];
            if (has_mito) {
                handle << "\t" << std::setprecision(6) << qc_res.subset_proportion[0][i];
            }
            handle << "\n";
        }
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/qc_thresholds.tsv");
        handle << std::setprecision(6);
        handle << "sums\t" << qc_filters.get_sum() << "\n";
        handle << "detected\t" << qc_filters.get_detected() << "\n";
        if (has_mito) {
            handle << "mito_proportions\t" << qc_filters.get_subset_proportion()[0] << "\n";
        }
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/qc_keep.tsv");
        for (auto x : keep_cells_index) {
            handle << x << "\n";
        }
        handle << std::flush;
    }

    // Mean-variance results.
    {
        std::ofstream handle(output + "/variances.tsv");
        handle << std::setprecision(6);
        handle << "mean\tvariance\tfitted\tresid\n";
        for (size_t i = 0, end = var_res.means.size(); i < end; ++i) {
            handle  
                << var_res.means[i] << "\t" 
                << var_res.variances[i] << "\t"
                << var_res.fitted[i] << "\t"
                << var_res.residuals[i] << "\n";
        }
        handle << std::flush;
    }

    // PCA-related results, stored in transposed form, effectively.
    {
        std::ofstream handle(output + "/pca.tsv");
        handle << std::setprecision(6);
        size_t NR = pca_res.components.rows(), NC = pca_res.components.cols();
        for (size_t i = 0; i < NC; ++i) {
            handle << pca_res.components(0, i); 
            for (size_t j = 1; j < NR; ++j) {
                handle << "\t" << pca_res.components(j, i);
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
        std::filesystem::create_directory(output + "/markers");
        for (size_t i = 0, nclusters = marker_res.mean.size(); i < nclusters; ++i) {
            std::ofstream handle(output + "/markers/" + std::to_string(i) + ".tsv");
            handle << "mean\tdetected\tlfc\tdelta_detected\tcohen\tauc\n";
            handle << std::setprecision(6);
            for (size_t j = 0, ngenes = mat->nrow(); j < ngenes; ++j) {
                handle << marker_res.mean[i][j] << "\t" 
                    << marker_res.detected[i][j] << "\t"
                    << marker_res.delta_mean[i].mean[j] << "\t"
                    << marker_res.delta_detected[i].mean[j] << "\t"
                    << marker_res.cohens_d[i].mean[j] << "\t"
                    << marker_res.auc[i].mean[j] << "\n";
            }
            handle << std::flush;
        }
    }

    // Dimensionality reduction results.
    {
        std::ofstream handle(output + "/tsne.tsv");
        handle << std::setprecision(6);
        for (size_t i = 0; i < nobs; ++i) {
            handle << tsne_output[2*i] << "\t" << tsne_output[2*i + 1] << "\n";
        }
        handle << std::flush;
    }

    {
        std::ofstream handle(output + "/umap.tsv");
        handle << std::setprecision(6);
        for (size_t i = 0; i < nobs; ++i) {
            handle << umap_output[2*i] << "\t" << umap_output[2*i + 1] << "\n";
        }
        handle << std::flush;
    }

    return 0;
}
