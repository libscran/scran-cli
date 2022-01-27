#include "scran/scran.hpp"

#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#include "tatami/ext/MatrixMarket_layered.hpp"

#include "qdtsne/qdtsne.hpp"
#include "umappp/Umap.hpp"

#include <map>
#include <memory>

int main(int argc, char* argv []) {
    /** Parsing the arguments. **/
    CLI::App app{"Single-cell RNA-seq analyses on the command-line"};
    
    std::string path;
    app.add_option("path", path, "Path to the Matrix Market file")->required();

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

    // Loading the data from an MatrixMarket file.
    size_t n = path.size();
    tatami::LayeredMatrixData<double, int> mat;
    if (n > 3 && path[n-3] == '.' && path[n-2] == 'g' && path[n-1] == 'z') {
        mat = tatami::MatrixMarket::load_layered_sparse_matrix_gzip(path.c_str());
    } else {
        mat = tatami::MatrixMarket::load_layered_sparse_matrix(path.c_str());
    }

    // Filtering out low-quality cells. 
    auto qc_res = scran::PerCellQCMetrics().run(mat.matrix.get(), { /* mito subset definitions go here */ });
    auto qc_filters = scran::PerCellQCFilters().set_nmads(nmads).run(qc_res);
    auto filtered = scran::FilterCells().run(mat.matrix, qc_filters.overall_filter.data());

    // Computing log-normalized expression values, re-using the total count from the QC step.
    auto size_factors = scran::subset_vector<false>(qc_res.sums, qc_filters.overall_filter.data());
    auto normalized = scran::LogNormCounts().run(filtered, std::move(size_factors));

    // Identifying highly variable genes.
    auto var_res = scran::ModelGeneVar().set_span(span).run(normalized.get());
    auto keep = scran::ChooseHVGs().set_top(nhvgs).run(var_res.residuals.size(), var_res.residuals.data());

    // Performing a PCA on the HVGs. We transpose the output so cells are columns again.
    auto pca_res = scran::RunPCA().set_rank(npcs).run(normalized.get(), keep.data());
    pca_res.pcs.adjointInPlace();

    // Building the nearest neighbor index.
    std::unique_ptr<knncolle::Base<> > ptr;
    if (approx) {
        ptr.reset(new knncolle::AnnoyEuclidean<>(npcs, pca_res.pcs.cols(), pca_res.pcs.data()));
    } else {
        ptr.reset(new knncolle::VpTreeEuclidean<>(npcs, pca_res.pcs.cols(), pca_res.pcs.data()));
    }

    // Performing clustering.
    auto clust_res = scran::ClusterSNNGraphMultiLevel().set_weighting_scheme(snn_scheme).set_neighbors(snn_neighbors).run(ptr.get());
    const auto& best_clustering = clust_res.membership[clust_res.max];

    // Throw in some marker detection.
    auto marker_res = scran::ScoreMarkers().run(normalized.get(), best_clustering.data());

    // Performing the t-SNE and UMAP.
    std::vector<double> tsne_output(ptr->nobs() * 2);
    qdtsne::Tsne<>().set_perplexity(tsne_perplexity).set_max_iter(tsne_iterations).run(ptr.get(), tsne_output.data());

    std::vector<double> umap_output(ptr->nobs() * 2);
    umappp::Umap<>().set_num_neighbors(umap_neighbors).set_min_dist(umap_mindist).set_num_epochs(umap_epochs).run(ptr.get(), 2, umap_output.data());

    return 0;
}
