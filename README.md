# CLI for single-cell analyses

![Unit tests](https://github.com/libscran/scran-cli/actions/workflows/run-tests.yaml/badge.svg)

## Overview

This repository implements a command-line tool for single-cell RNA-seq data analysis from a Matrix Market file.
It serves as an minimal example of how the various [**libscran**](https://github.com/libscran) components can be used together,
and is intended for developers wishing to use or extend **libscran** for their own applications.
End-users may prefer to use the [R](https://github.com/libscran/scrapper) or [Python](https://github.com/BiocPy/scranpy) bindings for actual analyses,
which are more ergonomic and have more features.

## Build instructions

First we need to build the **igraph** library, if it isn't already installed on the system.

```bash
cd extern/igraph
bash build.sh
```

Then we can proceed to building the tool itself:

```bash
cmake \
    -S . \
    -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=extern/igraph
cmake --build build
```

This will produce a `scran` executable in `build`, which can be run on the command line:

```bash
./build/scran --help
## Single-cell RNA-seq analyses on the command-line
## Usage: ./build/scran [OPTIONS] counts
## 
## Positionals:
##   counts TEXT REQUIRED        Path to the MatrixMarket file containing the counts.
## 
## Options:
##   -h,--help                   Print this help message and exit
##   -o,--output TEXT [output]   Path to the output directory. If empty, results are not saved.
##   --mito-list TEXT            Comma-separated list of the 0-based row indices of the mitochondrial genes. Closed intervals are also accepted as 'X-Y'.
##   --num-mads FLOAT [3]        Number of MADs to use for defining QC thresholds.
##   --fit-span FLOAT [0.3]      LOWESS span for fitting the mean-variance trend.
##   --num-hvgs INT [2500]       Number of HVGs to use for PCA.
##   --num-pcs INT [25]          Number of PCs to keep.
##   --nn-approx INT [1]         Whether to use an approximate neighbor search.
##   --snn-neighbors INT [10]    Number of neighbors to use for the SNN graph.
##   --snn-scheme TEXT:{ranked,number,jaccard} [ranked]
##                               Edge weighting scheme for SNN graph construction.
##   --snn-res FLOAT [0.5]       Resolution to use in multi-level community detection.
##   --tsne-perplexity FLOAT [30]
##                               Perplexity to use in t-SNE.
##   --tsne-iter FLOAT [500]     Number of iterations to use in t-SNE.
##   --umap-neighbors INT [15]   Number of neighbors to use in the UMAP.
##   --umap-mindist FLOAT [0.1]  Minimum distance to use in the UMAP.
##   --umap-epochs INT [500]     Number of epochs to use in the UMAP.
##   -t,--nthreads INT [1]       Number of threads to use (+2 for UMAP and t-SNE, which use their own threads).
```

## Usage instructions

To illustrate, we'll use the Bach mammary dataset (25k cells) [here](https://github.com/kanaverse/random-test-files).
This requires some knowledge of the row indices of the mitochondrial genes, obtained from the accompanying `features.tsv.gz` file.
With 8 threads, we can run:

```sh
time ../../build/scran \
    -t 8 \
    --mito-list 27908-27920 \
    matrix.mtx.gz

## Reading matrix from file... 4.542s
## Computing QC metrics... 0.205s
## Computing QC thresholds... 0.003s
## Filtering out low-quality cells... 0.001s
## Log-normalizing the counts... 0.001s
## Mean-variance modelling... 0.485s
## Principal components analysis... 7.249s
## Building the neighbor index... 1.122s
## Finding nearest neighbors... 1.341s
## SNN graph construction... 0.238s
## Multi-level graph clustering... 4.267s
## Marker detection... 3.121s
## t-SNE calculation... 18.011s
## UMAP calculation... 21.767s
## 
## real	0m38.073s
## user	1m36.467s
## sys	0m2.142s
```

This produces a directory at the specified `output` path, containing a variety of result files.

- `pca.tsv`: tab-separated principal components.
  This contains one row per cell after filtering and one column per component.
- `pca_varexp.tsv`: proportion of variance explained by each component.
  This has number of lines equal to the number of components.
- `qc_keep.tsv`: zero-based column indices of the cells retained after filtering.
  This has number of lines equal to the number of retained cells.
- `qc_metrics.tsv`: QC metrics computed from the input count matrix.
  This contains one row per cell in the input matrix.
  Columns correspond to different QC metrics.
- `qc_thresholds.tsv`: thresholds on the QC metrics to identify high-quality cells.
  This contains one line per metric.
- `snn_cluster_multilevel.tsv`: cluster assignments for each cell after filtering.
  This contains one line per cell after filtering, in the same order as listed in `qc_keep.tsv`.
- `tsne.tsv`: t-SNE coordinates for each cell after filtering.
  This contains one line per cell after filtering, in the same order as listed in `qc_keep.tsv`.
  Each column corresponds to one of the t-SNE dimensions.
- `umap.tsv`: UMAP coordinates for each cell after filtering.
  This contains one line per cell after filtering, in the same order as listed in `qc_keep.tsv`.
  Each column corresponds to one of the UMAP dimensions.
- `variances.tsv`: variance estimate for each gene.
  This contains one line per gene in the count matrix.
  Each column represents a different statistic.
- `markers/<CLUSTER>.tsv`: marker statistics for each gene in the specified `CLUSTER`.
  This contains one line per gene in the count matrix.
  Each column represents a different statistic.
  Cluster assignments correspond to those in `snn_cluster_multilevel.tsv`.
  
## Other comments

If I were to actually use this tool in real analyses, I would probably flesh it out to include:

- Batch correction with [**mnncorrect**](https://github.com/libscran/mnncorrect)
- Downsampling with [**nenesub**](https://github.com/libscran/nenesub)
- Gene set analyses with [**gsdecon**](https://github.com/libscran/gsdecon) and [**phyper**](https://github.com/libscran/phyper)
- Pseudo-bulk calculations with [**scran_aggregate**](https://github.com/libscran/scran_aggregate)
- ADT normalization with [**CLRm1**](https://github.com/libscran/clrm1)
- Multi-modal integration with [**mumosa**](https://github.com/libscran/mumosa)
- Count matrices from [HDF5](https://github.com/tatami-inc/tatami_hdf5) or [TileDB](https://github.com/tatami-inc/tatami_tiledb)

But I already did this for the R and Python bindings so I'm not going to do it again here. 
Any extensions are left as exercises for interested readers. 
