# CLI for single-cell analyses

This repository provides a no-frills command-line interface for single-cell RNA-seq data analysis from a Matrix Market file.
It is mostly intended for testing performance of the underlying C++ libraries without needing to fiddle with data analysis frameworks like R or Python.
To build, just:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

This will produce a `scran` executable in `build`, which can be run on the command line:

```bash
./build/scran --help
## Single-cell RNA-seq analyses on the command-line
## Usage: ./build/scran [OPTIONS] path
## 
## Positionals:
##   path TEXT REQUIRED          Path to the Matrix Market file
## 
## Options:
##   -h,--help                   Print this help message and exit
##   -t,--nthreads FLOAT=1       Number of threads to use (+2 for UMAP and t-SNE, which use their own threads)
##   -o,--output TEXT=output     Path to the output directory
##   --skip-output BOOLEAN=0     Run the analysis but do not save results
##   --qc-nmads FLOAT=3          Number of MADs to use for filtering
##   --hvg-span FLOAT=0.4        LOWESS span for variance modelling
##   --hvg-num INT=2500          Number of HVGs to use for PCA
##   --pca-num INT=25            Number of PCs to keep
##   --nn-approx BOOLEAN=1       Whether to use an approximate neighbor search
##   --snn-neighbors INT=10      Number of neighbors to use for the SNN graph
##   --snn-scheme ENUM:value in {jaccard->2,number->1,ranked->0} OR {2,1,0}=0
##                               Edge weighting scheme: ranked, number or jaccard
##   --snn-res FLOAT=1           Resolution to use in multi-level community detection
##   --tsne-perplexity FLOAT=30  Perplexity to use in t-SNE
##   --tsne-iter INT=500         Number of iterations to use in t-SNE
##   --umap-neighbors INT=15     Number of neighbors to use in the UMAP
##   --umap-mindist FLOAT=0.01   Minimum distance to use in the UMAP
##   --umap-epochs INT=500       Number of epochs to use in the UMAP
```

To illustrate, we'll use the Bach mammary dataset (25k cells) [here](https://github.com/clusterfork/random-test-files).
Running with 8 threads and omitting the output, we can do:

```sh
time ../../build/scran -t 8 --skip-output matrix.mtx.gz features.tsv.gz
## Initializing matrix... 9.084s
## Initializing gene annotation... 0.01s
## Computing QC metrics... 0.071s
## Computing QC thresholds... 0.003s
## Filtering cells... 0s
## Log-normalizing the counts... 0s
## Mean-variance modelling... 0.287s
## Principal components analysis... 7.805s
## Building the neighbor index... 1.137s
## Finding neighbors for clustering... 0.345s
## Finding neighbors for t-SNE... 1.317s
## Finding neighbors for UMAP... 0.413s
## SNN graph construction... 0.067s
## Multi-level clustering... 3.774s
## Marker detection... 3.443s
## UMAP calculation... 21.993s
## t-SNE calculation... 28.989s
## 
## real	0m49.508s
## user	1m53.869s
## sys	0m0.604s
```

If we were to keep the output (which is the default behavior), we would get a directory at the specified `output` path.
This contains QC metrics, variance modelling results, PCA coordinates, cluster assignments, UMAP/t-SNE values and marker statistics for each cluster.
