# CLI for single-cell analyses

This repository provides a no-frills command-line interface for single-cell RNA-seq data analysis from a Matrix Market file.
It is mostly intended for testing performance of the underlying C++ libraries without needing to fiddle with data analysis frameworks like R or Python.
To build, just:

```
cmake -S . -B build
cmake --build build
```

This will produce a `scran` executable in `build`, which can be executed to obtain a help page:

```
$ ./build/scran --help
Single-cell RNA-seq analyses on the command-line
Usage: ./build/scran [OPTIONS] path

Positionals:
  path TEXT REQUIRED          Path to the Matrix Market file

Options:
  -h,--help                   Print this help message and exit
  --qc-nmads FLOAT=3          Number of MADs to use for filtering
  --hvg-span FLOAT=0.4        LOWESS span for variance modelling
  --hvg-num INT=2500          Number of HVGs to use for PCA
  --pca-num INT=25            Number of PCs to keep
  --nn-approx BOOLEAN=1       Whether to use an approximate neighbor search
  --snn-neighbors INT=10      Number of neighbors to use for the SNN graph
  --snn-scheme ENUM:value in {jaccard->2,number->1,ranked->0} OR {2,1,0}=0
                              Edge weighting scheme: ranked, number or jaccard
  --tsne-perplexity FLOAT=30  Perplexity to use in t-SNE
  --tsne-iter INT=500         Number of iterations to use in t-SNE
  --umap-neighbors INT=15     Number of neighbors to use in the UMAP
  --umap-mindist FLOAT=0.01   Minimum distance to use in the UMAP
  --umap-epochs INT=500       Number of epochs to use in the UMAP
```
