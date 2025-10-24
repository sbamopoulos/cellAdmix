# cellAdmix

Evaluating and correcting cell admixtures in imaging-based spatial transcriptomics data.

- [Introduction](#introduction)
- [Installation](#installation)
- [Tutorials](#tutorials)
- [Citation](#citation)

## Introduction
Segmentation is a necessary part of analysing imaging-based spatial transcriptomics datasets. However, we recently found that even small errors in the segmentation can lead to apparent admixtures between adjacent cells, which creates technical artifacts in downstream analyses. This package is designed to help further clean up spatial datasets, by identifying and removing molecules believed to be assigned to the wrong cells. Our method, cellAdmix, works by identifying gene expression patterns in small regions of space. It effectively assigns these patterns to cell types and removes their molecules from cell types where they don't belong.

## Installation
Please install the latest version of the cellAdmix package from GitHub by running the following in R. It should only take a few seconds to install.
``` r
# Install the 'remotes' package if you don't have it yet
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Then install the package from GitHub
remotes::install_github("kharchenkolab/cellAdmix")
```

## Tutorials
The repository currently contains a [vignette](https://github.com/kharchenkolab/cellAdmix/blob/main/vignettes/NSCLC_tutorial_fulldata.ipynb) illustrating admixture impact on regional DE analysis, and the effect of factorization-based correction on the NSCLC dataset analyzed in our paper. We also have a second [vignette](https://github.com/kharchenkolab/cellAdmix/blob/main/vignettes/NSCLC_tutorial_sample_data.ipynb) for testing the package functionality using a provided downsampled version of this dataset (estimated runtime is ~25 minutes). We also developed a Bayesian metric to help determine how well the method worked on your data. We also include a
[vignette](https://github.com/kharchenkolab/cellAdmix/blob/main/vignettes/p_admix_tutorial.ipynb) demonstrating how to compute this metric and visualize the results.

## Citation
If you find cellAdmix useful for your publication, please cite:

[Mitchel, J., Gao, T., Cole, E., Petukhov, V. & Kharchenko, P. V. Impact of Segmentation Errors in Analysis of Spatial Transcriptomics Data. 2025.01.02.631135 Preprint at https://doi.org/10.1101/2025.01.02.631135 (2025).](https://www.biorxiv.org/content/10.1101/2025.01.02.631135v1)

