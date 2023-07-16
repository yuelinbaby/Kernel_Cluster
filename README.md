# Kernel_Cluster

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/yuelinbaby/Kernel_Cluster))
<!-- badges: end -->

## Overview

KernelCluster utilizes convolution kernels to aggregate gene expression
from adjacent spots, thereby increasing the robustness of clustering. 
By incorporating neighborhood information for clustering, KernelCluster is able to

-   decrease the sparsity of the expression matrix creaated in the spatial transcriptomics 
-   improve the accuracy to infer the spatial domains in spatial transcriptomics datasets
-   allows for the incorporation of various existing clustering
algorithms

KernelCluster is applicable to spatial transcriptomics datasets based in situ sequencing (ISS) technology (e.g. 10x
Visium) and scales well to large datasets. 

-   the
    [preprint](https://www.biorxiv.org/content/10.1101/2022.04.14.488259v1),


## Installation

The *Banksy* package can be installed via `remotes`:

``` r
remotes::install_github("prabhakarlab/Banksy", dependencies = TRUE)
```

Installation should take less than three minutes.

Installation should take less than three minutes.

**Known installation issues**

1.  Installation of `leidenAlg` has non-zero exit status

-   Refer to the [package
    website](https://github.com/kharchenkolab/leidenAlg#installation)
    for *leidenAlg* installation details. Otherwise, users may also
    install a separate branch of *Banksy* with

``` r
remotes::install_github("prabhakarlab/Banksy@feat-igraph-leiden")
```

## Documentation

Detailed description of *Banksy* functionality and example analyses are
available at the [package webpage]().

-   [DLPFC dataset: Human dorsolateral prefrontal cortex data](http://research.libd.org/spatial LIBD/)):
    Illustrates a grid search of parameters which best cluster cells.
(1) DLPFC dataset: Human dorsolateral prefrontal cortex data (http://research.libd.org/spatial LIBD/) 
-   [Human dorsolateral prefrontal cortex 10x Visium
    Illustrates analysis of multiple spatial transcriptomic datasets.
    
*KernelCluster* is also interoperable with
[Seurat](https://satijalab.org/seurat/) via *SeuratWrappers*.
Documentation on how to run KernelCluster on Seurat objects can be found
[here](/banksy.md).


