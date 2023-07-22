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


## Installation

The KernelCluster package can be used by copying the bash file  and chmod the bash file as executable file

The KernelCluster [bash file](https://github.com/yuelinbaby/Kernel_Cluster/blob/main/KernelKcluster.bash)


## Documentation

datasets

-   [DLPFC dataset: Human dorsolateral prefrontal cortex data](http://research.libd.org/spatial LIBD/)):
    Illustrates a grid search of parameters which best cluster cells.
(1) DLPFC dataset: Human dorsolateral prefrontal cortex data (http://research.libd.org/spatial LIBD/) 
-   [Human dorsolateral prefrontal cortex 10x Visium
    Illustrates analysis of multiple spatial transcriptomic datasets.
    
*KernelCluster* is also interoperable with
[Seurat](https://satijalab.org/seurat/) via *SeuratWrappers*.
Documentation on how to run KernelCluster on Seurat objects can be found
[here](/banksy.md).


