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

The KernelCluster package can be used by downloading the bash and script bin file.
- The KernelCluster base file [bash file](https://github.com/yuelinbaby/Kernel_Cluster/blob/main/KernelKcluster.bash)
- Scripts:[bin file](https://github.com/yuelinbaby/Kernel_Cluster/tree/master/KernelCluster) 

## Running Instructions
 ./KernelKcluster.bash  sample_dir/ bin_dir/ sample_id  project_id
 - sample_dir: the directory of sample, with two subdirectory, "input" and "ouput". The "input subdirectory" store SRT raw data.
 - bin_dir: the directory of KernelCluster
 - sample_id: the sample or slice ID of SRT dataset
 - project_id: the project id of SRT dataset


## Dependencies

KernelCluster relies on the following packages:

- [Seurat] R packages(v4.2.0)
- [dplyr] R packages(v1.1.2)
- [tidyverse] R packages(v1.3.2)
- [cluster] R packages(v2.1.4)
- [keras] Python packages(v2.12.0)
- [numpy] Python packages(v1.22.4)
- [pandas] Python packages(v1.4.3)

To install the required packages, you can use:
- for R packages: BiocManager::install('**')
- for python packages: pip install **

The spatial clustering packages listed below are optional in KernelCluster, serving solely for validating the performance of the method.
- [BayesSpace] R packages(v1.6.0)
- [spatialLIBD] R packages(v1.8.11)
- [SpaGCN] Python packages(v1.2.5)

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


