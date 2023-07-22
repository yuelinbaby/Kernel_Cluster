#!/usr/bin//bash

directory_path=$1
bin_path=$2
samid=$3
project=$4


echo "Directory path: $directory_path/location"
#generate seurat object for each spatial transcriptomics slide
Rscript   $bin_path/00.data_preprocessing.R  $project $samid $directory_path

#generate downsampling dataset based on the raw counts expression data
Rscript $bin_path/01.downsampling.R DLFPC 151673   $directory_path

#perform convolution evolution for the original expression profile
Rscript $bin_path/02.convolutional_kernel_IDW.R $project $samid $directory_path/output/

#perform autodecodor and calculated the MSE value for each kenel-optimized datasets
python $bin_path/03.autoencorder.py --indir  $directory_path/output/ --samid $samid

#perform non-spatail cluster for the convolution kernel optimized dataset 
Rscript $bin_path/04.Kernel_cluster.R  $samid  1 $directory_path/output/


