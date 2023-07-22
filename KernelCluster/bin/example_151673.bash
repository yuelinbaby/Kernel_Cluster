#!/usr/bin//bash
directory_path=/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/
bin_path=/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/bin/
samid=151673
project=DLFPC

echo "Directory path: $directory_path/location"
#generate seurat object for each spatial transcriptomics slide
Rscript   $bin_path/00.data_preprocessing.R  $project $samid $directory_path

# generate downsampling dataset based on the raw counts expression data
Rscript $bin_path/01.downsampling.R DLFPC 151673   $directory_path

#perform convolution evolution for the original expression profile
Rscript $bin_path/02.convolutional_kernel_IDW.R $project $samid $directory_path/output/

#perform convolution for each downsampling profile
step=0.1
for ((i=1; i<=9; i++))
do
	j=0$( echo "$i * $step" |bc)
    Rscript $bin_path/02.convolutional_kernel_IDW.R $project $samid $directory_path/output/simulation/downsample_$j/ 
    echo "Running: Rscript $bin_path/02.convolutional_kernel_IDW.R $project $samid $directory_path/output/simulation/downsample_$j/";
done

#perform autodecodor and calculated the MSE value for each kenel-optimized datasets
python $bin_path/03.autoencorder.py --indir  $directory_path/output/ --samid $samid

#perform non-spatail cluster for the orginal data 
step=0.1
for ((i=1; i<=9; i++))
do
    j=0$( echo "$i * $step" |bc)
   python  $bin_path/03.autoencorder.py --indir  $directory_path/output/simulation/downsample_$j/ --samid $samid
done



#perform non-spatail cluster for the convolution kernel optimized dataset 
Rscript $bin_path/04.cluster_nonspatial_originexp.R  $samid  1 $directory_path/output/
Rscript $bin_path/04.Kernel_cluster.R  $samid  1 $directory_path/output/
Rscript $bin_path/04.cluster_BayesSpace_SpatialLIBD.R  $directory_path/ $directory_path/output/$samid.rds  $samid 

step=0.1
for ((i=1; i<=9; i++))
do
    j=0$( echo "$i * $step" |bc)
    echo "Runing clustering for simulation data"
    Rscript $bin_path/04.cluster_nonspatial_originexp.R  $samid  1 $directory_path/output/simulation/downsample_$j/
    Rscript $bin_path/04.Kernel_cluster.R  $samid  1 $directory_path/output/simulation/downsample_$j/
    Rscript $bin_path/04.cluster_BayesSpace_SpatialLIBD.R  $directory_path/ $directory_path/output/$samid.rds $samid
done





