#!/usr/bin//bash
directory_path=/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/
bin_path=/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/bin/
samid=151673
project=DLFPC

step=0.1
for ((i=1; i<=9; i++))
do
    j=0$( echo "$i * $step" |bc)
    echo "Rscript $bin_path/04.cluster_nonspatial_originexp.R  $samid  $j $directory_path/output/simulation/downsample_$j/"
done
