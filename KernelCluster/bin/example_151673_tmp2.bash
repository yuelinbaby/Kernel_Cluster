#!/usr/bin//bash
directory_path=/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/
bin_path=/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/bin/
samid=151673
project=DLFPC


##run SpaGCN for the original dataset
#python /media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/bin/04.cluster_SpaGCN.py 
python $bin_path/04.cluster_SpaGCN.py   --indir  $directory_path/input/  --outdir $directory_path/output/  --samid 151673
#run SpaGCN for the  downsampling simulation dataset
step=0.1
for ((i=1; i<=9; i++))
do
    j=0$( echo "$i * $step" |bc)
    python $bin_path/04.cluster_SpaGCN_downsampling_data.py   --indir  $directory_path/input/  --outdir $directory_path/output/  --samid 151673  --iterid  $j
done

