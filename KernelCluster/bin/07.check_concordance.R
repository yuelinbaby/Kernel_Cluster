
indir=paste0(dir,"/input/")
outdir=paste0(dir,"/output/")



#check rds file
rds1 = "/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673/151673.rds"
rds2 = "/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output/151673.rds"

rds1="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673/simulation/downsample_0.4/151673.rds"
rds2="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output/simulation/downsample_0.4/151673.rds"
seurat1=readRDS(rds1)
seurat2=readRDS(rds2)

apply(seurat1@assays$Spatial@data,2,function(x) mean(x,na.rm=T))[1:10]
apply(seurat2@assays$Spatial@data,2,function(x) mean(x,na.rm=T))[1:10]

seurat1@reductions$pca@cell.embeddings[1,]
seurat2@reductions$pca@cell.embeddings[1,]

table(seurat1$origin_dataset_Louvain,seurat2$origin_dataset_Louvain) #一致
table(seurat1$SpatialLIBD,seurat2$SpatialLIBD) #一致
table(seurat1$BayesSpace,seurat2$BayesSpace) #不一致


rds1="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673/simulation/downsample_0.2/kernel_3_3.rds"
rds2="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output/simulation/downsample_0.2/151673_kernel_optimized.rds"

rds1="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/151673/kernel_3_3.rds"
rds2="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output/151673_kernel_optimized.rds"
#rds1="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output/kernel_3_3.rds"
seurat1=readRDS(rds1)
seurat2=readRDS(rds2)
table(seurat1$optimal_dataset_Louvain,seurat2$optimal_dataset_Louvain) #一致
table(seurat1$optimal_dataset_SLM,seurat2$optimal_dataset_SLM) #一致
table(seurat1$optimal_dataset_Leiden,seurat2$optimal_dataset_Leiden) #一致



