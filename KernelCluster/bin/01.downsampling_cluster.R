args=commandArgs(T)
in_dir= args[1]
out_dir= args[2]
samid=args[3]

library(ggplot2)
library(Seurat)
library(SpatialExperiment)
library(spatialLIBD)
library(scater)
library(ggspavis)
library(scran)
library(SCENIC)
library(SCopeLoomR)
library(foreach)
library(tidyverse)
library(cluster)
library(cowplot)
library(ggpubr)

#说明：这个是用利用grid方法，每个位点上下左右各扩展k个坐标，然后求得均值/stouffer——Z值用来代替当前这个位点的表达值。
#这样求出来的表达谱矩阵的大小是没有改变的。
#saminfo <- read.table("/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/data_back_bone",sep="\t",header=T)
#sams <- saminfo$sample_id
#indir="/media/disk1/SYSU/datasets/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/"
#outdir="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/"

#构建seurat分析对象

#methods=c("seurat_clusters","opt_clust","SpatialLIBD","SpaGCN","BayesSpace","STdeconvolve")

#for(sam in 1:length(sams)){

#samid=sams[sam]  #
#samid="151670"
matrix_dir=in_dir 
image_dir=paste0(in_dir,"/",samid)
out_path <- paste0(out_dir,"/",samid)
file_h5=paste0(matrix_dir,"/",samid,"_filtered_feature_bc_matrix.h5")
matrix.data <- Read10X_h5(file_h5)
rdsfile=paste0(out_dir,"/",samid,"/",samid,".rds")
seurat <- readRDS(rdsfile)


###ground_truth的细胞类别信息
ref_cell_type_file="/media/disk1/SYSU/datasets/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/cell_type_ground_truth.class"
ref_cell_type <- read.table(ref_cell_type_file,sep="\t",header=T)

ref_cell_type_sub=data.frame(cbind(ref_cell_type$barcode,ref_cell_type$sample_name,ref_cell_type$layer_guess))
colnames(ref_cell_type_sub)=c("barcode","sample","ground_truth")
ref_cell_type_sub[is.na(ref_cell_type_sub[,3]),]$ground_truth<-"unknown"
ref_cell_type_sub=ref_cell_type_sub[ref_cell_type_sub$sample==samid,] #dim(ref_cell_type_sub) #3592    3
ground_truth <- ref_cell_type_sub$ground_truth

#构建seurat分析对象
sample_seurat <- CreateSeuratObject(counts = matrix.data, project = "humanDLPFC",assay = "Spatial")
#加入空间转录组的图像信息
img <- Seurat::Read10X_Image(image.dir = image_dir,image.name = "tissue_lowres_image.png")
Seurat::DefaultAssay(object = img) <- 'Spatial'
#img <- img[colnames(x = anterior1)]
sample_seurat[['image']] <- img


#打算在matrix.data的基础上进行下采样

spot_count_sum <- apply(matrix.data,2,function(x) sum(x))
gene_count_sum <- apply(matrix.data,1,function(x) sum(x))


#下采样的比例 0.8,0.6,0.4,0.2,0.1
#先尝试0.1比例下采样的效果
out_path_sim <- paste0(out_path,"/simulation/")
if(!file.exists(out_path_sim)) {
  dir.create(out_path_sim)
}


samfracs <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
for (iter in 1:length(samfracs)){
  sample_frac=samfracs[iter]
out_path_down <- paste0(out_path,"/simulation/downsample_",sample_frac)

if(!file.exists(out_path_down)) {
  dir.create(out_path_down)
}

#以每个细胞为单位，抽取每个细胞中原来总的counts*dsr个转录结果
matrix.data.downsample <- matrix.data
matrix.data.downsample[,] = 0

#抽样结果存入矩阵matrix.data.downsample
for (i in 1:ncol(matrix.data)){
  s <- matrix.data[,i]
  pools <- rep(1:length(s),s)
  poolsub <- sample(pools,length(pools)*sample_frac,replace = FALSE)
  counts <- table(poolsub)
  samout <- cbind(as.numeric(names(counts)),as.numeric(counts))
  matrix.data.downsample[samout[,1],i] <- samout[,2] 
}

#累积概率分布图
spot_count_sum_down <- apply(matrix.data.downsample,2,function(x) sum(x))
gene_count_sum_down <- apply(matrix.data.downsample,1,function(x) sum(x))

pdf (paste0(out_path_down,"/spot_gene.pdf"),width=12,height=8)
par(mfrow=c(2,2))
plot(ecdf(gene_count_sum),xlim=c(0,6000),col="black",xlab="exp_count_sum_per_gene",ylab="ecdf",main="")
lines(ecdf(gene_count_sum_down),col="red")

plot(density(gene_count_sum_down),xlim=c(0,6000),col="red",xlab="exp_count_sum_per_gene",ylab="density",main="")
lines(density(gene_count_sum),col="black")

plot(ecdf(spot_count_sum),xlim=c(0,20000),col="black",xlab="exp_count_sum_per_spot",ylab="ecdf",main="")
lines(ecdf(spot_count_sum_down),col="red")

plot(density(spot_count_sum_down),xlim=c(0,20000),col="red",xlab="exp_count_sum_per_spot",ylab="density",main="")
lines(density(spot_count_sum),col="black")

dev.off()

#构建下采样seurat分析对象
sample_seurat <- CreateSeuratObject(counts = matrix.data.downsample, project = samid,assay = "Spatial")
#加入空间转录组的图像信息
img <- Seurat::Read10X_Image(image.dir = image_dir)
Seurat::DefaultAssay(object = img) <- 'Spatial'
#img <- img[colnames(x = anterior1)]
img@coordinates <- img@coordinates[rownames(img@coordinates) %in% rownames(sample_seurat@meta.data),] #增加这一步骤，因为会报错Error: All cells in the image must be present in assay.
sample_seurat[['image']] <- img

out_file = paste0(out_path_down,"/", samid, ".rds")

coldata <- GetTissueCoordinates(sample_seurat,
                                cols = c("row", "col", "tissue"),
                                scale = NULL)
#sample_seurat$tissue组织的代号，这里没有组织分类，所以都是1
sample_seurat$tissue <- coldata[rownames(sample_seurat@meta.data), "tissue"]

sample_seurat$tissue <- ifelse(sample_seurat$tissue == 1,
                               "on_tissue",
                               "not_on_tissue")
# Filter useful spots 
sample_seurat <- subset(sample_seurat, subset = tissue == "on_tissue")

# Continue with analysis
#给sample_seurat[["orig.ident"]] 操作给sample_seurat@meta.data增加"orig.ident"这一列信息，记录spot所属的样本ID。
sample_seurat[["orig.ident"]] <- samid

mt_genes <- row.names(sample_seurat)[grepl("^MT-", row.names(sample_seurat))] #线粒体基因
rps_genes <- row.names(sample_seurat)[grepl("^RPS", row.names(sample_seurat))] #核糖体蛋白基因
mrp_genes <- row.names(sample_seurat)[grepl("^MRP", row.names(sample_seurat))] #线粒体核小体蛋白？
rpl_genes <- row.names(sample_seurat)[grepl("^RPL", row.names(sample_seurat))] #RPL核糖体蛋白基因
rb_genes <- c(rps_genes, mrp_genes, rpl_genes)


sample_seurat <- sample_seurat[!rownames(sample_seurat) %in% c(rb_genes, mt_genes), ]
sample_seurat <- sample_seurat[rowSums(GetAssayData(sample_seurat, assay = "Spatial") > 0) > 10, ]
sample_seurat$nFeature_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial") > 0)
sample_seurat$nCount_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial"))

# SCT transform normalization ---------------------------------------------------------
sample_seurat <- SCTransform(sample_seurat,
                              assay = "Spatial",
                              verbose = FALSE)
DefaultAssay(sample_seurat) <- "SCT" #这一步做了设定，后续使用sample_seurat的数据都某人要用SCT转换后的结果。
# cpm normalization ---------------------------------------------------------
sample_seurat <- NormalizeData(sample_seurat,
                                normalization.method = 'LogNormalize',
                                scale.factor = 10000,
                                verbose = FALSE)

#sct 与log-归一化相比，结果如何?
#  为了探究规范化方法的不同，我们研究了sctransform和log规范化结果如何与UMIs的数量相关。
# also run standard log normalization for comparison
sample_seurat <- NormalizeData(sample_seurat, verbose = FALSE, assay = "Spatial")

DefaultAssay(sample_seurat) <- "Spatial"
#DefaultAssay(sample_seurat) <- "SCT" #silhouette_res不能正常运行会报错

sample_seurat<-FindVariableFeatures(sample_seurat) #增加这一步骤，因为报错Variable features haven't been set. Run FindVariableFeatures()

sample_seurat <- ScaleData(sample_seurat,
                            verbose = FALSE,
                            features = rownames(sample_seurat)) %>%
  Seurat::RunPCA() %>%
  Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)


sample_seurat <- FindNeighbors(sample_seurat, reduction = "pca", dims = 1:30)

seq_res <- seq(0.5, 1.5, 0.1)

sample_seurat <- FindClusters(sample_seurat,
                               resolution = seq_res,
                               verbose = F)
# Optimize ---------------------------------------------------------------------------------

cell_dists <- dist(sample_seurat@reductions$pca@cell.embeddings,
                   method = "euclidean")


cluster_info <- sample_seurat@meta.data[,grepl(paste0(DefaultAssay(sample_seurat), "_snn_res"),
                                                colnames(sample_seurat@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)


#根据给定的墨迹聚类计算轮廓信息。 
silhouette_res <- apply(cluster_info, 2, function(x){
  si <- silhouette(x, cell_dists)
  #if(!is.na(si)) {
  mean(si[, 'sil_width'])
  #} else {
  #  NA
  #}
})
optm_res <- names(which.max(silhouette_res))
sample_seurat[["opt_clust"]] <- sample_seurat[[optm_res]]
#mclust::adjustedRandIndex(ground_truth,sample_seurat@meta.data$opt_clust) 0.2657
saveRDS(sample_seurat,file=out_file)

library(DropletUtils)
# 将Seurat对象写入10X文件夹
write10xCounts(sample_seurat@assays$Spatial@counts,path=paste0(out_path_down,"/filtered_feature_bc_matrix"))


}

