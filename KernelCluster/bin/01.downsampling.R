library(optparse)
library(Seurat)
library(dplyr)

args=commandArgs(T)
project =args[1]
samid = args[2]
dir = args[3]

in_dir=paste0(dir,"/input/")
out_path=paste0(dir,"/output/")
file_h5=paste0(in_dir,"/",samid,"_filtered_feature_bc_matrix.h5")
matrix.data <- Read10X_h5(file_h5)


#downsampling data based on the matrix.data
spot_count_sum <- apply(matrix.data,2,function(x) sum(x))
gene_count_sum <- apply(matrix.data,1,function(x) sum(x))


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

  #sampling **percent counts of each spot
  matrix.data.downsample <- matrix.data
  matrix.data.downsample[,] = 0

  #resampling without replacing. matrix.data.downsample
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

  #object downsampling seurat data
  sample_seurat <- CreateSeuratObject(counts = matrix.data.downsample, project = project,assay = "Spatial")
  #imaage info
  img <- Seurat::Read10X_Image(image.dir = in_dir)
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  #img <- img[colnames(x = anterior1)]
  img@coordinates <- img@coordinates[rownames(img@coordinates) %in% rownames(sample_seurat@meta.data),] #增加这一步骤，因为会报错Error: All cells in the image must be present in assay.
  sample_seurat[['image']] <- img

  out_file = paste0(out_path_down,"/", samid, ".rds")

  coldata <- GetTissueCoordinates(sample_seurat,
                                cols = c("row", "col", "tissue"),
                                scale = NULL)
  #sample_seurat$tissue
  sample_seurat$tissue <- coldata[rownames(sample_seurat@meta.data), "tissue"]

  sample_seurat$tissue <- ifelse(sample_seurat$tissue == 1,
                               "on_tissue",
                               "not_on_tissue")
  # Filter useful spots 
  sample_seurat <- subset(sample_seurat, subset = tissue == "on_tissue")

  # Continue with analysis
  #给sample_seurat[["orig.ident"]] 操作给sample_seurat@meta.data增加"orig.ident"这一列信息，记录spot所属的样本ID。
  sample_seurat[["orig.ident"]] <- samid

  mt_genes <- row.names(sample_seurat)[grepl("^MT-", row.names(sample_seurat))] 
  rps_genes <- row.names(sample_seurat)[grepl("^RPS", row.names(sample_seurat))] 
  mrp_genes <- row.names(sample_seurat)[grepl("^MRP", row.names(sample_seurat))] 
  rpl_genes <- row.names(sample_seurat)[grepl("^RPL", row.names(sample_seurat))] 
  rb_genes <- c(rps_genes, mrp_genes, rpl_genes)


  sample_seurat <- sample_seurat[!rownames(sample_seurat) %in% c(rb_genes, mt_genes), ]
  sample_seurat <- sample_seurat[rowSums(GetAssayData(sample_seurat, assay = "Spatial") > 0) > 10, ]
  sample_seurat$nFeature_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial") > 0)
  sample_seurat$nCount_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial"))

  # SCT transform normalization ---------------------------------------------------------
  sample_seurat <- SCTransform(sample_seurat,
                              assay = "Spatial",
                              verbose = FALSE)
  # cpm normalization ---------------------------------------------------------
  sample_seurat <- NormalizeData(sample_seurat,
                                normalization.method = 'LogNormalize',
                                scale.factor = 10000,
                                verbose = FALSE)


  # also run standard log normalization for comparison
  sample_seurat <- NormalizeData(sample_seurat, verbose = FALSE, assay = "Spatial")

  DefaultAssay(sample_seurat) <- "Spatial"


  sample_seurat<-FindVariableFeatures(sample_seurat) #增加这一步骤，因为报错Variable features haven't been set. Run FindVariableFeatures()

  sample_seurat <- ScaleData(sample_seurat,
                            verbose = FALSE,
                            features = rownames(sample_seurat)) %>%
  Seurat::RunPCA() %>%
  Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  sample_seurat <- FindNeighbors(sample_seurat, reduction = "pca", dims = 1:30)
  library(DropletUtils)
  # 将Seurat对象写入10X文件夹
  outmatfile=paste0(out_path_down,"/filtered_feature_bc_matrix")
  if (file.exists(outmatfile)==FALSE){
 	 write10xCounts(sample_seurat@assays$Spatial@counts,path=paste0(out_path_down,"/filtered_feature_bc_matrix/"))
  }
  saveRDS(sample_seurat,out_file) 
}

