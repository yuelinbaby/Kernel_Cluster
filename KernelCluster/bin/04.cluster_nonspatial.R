

##############################################################################################################
#利用IDW矩阵对数据进行smooth，然后聚类的子程序
##############################################################################################################

smooth_cluster <- function(inseurat, Weight_Matrix,k,p){
  result=list()
  data <- inseurat@assays$Spatial@data
  outtag = paste0("kenel_",k*2+1,"_",p)
  data_smooth <- data %*% as.matrix(t(Weight_Matrix))
  seurat_smooth <- inseurat
  seurat_smooth@assays$Spatial@data <-  as.matrix(data_smooth)
  seurat_smooth@assays$SCT@data <-  as.matrix(data_smooth)
  seurat_smooth <-  FindVariableFeatures(seurat_smooth, selection.method = "vst", nfeatures = 2000)
  seurat_smooth  <- ScaleData(seurat_smooth,verbose = FALSE,features = rownames(seurat_smooth))
  seurat_smooth <- Seurat::RunPCA(seurat_smooth,npcs = 50) 
  seurat_smooth <- Seurat::RunUMAP(seurat_smooth,reduction = 'pca', dims = 1:30, verbose = FALSE)
  
  seq_res <- seq(0.5, 1.5, 0.1)
  seurat_smooth <- FindNeighbors(seurat_smooth, reduction = "pca", dims = 1:50)
  seurat_smooth<- FindClusters(seurat_smooth,
                               resolution = seq_res,
                               verbose = F)
  # Optimize ---------------------------------------------------------------------------------
  
  cell_dists <- dist(seurat_smooth@reductions$pca@cell.embeddings,
                     method = "euclidean")
  cluster_info <- seurat_smooth@meta.data[,grepl(paste0(DefaultAssay(seurat_smooth), "_snn_res"),colnames(seurat_smooth@meta.data))] 
  cluster_info <- dplyr::mutate_all(cluster_info,as.character) 
  cluster_info <- dplyr::mutate_all(cluster_info,as.numeric)
  #根据给定的墨迹聚类计算轮廓信息。 
  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    #if(!is.na(si)) {
    mean(si[, 'sil_width'])
    #} else {
    #  NA
    #}
  })
  optm_res <- names(which.max(silhouette_res)) #选择可以另silhouette_res值最大化的分辨率
  seurat_smooth[["opt_clust"]] <- seurat_smooth[[optm_res]]
  
  #画出top10个高变异基因的空间表达谱分布图
  genes=TopFeatures(seurat_smooth,nfeatures = 10)
  #SpatialFeaturePlot(seurat,features=genes,ncol=5) +ggtitle("seurat")
  result[[1]] <- SpatialFeaturePlot(seurat_smooth,features=genes,ncol=5) + ggtitle("seurat_smooth")
  #SpatialDimPlot(seurat,group.by = "opt_clust") + ggtitle("seurat")
  result[[2]] <- SpatialDimPlot(seurat_smooth,group.by = "opt_clust") + ggtitle(outtag)
  result[[3]]<-seurat_smooth@meta.data$opt_clust
  result[[4]] <- seurat_smooth
  result[[5]] <- DimPlot(object = seurat_smooth, reduction = 'umap') + ggtitle(outtag)
  return(result)
}



