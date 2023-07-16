args=commandArgs(T)
samid = args[1]
samfrac = as.numeric(as.character(args[2]))
dir = args[3]

#library(ggplot2)
library(Seurat)
#library(SpatialExperiment)
#library(spatialLIBD)
#library(scater)
#library(ggspavis)
#library(scran)
#library(SCENIC)
#library(SCopeLoomR)
#library(foreach)
library(tidyverse)
library(cluster)
#library(cowplot)
#library(ggpubr)


#dir="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github//test/output/"

seurat_cluster <- function(seurat,algorithmid,tag){
  seurat2=seurat
  seq_res <- seq(0.5, 1.5, 0.1)
  seurat2 <- FindNeighbors(seurat2, reduction = "pca", dims = 1:50)
  seurat2 <- FindClusters(seurat2,resolution = seq_res,verbose = F, algorithm=algorithmid)
  # Optimize ---------------------------------------------------------------------------------
  cell_dists <- dist(seurat2@reductions$pca@cell.embeddings,method = "euclidean")
  cluster_info <- seurat2@meta.data[,grepl(paste0(DefaultAssay(seurat2), "_snn_res"),colnames(seurat2@meta.data))]
  cluster_info <- dplyr::mutate_all(cluster_info,as.character)
  cluster_info <- dplyr::mutate_all(cluster_info,as.numeric)
  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    mean(si[, 'sil_width'])
  })
  optm_res <- names(which.max(silhouette_res)) #选择可以另silhouette_res值最大化的分辨率
  print (optm_res)
  seurat2[["opt_clust"]] <- seurat2[[optm_res]]
  clu<-seurat2[[optm_res]][,1]
  return(list(clu))
}

calDBI <- function(x=data,labels=labesls)
  ##DBI：任意两类别的类内样本到类中心平均距离之和除以两类中心点之间的距离，取最大值。DBI越小意味着类内距离越小，同时类间距离越大。
  ##data必须行为样本，列为特征 戴维森堡丁指数（DBI）
{
  labels2=labels
  labels=labels2[labels2 %in% names(table(labels2))[table(labels)>2]]
  clusters_n <- length(unique(labels))
  cluster_k <- list()
  for (i in c(1:clusters_n)) {
    cluster_k[[i]] <- x[which(labels==i),]
  }
  
  centroids <- list()
  for (i in c(1:clusters_n)) {
    centroids[[i]] <- apply(cluster_k[[i]],2,mean)
  }
  
  s <- list()
  for (i in c(1:clusters_n)) {
    a <- c()
    if (nrow(cluster_k[[i]])>5){
    for (j in c(1:nrow(cluster_k[[i]]))) {
      b <- dist(rbind(cluster_k[[i]][j,],centroids[[i]]),method = "euclidean")
      a <- c(a,b)
    }
    s[[i]] <- mean(a)
    }else{
      s[[i]] <- 0
    }
  }
  
  Ri <- list()
  for (i in c(1:clusters_n)){
    if (nrow(cluster_k[[i]])>5){
    r <- c()
    for (j in c(1:clusters_n)){
      if (j!=i){
        h <- (s[[i]]+s[[j]])/dist(rbind(centroids[[i]],centroids[[j]]),method = "euclidean")
        r <- c(r,h)
      }
    }
    Ri[[i]] <-  max(r,na.rm = T)
    }else{
      Ri[[i]] <- "NA"
    }
  }
  ris <- unlist(Ri);ris=ris[! ris=="NA"]
  dbi <-mean(as.numeric(ris[! ris=="NA"]))
  return(dbi)
}
#sample
#dbi <- calDBI(x,labels)#x为样本——特征矩阵（行为样本，列为特征），labels为聚类结果


calCH <- function(X,labels){ 
  ##X必须行为样本，列为特征
  labels2=labels
  labels=labels2[labels2 %in% names(table(labels2))[table(labels)>2]]
  labels_n <- length(unique(labels))
  samples_n <- nrow(X)
  X_mean <- apply(X,2,mean)
  ex_disp <- c()
  in_disp <- c()
  for (i in c(1:labels_n)) {
    cluster_k <- X[which(labels==i),]
    mean_k <- apply(cluster_k,2,mean)
    a1 <- nrow(cluster_k)*sum((mean_k-X_mean)^2)
    ex_disp <- c(ex_disp,a1)
    a2 <- sum((t(t(cluster_k)-mean_k))^2)
    in_disp <- c(in_disp,a2)
  }
  k1<- sum(ex_disp,na.rm=T)
  k2<- sum(in_disp,na.rm=T)
  if(k2==0)
  {
    return(1)
  }
  else
  {
    return((k1*(samples_n-labels_n))/(k2*(labels_n-1)))
  }
}
#sample
#ch<- calCH(X,labels)#X为样本——特征矩阵（行为样本，列为特征），labels为聚类结果
#Calinski-Harabaz（CH）

seurat_ori_cluster_sta_function <- function(seurat,tag){
  cell_dists <- dist(seurat@reductions$pca@cell.embeddings,method = "euclidean")
  x<-seurat@reductions$pca@cell.embeddings
  if (tag=="original"){
    bins=c("origin_dataset_Louvain","origin_dataset_SLM","origin_dataset_Leiden")
  }
  if (tag=="kernel"){
    bins=c("optimal_dataset_Louvain","optimal_dataset_SLM","optimal_dataset_Leiden")
  }
  tmpsta_result <- data.frame(matrix(NA,length(bins),4))
  for (binid in 1:length(bins)){
    bin=bins[binid]
    labels=as.numeric(seurat[[bin]][,1])
    SC <- mean(silhouette(labels, cell_dists)[,'sil_width']) # silhouette
    DBI <-  calDBI(x,labels) # 戴维森堡丁指数
    ch<- calCH(x,labels) #Calinski-Harabaz
    clusnum=length(unique(labels))
    tmpsta_result[binid,] <- c(clusnum,SC,DBI,ch);
  }
  return(tmpsta_result)
}


algorithms=c("Louvain","SLM","Leiden")
file_path=paste0(dir,"/",samid,".rds")
sign_change="NA";

if (file.exists(file_path)) {
    seurat_file <- paste0(dir,"/",samid,".rds")
    seurat <- readRDS(seurat_file)
    
    algorithms=c("Louvain","refine_Louvain","SLM","Leiden")
    #for(algid in 1:length(algorithms)){
    for(algid in c(1,3:4)){
    #for(algid in 1:length(algorithms)){
      algname=algorithms[algid];tag=paste0("origin_dataset_",algname);
      cluresult=seurat_cluster(seurat,algid,paste0("origin_dataset_",algname))
      clu=cluresult[[1]];
      seurat[[tag]] <- clu;
    }
    #接下来，统计每种方法得到的聚类结果的指标。
    Method=c("Louvain","SLM","Leiden")
    clusta  <- seurat_ori_cluster_sta_function(seurat,"original")
    clusta<-cbind(samid,samfrac,sign_change,Method,clusta)
    colnames(clusta)=c("sample","down.sampling.percent","sign_change","Method","num.cluster","SC","DBI","CH")
    rownames(clusta)=paste0(samid,"_",samfrac,"_",Method)
    saveRDS(seurat,seurat_file)
    write.table(clusta,paste0(dir,"/sample_original_cluster.sta"),row.names = T,col.names = T,quote=F,sep="\t")
}

