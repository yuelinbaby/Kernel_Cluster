args=commandArgs(T)
ref_cell_type_file=arg[1]
dir=args[2]
samid=args[3]

library(ggplot2)
library(Seurat)
library(tidyverse)
library(cluster)
library(cowplot)
library(ggpubr)

ref_cell_type_file="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/input/cell_type_ground_truth.class"
dir="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/output/"
samid="151673"


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

seurat_ori_cluster_sta_function <- function(seurat,ground_truth,tag){
  cell_dists <- dist(seurat@reductions$pca@cell.embeddings,method = "euclidean")
  x<-seurat@reductions$pca@cell.embeddings
  if (tag=="original"){
    bins=c("origin_dataset_Louvain","origin_dataset_SLM","origin_dataset_Leiden","SpatialLIBD","BayesSpace","SpaGCN")
  }
  if (tag=="kernel"){
    bins=c("optimal_dataset_Louvain","optimal_dataset_SLM","optimal_dataset_Leiden")
  }
  tmpsta_result <- data.frame(matrix(NA,length(bins),5))
  for (binid in 1:length(bins)){
    bin=bins[binid]
    labels=as.numeric(seurat[[bin]][,1])
    SC <- mean(silhouette(labels, cell_dists)[,'sil_width']) # silhouette
    DBI <-  calDBI(x,labels) # 戴维森堡丁指数
    ch<- calCH(x,labels) #Calinski-Harabaz
    clusnum=length(unique(labels))
    ari=mclust::adjustedRandIndex(ground_truth,labels)
    tmpsta_result[binid,] <- c(clusnum,ari,SC,DBI,ch);
  }
  return(tmpsta_result)
}



algorithms=c("Louvain","SLM","Leiden")
algorithm_spatials =c("SpatialLIBD", "BayesSpace", "SpaGCN")
ref_cell_type <- read.table(ref_cell_type_file,sep="\t",header=T)
ref_cell_type_sub=data.frame(cbind(ref_cell_type$barcode,ref_cell_type$sample_name,ref_cell_type$layer_guess))
colnames(ref_cell_type_sub)=c("barcode","sample","ground_truth")
ref_cell_type_sub[is.na(ref_cell_type_sub[,3]),]$ground_truth<-"unknown"
ref_cell_type_sub=ref_cell_type_sub[ref_cell_type_sub$sample==samid,] #dim(ref_cell_type_sub) #3592    3
ground_trutha <- ref_cell_type_sub$ground_truth


samfracs <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
result=data.frame(matrix(NA,10,5+length(algorithms)*2 ))
colnames(result) <- c("samid","downsamfrac",paste0("original_dataset_",algorithms),paste0("optimal_dataset_",algorithms),"SpatialLIBD", "BayesSpace", "SpaGCN")
idd=0
for (iter in 1:length(samfracs)){
  print(iter)
  idd=idd+1
  samfrac <- samfracs[iter]  
  out_path=paste0(dir,"/simulation/downsample_",samfrac,"/")
  if (iter==1){
    out_path=dir
  }
  seurat_ori_file <- paste0(out_path,"/",samid,".rds")
  seurat_kernel_file <- paste0(out_path,"/",samid,"_kernel_optimized.rds")
  seurat_ori <- readRDS(seurat_ori_file)
  seurat_kernel <- readRDS(seurat_kernel_file)
  
  SpaGCN_file <- paste0(dir,"/SpaGCN.result")
  SpaGCN_result = read.csv(SpaGCN_file)
  seurat_ori[["SpaGCN"]]=SpaGCN_result$refined_pred
  seurat_ori[["ground_truth"]]=ground_trutha
  seurat_kernel[["ground_truth"]]=ground_trutha
  
  result[iter,1:2]=c(samid,samfrac)
  for(algid in 1:length(algorithms)){
    algname=algorithms[algid];
    col_ori = paste0("origin_dataset_",algname); ori_label = t(seurat_ori[[as.character(col_ori)]]);
    ari_ori=mclust::adjustedRandIndex(ground_trutha,ori_label)
    col_ker =  paste0("optimal_dataset_",algname) ; ker_label = t(seurat_kernel[[as.character(col_ker)]]); 
    ari_ker=mclust::adjustedRandIndex(ground_trutha,ker_label) 
    result[iter,2+algid]=ari_ori;result[iter,5+algid]=ari_ker
  }
  for (algid in 1:length(algorithm_spatials)){
    algname=algorithm_spatials[algid]
    ori_label = t(seurat_ori[[as.character(algname)]]);
    ari_ori=mclust::adjustedRandIndex(ground_trutha,ori_label)
    result[iter,8+algid]=ari_ori;
  }
  
  #接下来，统计每种方法得到的聚类结果的指标。
  Method=c("Louvain","SLM","Leiden","kernel_Louvain","kernel_SLM","kernel_Leiden","SpatialLIBD","BayesSpace","SpaGCN")
  origin_clu_sta <- seurat_ori_cluster_sta_function(seurat_ori,ground_trutha,"original")
  kernel_clu_sta <- seurat_ori_cluster_sta_function(seurat_kernel,ground_trutha,"kernel")
  clusta <- rbind(origin_clu_sta[1:3,],kernel_clu_sta,origin_clu_sta[4:6,])
  clusta<-cbind(samid,samfrac,Method,clusta)
  colnames(clusta)=c("sample","down.sampling.percent","Method","num.cluster","ARI","SC","DBI","CH")
  rownames(clusta)=paste0(samid,"_",samfrac,"_",Method)
  saveRDS(seurat_ori,seurat_ori_file)
  saveRDS(seurat_kernel,seurat_kernel_file)
  if (idd==1){
    clustaall=clusta
  }else{
    clustaall=rbind(clustaall,clusta)
  }
}


write.table(clustaall,paste0(dir,"/",samid,".ARI.sta"),row.names = T,col.names = T,quote=F)
write.table(result,paste0(dir,"/",samid,".ARI_performance.sta"),row.names = T,col.names = T,quote=F)



