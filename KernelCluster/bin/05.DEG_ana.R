args=commandArgs(T)
dir=args[1]
samid=args[2]

library(ggplot2)
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
library(cowplot)
library(ggpubr)

indir=paste0(dir,"/input/")
outdir=paste0(dir,"/output/")


##DEG analysis 
  samoutdir=paste0(outdir,"/")
  setwd(samoutdir)
  samfracs <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
  
  #spot类别的ground-truth
  ref_cell_type_file=paste0(indir,"/cell_type_ground_truth.class")
  ref_cell_type <- read.table(ref_cell_type_file,sep="\t",header=T)
  ref_cell_type_sub=data.frame(cbind(ref_cell_type$barcode,ref_cell_type$sample_name,ref_cell_type$layer_guess))
  colnames(ref_cell_type_sub)=c("barcode","sample","ground_truth")
  ref_cell_type_sub[is.na(ref_cell_type_sub[,3]),]$ground_truth<-"unknown"
  ref_cell_type_sub=ref_cell_type_sub[ref_cell_type_sub$sample==samid,] #dim(ref_cell_type_sub) #3592    3
  ground_truth <- ref_cell_type_sub$ground_truth

  #stafile=paste0(outdir,"/",samid,"/kenel_louvain_clustering.sta")
  #sta_result <- data.frame(matrix(NA,length(samfracs),24))
  for (iter in 1:length(samfracs)){
    idd=idd+1
    print(iter)
    samfrac <- samfracs[iter]  
    out_path=paste0(samoutdir,"/simulation/downsample_",samfrac,"/")
    if (iter==1){
      out_path=samoutdir
    }
    encorfile=paste0(out_path,"/",samid,".encorder.mse.csv")
    arifile=paste0(out_path,"/IDW_smooth_clustering_ARI.sta")
    encor <- read.table(encorfile,sep=" ",header=F)
    ari <- read.table(arifile,header=T,sep="\t");ari=t(ari)[,1:10]
    rownames(encor)= seq(1,5,1)*2+1
    colnames(encor)= paste0("p",seq(0,9,1))
    #每个winsize都找到导数值最小的power
    #for (i in 1:nrow(encor)){
    i=3 #只看k=3的结果
    dy=round(diff(as.numeric(encor[i,])),5)
    sign_change <- min(which(diff(sign(dy)) != 0))
    print(sign_change) #曲线第一个最低点（encorder mse的第一个局部最优解）对应的power值
    loc=sign_change+1 
    ariv=ari[i,loc]
    seurat_ori_file=paste0(out_path,"/",samid,".rds")
    seurat_kernel_file=paste0(out_path,"/kernel_3_",sign_change,".rds")
    seurat_ori=readRDS(seurat_ori_file)
    seurat_ori@meta.data$ground_truth = ground_truth
    seurat_kernel=readRDS(seurat_kernel_file)
    #groundtruth
    ari1=mclust::adjustedRandIndex(ground_truth,seurat_ori@meta.data$opt_clust)
    ari2=mclust::adjustedRandIndex(ground_truth,seurat_kernel@meta.data$opt_clust)
    table(ground_truth,seurat_ori@meta.data$opt_clust)
    table(ground_truth,seurat_kernel@meta.data$opt_clust)
    
    DEG_ana <- function(seurat_smooth){
    Idents(seurat_smooth) <- seurat_smooth@meta.data$ground_truth
    cluster_tag <- "ground_truth"
    clusters <- unique(seurat_smooth@meta.data[,cluster_tag]);clusters=clusters[order(clusters)]
    clusters <- clusters[table(seurat_smooth@meta.data$ground_truth)>10]
    clusters_all <- c("Layer1","Layer2","Layer3","Layer4","Layer5","Layer6","WM","unknown")
    clusters <- clusters_all[clusters_all %in% clusters]
    DEG_result <- data.frame(matrix(0,nrow(seurat_smooth@assays$Spatial@data),length(clusters)))
    row.names(DEG_result)<- row.names(seurat_smooth@assays$Spatial@data)
    colnames(DEG_result) <- paste0(clusters)
    DEG_result_value <- DEG_result
    for (i in 1:length(clusters)){
      cid=clusters[i]
      DEG=FindMarkers(seurat_smooth,ident.1 = cid,slot="data",min.cells.feature =3,logfc.threshold = 0.25,min.pct = 0.25,test.use = "wilcox",min.cells.group = 3)
      upgene <- row.names(DEG[DEG$avg_log2FC>0.25 & DEG$p_val<1e-10,])
      downgene <- row.names(DEG[DEG$avg_log2FC< (-0.25) & DEG$p_val<1e-10,])
      DEG_result[upgene,i]=1
      DEG_result[downgene,i]=-1
      DEG_result_value[upgene,i]<-DEG[DEG$avg_log2FC>0.25 & DEG$p_val<1e-10,]$avg_log2FC
      DEG_result_value[downgene,i]<-DEG[DEG$avg_log2FC< (-0.25) & DEG$p_val<1e-10,]$avg_log2FC
    }
    ##画出heatmap图观察一下给DEG矩阵能否给当前的细胞进行聚类，看看利用离散值聚类（-1,0,1）和利用foldchange聚类结果差异大吗
    DEG_result_sub<-DEG_result[apply(DEG_result,1,function(x) sum(x[x>0]))>0 | apply(DEG_result,1,function(x) sum(x[x<0]))>0,]
    DEG_result_value_sub<-DEG_result_value[apply(DEG_result,1,function(x) sum(x[x>0]))>0 | apply(DEG_result,1,function(x) sum(x[x<0]))>0,]
    DEG_out=data.frame(matrix(NA,nrow(DEG_result_sub),2*length(clusters_all)))
    for (i in 1:length(clusters_all)){
        clu=clusters_all[i]
        if ( length(colnames(DEG_result_sub)[colnames(DEG_result_sub)==clu])>0){
          DEG_out[,i]=DEG_result_sub[,colnames(DEG_result_sub)==clu]
          DEG_out[,(length(clusters_all)+i)]=DEG_result_value_sub[,colnames(DEG_result_sub)==clu]
          DEG_out[1:4,]
        }
    }
    return(DEG_out)
    }
    ##收集原始数据vskernel优化后的各类别的差异表达基因集合
    DEG_a <- DEG_ana(seurat_ori); DEG_a <- data.frame(cbind(samid,samfrac,"origin",DEG_a)); colnames(DEG_a) <- c("samid","down_sam_frac","dataset",paste0(rep(clusters_all,2),"_",c(rep("tag",8),rep("value",8))))
    DEG_b <- DEG_ana(seurat_kernel); DEG_b <- data.frame(cbind(samid,samfrac,"kernel",DEG_b)) ; colnames(DEG_b) <- c("samid","down_sam_frac","dataset",paste0(rep(clusters_all,2),"_",c(rep("tag",8),rep("value",8))))
    DEG_M <- rbind(DEG_a, DEG_b)
    if (idd==1){
      DEG_A <- DEG_M
    }else{
      DEG_A <- rbind(DEG_A,DEG_M)
    }
  }  

write.table(DEG_A,file=paste0(outdir,"/encorder_optim_ari_DEGgene.sta"),col.names = T,row.names = F,quote=F,sep="\t")




