args=commandArgs(T)
project =args[1]
samid = args[2]
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
#library(cluster)
#library(cowplot)
#library(ggpubr)

setwd(dir)

rdsfile=paste0(dir,"/",samid,".rds")
seurat <- readRDS(rdsfile)

##############################################################################################################################
data <- seurat@assays$Spatial@data

##############################################################################################################
#obtain the value of each elements in kernel, based on IDW algorithm
##############################################################################################################

coordinates <- seurat@images$image@coordinates
x <- seurat@images$image@coordinates$row
y <- seurat@images$image@coordinates$col

creat_IDW_W_matrix <- function(k,dp){
  W=data.frame(matrix(0,nrow(coordinates),nrow(coordinates)))
  rownames(W)=rownames(coordinates)
  colnames(W)=rownames(coordinates)
  for(i in 1:nrow(coordinates)){
    #  print(i)
    x <- coordinates$row[i]
    y <- coordinates$col[i]
    spot=rownames(coordinates[i,])
    xs<-x-k;xe<-x+k
    ys<-y-k;ye<-y+k
    block <-  coordinates[coordinates$row>=xs & coordinates$row<=xe & coordinates$col>=ys & coordinates$col<=ye ,] 
    dis <- sqrt((block$row-x)^2 + (block$col-y)^2);dis[dis==0]=1
    #set the self distance to be 1
    w <- 1/dis^dp;  w <- w/sum(w) # normalization the sum 
    W[i, rownames(block)]=w #
  }
  return(W)
}

convolution_kernel_perform<- function(inseurat, Weight_Matrix,k,p){
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
  seurat_smooth <- FindNeighbors(seurat_smooth, reduction = "pca", dims = 1:50)
  result[[1]] <- seurat_smooth
  return(result)
}

##############################################################################################################
#convolution kernel data optimization
##############################################################################################################


pbins <- seq(0,9,1)  ## test the power
  k=1; kernel_size=k*2+1; ##setting the size of kernel as three
  for (j in 1:length(pbins)){
      p<-pbins[j]
      print(c(k,p));
      IDW <-  creat_IDW_W_matrix(k,p) 
      seurat_smooth <- convolution_kernel_perform(seurat, IDW, k,p);
      seurat_out <- seurat_smooth[[1]]
      outseuratfile <- paste0(dir,"/kernel_",kernel_size,"_",p,".rds");
      saveRDS(seurat_out,outseuratfile)
      exp <- seurat_out@assays$Spatial@data
      outfile=paste0(dir,"/kernel_",kernel_size,"_",p,".exp_data.csv")
      write.csv(exp,outfile)
   }



