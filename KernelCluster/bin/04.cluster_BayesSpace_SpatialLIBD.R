args=commandArgs(T)
dir=args[1]
seurat_file =args[2]
samid=args[3]
samid=as.character(samid)

library(Seurat)
library(SpatialExperiment)
library(BayesSpace)
library(scater)
library(scran)

#######################################################################
indir=paste0(dir,"/input/")
outdir=paste0(dir,"/output/")
imgpath = indir
#make sure that these following files in the input director
fnm <- paste0(indir,"/",samid,"_filtered_feature_bc_matrix.h5");
fgene <- paste0(indir,"/",samid,"_filtered_feature_bc_matrix__features.tsv.gz")
#ENSG00000243485 MIR1302-2HG     Gene Expression
#ENSG00000237613 FAM138A Gene Expression
floc <- paste0(indir,"/tissue_positions_list.csv")
#ACGCCTGACACGCGCT-1,0,0,0,2427,2811
#TACCGATCCAACACTT-1,0,1,1,2547,2879
#/media/disk1/SYSU/project/bin/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/SpatialLIBD_BayesSpace.R

#help(SpatialExperiment))
sce <- Read10X_h5(fnm)

# read in image data
img <- readImgData(
  path = imgpath,
  sample_id=samid)

xyz <- read.csv(floc, header = FALSE)
colnames(xyz)<-c("barcode", "in_tissue", "array_row", "array_col","pxl_row_in_fullres", "pxl_col_in_fullres")
spots<-colnames(sce)
xyz_1 <- xyz[xyz$barcode %in% spots,]  #删除掉表达谱矩阵中没有的点，否则会报错
xyz<-xyz_1
for(i in 1:length(spots)){
  spot<-spots[i]
  xyz[i,]=xyz_1[xyz_1$barcode==spot,]
}

rd <- data.frame(read.table(fgene,sep="\t",header=F))
colnames(rd) <- c("ensemblid","gene_name","type")

####################################################
#SpatialLIBD
####################################################
library(spatialLIBD)

spe <- SpatialExperiment(
  assays = list(counts = sce),
  rowData = rd,
  colData = DataFrame(xyz),
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  imgData = img,
  sample_id = samid)


#17.4 Quality control (QC)
spe <- spe[, colData(spe)$in_tissue == 1]
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 600
qc_detected <- colData(spe)$detected < 400
qc_mito <- colData(spe)$subsets_mito_percent > 28
qc_cell_count <- colData(spe)$cell_count > 10

# number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito, qc_cell_count), 2, sum)
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
table(discard)
#17.5 Normalization

# calculate library size factors
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))

# calculate logcounts and store in object
spe <- logNormCounts(spe)

assayNames(spe)

#17.6 Feature selection
# remove mitochondrial genes
spe <- spe[!is_mito, ]
# fit mean-variance relationship
dec <- modelGeneVar(spe)
# visualize mean-variance relationship
fit <- metadata(dec)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
#17.8 Dimensionality reduction
# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")

reducedDimNames(spe)

dim(reducedDim(spe, "UMAP"))
# update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
SpatialLIBD_clus <- g_walk$membership
table(SpatialLIBD_clus)

####################################################
#BayesSpace
####################################################
library(BayesSpace)

sqe <- SingleCellExperiment(assays=list(counts = as(sce,"dgCMatrix")),
                            rowData=rd,
                            colData=DataFrame(xyz))

metadata(sqe)$BayesSpace.data <- list()
metadata(sqe)$BayesSpace.data$platform <- "Visium"
metadata(sqe)$BayesSpace.data$is.enhanced <- FALSE

set.seed(1314)
?spatialPreprocess
sqe <- scater::logNormCounts(sqe)
sqe <- spatialPreprocess(sqe, platform="Visium", 
                         n.PCs=7, n.HVGs=2000, log.normalize=FALSE)
#注意，这里colnames需要修正成软件所需要的正确的列名，否则会报错
colnames(colData(sqe)) <- c("barcode","in_tissue","row","col","imagerow","imagecol","sizeFactor")
colData(sqe)$col = as.numeric(colData(sqe)$col)
colData(sqe)$row = as.numeric(colData(sqe)$row)
#sqe <- qTune(sqe, qs=seq(2, 10), platform="Visium", d=7)
#qPlot(sqe)
##setting the q value to be 7 by default
sqe <- spatialCluster(sqe, q=7, platform="Visium", d=7,
                      init.method="mclust", model="t", gamma=2,
                      nrep=10000, burn.in=100,   
                      save.chain=TRUE)


seurat <- readRDS(seurat_file)
seurat[["SpatialLIBD"]] <- SpatialLIBD_clus;
seurat[["BayesSpace"]] <- data.frame(colData(sqe))$spatial.cluster
saveRDS(seurat,seurat_file)








