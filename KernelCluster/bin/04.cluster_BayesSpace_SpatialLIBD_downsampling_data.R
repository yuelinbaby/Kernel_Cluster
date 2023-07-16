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
#data(dlpfc151673)
#banksy_data <- dlpfc151673
#dlpfc <- list(dlpfc151673)


####################################################
#SpatialLIBD
####################################################
#samid="Adult_Mouse_OlfactoryBulb"
#indir="/media/disk1/SYSU/datasets/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/"
#outdir="/media/disk1/SYSU/project/2021_human_dorsolateral_prefrontal_cortex_PMID33558695/output/processed_visium/"
#samid="151507"
indir=paste0(dir,"/input/")
outdir=paste0(dir,"/output/")

imgpath=indir
floc <- paste0(indir,"/tissue_positions_list.csv")


samfracs <- c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)

for (iter in 1:length(samfracs)){
  
sample_frac=samfracs[iter]
out_path_down <- paste0(outdir,"/simulation/downsample_",sample_frac)

seurat_file <- paste0(out_path_down,"/", samid, ".rds")
seurat <-readRDS(seurat_file)

ReadMtx_file <- paste0(out_path_down,"/filtered_feature_bc_matrix/matrix.mtx")
features_file <- paste0(out_path_down,"/filtered_feature_bc_matrix/genes.tsv")
cell_file <- paste0(out_path_down,"/filtered_feature_bc_matrix/barcodes.tsv")
sce <- ReadMtx(mtx=ReadMtx_file,features=features_file,cells=cell_file)
xyz <- read.csv(floc, header = FALSE)
colnames(xyz)<-c("barcode", "in_tissue", "array_row", "array_col","pxl_row_in_fullres", "pxl_col_in_fullres")
#xyz$barcode <- gsub("-", ".", xyz$barcode ,)
spots<-colnames(sce)
xyz_1 <- xyz[xyz$barcode %in% spots,]  #删除掉表达谱矩阵中没有的点，否则会报错
xyz<-xyz_1
for(i in 1:length(spots)){
  spot<-spots[i]
  xyz[i,]=xyz_1[xyz_1$barcode==spot,]
}
rownames(xyz) <- xyz$barcode


# read in image data
img <- readImgData(
  path = imgpath,
  sample_id=samid)


rowData <- read.table(file.path(features_file), header=FALSE, sep = "\t")

#genes <- rownames(sce)
#rowData <- rowDataall[rowDataall$gene_name %in% genes,]  #删除掉表达谱矩阵中没有的点，否则会报错
#for(i in 1:length(genes)){
#gene<-genes[i]
#rowData[i,]=rowDataall[rowDataall$gene_name==gene,]
#}
#rownames(xyz) <- xyz$barcode

colnames(rowData) <- c("gene_1","gene_name")

sce <- sce[,colnames(sce) %in% xyz$barcode]

xyz$pxl_row_in_fullres <- as.numeric(xyz$pxl_row_in_fullres)
xyz$pxl_col_in_fullres <- as.numeric(xyz$pxl_col_in_fullres)

spe <- SpatialExperiment(
  assays = list(counts=sce),
  rowData = rowData,
  colData = DataFrame(xyz),
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  #spatialCoordsNames = c("array_col", "array_row"),
  imgData = img,
  sample_id = samid)


##########################################################
#参考spatialLIBD流程对数据进行处理# subset to keep only spots over tissue
#17.4 Quality control (QC)
spe <- spe[, colData(spe)$in_tissue == 1]
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
rowData(spe)$gene_name[is_mito]
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

# filter low-quality spots
#spe <- spe[, !colData(spe)$discard]
#dim(spe)

#17.5 Normalization

# calculate library size factors
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))


hist(sizeFactors(spe), breaks = 20)

# calculate logcounts and store in object
spe <- logNormCounts(spe)

assayNames(spe)

#17.6 Feature selection
# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

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

#17.9 Clustering
# graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
SpatialLIBD_clus <- g_walk$membership


####################################################
#BayesSpace
####################################################
library(BayesSpace)

sqe <- SingleCellExperiment(assays=list(counts = as(sce,"dgCMatrix")),
                            rowData=rowData,
                            colData=DataFrame(xyz))

metadata(sqe)$BayesSpace.data <- list()
metadata(sqe)$BayesSpace.data$platform <- "Visium"
metadata(sqe)$BayesSpace.data$is.enhanced <- FALSE

set.seed(1314)
?spatialPreprocess
sqe <- scater::logNormCounts(sqe)
sqe <- spatialPreprocess(sqe, platform="Visium", 
                         n.PCs=7, n.HVGs=2000, log.normalize=FALSE)
#
colnames(colData(sqe)) <- c("barcode","in_tissue","row","col","imagerow","imagecol","sizeFactor")
colData(sqe)$col = as.numeric(colData(sqe)$col)
colData(sqe)$row = as.numeric(colData(sqe)$row)
#sqe <- qTune(sqe, qs=seq(2, 10), platform="Visium", d=8)
#qPlot(sqe)

#
sqe <- spatialCluster(sqe, q=7, platform="Visium", d=7,
                      init.method="mclust", model="t", gamma=2,
                      nrep=10000, burn.in=100,    #  nrep建议大于10000 ，奈何我们为了演示只能牺牲点了
                      save.chain=TRUE)

seurat <- readRDS(seurat_file)
seurat[["SpatialLIBD"]] <- SpatialLIBD_clus;
seurat[["BayesSpace"]] <- data.frame(colData(sqe))$spatial.cluster
saveRDS(seurat,seurat_file)
}
