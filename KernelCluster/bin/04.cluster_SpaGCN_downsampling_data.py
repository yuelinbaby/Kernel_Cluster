import os,csv,re
from numpyro import sample
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import pandas as pd
import sys
print (sys.argv)



import argparse
parser = argparse.ArgumentParser(description='para transfer')
parser.add_argument('--indir', type=str, help='indir -> inputdir.')
parser.add_argument('--outdir', type=str, help='outdir -> output dir.')
parser.add_argument('--samid', type=str, help='samid -> sample id.')
parser.add_argument('--iterid', type=str, help='iterid -> iterid')

args = parser.parse_args()

indir = args.indir
outdir = args.outdir
samid = args.samid
sample_frac =args.iterid

import cv2

spg.__version__


out_path=f'{outdir}/'
out_path_down=f'{outdir}/simulation/downsample_{sample_frac}/'
max_dir=f'{outdir}/simulation/downsample_{sample_frac}/filtered_feature_bc_matrix/'

adata=sc.read_10x_mtx(f'{max_dir}')
#Read in hitology image

spafile = f'{indir}/tissue_positions_list.csv'
if not os.path.exists(f'{spafile}'):
    spafile = f'{indir}/tissue_positions.csv'

spatial=pd.read_csv(spafile,sep=",",header=None,na_filter=False,index_col=0)

adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5]
adata.obs["x_array"]=adata.obs["x2"]
adata.obs["y_array"]=adata.obs["x3"]
#adata.obs["x_pixel"]=adata.obs["x4"]
#adata.obs["y_pixel"]=adata.obs["x5"]
adata.obs["x_pixel"]=adata.obs["x2"] #这两列才是pixel
adata.obs["y_pixel"]=adata.obs["x3"]
#Select captured samples
adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.write_h5ad(f'{out_path_down}/sample_data.h5ad')
#Read in hitology image
#img=cv2.imread(f'{indir}/{samid}/spatial/detected_tissue_image.jpg')
img=cv2.imread(f'{indir}/tissue_lowres_image.png')
#img=cv2.imread(f'{indir}/{samid}/spatial/detected_tissue_image.jpg')
#img=cv2.imread(f'{indir}/{samid}/{samid}_full_image.tif')

### 4. Integrate gene expression and histology into a Graph

#Set coordinates
x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

#Test coordinates on the image
img_new=img.copy()
for i in range(len(x_pixel)):
    x=x_pixel[i]
    y=y_pixel[i]
    img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0

cv2.imwrite(f'{out_path_down}/{samid}_map.jpg', img_new)


#- The ‘s’ parameter determines the weight given to histology when calculating Euclidean distance between every two spots. ‘s = 1’ means that the histology pixel intensity value has the same scale variance as the (x,y) coordinates, whereas higher value of ‘s’ indicates higher scale variance, hence, higher weight to histology, when calculating the Euclidean distance. 

#- The "b"parameter determines the area of each spot when extracting color intensity.

 #Calculate adjacent matrix, 这里adj是所有细胞两两之间边的权重
s=1
b=49
adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
#If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
#adj=calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
np.savetxt(f'{out_path_down}/adj.csv', adj, delimiter=',')

#### 5. Spatial domain detection using SpaGCN
##### 5.1 Expression data preprocessing
#adata=sc.read("./data/151673/sample_data.h5ad")
#adj=np.loadtxt('./data/151673/adj.csv', delimiter=',')

adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

#### 5.2 Set hyper-parameters
#- p: Percentage of total expression contributed by neighborhoods.
#- l: Parameter to control p.

p=0.5 
#Find the l value given p
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

#- n_clusters: Number of spatial domains wanted.
#- res: Resolution in the initial Louvain's Clustering methods. If the number of clusters is known, we can use the spg.search_res() fnction to search for suitable resolution(optional)


#If the number of clusters known, we can use the spg.search_res() fnction to search for suitable resolution(optional)
#For this toy data, we set the number of clusters=7 since this tissue has 7 layers
n_clusters=10
#Set seed
r_seed=t_seed=n_seed=100
#Seaech for suitable resolution
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

#### 5.3 Run SpaGCN
clf=spg.SpaGCN() #建立一个空的spacgn对象
clf.set_l(l) #设定l，即超参数
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
#Save results
adata.write_h5ad(f'{out_path_down}/results.h5ad')
adata.obs.to_csv(f'{out_path_down}/SpaGCN.result')
