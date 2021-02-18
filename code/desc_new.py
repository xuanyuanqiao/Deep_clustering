import numpy as np
import pandas as pd
import scanpy as sc
import platform
import desc
from time import time                                                       
import sys
import matplotlib
import matplotlib.pyplot as plt
import tensorflow as tf
import keras 
sc.settings.set_figure_params(dpi=150)
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=200, facecolor='white')


adata = sc.read_10x_mtx(
    'DATA_DIR',
    var_names='gene_symbols',                
    cache=True)                         

adata.var_names_make_unique()
###filter
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
##remove cells with a high proportion of mitochondria genes expression.
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


#adata = adata[adata.obs['n_genes'] < 2500, :]
adata = adata[adata.obs['n_counts']<1500,:]
adata = adata[adata.obs['percent_mito'] < 0.5, :]
##normalization
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw=adata
##Selection of highly variable genes

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True)
adata = adata[:, adata.var['highly_variable']]

desc.scale(adata, zero_center=True, max_value=3)# if the the dataset has two or more batches you can use `adata=desc.scale(adata,groupby="BatchID")`
save_dir="h5_result"

adata=desc.train(adata,
        dims=[adata.shape[1],128,32],
        tol=0.005,
        n_neighbors=10,
        batch_size=256,
        louvain_resolution=[1.0],
        save_dir=str(save_dir),
        do_tsne=True,
        learning_rate=200, 
        use_GPU=False,
        num_Cores=2, 
        num_Cores_tsne=4,
        save_encoder_weights=False,
        save_encoder_step=3,
        use_ae_weights=False,
        do_umap=True,
        max_iter=100,
        pretrain_epochs=5) 
adata.obs['max.prob']=adata.uns["prob_matrix1.0"].max(1)
sc.pl.scatter(adata,basis="tsne1.0",color=['desc_1.0','max.prob'])
sc.pl.scatter(adata,basis="umap1.0",color=['desc_1.0','max.prob'])
