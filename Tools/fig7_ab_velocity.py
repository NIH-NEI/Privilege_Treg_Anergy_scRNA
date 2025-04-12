#Install the following
#python=3.9 
#scvelo=0.3.3
#scanpy=1.10.3
#scvelo=0.3.3
#numpy=1.23.5

# Import packages
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad

# Load the velocity results data (provided in the same folder that contains this script)
adata = sc.read_h5ad('my_data_lantenteye.h5ad')

# Plot velocity stream (outputs to 'figures' folder in current working directory)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='RNA_snn_res.0.2', save='embedding_stream0731.png', title='',dpi = 300,legend_loc='none',palette={'0':'#bfa708','2':'#48c785','3':'#33a7f5','4':'#d272ed'})

# Plot lantent time on the umap (outputs to 'figures' folder in current working directory)
scv.pl.scatter(adata, color='latent_time', color_map='viridis', size=40, save='latenteye.png')
