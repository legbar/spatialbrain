#!/usr/bin/env python
# coding: utf-8

# # Stage 2

# DA neuron identification

# In[39]:


import pandas as pd
import numpy as np
import scanpy as sc
import os
import glob
import re
from joblib import Parallel, delayed
from functools import partial
import scipy.sparse
import seaborn as sns
import scanpy_gpu_funcs as rsf
import cudf
import cupy as cp
from cuml.decomposition import PCA
from scipy.sparse import issparse
from SCTransform import SCTransform
from tqdm import tqdm
import pickle
# import bbknn

from sklearn.neighbors import LocalOutlierFactor

# import matplotlib.pyplot as plt
# from matplotlib import rcParams
# sc.set_figure_params(dpi= 100, dpi_save = 300)
# rcParams['figure.figsize'] = 5,5

# import rmm
# rmm.reinitialize(
#     managed_memory=True, # Allows oversubscription
#     pool_allocator=False, # default is False
#     devices=0, # GPU device IDs to register. By default registers only GPU 0.
# )
# cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

os.chdir('/active/paper/')


# # Logreg counts function

# In[40]:


def log_counts_for_logreg(adata):
    # log transform the counts
    adata.X = adata.layers['counts'].copy()
    if 'log1p' in adata.uns.keys(): 
        del adata.uns['log1p']
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # create lognorm_counts
    lognorm_counts = pd.DataFrame(adata.X.A, index=adata.obs_names, columns=adata.var_names)
    # restore original counts
    adata.X = adata.layers['counts'].copy()
    return(lognorm_counts)


# # Log reg

# In[41]:


from sklearn import metrics
from cuml.linear_model import LogisticRegression
# from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from matplotlib import colors
from matplotlib import cm
from sklearn.metrics import precision_recall_curve    

def logistic_model(data, cell_types, sparsity=0.2, fraction=0.5):
    X = data
    X_train, X_test, y_train, y_test = train_test_split(X, cell_types, test_size=fraction, random_state=1)
    lr = LogisticRegression(penalty='l1', C=sparsity, solver='qn', max_iter=10000)
    lr.fit(X_train, y_train)
    y_prob = lr.predict_proba(X_test)
    lr.coef_ = lr.coef_.transpose()
    lr_res = pd.DataFrame.from_records(lr.coef_, columns=X.columns)

    return(y_prob, y_test, lr_res, lr)


# # SCT and cluster function

# In[42]:


def sct_cluster(input_adata, res = 0.1, n_HVG = 5000):
    
    adata = input_adata.copy()
    
    # find var genes
    # sc.pp.highly_variable_genes(adata, 
    #                             flavor='seurat_v3', 
    #                             n_top_genes=n_HVG, 
    #                             batch_key='sample_id', 
    #                             span=1)
    
    # SCTransform

    # subset to HVG
    # adata = adata[:,adata.var['highly_variable']].copy()
    # run SCT
    SCTransform(adata, 
            min_cells=100, 
            gmean_eps=1, 
            n_genes=None, 
            n_cells=None, # use all cells
            bin_size=100000, 
            bw_adjust=3, 
            inplace=True)
        # store SCT layer
        # adata.layers['SCT'] = adata.X.copy()
    # else:
    #     adata.X = adata.layers['SCT'].copy()

    # delete any 'leiden_' columns
    adata.obs = adata.obs[adata.obs.columns.drop(list(adata.obs.filter(regex='^leiden_')))]

    sc.pp.pca(adata, random_state=1, use_highly_variable=False) # they are all variable
    sc.pp.neighbors(adata, method='rapids', random_state=1, n_neighbors=100)
    sc.tl.umap(adata, method='rapids', random_state=1)
    rsf.leiden(adata, resolution = res)
    
    n_clusters = len(adata.obs['leiden'].unique())
    
    # keep increasing the resolution until more than 1 cluster identified
    while n_clusters == 1:
        res = res + 0.1
        rsf.leiden(adata, resolution = res)
        n_clusters = len(adata.obs['leiden'].unique())
        
    leiden = adata.obs['leiden'].astype('str')
    
    # restore raw counts
    # adata.X = adata.layers['counts'].copy()
    
    return(leiden, res)


# # Decide cluster outcome function

# In[43]:


def decide_cluster_outcome(input_adata, leiden, parent_cluster):
    y_prob, y_test, lr_res, lr = logistic_model(log_counts_for_logreg(input_adata), 
                                                leiden.astype('int'))

    f1_max = {}
    for i, cell_type in enumerate(lr.classes_):
        
        if (y_test == lr.classes_[i]).sum() < 200:
            f1_max[cell_type] = 0
        else:
            precision, recall, thresholds = precision_recall_curve(y_test == cell_type, y_prob[:, i])
            f1_scores = [metrics.fbeta_score(y_test == cell_type, y_prob[:, i] > t, beta=1, zero_division=0) for t in np.random.choice(thresholds, 100)]
            f1_max[cell_type] = f1_scores[np.argmax(f1_scores)]

    print(f1_max)
    
    leiden = pd.DataFrame(leiden)

    for key, f1 in f1_max.items():
        if f1 >= 0.8:
            # if leiden.loc[leiden['leiden'] == str(key), 'leiden'].shape[0] < 200:
                # leiden.loc[leiden['leiden'] == str(key), 'leiden'] = parent_cluster + '_' + leiden.loc[leiden['leiden'] == str(key), 'leiden'] + '_FINAL_SMALL'
                # print(str(key) + ' finished (small)')
            # else:
            leiden.loc[leiden['leiden'] == str(key), 'leiden'] = parent_cluster + '_' + leiden.loc[leiden['leiden'] == str(key), 'leiden']
            print(str(key) + ' continuing')
        else:
            # may need to edit this to include last cluster ID before '_FINAL', as if there is more than one '_FINAL' cluster in a given round, they end up grouped together.
            leiden.loc[leiden['leiden'] == str(key), 'leiden'] = parent_cluster + '_FINAL'
            print(str(key) + ' finished')

    return(leiden)


# In[44]:


# # if first round, assign cluster to parent column
#     if first_round == False:
#         adata.obs['parent'] = adata.obs['leiden'].copy()


# # Load full adata

# In[45]:


with open('input/adata/midbrain/adata_spatial.pickle', 'rb') as f:
    adata = pickle.load(f)
    
# make float 64 for reproducibility
# https://github.com/theislab/scanpy/issues/313
#adata.layers['counts'] = adata.layers['counts'].astype('float64')
#adata.X = adata.X.astype('float64')


# In[46]:


import random

def seed_everything(seed=42):
    """"
    Seed everything.
    """   
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    
seed_everything(12)


# # SCT decision

# In[47]:



#adata = adata[adata.obs_names.isin(random.sample(adata.obs_names.to_list(), 1000))].copy()
sc.pp.filter_genes(adata, min_cells=50)
adata


# In[48]:


leiden, res = sct_cluster(adata, n_HVG=adata.shape[1])
leiden_df = decide_cluster_outcome(adata, leiden, 's')
adata.obs['current_leiden'] = leiden_df['leiden'].copy()


# In[49]:


adata.obs.groupby('current_leiden').size()


# In[50]:



while sum(['_FINAL' not in x for x in adata.obs['current_leiden'].unique()]):
    
    adatas = {}
    for cl in adata.obs['current_leiden'].unique():
        print(cl)
        if '_FINAL' not in cl:
            # try:
            adatas[cl] = adata[adata.obs['current_leiden'] == cl].copy()
            sc.pp.filter_genes(adatas[cl], min_cells=50)
            print(adatas[cl].shape)
            leiden, res = sct_cluster(adatas[cl], n_HVG=adata.shape[1])
            leiden_df = decide_cluster_outcome(adatas[cl], leiden, cl)
            adata.obs.loc[adata.obs['current_leiden'] == cl, 'current_leiden'] = leiden_df['leiden']
            with open('input/adata/midbrain/adata.pickle', 'wb') as f:
                pickle.dump(adata, f)
            # except:
                # pass
    


# In[51]:


adata.obs.groupby(['current_leiden']).size()


# In[52]:


# with open('input/adata/midbrain/adata.pickle', 'wb') as f:
#     pickle.dump(adata, f)


# In[ ]:




