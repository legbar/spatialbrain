# import scanpy as sc
import anndata as ad
import rapids_singlecell as rsc
from rapids_singlecell.cunnData import cunnData

import cupy as cp
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=False, # Allows oversubscription
    pool_allocator=False, # default is False
    devices=0, # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from joblib import Parallel, delayed
from functools import partial
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix, issparse
import seaborn as sns

import sklearn.neighbors

import scvi
scvi.settings.seed = 1 

from joblib import Parallel, delayed

import os
os.chdir('/f_active/paper_23/') # laune

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

from sklearn import metrics
from cuml.linear_model import LogisticRegression
# from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from matplotlib import colors
from matplotlib import cm
from sklearn.metrics import precision_recall_curve    

import scanpy as sc

def logistic_model(data, cell_types, sparsity=0.2, fraction=0.5):
    X = data
    X_train, X_test, y_train, y_test = train_test_split(X, cell_types, test_size=fraction, random_state=1)
    lr = LogisticRegression(penalty='l1', C = sparsity, max_iter=10000)
    lr.fit(X_train, y_train)
    y_prob = lr.predict_proba(X_test)
    lr_coef = lr.coef_.transpose()
    # lr_res = pd.DataFrame.from_records(lr_coef, columns=X.columns)
    # return(y_prob)
    return(y_prob, y_test, lr)

def norm_cluster(input_adata, res = 0.1, n_HVG = 1000, prime = False):

    adata = input_adata
    adata.X = adata.layers['counts'].copy()
    adata_full = adata.copy()
    # print('filtering genes and cells')
    # sc.pp.filter_genes(adata_temp, min_cells=50)
    # sc.pp.filter_cells(adata, min_genes = 200)
    # sc.pp.filter_genes(adata, min_cells = 20)
    print('adata shape...')
    print(adata.shape)
    print('adata shape after cell filtering...')
    sc.pp.filter_genes(adata, min_cells=100)
    print(adata.shape)
    # print('normalizing...')
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # print('log1p...')
    # sc.pp.log1p(adata)
    # adata.raw = adata  # freeze the state in `.raw`
    # print(adata)
    
    print('filtering highly variable...')
    # try:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_HVG,
        subset=True,
        layer="counts",
        flavor="seurat_v3",
        batch_key="sample_name",
        span = 1
    )
    # except:
    #     print('!!! HVG ERROR: BATCH ISSUE LIKELY !!!')
    #     print(adata.obs.groupby('sample_name').size())
    #     sc.pp.highly_variable_genes(
    #         adata,
    #         n_top_genes=n_HVG,
    #         subset=True,
    #         layer="counts",
    #         flavor="seurat_v3"
    #     )

    print('adata shape after HVG selection...')
    print(adata.shape)

    print('number of samples...')
    print(len(adata.obs['sample_id'].unique()))

    min_batch_size = adata.obs.groupby('sample_id').size().min()
    print('smallest group size: ' + str(min_batch_size))
    # batch_size = 128
    # if min_batch_size < 200:
    #     batch_size = min_batch_size
    #     scvi.settings.batch_size = min_batch_size
    # while min_batch_size % batch_size == 1:
    #     batch_size += 1
    #     scvi.settings.batch_size += 1
    # print('batch size: ' + str(batch_size))
    print('starting SCVI...')
    # print(adata.obs.columns)
    scvi.model.SCVI.setup_anndata(adata, 
                                  layer="counts",
                                  categorical_covariate_keys=['batch'],
                                  continuous_covariate_keys=['n_genes_by_counts']
                                 )
    print('preparing SCVI model...')
    model_scvi = scvi.model.SCVI(adata)
    print('training model...')
    batch_size = 128
    for i in range(20):
    # while True:
        mult = i + 1
        try:
            bs = batch_size*mult
            print("Training with batch size: " + str(bs))
            model_scvi.train(early_stopping=True, batch_size = bs)
        except:
            print("NEED A LARGER BATCH SIZE")
        else:
            break
    adata.obsm["X_scVI"] = model_scvi.get_latent_representation()

    print('nearest neighbours...')
    rsc.tl.neighbors(adata, use_rep="X_scVI")
    print('leiden res ' + str(res) + '...')
    sc.tl.leiden(adata, resolution = res)
    n_clusters = len(adata.obs['leiden'].unique())
    print('N clusters: ' + str(n_clusters))
    
    # keep increasing the resolution until more than 1 cluster identified
    while n_clusters == 1:
        res = res + 0.1
        print('leiden res ' + str(res) + '...')
        sc.tl.leiden(adata, resolution = res)
        n_clusters = len(adata.obs['leiden'].unique())
        print('N clusters: ' + str(n_clusters))

    leiden = adata.obs['leiden'].astype('str')
    
    if prime:
        leiden_df = decide_cluster_outcome(adata, leiden, 's')
        adata_full.obs['current_leiden'] = leiden_df['leiden'].astype('str').copy()
        # adata_full.write_h5ad('input/03_norm/adata/adata_clustered_231007_full_primed.h5ad')
        return(adata_full)

    else:
        return(leiden)

def decide_cluster_outcome(input_adata, leiden, parent_cluster):
    y_prob, y_test, lr = logistic_model(log_counts_for_logreg(input_adata), 
                                                leiden.astype('int'))

    f1_max = {}
    for i, cell_type in enumerate(lr.classes_):
        
        if (y_test == lr.classes_[i]).sum() < MIN_CLUSTER_SIZE:
            f1_max[cell_type] = 0
        else:
            precision, recall, thresholds = precision_recall_curve(y_test == cell_type, y_prob.loc[:, i])
            f1_scores = [metrics.fbeta_score(y_test == cell_type, y_prob.loc[:, i] > t, beta=1, zero_division=0) for t in np.random.choice(thresholds, 100)]
            f1_max[cell_type] = f1_scores[np.argmax(f1_scores)]

    print(f1_max)
    
    leiden = pd.DataFrame(leiden)

    for key, f1 in f1_max.items():
        if f1 >= 0.75:
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

def meat_and_potatoes(adata, cl, n_HVG):
    # if '_FINAL' not in cl:
    print('making adata_temp...')
    adata_temp = adata[adata.obs['current_leiden'] == cl].copy()
    leiden = norm_cluster(adata_temp, n_HVG = n_HVG)
    print('deciding cluster outcome...')
    leiden_df = decide_cluster_outcome(adata_temp, leiden, cl)
    print('Updating leiden labeling...')
    leiden_df['leiden'] = leiden_df['leiden'].astype('str')
    leiden_d = leiden_df.loc[leiden_df['leiden'].str.contains(cl), 'leiden'].to_dict()
    # print(leiden_d)
    # leiden_d = leiden_df['leiden'].to_dict()
    # current_leiden_d = adata.obs['current_leiden'].to_dict()
    # print(current_leiden_d)
    # for k in current_leiden_d.keys():
    #     if k in leiden_d.keys():
    #         current_leiden_d[k] = leiden_d[k]
    # df = pd.DataFrame.from_dict(current_leiden_d, orient = 'index', columns = ['current_leiden'])
    return(leiden_d)

# def update_adata_obs(adata, cl, d):

#     adata.obs['current_leiden'] = adata.obs['current_leiden'].astype('str')
#     adata.obs['current_leiden'].update(d)
    
#     # adata.obs = adata.obs.drop('current_leiden', axis = 1).join(df)
#     # adata.obs['current_leiden'] = adata.obs['current_leiden'].astype('str')
   
#     # print("dropping singletons...")
#     # print(adata.shape)
#     # adata = adata[adata.obs['current_leiden'].str.startswith('s_').astype('bool'), :].copy() # remove singleton clusters
#     # print(adata.shape)
    
#     if sum(adata.obs['current_leiden'].isna()) > 0:
#         print("removing nan clusters...")
#         adata = adata[~adata.obs['current_leiden'].isna(), :].copy() # remove nan clusters
#     return(adata)
    
n_HVG = 1000
MIN_CLUSTER_SIZE = 200

adata = sc.read('input/02_stereoseq_qc/adata/adata.h5ad')
# adata = sc.read('input/03_norm/adata/adata_clustered_231003_full_26.h5ad')
# adata = sc.read('input/03_norm/adata/adata_clustered_231003_full_58restart.h5ad')

print(adata)

print('subsetting data for testing...')
# adata = adata[adata.obs['n_genes_by_counts'] > 500].copy()
# adata = adata[adata.obs['batch'] == 'batch1', :].copy()
# sc.pp.subsample(adata, n_obs=5000)
print('data shape before filtering genes...')
print(adata.shape)
sc.pp.filter_genes(adata, min_cells=50)
print('data shape after filtering genes...')
print(adata.shape)

print('storing raw counts...')
adata.layers["counts"] = adata.X.copy()  # preserve counts

if 'current_leiden' not in adata.obs.columns:
    adata = norm_cluster(adata, n_HVG = n_HVG, prime = True)
    # adata.write_h5ad('input/03_norm/adata/temp.h5ad')
    
while sum(['_FINAL' not in x for x in adata.obs['current_leiden'].unique()]):

    print(adata.obs.groupby('current_leiden').size())

    clusters_to_process = [x for x in list(adata.obs['current_leiden'].unique()) if '_FINAL' not in x]
    ds = Parallel(n_jobs=len(clusters_to_process))(delayed(meat_and_potatoes)(adata, cl, n_HVG) for cl in clusters_to_process)
    # d_ds = dict(zip(clusters_to_process, ds))
    for d in ds:
        adata.obs['current_leiden'] = adata.obs['current_leiden'].astype('str')
        adata.obs['current_leiden'].update(d)
        # adata = update_adata_obs(adata, i)

                   
    # for count, cl in enumerate(list(adata.obs['current_leiden'].unique())):
    #     if '_FINAL' not in cl:
    #         print('making adata_temp...')
    #         adata_temp = adata[adata.obs['current_leiden'] == cl].copy()
    #         leiden = norm_cluster(adata_temp, n_HVG = n_HVG)
    #         print('deciding cluster outcome...')
    #         leiden_df = decide_cluster_outcome(adata_temp, leiden, cl)
    #         print('Updating leiden labeling...')
    #         leiden_df['leiden'] = leiden_df['leiden'].astype('str')
    #         leiden_d = leiden_df['leiden'].to_dict()
    #         current_leiden_d = adata.obs['current_leiden'].to_dict()
    #         for k in current_leiden_d.keys():
    #             if k in leiden_d.keys():
    #                 current_leiden_d[k] = leiden_d[k]
    #         df = pd.DataFrame.from_dict(current_leiden_d, orient = 'index', columns = ['current_leiden'])
    #         adata.obs = adata.obs.drop('current_leiden', axis = 1).join(df)
    #         adata.obs['current_leiden'] = adata.obs['current_leiden'].astype('str')
    #         # print("dropping singletons...")
    #         # print(adata.shape)
    #         # adata = adata[adata.obs['current_leiden'].str.startswith('s_').astype('bool'), :].copy() # remove singleton clusters
    #         # print(adata.shape)
    #         if sum(adata.obs['current_leiden'].isna()) > 0:
    #             print("removing nan clusters...")
    #             adata = adata[~adata.obs['current_leiden'].isna(), :].copy() # remove nan clusters

    print('Saving checkpoint...')
    adata.write_h5ad('input/03_norm/adata/adata_clustered_full_HVG1000_groupSize200_f0.75_categoricalBatch_continuousNGenesByCounts.h5ad')

    # adata.write_h5ad('input/03_norm/adata/adata_clustered_full_' + str(count) + '_231009_HVG_groupSize200_f0.75_categoricalSampleName_continuousNone.h5ad')
