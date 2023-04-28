import pandas as pd
import numpy as np
import scanpy as sc
import os
import glob
import re
from joblib import Parallel, delayed
from functools import partial
import scipy.sparse

os.chdir('/zfs/analysis/paper')

sample_metadata = pd.read_csv('input/metadata/sample_metadata_midbrain.txt', sep='\t', index_col=['sample'])

# # set names
# names = list(sample_metadata.index)

# # load gems
# gems = Parallel(n_jobs=-2)(delayed(partial(pd.read_csv, sep='\t', index_col=['x', 'y']))(x) for x in glob.glob('input/gems/midbrain/*'))
# names = [re.findall('(?<=midbrain\/)FP.*(?=.bin1)', x)[0] for x in glob.glob('input/gems/midbrain/*')]
# gems = dict(zip(names, gems))

# # load masks
# masks = Parallel(n_jobs=-2)(delayed(partial(pd.read_csv, index_col=['x', 'y']))(x) for x in glob.glob('input/masks_dfs/midbrain/20211214/*'))
# names = [re.findall('(?<=20211214\/)FP.*(?=_mask)', x)[0] for x in glob.glob('input/masks_dfs/midbrain/20211214/*')]
# masks = dict(zip(names, masks))

# # join gems and masks: do not parallelise
# masked_gems = {}
# for name in names:
#     masked_gems[name] = gems[name].join(masks[name], how='inner').to_csv('input/masked_gems/'+name+'.csv.gz')

# load masked gems
print('Loading masked gems...')
masked_gems = Parallel(n_jobs=-2)(delayed(partial(pd.read_csv))(x) for x in glob.glob('input/masked_gems/*'))
names = [re.findall('(?<=masked_gems\/)FP.*(?=.csv.gz)', x)[0] for x in glob.glob('input/masked_gems/*')]
masked_gems = dict(zip(names, masked_gems))

# reset cell ID
def reset_cell_id(name):
    df = masked_gems[name]
    df['cell_ID'] = masked_gems[name].groupby('label').ngroup()
    df.drop('label', axis=1, inplace=True)
    return(df)
masked_gems = Parallel(n_jobs=-2)(delayed(reset_cell_id)(x) for x in names)
masked_gems = dict(zip(names, masked_gems))

# create xy df
print('Exporting xy information...')
for name in names:
    xy = masked_gems[name][['x', 'y', 'cell_ID']]
    xy.groupby('cell_ID')[['x', 'y']].mean().round().astype('int').to_csv('input/xy/'+name+'.csv.gz')

# produce long counts
def masked_gems_to_long_counts(name):
    counts = masked_gems[name].reset_index()[['geneID', 'MIDCounts', 'cell_ID']]
    # counts['geneID_int'] = counts['geneID'].rank(method='dense').astype(int)
    counts['geneID_int'] = counts.groupby('geneID').ngroup()
    counts.to_csv('input/counts_long/'+name+'.csv.gz')
print('Exporting long counts...')
Parallel(n_jobs=4)(delayed(masked_gems_to_long_counts)(x) for x in names)

# clean up
del masked_gems

# load long counts
counts = Parallel(n_jobs=-2)(delayed(pd.read_csv)(x) for x in glob.glob('input/counts_long/*'))
names = [re.findall('(?<=counts_long\/)FP.*(?=.csv.gz)', x)[0] for x in glob.glob('input/counts_long/*')]
counts = dict(zip(names, counts))

# produce gene_id dfs
print('Exporting gene ids...')
def create_gene_ids(name):
    counts[name][['geneID', 'geneID_int']].set_index('geneID_int').sort_index().drop_duplicates().to_csv('input/gene_ids/'+name+'.csv.gz')
Parallel(n_jobs=-2)(delayed(create_gene_ids)(x) for x in names)

# create sparse matrices
print('Exporting sparse matrices...')
for name in names:
    sparse = scipy.sparse.csr_matrix((counts[name].MIDCounts, (counts[name].cell_ID, counts[name].geneID_int)))
    scipy.sparse.save_npz('input/sparses/'+name+'.npz', sparse)