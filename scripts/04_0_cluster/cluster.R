
# renv::install("bioc::scry")
# renv::install("bioc::BiocNeighbors")
# renv::install("igrabski/sc-SHC")

library(scSHC)
library(reticulate)
library(tidyverse)
# library(MatrixExtra)

setwd("/f_active/paper_23/")

scipy <- import("scipy.sparse")
sparse_mat <- scipy$load_npz('input/03_norm_bc_dim_red/sparse.npz')
# sparse_mat <- scipy$load_npz('input/clustering/sparse.npz')
data(counts)
dim(counts)

sparse_mat <- MatrixExtra::t_deep(sparse_mat)
sparse_mat <- as(sparse_mat, "dgCMatrix")
dim(sparse_mat)
var <- read_csv("input/03_norm_bc_dim_red/var.csv")$var
obs <- read_csv("input/03_norm_bc_dim_red/obs.csv")$obs
batch <- read_csv("input/03_norm_bc_dim_red/batch.csv")$batch
rownames(sparse_mat) <- var
colnames(sparse_mat) <- obs
sub <- sparse_mat[,sample(seq(1, ncol(sparse_mat)), size = round(ncol(sparse_mat)/100))]
dim(sub)

clusters <- scSHC(sub, num_features = 4000, batch = batch, parallel = T, cores = 16)
clusters
