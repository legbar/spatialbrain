# setup ----

library(tidyverse)
library(ggsci)
library(cowplot)
library(parallel)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(MAST)

setwd("/active/paper_shiny/")

# Interquartile range outlier detection function (IROF) ----
outlier_iqr <- function(column, low_high = "both") {
  if (low_high == "both") {
    ifelse(
      column < quantile(column, 0.25) - 1.5 * IQR(column) |
        column > quantile(column, 0.75) + 1.5 * IQR(column),
      TRUE,
      FALSE
    )
  } else if (low_high == "low") {
    ifelse(column < quantile(column, 0.25) - 1.5 * IQR(column),
           TRUE,
           FALSE)
  } else if (low_high == "high") {
    ifelse(column > quantile(column, 0.75) + 1.5 * IQR(column),
           TRUE,
           FALSE)
  }
}

# spatial metadata ----
# load all cells spatial metadata
metadata_all_cells <- read_csv("input/setup/all_cells_metadata.csv")
colnames(metadata_all_cells)[1] <- "cell_id"

metadata_all_cells <- metadata_all_cells %>%
  left_join({
    tibble(
      mouse_id = unique(metadata_all_cells$mouse_id),
      mouse_id_presentation_no_genotype = c("Young 1",
                                            "Young 2",
                                            "Young 3",
                                            "Old 1",
                                            "Old 2",
                                            "Old 3",
                                            "Old 4")
    )
  })

# add public cell type names
# files <- list.files("input/setup/markers_spatial",
#                     pattern = ".csv",
#                     full.names = T)
# names(files) <-
#   str_extract(files, "(?<=markers_spatial/)[:graph:]+(?=.csv)")
# tibble(orig = names(files)) %>%
#   write_csv("input/setup/cell_type_names.csv")
cell_type_names <- read_csv("input/setup/cell_type_names.csv", 
                            col_names = c("cell_type_internal", "cell_type_public"), 
                            skip = 1)

cell_type_names_list <- cell_type_names$cell_type_internal
names(cell_type_names_list) <- cell_type_names$cell_type_public
cell_type_names_list %>% saveRDS("input/startup/cell_type_names.rds")

metadata_all_cells <- metadata_all_cells %>%
  select(everything(), 
         "cell_type_internal" = cell_type_publish, 
         -cell_type) %>%
  left_join(cell_type_names)
metadata_all_cells %>%
  saveRDS("input/startup/metadata_all_cells.rds")

# extract DA metadata and recode mouse names
metadata_da <- metadata_all_cells %>%
  filter(str_detect(cell_type_internal, "^DA_")) %>%
  mutate(region = ifelse(str_detect(cell_type_internal, "SN"), "SN", "VTA")) %>%
  select(cell_id, mouse_id_presentation_no_genotype, cell_type_public, mouse_id, age, genotype, region, cell_type_internal, x, y)   

metadata_da %>% saveRDS("input/startup/da_metadata.rds")

# filter spatial outliers for better plot presentation
metadata_da %>%
  group_by(mouse_id) %>%
  filter(!outlier_iqr(x)) %>%
  filter(!outlier_iqr(y)) %>%
  write_csv("input/sn_vta/da_metadata_spatial_outliers_removed.csv")

# spatial markers ----
files <- list.files("input/setup/markers_spatial",
                    pattern = ".csv",
                    full.names = T)
names(files) <-
  str_extract(files, "(?<=markers_spatial/)[:graph:]+(?=.csv)")
tibble(orig = names(files)) %>%
  write_csv("input/setup/cell_type_names.csv")

 mclapply(names(files), function(n) {
  read_csv(files[[n]]) %>%
    select("gene" = names,
           "lfc" = logfoldchanges,
           "padj" = pvals_adj) %>%
    filter(padj < 0.01) %>%
    saveRDS(paste0("input/markers/cell_types/", n, ".rds"))
}, mc.cores = detectCores()-2)

# create counts per gene for all cells
# load matrix for all cells
m <- readMM("input/setup/adata.mtx")
# m <- as(m, "dgCMatrix")
dimnames(m) <-
  list(metadata_all_cells$cell_id,
       read_csv("input/setup/var_names.csv")$geneID)
sce <- SingleCellExperiment(list(counts = t(m)),
                            colData = metadata_all_cells)
mat <- t(as.matrix(counts(sce)))
rm(sce)

mclapply(colnames(mat), function(gene){
  enframe(mat[,gene], name = "cell_id", value = "count") %>%
    # left_join(metadata_all_cells) %>% 
    saveRDS(paste0("input/markers/counts_per_gene/", gene, ".rds"))
}, mc.cores = detectCores()-2)

all(readRDS(paste0("input/markers/counts_per_gene/", selected_gene, ".rds"))$cell_id == metadata_all_cells$cell_id)

cell_types <- metadata_all_cells$cell_type_internal
names(cell_types) <- metadata_all_cells$cell_type_public

selected_gene <- "Tagln3"
# selected_cell_type <- "DA_VTA"
counts <- readRDS(paste0("/active/paper_shiny/input/markers/counts_per_gene/", selected_gene, ".rds")) %>%
  cbind(metadata_all_cells[,-1])

rbind(
  slice_sample(counts[counts$count == 0,], prop = 0.1),
  counts[counts$count > 0,]
) %>%
  # filter(count > 0) %>%
  arrange(count) %>%
  # mutate(alpha = ifelse(count == 0, 0.01, 1)) %>%

  # slice_sample(prop = 0.1) %>%
  ggplot(aes(x = x, 
             y = y, 
             colour = count, 
             alpha = count,
             size = count)) +
  geom_point() +
  # geom_point(data = counts[counts$cell_type_publish == selected_cell_type,], 
  #            colour = "orange") +
  facet_wrap(vars(mouse_id_presentation_no_genotype), 
             scales = "free") +
  theme_cowplot() +
  panel_border() +
  scale_y_reverse() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(), 
    # legend.position = "none"
  ) +
  scale_colour_viridis_c(option = "A") +
  # scale_color_gradientn(colours = c("white", "white", "red"), 
  #                       breaks = c(0, 3, 12),
  #                       limits = c(0, 12), 
  #                       oob = scales::squish) +
  scale_size(range = c(0.01, 1.5)) +
  scale_alpha(range = c(0, 1))

# sn_vta ----

# load da metadata
metadata_da <- read_csv("input/setup/da_metadata.csv")

# create da sce object
# load matrix for all cells
m <- readMM("input/setup/adata.mtx")
# m <- as(m, "dgCMatrix")
dimnames(m) <-
  list(metadata_all_cells$cell_id,
       read_csv("input/setup/var_names.csv")$geneID)
sce <- SingleCellExperiment(list(counts = t(m)),
                            colData = metadata_all_cells)
saveRDS(sce, "input/setup/sce.rds")
sce_da <- sce[, str_detect(colData(sce)$cell_type_publish, "^DA_")]
sce_da <- logNormCounts(sce_da)
colData(sce_da)$cell_type_publish <-
  factor(colData(sce_da)$cell_type_publish,
         levels = c("DA_VTA", "DA_SN"))
colData(sce_da)$age <-
  factor(colData(sce_da)$age, levels = c("YOUNG", "OLD"))
colData(sce_da)$genotype <-
  factor(colData(sce_da)$genotype, levels = c("WT", "OVX"))
sca_da <- SceToSingleCellAssay(sce_da)
sca_da <-
  sca_da[freq(sca_da) >= 0.1,] # only assess genes detected in > 0.1 of cells
# cell detection rate
cdr <- colSums(assay(sca_da) > 0)
colData(sca_da)$cdr <- scale(cdr)
saveRDS(sce_da, "input/setup/sce_da.rds")
saveRDS(sca_da, "input/setup/sca_da.rds")

# save counts per gene
dir.create("input/sn_vta/da_counts_per_gene", recursive = T)
sca_da_counts <-
  as_tibble(t(as.matrix(counts(sca_da))), rownames = "cell_id")
mclapply(colnames(sca_da_counts)[-1], function(gene) {
  c <- sca_da_counts[, c("cell_id", gene)]
  colnames(c)[2] <- "count"
  write_csv(c, paste0("input/sn_vta/da_counts_per_gene/", gene, ".csv"))
}, mc.cores = detectCores() - 2)

# sn vs vta MAST
z <- zlm(
  ~ cdr + cell_type_publish + (1 | mouse_id),
  sca = sca_da,
  method = "glmer",
  ebayes = F,
  strictConvergence = F,
  exprs_values = 'logcounts'
)
summary(z)
# LRT
summaryCond <- summary(z, doLRT = 'cell_type_publishDA_SN')
summaryDt <- summaryCond$datatable
fcHurdle <-
  merge(summaryDt[contrast == 'cell_type_publishDA_SN' &
                    component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryDt[contrast == 'cell_type_publishDA_SN' &
                    component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
# FCTHRESHOLD <- log2(1.1)
# [fdr<.05 & abs(coef)>FCTHRESHOLD]
fcHurdle <-
  merge(fcHurdle, as_tibble(mcols(sca_da)), by = 'primerid') %>%
  as_tibble %>%
  arrange(fdr)
# write out in simple format
fcHurdle %>%
  select("gene" = primerid,
         "lfc" = coef,
         fdr) %>%
  # filter(fdr < 0.1) %>%
  write_csv("input/sn_vta/sn_vta_mast.csv")


# da markers MAST: MEMORY CONSTRAINTS ----
# 
# # load da metadata
# metadata_da <- read_csv("input/setup/da_metadata.csv")
# 
# # create da sce object
# # load matrix for all cells
# m <- readMM("input/setup/adata.mtx")
# # m <- as(m, "dgCMatrix")
# dimnames(m) <-
#   list(metadata_all_cells$cell_id,
#        read_csv("input/setup/var_names.csv")$geneID)
# sce <- SingleCellExperiment(list(counts = t(m)),
#                             colData = metadata_all_cells)
# saveRDS(sce, "input/setup/sce.rds")
# sce <- logNormCounts(sce)
# colData(sce)$DA <-
#   factor(ifelse(str_detect(colData(sce)$cell_type_internal, "DA_"), 
#                 "DA", "Other"), 
#          levels = c("Other", "DA"))
# sum(colData(sce)$DA == "DA")
# sca <- SceToSingleCellAssay(sce)
# sca <-
#   sca[freq(sca) >= 0.01,] # only assess genes detected in > 0.1 of cells
# # cell detection rate
# cdr <- colSums(assay(sca) > 0)
# colData(sca)$cdr <- scale(cdr)
# saveRDS(sce, "input/setup/sce.rds")
# saveRDS(sca, "input/setup/sca.rds")
# 
# # save counts per gene
# dir.create("input/sn_vta/da_counts_per_gene", recursive = T)
# sca_da_counts <-
#   as_tibble(t(as.matrix(counts(sca_da))), rownames = "cell_id")
# mclapply(colnames(sca_da_counts)[-1], function(gene) {
#   c <- sca_da_counts[, c("cell_id", gene)]
#   colnames(c)[2] <- "count"
#   write_csv(c, paste0("input/sn_vta/da_counts_per_gene/", gene, ".csv"))
# }, mc.cores = detectCores() - 2)
# 
# # sn vs vta MAST
# z <- zlm(
#   ~ cdr + cell_type_publish + (1 | mouse_id),
#   sca = sca_da,
#   method = "glmer",
#   ebayes = F,
#   strictConvergence = F,
#   exprs_values = 'logcounts'
# )
# summary(z)
# # LRT
# summaryCond <- summary(z, doLRT = 'cell_type_publishDA_SN')
# summaryDt <- summaryCond$datatable
# fcHurdle <-
#   merge(summaryDt[contrast == 'cell_type_publishDA_SN' &
#                     component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
#         summaryDt[contrast == 'cell_type_publishDA_SN' &
#                     component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
# fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
# # FCTHRESHOLD <- log2(1.1)
# # [fdr<.05 & abs(coef)>FCTHRESHOLD]
# fcHurdle <-
#   merge(fcHurdle, as_tibble(mcols(sca_da)), by = 'primerid') %>%
#   as_tibble %>%
#   arrange(fdr)
# # write out in simple format
# fcHurdle %>%
#   select("gene" = primerid,
#          "lfc" = coef,
#          fdr) %>%
#   # filter(fdr < 0.1) %>%
#   write_csv("input/sn_vta/sn_vta_mast.csv")



# TRAP ----
# TRAP enrichment ----

readRDS("input/trap/MB_FRACTION_META.rds") %>%
  select("Gene" = external_gene_name, 
         "Log2 Fold Enrichment" = log2FoldChange, 
         "FDR-P" = sumz_adj) %>%
  mutate(across(where(is.numeric), signif, 3)) %>%
  saveRDS("input/startup/MB_FRACTION_META.rds")

readRDS("input/trap/MB_FRACTION_META.rds") %>%
  mutate(enrichment = ifelse(sumz_adj > 0.01, "Unchanged",
                             ifelse(log2FoldChange > 0, "Enriched", "Depleted"))) %>%
  select(baseMean_C1, 
         log2FoldChange, 
         enrichment, 
         external_gene_name) %>%
  saveRDS("input/startup/MA_plot_data.rds")
