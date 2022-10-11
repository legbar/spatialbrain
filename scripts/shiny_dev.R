# setup ----

library(ggsci)
library(cowplot)
library(parallel)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(MAST)
library(tidyverse)

setwd("/active/paper_shiny/")

anno <- readRDS("input/setup/anno.rds")

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
names(files) <- str_extract(files, "(?<=spatial/)[:graph:]+(?=.csv)")

mclapply(names(files), function(n) {
  # print(n)
  read_csv(files[[n]]) %>%
    filter(pvals_adj < 0.01) %>%
    filter(logfoldchanges > 0) %>%
     mutate(score = logfoldchanges * -log10(pvals_adj)) %>%
     arrange(desc(score)) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
    select("Gene" = names,
           "LFC" = logfoldchanges,
           "FDR-P" = pvals_adj) %>%
    saveRDS(paste0("input/markers/cell_types/", n, ".rds"))
}, mc.cores = parallel::detectCores()-1)
 
# readRDS("input/markers/counts_per_gene/Th.rds") %>% 
#   full_join(metadata_all_cells) %>%
#   ggplot(aes(x = cell_type_public, 
#              y = count)) +
#   geom_violin(scale = "width") +
#   theme_cowplot() +
#   coord_flip() +
#   labs(y = "Count") +
#   theme(axis.title.y = element_blank())
 
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

# TRAP splicing ----

# load the drimseq processed object: MB FRACTION (Cohort 1)
d_MB_FRACTION <- readRDS("input/setup/d_MB_FRACTION.rds")
# plot the number of transcripts per gene
counts_MB_FRACTION_DTU <- readRDS("input/setup/counts_MB_FRACTION_DTU.rds")

single_isoform_genes <- substr(names(table(counts_MB_FRACTION_DTU$gene_id)[table(counts_MB_FRACTION_DTU$gene_id) == 1]), 1, 18)
# filter zeros function
filter_zeros <- function(dds) {
  dds <- dds[rowSums(counts(dds)) > 0,]
}
filter_genes <- function(dds,
                         grouping,
                         prop = 3/5,
                         min = 0) {
  counts(dds) %>%
    as_tibble(rownames = "gene_id") %>%
    pivot_longer(-gene_id,
                 names_to = "sample_name",
                 values_to = "count") %>%
    inner_join(colData(dds), copy = TRUE) %>%
    mutate(.after = count,
           detected = count > min) %>%
    group_by(across(all_of(!!grouping))) %>%
    mutate(.after = count,
           proportion = sum(detected) / length(gene_id)) %>%
    group_by(gene_id) %>%
    filter(proportion >= prop) %>%
    pull(gene_id) %>%
    unique
}
dds <- readRDS("input/setup/dds.rds")
dds_C1 <- dds[, colData(dds)$cohort == "C1"] %>% filter_zeros()
dds_C1_MB <- dds_C1[, colData(dds_C1)$compartment == "MB"] %>% filter_zeros()
dds_C1_MB_FILTER <- filter_genes(dds_C1_MB, grouping = c("fraction", "gene_id"))
single_isoform_genes <- single_isoform_genes[single_isoform_genes %in% dds_C1_MB_FILTER]

# get results
res_MB_FRACTION <- DRIMSeq::results(d_MB_FRACTION)
res.txp_MB_FRACTION <- DRIMSeq::results(d_MB_FRACTION, level = "feature")

# replace NA values with 1
no.na <- function(x){
  ifelse(is.na(x), 1, x)
}

res_MB_FRACTION$pvalue <- no.na(res_MB_FRACTION$pvalue)
res.txp_MB_FRACTION$pvalue <- no.na(res.txp_MB_FRACTION$pvalue)

# remove gene id version number
res_MB_FRACTION$gene_id_simple <- str_extract(res_MB_FRACTION$gene_id,
            "[:alnum:]+(?=\\.[:digit:]+)")

# load anno_tx
anno_tx <- readRDS("/active/paper/input/trap/anno_tx.rds")

res_MB_FRACTION %>%
  select(ensembl_gene_id = gene_id_simple,
         -gene_id,
         everything()) %>%
  left_join(anno) %>%
  arrange(adj_pvalue) %>%
  # filter(external_gene_name %in% homologs$mmusculus_homolog_associated_gene_name) %>% # to focus on GWAS candidate genes
  mutate(adj_pvalue = ifelse(pvalue == 1, 1, adj_pvalue)) %>% 
  select("Gene" = external_gene_name, 
         "FDR-P" = adj_pvalue) %>%
  mutate(across(c(`FDR-P`), ~ signif(.x, 3))) %>%
  saveRDS("input/startup/splicing_meta.rds")

genes <- counts(d_MB_FRACTION) %>%
  as_tibble %>% 
  mutate(ensembl_gene_id = str_extract(gene_id,
                                       "[:alnum:]+(?=\\.[:digit:]+)"),
         ensembl_transcript_id = str_extract(feature_id,
                                             "[:alnum:]+(?=\\.[:digit:]+)")) %>%
  left_join(anno, 
            by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
  select(ensembl_gene_id, 
         external_gene_name) %>%
  distinct

d <- counts(d_MB_FRACTION) %>%
  as_tibble %>%
  select(ensembl_gene_id = gene_id,
         ensembl_transcript_id = feature_id,
         everything()) %>%
  # filter(ensembl_gene_id == g) %>%
  mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                       "[:alnum:]+(?=\\.[:digit:]+)"),
         ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                             "[:alnum:]+(?=\\.[:digit:]+)")) %>%
  left_join(anno_tx) %>%
  # mutate(feature_id = paste0("Transcript ", seq(1, nrow(.), 1))) %>%
  pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
               names_to = "sample_name",
               values_to = "count") %>%
  group_by(ensembl_gene_id, sample_name) %>%
  filter(median(count) > 0) %>%
  mutate(count = count / sum(count)) %>%
  inner_join(colData(dds), copy = TRUE) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
  left_join(anno)

mclapply(genes$ensembl_gene_id, function(g){
  m <- d %>%
    filter(ensembl_gene_id == g) %>%
    ungroup %>%
    select(external_gene_name, external_transcript_name, 
           sample_name, count) %>%
    mutate(fraction = ifelse(str_detect(sample_name, "TOTAL"), "TOTAL", "TRAP")) %>%
    select("Gene" = external_gene_name, 
           "Transcript" = external_transcript_name, 
           "Fraction" = fraction,
           "Count" = count)
  
  t <- m$external_gene_name %>% unique
  
  saveRDS(m, paste0("input/splicing/", t, ".rds"))
}, mc.cores = parallel::detectCores()-1)



# ageing ----
MB_AGE_META <-
  readRDS("input/setup/MB_AGE_META.rds")

MB_AGE_META %>%
  arrange(padj) %>%
  mutate(across(where(is.numeric), signif, 3)) %>%
  select("Gene" = external_gene_name, 
         "LFC" = log2FoldChange, 
         "FDR-P" = padj) %>%
  saveRDS("input/startup/MB_AGE_META.rds")

# # ONT object 
# # Load ONT Data & QC ----
# 
# # Load counts
# ont_tx_counts <-
#   read_tsv(
#     "input/setup/merged_counts.tsv") %>%
#   separate(
#     Reference,
#     into = c(
#       "ensembl_transcript_id",
#       "ensembl_gene_id",
#       NA,
#       NA,
#       NA,
#       NA,
#       NA,
#       NA
#     ),
#     sep = "([|])"
#   )
# 
# # The reference column has 8 components
# # all(unlist(lapply(strsplit(ont_tx_counts$Reference, "|", fixed = T), length)) == 8)
# # strsplit(ont_tx_counts$Reference[1], "|", fixed = T)
# ont_metadata <- tibble(flowcell_barcode = colnames(ont_tx_counts)[-c(1, 2)],
#                        barcode = str_extract(flowcell_barcode, "(?<=_)barcode[:digit:]{2}"),
#                        flowcell = str_extract(flowcell_barcode, "[:alnum:]+(?=_)"))
# 
# ont_metadata0 <- read_csv("input/setup/trap_2020_ont_metadata_mb.csv") %>%
#   select("sample_name0" = sample_id,
#          barcode)
# 
# ont_metadata <- ont_metadata %>%
#   left_join(ont_metadata0,
#             by  = "barcode") %>%
#   left_join(
#     colData(dds),
#     copy = T
#   ) %>%
#   select(sample_name,
#          barcode,
#          flowcell,
#          flowcell_barcode,
#          cohort,
#          region,
#          age,
#          genotype,
#          sex,
#          collection,
#          fraction,
#          compartment,
#          sample_code)
# 
# # collapse transcripts
# gene_id <- as.factor(ont_tx_counts$ensembl_gene_id)
# length(gene_id) == nrow(ont_tx_counts)
# sp <- split(seq(along = gene_id), gene_id)
# countData <- as.matrix(ont_tx_counts[,-c(1,2)])
# rownames(countData) <- ont_tx_counts$ensembl_transcript_id
# countData <- sapply(sp, function(i) colSums(countData[i,, drop = FALSE]))
# dim(countData)
# countData <- t(countData)
# 
# dds_ONT <- DESeqDataSetFromMatrix(countData = countData,
#                                   colData = ont_metadata,
#                                   design = ~1)
# # remove flongle data
# dds_ONT <- dds_ONT[,colData(dds_ONT)$flowcell != "AEL855"]
# 
# # collapse flowcell replicates
# dds_ONT <- collapseReplicates(dds_ONT,
#                               groupby = dds_ONT$sample_name,
#                               run = dds_ONT$flowcell)
# # remove version number from dds_ONT rownames
# rownames(dds_ONT) <- substr(rownames(dds_ONT), 1, 18)
# 
# counts(dds_ONT) %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(ensembl_gene_id %in% MB_AGE_META$ensembl_gene_id) %>%
#   pivot_longer(-ensembl_gene_id, 
#                names_to = "sample_name", 
#                values_to = "count") %>%
#   left_join(anno)

dds_C1_MB <- DESeq(dds_C1_MB, minReplicatesForReplace = Inf)

counts <- counts(dds_C1_MB, normalized = T) %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  filter(ensembl_gene_id %in% MB_AGE_META$ensembl_gene_id) %>%
  pivot_longer(-ensembl_gene_id,
               names_to = "sample_name",
               values_to = "count") %>%
  left_join(anno) %>%
  mutate(fraction = ifelse(str_detect(sample_name, "TOTAL"), 
                           "TOTAL", "TRAP")) %>%
  mutate(count = signif(count, 3)) %>%
  mutate(age = ifelse(str_detect(sample_name, "OLD"), "OLD", "YOUNG")) %>%
  select("Gene" = external_gene_name, 
         "Age" = age,
         "Fraction" = fraction,
         "Count" = count)

dir.create("input/ageing", showWarnings = F)
mclapply(MB_AGE_META$external_gene_name, function(g){
  counts[counts$Gene == g,] %>%
    select(-Gene) %>%
    saveRDS(paste0("input/ageing/", g, ".rds"))
}, mc.cores = parallel::detectCores()-1)

  

# xy for ageing cell type numbers ----

all_cells_metadata <- read_csv("/active/paper_shiny/input/setup/all_cells_metadata.csv")
all_cells_metadata %>%
  group_by(mouse_id, cell_type_publish, age) %>%
  summarise(n = n()) %>%
  ungroup %>%
  select(-mouse_id) %>%
  mutate(age = factor(age, levels = c("YOUNG", "OLD"))) %>%
  saveRDS("input/ageing_cell_type_numbers/actual_numbers.rds")

all_cells_metadata %>%
  group_by(mouse_id, cell_type_publish, age) %>%
  summarise(n = n()) %>%
  ungroup %>%
  select(-mouse_id) %>%
  mutate(age = factor(age, levels = c("YOUNG", "OLD"))) %>%
  left_join(read_csv("input/setup/cell_type_names.csv", col_names = c("cell_type_publish", "Cell Type"))) %>%
  select(`Cell Type`,
         "Age" = age, 
         "N" = n) %>%
  group_by(`Cell Type`, Age) %>%
  mutate("Mouse ID" = paste0(Age, ": ", row_number())) %>% 
  arrange(`Cell Type`) %>%
  ungroup %>%
  saveRDS("input/ageing_cell_type_numbers/actual_numbers_with_ids.rds")

actual_numbers <- readRDS("input/ageing_cell_type_numbers/actual_numbers.rds")

actual_numbers %>%
  filter(cell_type_publish == "AGE_ASTRO") %>%
  ggplot(aes(x = age, 
             y = n)) +
  geom_point() +
  theme_cowplot() +
  labs(x = "Age", y = "Number of Cells")

data <- readRDS("input/ageing_cell_type_numbers/data.rds") %>%
  arrange(fdr)

data %>%
  ggplot(aes(x = estimate, 
             y = -log10(fdr), 
             label = ifelse(fdr < 0.005, 
                            cell_type_full, 
                            ""), 
             color = ifelse(cell_type_full == "DA SN", 
                            "red", "black"))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.01), linetype = "dotted") +
  geom_label_repel(size = 1.5) +
  scale_x_continuous(limits = c(-1.5, 1.5), 
                     oob = scales::squish) +
  labs(x = expression(Log[2] ~ OR), 
       y = expression(-Log[10] ~ FDR)) +
  scale_color_manual(values = c("black", "red")) +
  theme_cowplot(6) +
  theme(legend.position = "none")


all_cells_metadata %>%
  group_by(age) %>%
  mutate(mouse_id_no = as.numeric(factor(mouse_id)),
         mouse_id_label = paste0(age, ": ", mouse_id_no)) %>%
  mutate(mouse_id_label = factor(
    mouse_id_label,
    levels = c("YOUNG: 1",
               "YOUNG: 2",
               "YOUNG: 3",
               "OLD: 1",
               "OLD: 2",
               "OLD: 3", 
               "OLD: 4")
  )) %>%
  ungroup %>%
  select(mouse_id_label, cell_type_publish, x, y) %>%
  saveRDS("input/ageing_cell_type_numbers/xy.rds")
