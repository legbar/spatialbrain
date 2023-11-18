library(tidyverse)
Sys.setenv("MC_CORES"=18L)
library(parallel)
options(mc.cores = 18)
library(cowplot)
library(ggsci)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(MAST)
# library(Libra)
library(Seurat)
library(ggh4x)
library(ggrepel)

spatial_all_cells <-
  read_csv("output/spatial_coordinates/all_cells.csv")

# Sox6 ----

"Agtr1a" %in% rownames(sce)

sca_split_no_filter <-
  mclapply(unique(colData(sce)$cell_type_publish), function(x) {
    # sce <- sce[, colData(sce)$cell_type_DA_grouped == x]
    sce <- sce[, colData(sce)$cell_type_publish == x]
    sce <- logNormCounts(sce)
    colData(sce)$age <-
      factor(colData(sce)$age, levels = c("YOUNG", "OLD"))
    colData(sce)$genotype <-
      factor(colData(sce)$genotype, levels = c("WT", "OVX"))
    sca <- SceToSingleCellAssay(sce)
    # sca <- sca[freq(sca) >= 0.1,] # only assess genes detected in > 0.1 of cells
  }, mc.cores = 18)
names(sca_split_no_filter) <- unique(colData(sce)$cell_type_publish)
gene_id_conversion <- MB_FRACTION_META %>%
  select(ensembl_gene_id, external_gene_name)
MB_FRACTION_TRANSLATED_GENES_external_gene_name <-
  gene_id_conversion$external_gene_name[match(MB_FRACTION_TRANSLATED_GENES,
                                              gene_id_conversion$ensembl_gene_id)]
# sce_DA_no_filter <- sce[rownames(sce) %in% MB_FRACTION_TRANSLATED_GENES_external_gene_name,
#               colData(sce)$cell_type_publish %in% c("DA_SN", "DA_VTA")]
sce_DA_no_filter <-
  sce[, colData(sce)$cell_type_publish %in% c("DA_SN", "DA_VTA")]
sce_DA_no_filter <- logNormCounts(sce_DA_no_filter)
colData(sce_DA_no_filter)$age <-
  factor(colData(sce_DA_no_filter)$age, levels = c("YOUNG", "OLD"))
colData(sce_DA_no_filter)$genotype <-
  factor(colData(sce_DA_no_filter)$genotype, levels = c("WT", "OVX"))
sca_DA_no_filter <- SceToSingleCellAssay(sce_DA_no_filter)

"Agtr1a" %in% rownames(sca_DA_no_filter)
enframe(assay(sca_DA_no_filter)[rownames(sca_DA_no_filter) == 'Sox6',],
        name = 'cell_id', value = 'count') %>%
  left_join(meta) %>%
  # mutate(Age = factor(age, levels = c("YOUNG", "OLD"))) %>%
  filter(count > 0) %>%
  ggplot(aes(y = cell_type_publish,
             x = count,
             fill = cell_type_publish)) +
  geom_violin() +
  geom_jitter(alpha = 0.1) +
  # facet_wrap(vars(cell_type_publish),
  #            nrow = 1) +
  theme_cowplot() +
  panel_border() +
  scale_fill_d3() +
  labs(y = "Count",
       title = "Cplx1: DA Neurons")

#Interquartile range outlier detection function (IROF) ----
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

# FIG 5 CELL TYPES WITH AGE ----

spatial_all_cells$cell_type_publish %>% unique

ct_with_age <- spatial_all_cells %>%
  filter(mouse_id != "OLD_WT_REL121.1a") %>%
  group_by(mouse_id, age, cell_type_publish) %>%
  tally %>%
  group_by(mouse_id) %>%
  mutate(frac_of_total = n / sum(n)) %>%
  group_by(age, cell_type_publish) %>%
  mutate(age_ct_median = median(n)) %>%
  ungroup

# get young ct count median
ct_young_medians <- ct_with_age %>%
  filter(age == "YOUNG") %>%
  group_by(cell_type_publish) %>%
  mutate(median_n = median(n),
         median_frac_of_total = median(frac_of_total)) %>%
  arrange(cell_type_publish) %>%
  ungroup() %>%
  select(cell_type_publish,
         median_n,
         median_frac_of_total) %>%
  distinct()

# calculate log2 ratio and remove small cell types
ct_with_age <- ct_with_age %>%
  full_join(ct_young_medians) %>%
  mutate(
    log2_ratio_n = log2(n / median_n),
    log2_ratio_frac_of_total = log2(frac_of_total / median_frac_of_total)
  ) %>%
  arrange(cell_type_publish) %>%
  group_by(cell_type_publish) %>%
  filter(sum(n) > 100)

# set cell type plot order based on ratio
ct_plot_order <- ct_with_age %>%
  filter(age == "OLD") %>%
  group_by(cell_type_publish) %>%
  summarise(ratio = median(log2_ratio_frac_of_total)) %>%
  arrange(ratio) %>%
  mutate(cell_type_publish = str_remove(cell_type_publish, "_NEU")) %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  pull(cell_type_publish)

# plot
ct_with_age %>%
  mutate(cell_type_publish = str_remove(cell_type_publish, "_NEU")) %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  mutate(
    cell_type_publish = factor(cell_type_publish,
                               levels = ct_plot_order),
    age = factor(age, levels = c("YOUNG", "OLD"))
  ) %>%
  arrange(cell_type_publish, age) %>%
  ungroup %>%
  group_by(cell_type_publish) %>%
  # filter(median(n) > 100) %>%
  # filter(!outlier_iqr(n)) %>%
  
  # filter(!cell_type_publish %in% c("CORTICAL_THALAMIC_IL31RA_NEU",
  #                                  ))
  arrange(cell_type_publish, age) %>%
  
  mutate(
    cell_type_publish = factor(cell_type_publish,
                               levels = ct_plot_order),
    age = factor(age, levels = c("YOUNG", "OLD"))
  ) %>%
  group_by(cell_type_publish, age) %>%
  mutate(var = sd(n) / age_ct_median) %>%
  group_by(cell_type_publish) %>%
  filter(max(var) < 1) %>%
  # filter(n() == 6) %>%
  ggplot(aes(x = cell_type_publish,
             y = log2_ratio_frac_of_total)) +
  geom_boxplot(aes(color = age,
                   fill = age)) +
  geom_point(aes(group = age),
             size = 0.025,
             position = position_dodge(width = 0.75)) +
  # geom_jitter(aes(fill = age),
  #             width = 0.5,
  #             size = 0.5) +
  theme_cowplot(4) +
  geom_hline(yintercept = 0) +
  scale_fill_d3() +
  scale_color_d3() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    plot.margin = margin(7, 7, 7, 7, unit = "pt")
  ) +
  # scale_y_continuous(limits = c(-2.5, 2.5)) +
  labs(
    y = expression(Log[2] ~ Normalized ~ Ratio),
    x = "",
    fill = "Age",
    color = "Age"
  )

ggsave(
  "output/plots/figure_4/cell_types_ageing.png",
  bg = "transparent",
  width = 110,
  height = 70,
  units = "mm",
  dpi = 300
)

# FIG 5 MICROGLIA AGE ----

spatial_all_cells %>%
  # filter(cell_type == "Plp1_Plp1+++_Cathepsin-Tyrobp") %>%
  filter(mouse_id != "OLD_WT_REL121.1a") %>%
  group_by(age) %>%
  mutate(mouse_id_no = as.numeric(factor(mouse_id)),
         mouse_id_label = paste0(age, ": ", mouse_id_no)) %>%
  mutate(
    interest = ifelse(
      cell_type == "Plp1_Plp1+++_Cathepsin-Tyrobp",
      cell_type,
      "OTHER"
    ),
    size_interest = cell_type == "Plp1_Plp1+++_Cathepsin-Tyrobp"
  ) %>%
  mutate(mouse_id_label = factor(
    mouse_id_label,
    levels = c("YOUNG: 1",
               "OLD: 1",
               "YOUNG: 2",
               "OLD: 2",
               "YOUNG: 3",
               "OLD: 3")
  )) %>%
  mutate(interest = ifelse(str_detect(interest, "OTHER"), "OTHER", "MICROGLIA")) %>%
  arrange(desc(interest)) %>%
  ggplot(aes(
    x = x,
    y = -y,
    colour = interest,
    size = size_interest,
    alpha = size_interest
  )) +
  geom_point() +
  scale_color_manual(values = c(pal_d3()(1), "grey")) +
  scale_size_manual(values = c(0.02, 1), guide = 'none') +
  scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
  facet_wrap(vars(mouse_id_label),
             ncol = 2,
             scales = "free") +
  theme_minimal() +
  theme(
    rect = element_rect(fill = "transparent"),
    # panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.text = element_text(face = "bold", size = 14), 
    legend.position = "none"
  ) +
  labs(title = "AGEING MICROGLIA",
       color = "")
ggsave(
  "output/plots/figure_4/ageing_microglia.png",
  width = 160,
  height = 180,
  dpi = 300,
  units = "mm",
  bg = "transparent"
)

# FIG 5 AGE DE ----

# get a list of glial markers
glial_markers <- read_csv("output/markers/GLIAL.csv") %>%
  filter(pvals_adj < 0.01) %>%
  filter(logfoldchanges > 0) %>%
  pull(names)

# load matrix
m <- readMM("output/counts/adata.mtx")
# m <- as(m, "dgCMatrix")
meta <- spatial_all_cells
colnames(meta)[1] <- "cell_id"
var_names <- read_delim("output/counts/var_names.csv", delim = ',')
dimnames(m) <- list(meta$cell_id, var_names$geneID)

sce <- SingleCellExperiment(list(counts = t(m)),
                            colData = meta)
# detection <- tibble(count = colSums(counts(sce) > 0))
# detection <- tibble(count = rowSums(counts(sce) > 0))
# summary(detection)
colData(sce)$neu_glia <-
  ifelse(
    str_detect(colData(sce)$cell_type_publish, "GLIA|RBC|ASTRO|MICRO|OLIG"),
    "GLIA",
    "NEU"
  )

# mast ----

# combine SN and VTA DA
colData(sce)$cell_type_publish <- ifelse(colData(sce)$cell_type_publish == "DA_SN", 
                                         "DA", 
                                         ifelse(colData(sce)$cell_type_publish == "DA_VTA", 
                                                "DA", colData(sce)$cell_type_publish))

sca_split <-
  mclapply(unique(colData(sce)$cell_type_publish), function(x) {
    # sce <- sce[, colData(sce)$cell_type_DA_grouped == x]
    sce <- sce[, colData(sce)$cell_type_publish == x]
    sce <- logNormCounts(sce)
    colData(sce)$age <-
      factor(colData(sce)$age, levels = c("YOUNG", "OLD"))
    colData(sce)$genotype <-
      factor(colData(sce)$genotype, levels = c("WT", "OVX"))
    sca <- SceToSingleCellAssay(sce)
    sca <-
      sca[freq(sca) >= 0.1,] # only assess genes detected in > 0.1 of cells
  }, mc.cores = 18)
names(sca_split) <- unique(colData(sce)$cell_type_publish)

age_mast <- function(sca){
  # cell detection rate
  cdr <- colSums(assay(sca) > 0)
  colData(sca)$cdr <- scale(cdr)
  # model
  print("running zlm")
  z <- zlm(
    ~ cdr + age + (1 | mouse_id),
    sca = sca,
    method = "glmer",
    ebayes = F,
    strictConvergence = F,
    exprs_values = 'logcounts', 
    parallel = T
  )
  # LRT age
  print("running LRT")
  summaryCond <- summary(z, doLRT = 'ageOLD')
  summaryDt <- summaryCond$datatable
  print("running fcHurdle")
  fcHurdle <-
    merge(summaryDt[contrast == 'ageOLD' &
                      component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
          summaryDt[contrast == 'ageOLD' &
                      component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
  fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
  # FCTHRESHOLD <- log2(1.1)
  # [fdr<.05 & abs(coef)>FCTHRESHOLD]
  fcHurdle <-
    merge(fcHurdle, as_tibble(mcols(sca)), by = 'primerid') %>%
    as_tibble %>%
    arrange(fdr)
  
  return(fcHurdle)
}
MAST::summary()
age_de_DA <- age_mast(sca_split[["DA"]])

# age DE per cell type
age_de_per_ct <- mclapply(sca_split, age_mast, 
                          mc.cores = parallel::detectCores() - 2)

names(age_de_per_ct) <- names(sca_split)

bind_rows(age_de_per_ct, .id = 'cell_type') %>%
  write_csv("output/de/age.csv")

age_de_per_ct <- read_csv("output/de/age.csv")

age_de_per_ct %>%
  filter(fdr < 0.05) %>% 
  mutate(glial = str_detect(cell_type, "GLIA|OLIG|MICRO")) %>%
  group_by(glial) %>%
  summarise(median_lfc = median(abs(coef)))
  pull(cell_type) %>% unique

# age in DA neurons, SN-VTA combined

# sce_DA <- sce[, colData(sce)$cell_type_publish %in% c("DA_SN", "DA_VTA")]
# sce_DA <- logNormCounts(sce_DA)
# colData(sce_DA)$age <- factor(colData(sce_DA)$age, levels = c("YOUNG", "OLD"))
# colData(sce_DA)$genotype <- factor(colData(sce_DA)$genotype, levels = c("WT", "OVX"))
# sca_DA <- SceToSingleCellAssay(sce_DA)
# sca_DA <- sca_DA[freq(sca_DA) >= 0.1,] # only assess genes detected in > 0.1 of cells
# # cell detection rate
# cdr <- colSums(assay(sca_DA)>0)
# colData(sca_DA)$cdr <- scale(cdr)
# # model
# z <- zlm(~ cdr + age + (1 | mouse_id),
#          sca = sca_DA,
#          method = "glmer",
#          ebayes = F,
#          strictConvergence = F,
#          exprs_values = 'logcounts')
# # LRT age
# summaryCond <- summary(z, doLRT='ageOLD')
# summaryDt <- summaryCond$datatable
# fcHurdle <- merge(summaryDt[contrast=='ageOLD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
#                   summaryDt[contrast=='ageOLD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
# fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
# # FCTHRESHOLD <- log2(1.1)
# # [fdr<.05 & abs(coef)>FCTHRESHOLD]
# fcHurdle <- merge(fcHurdle, as_tibble(mcols(sca_DA)), by='primerid') %>%
#   as_tibble %>%
#   arrange(fdr)

# age in DA neurons, SN-VTA combined, TRAP EXPRESSED genes
# MB_FRACTION_META <- readRDS("input/trap/MB_FRACTION_META.rds")
# MB_FRACTION_TRANSLATED_GENES <- readRDS("input/trap/MB_FRACTION_TRANSLATED_GENES.rds")
#
# gene_id_conversion <- MB_FRACTION_META %>%
#   select(ensembl_gene_id, external_gene_name)
# MB_FRACTION_TRANSLATED_GENES_external_gene_name <- gene_id_conversion$external_gene_name[match(MB_FRACTION_TRANSLATED_GENES, gene_id_conversion$ensembl_gene_id)]
# sce_DA <- sce[rownames(sce) %in% MB_FRACTION_TRANSLATED_GENES_external_gene_name,
#               colData(sce)$cell_type_publish %in% c("DA_SN", "DA_VTA")]
# sce_DA <- logNormCounts(sce_DA)
# colData(sce_DA)$age <- factor(colData(sce_DA)$age, levels = c("YOUNG", "OLD"))
# colData(sce_DA)$genotype <- factor(colData(sce_DA)$genotype, levels = c("WT", "OVX"))
# sca_DA <- SceToSingleCellAssay(sce_DA)
# sca_DA <- sca_DA[freq(sca_DA) >= 0.1,] # only assess genes detected in > 0.1 of cells
# # cell detection rate
# cdr <- colSums(assay(sca_DA)>0)
# colData(sca_DA)$cdr <- scale(cdr)
# # model
# z <- zlm(~ cdr + age + (1 | mouse_id),
#          sca = sca_DA,
#          method = "glmer",
#          ebayes = F,
#          strictConvergence = F,
#          exprs_values = 'logcounts')
# # LRT age
# summaryCond <- summary(z, doLRT='ageOLD')
# summaryDt <- summaryCond$datatable
# fcHurdle <- merge(summaryDt[contrast=='ageOLD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
#                   summaryDt[contrast=='ageOLD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
# fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
# # FCTHRESHOLD <- log2(1.1)
# # [fdr<.05 & abs(coef)>FCTHRESHOLD]
# fcHurdle <- merge(fcHurdle, as_tibble(mcols(sca_DA)), by='primerid') %>%
#   as_tibble %>%
#   arrange(fdr)
#
# age_de_DA_TRAP_translated <- fcHurdle
# age_de_DA_TRAP_translated %>% write_csv("output/de/age_DA_TRAP_translated.csv")
age_de_DA_TRAP_translated <-
  read_csv("output/de/age_DA_TRAP_translated.csv")

# enframe(assay(sca_DA)[rownames(sca_DA) == 'Cplx1',],
#         name = 'cell_id', value = 'count') %>%
#   left_join(meta) %>%
#   mutate(Age = factor(age, levels = c("YOUNG", "OLD"))) %>%
#   filter(count > 0) %>%
#   ggplot(aes(x = Age,
#              y = count,
#              fill = Age)) +
#   geom_violin() +
#   geom_jitter(alpha = 0.1) +
#   # facet_wrap(vars(cell_type_publish),
#   #            nrow = 1) +
#   theme_cowplot() +
#   panel_border() +
#   scale_fill_d3() +
#   labs(y = "Count",
#        title = "Cplx1: DA Neurons")

# genotype DE per cell type
genotype_mast <- function(sca){
  # cell detection rate
  cdr <- colSums(assay(sca) > 0)
  colData(sca)$cdr <- scale(cdr)
  # model
  print("running zlm")
  z <- zlm(
    ~ cdr + genotype + (1 | mouse_id),
    sca = sca,
    method = "glmer",
    ebayes = F,
    strictConvergence = F,
    exprs_values = 'logcounts', 
    parallel = T
  )
  # LRT genotype
  print("running LRT")
  summaryCond <- summary(z, doLRT = 'genotypeOVX')
  summaryDt <- summaryCond$datatable
  print("running fcHurdle")
  fcHurdle <-
    merge(summaryDt[contrast == 'genotypeOVX' &
                      component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
          summaryDt[contrast == 'genotypeOVX' &
                      component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
  fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
  # FCTHRESHOLD <- log2(1.1)
  # [fdr<.05 & abs(coef)>FCTHRESHOLD]
  fcHurdle <-
    merge(fcHurdle, as_tibble(mcols(sca)), by = 'primerid') %>%
    as_tibble %>%
    arrange(fdr)
  
  return(fcHurdle)
}

# genotype DE per cell type
genotype_de_per_ct <- mclapply(sca_split, genotype_mast, 
                          mc.cores = parallel::detectCores() - 2)

names(genotype_de_per_ct) <- names(sca_split)

bind_rows(genotype_de_per_ct, .id = 'cell_type') %>%
  write_csv("output/de/genotype.csv")

genotype_de_per_ct <- read_csv("output/de/genotype.csv")

genotype_de_per_ct %>%
  filter(fdr < 0.05) %>% View
  # pull(primerid) %>% unique
  pull(cell_type) %>% unique

# genotype_de_per_ct <- mclapply(names(sca_split), function(n) {
#   sca <- sca_split[[n]]
#   # cell detection rate
#   cdr <- colSums(assay(sca) > 0)
#   colData(sca)$cdr <- scale(cdr)
#   # model
#   z <- zlm(
#     ~ cdr + genotype + (1 | mouse_id),
#     sca = sca,
#     method = "glmer",
#     ebayes = F,
#     strictConvergence = F,
#     exprs_values = 'logcounts'
#   )
#   # LRT age
#   summaryCond <- summary(z, doLRT = 'genotypeOVX')
#   summaryDt <- summaryCond$datatable
#   fcHurdle <-
#     merge(summaryDt[contrast == 'genotypeOVX' &
#                       component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
#           summaryDt[contrast == 'genotypeOVX' &
#                       component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
#   fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
#   # FCTHRESHOLD <- log2(1.1)
#   # [fdr<.05 & abs(coef)>FCTHRESHOLD]
#   fcHurdle <-
#     merge(fcHurdle, as_tibble(mcols(sca)), by = 'primerid') %>%
#     as_tibble %>%
#     arrange(fdr)
#   
#   return(fcHurdle)
# }, mc.cores = parallel::detectCores() - 2)
# 
# names(genotype_de_per_ct) <- names(sca_split)
# 
# bind_rows(genotype_de_per_ct, .id = 'cell_type') %>% write_csv("output/de/genotype.csv")
# 
# genotype_de_per_ct <- read_csv("output/de/genotype.csv")

# genotype in DA neurons, SN-VTA combined, TRAP EXPRESSED genes
gene_id_conversion <- MB_FRACTION_META %>%
  select(ensembl_gene_id, external_gene_name)
MB_FRACTION_TRANSLATED_GENES_external_gene_name <-
  gene_id_conversion$external_gene_name[match(MB_FRACTION_TRANSLATED_GENES,
                                              gene_id_conversion$ensembl_gene_id)]
sce_DA <-
  sce[rownames(sce) %in% MB_FRACTION_TRANSLATED_GENES_external_gene_name,
      colData(sce)$cell_type_publish %in% c("DA_SN", "DA_VTA")]
sce_DA <- logNormCounts(sce_DA)
colData(sce_DA)$age <-
  factor(colData(sce_DA)$genotype, levels = c("YOUNG", "OLD"))
colData(sce_DA)$genotype <-
  factor(colData(sce_DA)$genotype, levels = c("WT", "OVX"))
sca_DA <- SceToSingleCellAssay(sce_DA)
sca_DA <-
  sca_DA[freq(sca_DA) >= 0.1,] # only assess genes detected in > 0.1 of cells
# cell detection rate
cdr <- colSums(assay(sca_DA) > 0)
colData(sca_DA)$cdr <- scale(cdr)
# model
z <- zlm(
  ~ cdr + genotype + (1 | mouse_id),
  sca = sca_DA,
  method = "glmer",
  ebayes = F,
  # strictConvergence = F,
  exprs_values = 'logcounts'
)
# LRT genotype
summaryCond <- summary(z, doLRT = 'genotypeOVX')
summaryDt <- summaryCond$datatable
fcHurdle <-
  merge(summaryDt[contrast == 'genotypeOVX' &
                    component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryDt[contrast == 'genotypeOVX' &
                    component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
# FCTHRESHOLD <- log2(1.1)
# [fdr<.05 & abs(coef)>FCTHRESHOLD]
fcHurdle <-
  merge(fcHurdle, as_tibble(mcols(sca_DA)), by = 'primerid') %>%
  as_tibble %>%
  arrange(fdr)

fcHurdle %>% write_csv("output/de/genotype_DA_TRAP_translated.csv")

# check for GWAS overlap - genotype
genes <- genotype_de_per_ct %>% filter(fdr < 0.1) %>% pull(primerid)
genes[genes %in% {
  nalls_trap_combined$mmusculus_homolog_associated_gene_name %>% unique
}]

# check for GWAS overlap - age
age_de_per_ct %>%
  # bind_rows(.id = "cell_type") %>%
  filter(fdr < 0.1) %>% View
genes <- age_de_per_ct %>%
  # bind_rows(.id = "cell_type") %>%
  filter(fdr < 0.1) %>%
  pull(primerid) %>%
  unique
age_gwas_genes <-
  genes[genes %in% {
    nalls_trap_combined$mmusculus_homolog_associated_gene_name %>% unique
  }]
age_de_per_ct %>%
  filter(fdr < 0.1) %>%
  filter(primerid %in% age_gwas_genes)

nalls_trap_combined %>%
  filter(mmusculus_homolog_associated_gene_name %in% age_gwas_genes) %>% View

# ageing DE gprofiler ----

glial_age_gprofiler <- age_de_per_ct %>%
  filter(fdr < 0.05) %>%
  filter(!str_detect(cell_type, "NEU")) %>% 
  arrange(fdr) %>%
  pull(primerid) %>%
  gost(organism = "mmusculus", 
       ordered_query = T, 
       exclude_iea = T, 
       evcodes = T)

neuronal_age_gprofiler <- age_de_per_ct %>%
  filter(fdr < 0.05) %>%
  filter(str_detect(cell_type, "NEU")) %>% 
  arrange(fdr) %>%
  pull(primerid) %>%
  gost(organism = "mmusculus", 
       ordered_query = T, 
       exclude_iea = T, 
       evcodes = T)

neuronal_age_gprofiler$result %>% View

# # glmgampoi ----
# 
# library(glmGamPoi)
# 
# fit <- glm_gp(
#   sce,
#   design = ~ cell_type_publish + age +  age:cell_type_publish - 1,
#   reference_level = "OLIG",
#   on_disk = T
# )
# summary(fit)
# 
# # pseudobulk preparation ----
# 
# library(Matrix.utils)
# 
# # Subset metadata to only include the cluster and sample IDs to aggregate across
# groups <- colData(sce)[, c("cell_type_publish", "mouse_id")]
# 
# # Aggregate across cluster-sample groups
# pb <- aggregate.Matrix(t(counts(sce)),
#                        groupings = groups, fun = "sum")
# 
# pb <-
#   mclapply(unique(colData(sce)$cell_type_publish), function(ct) {
#     pb_split <- pb[str_detect(rownames(pb), paste0("^", ct)), ]
#     rownames(pb_split) <-
#       str_remove(rownames(pb_split), "[:graph:]+(?=YOUNG|OLD)")
#     pb_split <- pb_split[, colMins(pb_split) > 0]
#     return(pb_split)
#   }, mc.cores = parallel::detectCores() - 2)
# names(pb) <- unique(colData(sce)$cell_type_publish)
# 
# # n cells per mouse per cell type
# table(sce$cell_type_publish, sce$mouse_id)
# 
# pb_meta <- lapply(names(pb), function(ct) {
#   tibble(
#     cell_type = ct,
#     mouse = rownames(pb[[ct]]),
#     age = factor(ifelse(
#       str_detect(mouse, "YOUNG"), "YOUNG", "OLD"
#     ), levels = c("YOUNG", "OLD")),
#     genotype = factor(ifelse(str_detect(mouse, "WT"), "WT", "OVX"), levels = c("WT", "OVX"))
#   ) %>%
#     as.data.frame
# })
# names(pb_meta) <- unique(colData(sce)$cell_type_publish)
# 
# library(DESeq2)
# 
# pb_res <- mclapply(names(pb_meta), function(n) {
#   dds <- DESeqDataSetFromMatrix(t(pb[[n]]),
#                                 colData = pb_meta[[n]],
#                                 design = ~ age)
#   # vsd <- vst(dds)
#   # DESeq2::plotPCA(vsd, intgroup = "age")
#   dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
#   # plotDispEsts(dds)
#   # resultsNames(dds)
#   res <- results(dds,
#                  alpha = 0.05)
#   res <- lfcShrink(dds,
#                    coef = resultsNames(dds)[2],
#                    res = res)
#   res <- res %>%
#     data.frame() %>%
#     rownames_to_column(var = 'gene') %>%
#     as_tibble %>%
#     arrange(padj) %>%
#     filter(padj < 0.05)
# }, mc.cores = parallel::detectCores() - 2)
# names(pb_res) <- unique(colData(sce)$cell_type_publish)


# # group dopa
# colData(sce)$cell_type <- ifelse(str_detect(colData(sce)$cell_type, 'VTA$'),
#                                               "Dopaminergic_VTA",
#        ifelse(str_detect(colData(sce)$cell_type, 'SN$'),
#               "Dopaminergic_SN",
#               colData(sce)$cell_type))
# colData(sce)$cell_type_DA_grouped <- ifelse(str_detect(colData(sce)$cell_type, 'opamin'),
#                                             "Dopaminergic",
#                                             colData(sce)$cell_type)

# create a subset object of cells within 1 mm of DA neurons ----
# meta_da <- meta %>%
#   filter(str_detect(cell_type_publish, "^DA_"))
#
# cells_da <- meta_da$cell_id
# str_extract(cells_da[1], "[:graph:]+(?=_[:digit:]+)")
# cells_proximal_to_da <- mclapply(cells_da, function(cell){
#   x_root <- meta_da %>%
#     filter(cell_id == cell) %>%
#     pull(x)
#   y_root <- meta_da %>%
#     filter(cell_id == cell) %>%
#     pull(y)
#   cells <- meta %>%
#     filter(str_detect(cell_id, str_extract(cell, "[:graph:]+(?=_[:digit:]+)"))) %>%
#     filter(abs(x-x_root) <= 50) %>%
#     filter(abs(y-y_root) <= 50) %>%
#     pull(cell_id)
# }, mc.cores = 18) %>%
#   unlist %>%
#   unique
#
# meta %>%
#   filter(cell_id %in% cells_proximal_to_da) %>%
#   group_by(cell_type_publish) %>%
#   filter(n() > 50) %>%
#   ggplot(aes(x, -y,
#              colour = cell_type_publish)) +
#   geom_point() +
#   facet_wrap(vars(mouse_id), scales = "free")





# the number of cell per cell type ----
# n_cells_per_ct <- meta %>%
#   group_by(cell_type_publish) %>%
#   summarise(n_cells = n()) %>%
#   mutate(prop_cells = n_cells/sum(n_cells))

# LIBRA ----------
# ALL CELL TYPES

# make a seurat obj
# seu <- as.Seurat(sce, data = NULL)

# seu@meta.data$neu_glia <- ifelse(str_detect(seu@meta.data$cell_type_publish, "GLIA|RBC|ASTRO|MICRO|OLIG"), "GLIA", "NEU")
# seu@meta.data$group <- 'group1'
#
# DE_LIBRA_GLIAL <- run_de(seu,
#                          cell_type_col = "group",
#                          label_col = "neu_glia",
#                          replicate_col = "mouse_id",
#                          de_method = "edgeR",
#                          de_family = "pseudobulk",
#                          de_type = "LRT",
#                          n_threads = 18)
#
# DE_LIBRA_AGE <- run_de(seu,
#                        cell_type_col = "cell_type_publish",
#                        label_col = "age",
#                        replicate_col = "mouse_id",
#                        de_method = "edgeR",
#                        de_family = "pseudobulk",
#                        de_type = "LRT",
#                        n_threads = 18)
#
# DE_LIBRA_GENOTYPE <- Libra::run_de(seu,
#                             cell_type_col = "cell_type_publish",
#                             label_col = "genotype",
#                             replicate_col = "mouse_id",
#                             de_method = "edgeR",
#                             de_family = "pseudobulk",
#                             de_type = "QLF",
#                             n_threads = 18)
#
# seu_DA <- seu[,str_detect(seu@meta.data$cell_type_publish, 'DA_')]
# seu_DA <- seu_DA[rowSums(seu_DA@assays$originalexp@counts > 0) > dim(seu_DA)[2] * 0.2,]
# 'Snca' %in% {seu_DA@assays$originalexp@counts %>% rownames}
#
#
#
# DE_LIBRA_GENOTYPE <- Libra::run_de(seu_DA,
#                                    cell_type_col = "cell_type_publish",
#                                    label_col = "genotype",
#                                    replicate_col = "mouse_id",
#                                    de_method = "MAST",
#                                    de_family = "singlecell",
#                                    # de_type = "LRT",
#                                    n_threads = 18)
#
# # DA genotyping genes
# DE_LIBRA_GENOTYPE %>%
#   filter(gene == "Snca") %>% View
#   filter(p_val_adj < 0.05) %>%
#   group_by(gene) %>%
#   mutate(n = n()) %>%
#   filter(str_detect(cell_type, "DA_"))
#
# # DA ageing genes
# DE_LIBRA_AGE %>%
#   filter(p_val_adj < 0.01) %>%
#   group_by(gene) %>%
#   mutate(n = n()) %>%
#   filter(str_detect(cell_type, "DA_"))

# DE ageing summary stats ----
DE_ageing_summary <- age_de_per_ct %>%
  # bind_rows({
  #   age_de_DA_TRAP_translated %>%
  #     mutate(.before = everything(),
  #            cell_type = "DA")
  # }) %>%
  # bind_rows(.id = 'cell_type') %>%
  filter(fdr < 0.05) %>%
  group_by(primerid) %>%
  mutate(n_ct = n()) %>%
  group_by(cell_type) %>%
  mutate(n_genes = n()) %>%
  group_by(primerid) %>%
  mutate(discordant = n_distinct(sign(`Pr(>Chisq)`)) == 2,
         incl_glia = sum(str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO")) != 0)

# histogram of number of cell types in which genes are DE
DE_ageing_summary %>%
  select(primerid, n_ct) %>%
  distinct %>%
  # mutate(neu_glia = ifelse(str_detect(cell_type, "GLIA|OLIG|RBC"),
  #                          "GLIAL", "NEURONAL")) %>%
  ggplot(aes(x = n_ct)) +
  geom_histogram(binwidth = 1,
                 colour = "black",
                 fill = pal_d3()(1)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = seq(1, 25, 1),
                     limits = c(0.5, 4.5)) +
  scale_fill_d3() +
  labs(x = "Number of Cell Types in which DE",
       y = "Number of DE Genes",
       fill = "") +
  theme(legend.position = c(0.75, 0.85))
ggsave(
  "output/plots/figure_5/age_DE_histogram_cell_types.png",
  width = 6,
  height = 4,
  dpi = 300
)

# breakdown of cell type for single cell type DE genes
DE_ageing_summary %>%
  mutate(neu_glia = ifelse(
    str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO"),
    "GLIAL",
    "NEURONAL"
  )) %>%
  # filter(n_ct == 1) %>%
  ggplot(aes(x = neu_glia,
             fill = neu_glia)) +
  geom_bar(colour = "black") +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_d3(palette = "category20") +
  labs(x = "Cell Type Category",
       y = "Number of DE Genes",
       fill = "") +
  theme(legend.position = "none")
ggsave(
  "output/plots/figure_5/age_DE_single_cell_type_breakdown.png",
  width = 5,
  height = 4,
  dpi = 300
)
# Glial expression of DE genes, relative to neuronal
# ageing genes that are commonly glia and neuronal DE
age_genes_incl_glia <- DE_ageing_summary %>%
  mutate(neu_glia = ifelse(
    str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO"),
    "GLIAL",
    "NEURONAL"
  )) %>%
  group_by(primerid) %>%
  filter(sum(neu_glia == "NEURONAL") > 0 &
           sum(neu_glia == "GLIAL") > 0) %>%
  pull(primerid) %>%
  unique %>%
  sort

age_genes_incl_glia_counts_meta <-
  m[, colnames(m) %in% age_genes_incl_glia] %>%
  as.matrix() %>%
  as_tibble(rownames = "cell_id") %>%
  pivot_longer(-cell_id,
               names_to = 'gene',
               values_to = 'count') %>%
  left_join(meta)

age_glial_ratios <- mclapply(age_genes_incl_glia, function(x) {
  print(x)
  filter_gene <- x
  age_genes_incl_glia_counts_meta %>%
    filter(gene == filter_gene) %>%
    filter(cell_type_publish %in% {
      DE_ageing_summary %>%
        filter(primerid == filter_gene) %>%
        pull(cell_type)
    }) %>%
    mutate(neu_glia = ifelse(
      str_detect(cell_type_publish, "GLIA|OLIG|RBC|MICRO|ASTRO"),
      "GLIAL",
      "NEURONAL"
    )) %>%
    group_by(neu_glia) %>%
    summarise(count = mean(count)) %>%
    pivot_wider(names_from = "neu_glia",
                values_from = "count") %>%
    mutate(ratio = log2(GLIAL / NEURONAL)) %>%
    pull(ratio) %>%
    unlist
}, mc.cores = 18)
names(age_glial_ratios) <- age_genes_incl_glia

enframe(unlist(age_glial_ratios),
        name = "gene",
        value = "ratio") %>%
  ggplot(aes(x = ratio)) +
  geom_histogram(colour = "black",
                 binwidth = 0.2,
                 fill = pal_d3()(1)) +
  theme_cowplot() +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  scale_y_continuous(breaks = seq(0, 8, 2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(Log[2] ~ Ratio ~ of ~ Glial:Neuronal ~ Expression),
       y = "Number of DE Genes")
ggsave(
  "output/plots/figure_5/age_DE_shared_glial_expression_ratio.png",
  width = 6,
  height = 4,
  dpi = 300
)

# tally of de genes per cell type, with and without shared glial genes
# show which genes were removed
DE_ageing_summary %>%
  mutate(remove = ifelse(primerid %in% age_genes_incl_glia,
                         "REMOVE", "KEEP")) %>%
  group_by(cell_type, remove) %>%
  tally %>%
  mutate(neu_glia = ifelse(
    str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO"),
    "GLIAL",
    "NEURONAL"
  )) %>%
  mutate(cell_type = str_remove(cell_type, "_NEU")) %>%
  mutate(cell_type = str_replace_all(cell_type, "_", " ")) %>%
  pivot_wider(names_from = "remove",
              values_from = "n") %>%
  mutate(
    KEEP = replace_na(KEEP, 0),
    KEEP = ifelse(neu_glia == "GLIAL",
                  KEEP + REMOVE,
                  KEEP),
    REMOVE = ifelse(neu_glia == "GLIAL",
                    0,
                    REMOVE)
  ) %>%
  mutate(sum = KEEP + REMOVE) %>%
  pivot_longer(c(KEEP, REMOVE),
               names_to = "status",
               values_to = "n") %>%
  mutate(status = ifelse(status == "KEEP",
                         "Exclusively neuronal",
                         "Shared glial")) %>%
  mutate(status = factor(status, levels = c("Shared glial", "Exclusively neuronal"))) %>%
  arrange(status) %>%
  ggplot(aes(x = n,
             y = cell_type,
             fill = status)) +
  geom_bar(stat = "identity",
           colour = "black") +
  scale_fill_manual(values = c("grey", pal_d3()(1))) +
  theme_cowplot() +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme(legend.position = c(0.5, 0.75)) +
  labs(x = "Number of DE Genes",
       y = "Cell Type",
       fill = "")

# show final ageing DE tally
DE_ageing_summary %>%
  # mutate(remove = ifelse(primerid %in% age_genes_incl_glia,
  #                        "REMOVE", "KEEP")) %>%
  group_by(cell_type) %>%
  tally %>%
  mutate(neu_glia = ifelse(
    str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO"),
    "GLIAL",
    "NEURONAL"
  )) %>%
  mutate(cell_type = str_remove(cell_type, "_NEU")) %>%
  mutate(cell_type = str_replace_all(cell_type, "_", " ")) %>%
  # pivot_wider(names_from = "remove",
  #             values_from = "n") %>%
  # mutate(KEEP = replace_na(KEEP, 0),
  #        KEEP = ifelse(neu_glia == "GLIAL",
  #                      KEEP + REMOVE,
  #                      KEEP),
  #        REMOVE = ifelse(neu_glia == "GLIAL",
  #                        0,
  #                        REMOVE)) %>%
  # mutate(sum = KEEP + REMOVE) %>%
  # pivot_longer(c(KEEP, REMOVE),
#              names_to = "status",
#              values_to = "n") %>%
# mutate(status = ifelse(status == "KEEP",
#                        "Exclusively neuronal",
#                        "Shared glial")) %>%
# mutate(status = factor(status, levels = c("Shared glial", "Exclusively neuronal"))) %>%
# arrange(status) %>%
# filter(status != "Shared glial" | str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO")) %>%
ungroup %>%
  filter(n != 0) %>%
  arrange(n) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type)) %>%
  ggplot(aes(x = n,
             y = cell_type,
             fill = neu_glia)) +
  geom_bar(stat = "identity",
           colour = "black") +
  scale_fill_d3() +
  theme_cowplot(6) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = "Number of DE Genes",
       y = "",
       fill = "") +
  facet_grid2(vars(neu_glia),
              scales = "free",
              space = "free") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "black")
  ) +
  panel_border() +
  geom_text(aes(label = n), hjust = -0.5, size = 2)

ggsave(
  "output/plots/figure_4/age_DE_tally.png",
  width = 70,
  height = 70,
  units = "mm",
  dpi = 300
)

# genes for ageing heatmap that don't include glia
DE_ageing_summary %>%
  filter(!incl_glia) %>%
  pull(primerid) %>% unique

DE_ageing_summary %>%
  filter(!incl_glia) %>%
  pull(primerid) %>%
  unique

DE_ageing_summary %>%
  filter(!incl_glia) %>%
  pull(primerid) %>%
  unique

plot_data <- DE_ageing_summary %>%
  # filter(!incl_glia) %>%
  ungroup %>%
  tidyr::expand(primerid, cell_type) %>%
  left_join({
    DE_ageing_summary %>%
      # filter(!incl_glia) %>%
      ungroup
  }) %>%
  mutate(neu_glia = ifelse(
    str_detect(cell_type, "GLIA|OLIG|RBC|MICRO|ASTRO"),
    "GLIAL",
    "NEURONAL"
  )) %>%
  filter(!str_detect(primerid, '^Gm[:digit:]+|Rik$')) %>%
  # mutate(coef = replace_na(coef, 0)) %>%
  mutate(cell_type = str_remove(cell_type, "_NEU")) %>%
  mutate(cell_type = str_replace_all(cell_type, "_", " ")) %>%
  mutate(coef = ifelse(coef > 1, 1,
                       ifelse(coef < -1, -1,
                              coef))) %>%
  # mutate(change = ifelse(coef > 0, "Upregulated", "Downregulated"), 
  #        change = factor(change, levels = c("Downregulated", "Upregulated"))) %>%
  group_by(cell_type) %>%
  mutate(total = sum(!is.na(coef)))

plot_data %>%
  ggplot() +
  geom_tile(aes(x = primerid,
                y = cell_type,
                fill = coef),
            colour = "black",
            size = 0.05) +
  # scale_fill_manual(
  #   values = c("blue", "red1"),
  #   na.value = "gray70",
  #   na.translate = F
  # ) +
  scale_fill_gradientn(colours = c("blue",
                                   "blue", 
                                   "lightgrey", 
                                   "red", "red"), 
                       values = scales::rescale(c(-1, -0.5, 0, 0.5, 1),
                                        from = c(-1, 1)),
                       # breaks = c(-1, -0.05, 0, 0.05, 1),
                       # limits = c(-0.2, 0.2),
                       na.value = "white") +
  # scale_y_discrete(position = "top") +
  theme_cowplot(6) +
  panel_border() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 5
    ),
    # axis.text.x = element_text(size = 4),
    axis.ticks.x = element_line(size = 0.05),
    strip.background = element_rect(fill = "white", colour = "black"),
    plot.margin = margin(7, 35, 14, 7, unit = "pt"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.margin = margin(-10, -10, -10, -10)
  ) +
  labs(x = "",
       y = "",
       fill = expression(Log[2] ~ Fold ~ Change ~ `in` ~ Age)) +
  facet_grid2(rows = vars(neu_glia),
              scales = "free",
              space = "free") +
  panel_border(colour = "black", size = 0.5) +
  geom_text(aes(y = cell_type, 
                label = total
                ), x = 60, 
            size = 2, 
            check_overlap = T) +
  coord_cartesian(xlim = c(0, 61), clip = "off")

ggsave(
  "output/plots/figure_4/age_DE_heatmap.png",
  width = 180,
  height = 45,
  units = "mm",
  dpi = 300
)

# olig ageing ----
library(ggrepel)

  age_de_per_ct %>%
  filter(str_detect(cell_type, "OLIG")) %>%
  ggplot(aes(
    x = coef,
    y = -log10(fdr),
    label = ifelse(fdr < 0.01, primerid, "")
  )) +
  geom_point(aes(
    colour = fdr < 0.05,
    alpha = fdr < 0.05,
    size = fdr < 0.05
  )) +
  scale_x_continuous(limits = c(-0.5, 0.5), oob = scales::squish) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_size_manual(values = c(1, 2)) +
  scale_color_d3() +
  theme_cowplot() +
  geom_label_repel() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  labs(
    x = expression(Log[2] ~ FC),
    y = expression(-Log[10] ~ FDR),
    colour = "FDR < 0.05",
    size = "FDR < 0.05",
    alpha = "FDR < 0.05"
  ) +
  theme(legend.position = "top")
ggsave(
  "output/plots/figure_5/age_DE_oligo_volcano.png",
  width = 5,
  height = 6,
  dpi = 300
)

# gprofiler2 ageing ----

library(gprofiler2)
genes <-
  bind_rows(age_de_per_ct, .id = "cell_type") %>%
  filter(fdr < 0.05) %>%
  filter(str_detect(cell_type, "OLIG")) %>%
  filter(coef > 0) %>%
  arrange(fdr) %>%
  pull(primerid)

gost(
  genes,
  organism = 'mmusculus',
  ordered_query = T,
  as_short_link = T
)

# DA ageing with TRAP ----
MB_AGE_META <-
  readRDS("input/trap/MB_AGE_META.rds")
# bind_rows(age_de_per_ct, .id = "cell_type") %>% View

MB_AGE_META %>%
  filter(padj < 0.05) %>%
  group_by(specific) %>%
  tally

# load GWAS homologs
homologs <- readRDS("input/trap/homologs_gwas.rds")

library(ggrepel)
MB_AGE_META %>%
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    size = -log10(padj)
  )) +
  geom_point(aes(colour = -log10(padj))) +
  geom_point(data = {
    MB_AGE_META %>%
      filter(padj < 0.05) %>%
      filter(external_gene_name %in% homologs$mmusculus_homolog_associated_gene_name)
  },
  shape = 21,
  colour = "black") +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dotted") +
  geom_label_repel(
    data = {
      MB_AGE_META %>%
        filter(external_gene_name %in% homologs$mmusculus_homolog_associated_gene_name) %>%
        filter(external_gene_name != "Gm20186")
    },
    aes(label = ifelse(padj < 0.01, external_gene_name, "")),
    size = 1.6,
    max.overlaps = 6000
  ) +
  annotate(
    geom = "text",
    x = -1.4,
    y = 2.2,
    label = "FDR = 0.01",
    color = "black",
    size = 2
  ) +

  theme_cowplot(8) +
  
  scale_x_continuous(limits = c(-1.5, 1.5), oob = scales::squish) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    limits = c(0, 7),
    oob = scales::squish
  ) +
  # coord_cartesian(ylim = c(0, 7.5), clip = "off") +
  # annotate(
  #   geom = "text",
  #   x = -1,
  #   y = 7.25,
  #   label = paste0({MB_AGE_META %>%
  #       filter(padj < 0.01, 
  #              log2FoldChange < 0) %>%
  #       nrow()}, " genes"),
  #   color = "black",
  #   size = 2
  # ) +
  # annotate(
  #   geom = "text",
  #   x = 1,
  #   y = 7.25,
  #   label = paste0({MB_AGE_META %>%
  #       filter(padj < 0.01, 
  #              log2FoldChange > 0) %>%
  #       nrow()}, " genes"),
  #   color = "black",
  #   size = 2
  # ) +
  # scale_alpha_manual(values = c(0.2, 1)) +
  scale_size_continuous(range = c(0.1, 0.75),
                        guide = "none") +
  # scale_color_d3() +
  scale_color_gradientn(
    colours = c(pal_d3()(1),
                pal_d3()(2)[2],
                pal_d3()(2)[2]),
    values = c(min(-log10(MB_AGE_META$padj)),
               5,
               max(-log10(MB_AGE_META$padj))) / max(-log10(MB_AGE_META$padj))
  ) +
  
  labs(x = expression(Log[2] ~ FC),
       y = expression(-Log[10] ~ FDR)) +
  theme(legend.position = "none")

ggsave(
  "output/plots/figure_4/age_DE_DA_TRAP_volcano.png",
  width = 85,
  height = 70,
  dpi = 300,
  units = 'mm'
)

paste0({MB_AGE_META %>%
          filter(padj < 0.01,
                 log2FoldChange < 0) %>%
          nrow()}, " genes")
paste0({MB_AGE_META %>%
    filter(padj < 0.01,
           log2FoldChange > 0) %>%
    nrow()}, " genes")

# TRAP age fgsea ----

library(fgsea)
anno_fgsea <- readRDS("input/trap/anno_fgsea.rds")
anno_human <- readRDS("input/trap/anno_human.rds")

plotEnrichment <-
  function (pathway,
            stats,
            gseaParam = 1,
            ticksSize = 0.2){
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))
    pathway <-
      unname(as.vector(na.omit(match(
        pathway, names(statsAdj)
      ))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj,
                            selectedStats = pathway,
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms)) / 8
    x = y = NULL
    print(max(tops))
    g <-
      ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = "black") + 
      # geom_hline(yintercept = max(tops),
      #                                                                         colour = "red",
      #                                                                         linetype = "dashed") + geom_hline(
      #                                                                           yintercept = min(bottoms),
      #                                                                           colour = "red",
      #                                                                           linetype = "dashed"
      #                                                                         ) + 
      geom_hline(yintercept = 0,
                                                                                             colour = "black") + 
      geom_line(color = ifelse(max(abs(tops)) > max(abs(bottoms)), 
                                                                                                                                          "red", "blue")) + theme_bw() +
      geom_segment(
        data = data.frame(x = pathway),
        mapping = aes(
          x = x,
          y = -diff /
            2,
          xend = x,
          yend = diff / 2
        ),
        size = ticksSize
      ) +
      theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "rank", y = "enrichment score")
    g
  }

# KEGG

# load gmt
# pathways <- fgsea::gmtPathways("R/objects/msigdb.v7.4.symbols.gmt.txt")
pathways <- fgsea::gmtPathways("input/trap/c2.cp.kegg.v7.5.1.symbols.gmt.txt")
# pathways <- fgsea::gmtPathways("input/trap/c5.go.v7.4.symbols.gmt.txt")
# pathways <- fgsea::gmtPathways("R/objects/c5.go.mf.v7.4.symbols.gmt.txt")

ranks <- MB_AGE_META %>%
  # filter(!conflict
  #          # ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES
  #        ) %>%
  mutate(score = -log10(padj) * log2FoldChange) %>%
  mutate(rank = score) %>%
  inner_join(anno_human) %>%
  select(hsapiens_homolog_associated_gene_name,
         rank) %>%
  filter(hsapiens_homolog_associated_gene_name != "") %>%
  distinct(across(-rank), .keep_all = T) %>%
  deframe

set.seed(1234)
MB_AGE_FGSEA_KEGG <- fgsea(pathways,
                      ranks,
                      eps = 0,
                      minSize = 10,
                      maxSize = 500,
                      nPermSimple = 10000,
                      BPPARAM = BiocParallel::MulticoreParam(),
                      nproc = 18)
library(ggpubr)
fgsea_size = 7

p1 <- plotEnrichment(pathways[["KEGG_PARKINSONS_DISEASE"]],
                ranks) + labs(title=paste0("KEGG Parkinson's Disease:\nP = ", {MB_AGE_FGSEA_KEGG %>%
                    filter(pathway == "KEGG_PARKINSONS_DISEASE") %>%
                    pull(padj) %>%
                    signif(3)}),
                    y = "Enrichment Score", 
                    x = "Rank") +
  theme_cowplot(fgsea_size) +
  theme(title = element_text(size = 6))
    # theme(axis.title.x = element_blank())

p3 <- plotEnrichment(pathways[["KEGG_LYSOSOME"]],
                      ranks) + labs(title=paste0("KEGG Lysosome:\nP = ", {MB_AGE_FGSEA_KEGG %>%
                          filter(pathway == "KEGG_LYSOSOME") %>%
                          pull(padj) %>%
                          signif(3)}),
                          y = "Enrichment Score", 
                          x = "Rank") +
  theme_cowplot(fgsea_size) +
  theme(title = element_text(size = 6))

# GO

pathways <- fgsea::gmtPathways("input/trap/c5.go.v7.4.symbols.gmt.txt")
set.seed(1234)
MB_AGE_FGSEA_GO <- fgsea(pathways,
                           ranks,
                           eps = 0,
                           minSize = 10,
                           maxSize = 500,
                           nPermSimple = 10000,
                           BPPARAM = BiocParallel::MulticoreParam(),
                           nproc = 18)

p2 <- plotEnrichment(pathways[["GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX"]],
                     ranks) + labs(title=paste0("GO:CC Mitochondrial Protein-\ncontaining Complex:\nP = ", {MB_AGE_FGSEA_GO %>%
                         filter(pathway == "GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX") %>%
                         pull(padj) %>%
                         signif(3)}),
                         y = "Enrichment Score", 
                         x = "Rank") +
  theme_cowplot(fgsea_size) +
  theme(title = element_text(size = 4.5))
# theme(axis.title.x = element_blank(), 
#       axis.title.y = element_blank()) 

p4 <- plotEnrichment(pathways[["GOBP_SYNAPSE_ORGANIZATION"]],
                     ranks) + labs(title=paste0("GO:BP Synapse Organisation:\nP = ", {MB_AGE_FGSEA_GO %>%
                         filter(pathway == "GOBP_SYNAPSE_ORGANIZATION") %>%
                         pull(padj) %>%
                         signif(3)}),
                         y = "Enrichment Score", 
                         x = "Rank") +
  theme_cowplot(fgsea_size) +
  theme(title = element_text(size = 6))
  # theme(axis.title.y = element_blank())

p_age_fgsea <- ggarrange(
  p4, p3, p1, p2,
  nrow = 1, 
  ncol = 4)
  
ggsave(
  plot = p_age_fgsea,
  filename = "output/plots/figure_4/age_DE_DA_TRAP_FGSEA.png",
  width = 160,
  height = 40,
  dpi = 300, 
  units = "mm"
)

# DA ageing gprofiler ----
# # upreg
# MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   filter(log2FoldChange > 0) %>%
#   arrange(padj) %>%
#   pull(external_gene_name) %>%
#   gprofiler2::gost(
#     organism = 'mmusculus',
#     ordered_query = T,
#     as_short_link = T
#   )
# # downreg
# MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   filter(log2FoldChange < 0) %>%
#   arrange(padj) %>%
#   pull(external_gene_name) %>%
#   gprofiler2::gost(
#     organism = 'mmusculus',
#     ordered_query = T,
#     as_short_link = T
#   )

# DA ageing heatmap ----
# 
# # get genes with highest degree
# gene_labels <- sapply(split(degree(igr), membership(clusters)), function(x) names(sort(x, decreasing = T)[1]) )
# # Convert to readable names
# gene_labels <- sapply(gene_labels, get_external_gene_name)
# 
# # expand to the number of genes
# gene_labels <- gene_labels[membership(clusters)]
# 
# # plot an ageing heatmap
# plot_data <- t(scale(t(assay(vst(dds_C1_MB_IP))[match(names(membership(clusters)), rownames(dds_C1_MB_IP)),order(colData(dds_C1_MB_IP)$age)])))
# plot_data_quantiles <- quantile(plot_data, probs = seq(0, 1, length.out = 10))
# 
# # make an annotation column df
# annot_cols <- data.frame(Age = str_extract(colnames(plot_data), pattern = "(?<=MB_)[:alpha:]+"))
# rownames(annot_cols) <- colnames(plot_data)
# # make an annotation row df
# annot_rows <- data.frame(Group = gene_labels)
# # correct non-unique rownames
# rownames(plot_data) <- make.unique(rownames(plot_data))
# rownames(annot_rows) <- rownames(plot_data)
# plot_data %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   pivot_longer(-ensembl_gene_id,
#                names_to = "sample",
#                values_to = "lfc") %>%
#   left_join(as_tibble(annot_cols, rownames = "sample")) %>%
#   left_join(as_tibble(annot_rows, rownames = "ensembl_gene_id")) %>%
#   left_join({
#     enframe(mat_colors$Age, name = "Age", value = "age_colour")
#   }) %>%
#   left_join({
#     enframe(mat_colors$Group, name = "Group", value = "group_colour")
#   }) %>%
#   ggplot(aes(y = ensembl_gene_id,
#              x = sample,
#              fill = lfc)) +
#   geom_tile(colour = "black") +
#   theme_cowplot() +
#   scale_fill_gradient2(low = "blue",
#                        mid = "grey",
#                        high = "red") +
#   facet_wrap(vars(Group), scales = "free_y") +
#   theme(
#     axis.text.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks = element_blank()
#   ) +
#   panel_border()
# 
# 
# 
# pheatmap(
#   mat = plot_data[order(names(gene_labels)),],
#   scale = "none",
#   treeheight_row = 0,
#   show_rownames = F,
#   show_colnames = F,
#   annotation_col = annot_cols,
#   annotation_row = annot_rows,
#   annotation_colors = mat_colors,
#   # color = brewer.pal(n = length(plot_data_quantiles) - 1,
#   #                    name = "PiYG"),
#   color = colorspace::diverging_hcl(
#     11,
#     h = c(250, 10),
#     c = 100,
#     l = c(37, 88),
#     power = c(0.7, 1.7)
#   ),
#   # color = colorspace::diverging_hcl(11,
#   #                                   h = c(180, 50), c = 80, l = c(20, 95), power = c(0.7, 1.3)),
#   # color = brewer.pal(n = 8,
#   #                    name = "RdBu"),
#   breaks = plot_data_quantiles,
#   border_color = NA,
#   cluster_rows = FALSE,
#   cluster_cols = FALSE,
#   gaps_row = which(unlist(
#     lapply(seq(2:length(names(gene_labels)[order(names(gene_labels))])), function(n)
#       names(gene_labels)[order(names(gene_labels))][n] > names(gene_labels)[order(names(gene_labels))][n - 1])
#   )),
#   gaps_col = which(unlist(
#     lapply(seq(2:ncol(plot_data)), function(n)
#       colData(dds_C1_MB_IP)$age[n] != colData(dds_C1_MB_IP)$age[n - 1])
#   ))
#   # cellwidth = 6,
#   # cellheight = 2
# )


# # calculate the number of cell types in which each gene is DEG ----
# # the number of up cases
# # the number of down cases
# # the number of glial cases
# 
# age_deg_meta <- age_de_per_ct %>%
#   bind_rows(.id = "cell_type_publish") %>%
#   left_join(n_cells_per_ct) %>%
#   group_by(primerid) %>%
#   summarise(
#     n = n(),
#     n_up = sum(coef > 0),
#     n_down = sum(coef < 0),
#     de_glia = sum(str_detect(cell_type_publish, "GLIA|OLIG|RBC")),
#     prop_cells = sum(prop_cells),
#     glial_marker = primerid %in% glial_markers
#   ) %>%
#   distinct()
# 
# # get DA age degs ----
# age_da_degs <-
#   unique(c(age_de_per_ct$DA_SN$primerid,
#            age_de_per_ct$DA_VTA$primerid))
# 
# age_de_per_ct %>%
#   bind_rows(.id = "cell_type") %>%
#   filter(str_detect(cell_type, "DA_")) %>%
#   left_join(age_deg_meta) %>%
#   arrange(prop_cells) %>% View
# 
# # exclusively DAergic
# age_de_da_exclusive <- age_de_per_ct %>%
#   bind_rows(.id = "cell_type") %>%
#   left_join(age_deg_meta) %>%
#   group_by(primerid) %>%
#   filter(n() <= 2 &
#            all(str_detect(cell_type, "^DA_"))) %>%
#   arrange(primerid)
# 
# # non-exclusively DAergic
# age_de_median_coef_da_nonexclusive <-
#   age_de_per_ct %>%
#   bind_rows(.id = "cell_type") %>%
#   left_join(age_deg_meta) %>%
#   filter(!str_detect(cell_type, "DA_")) %>%
#   filter(primerid %in% age_da_degs) %>%
#   filter(!primerid %in% age_de_da_exclusive$primerid) %>%
#   group_by(primerid) %>%
#   summarise(median_coef = median(coef))
# 
# age_de_per_ct %>%
#   bind_rows(.id = "cell_type") %>%
#   left_join(age_deg_meta) %>%
#   filter(str_detect(cell_type, "DA_")) %>%
#   filter(!primerid %in% age_de_da_exclusive$primerid) %>%
#   left_join(age_de_median_coef_da_nonexclusive) %>%
#   filter(sign(coef) == sign(median_coef)) %>%
#   mutate(log2_ratio_coef = log2(coef / median_coef)) %>%
#   arrange(desc(log2_ratio_coef)) %>%
#   # filter(de_glia == 0) %>%
#   mutate(rank = row_number()) %>%
#   ggplot(aes(x = prop_cells,
#              y = log2_ratio_coef)) +
#   geom_point()
# 
# age_de_per_ct %>%
#   bind_rows(.id = "cell_type") %>%
#   left_join(age_deg_meta) %>%
#   # filter(!str_detect(cell_type, "DA_")) %>%
#   filter(primerid == "Tmtc2") %>%
#   mutate(label = ifelse(str_detect(cell_type, "^DA_"),
#                         "DOPAMINERGIC", "OTHER")) %>%
#   ggplot(aes(x = cell_type,
#              y = coef,
#              fill = label)) +
#   geom_col()
# 
# 


# Ageing STRING plot ----

library(igraph)
library(STRINGdb)
library(ggsci)
library(ggplotify)
library(ggpubr)

anno_human <- readRDS("input/trap/anno_human.rds")

# build STRING Interaction database ---------------------------------------

# setting to human database, for human interactions
string_db <- STRINGdb$new(version = "11",
                          # species = 10090, # mus musculus
                          species = 9606,
                          score_threshold = 400, # medium confidence
                          input_directory = "")

# setting to human database, for human interactions
string_db_mouse <- STRINGdb$new(version = "11",
                                species = 10090, # mus musculus
                                score_threshold = 400, # medium confidence
                                input_directory = "")

string_genes_df <- data.frame(gene = anno_human %>%
                                filter(hsapiens_homolog_ensembl_gene != "") %>%
                                pull(unique(hsapiens_homolog_ensembl_gene))) %>%
  left_join(anno_human, by = c("gene" = "hsapiens_homolog_ensembl_gene")) %>%
  distinct()

#String can't map non-protein-coding genes (and many others!)
string <- string_db$map(string_genes_df,
                        "hsapiens_homolog_associated_gene_name",
                        removeUnmappedRows = TRUE)

# STRING neighbour plots --------------------------------------

get_string_neighbours <- function(external_gene_name,
                                  top_n = 10) {
  # gene <- string_genes_df[string_genes_df$external_gene_name == external_gene_name, ]$hsapiens_homolog_associated_gene_name
  gene <- string_db$mp(external_gene_name)
  
  top_n <- tibble(to = string_db$get_neighbors(c(gene))) %>%
    mutate(.before = to,
           from = gene) %>%
    mutate(combined_score = sapply(to, function(x) {
      string_db$get_interactions(c(from, x))[1, 3]})) %>%
    # mutate(from = external_gene_name,
    #        to = string$external_gene_name[match(to, string$STRING_id)]) %>%
    drop_na() %>%
    slice_max(order_by = combined_score,
              n = top_n) %>%
    pull(to)
  
  neighbours <- string_db$get_interactions(c(gene, top_n)) %>%
    mutate(from = string$hsapiens_homolog_associated_gene_name[match(from, string$STRING_id)],
           to = string$hsapiens_homolog_associated_gene_name[match(to, string$STRING_id)]) %>%
    drop_na() %>%
    distinct()
}

# plot neighbours
plot_neighbours <- function(gene, neighbours){
  igr <- graph_from_data_frame(neighbours,
                               directed = F)
  
  V(igr)$size <- 5
  # V(igr)$frame.color <- "white"
  # V(igr)$color <- "orange"
  # V(igr)$label <- ""
  # E(igr)$arrow.mode <- 0
  E(igr)$width <- 2
  
  inc.edges <- incident(igr,  gene, mode="all")
  # Set colors to plot the selected edges.
  ecol <- rep("gray80", ecount(igr))
  ecol[inc.edges] <- "orange"
  E(igr)[inc.edges]$width <- 3
  vcol <- rep("lightblue", vcount(igr))
  vcol[names(V(igr))==gene] <- "gold"
  vfont <- rep(1, vcount(igr))
  vfont[names(V(igr)) == gene] <- 2
  vcol[names(V(igr)) %in% sapply(MB_FRACTION_ENRICHED_GENES, get_human_gene_name)] <- "green"
  
  vframecol <- rep("grey10", vcount(igr))
  # vframecol[names(V(igr)) %in% sapply(MB_FRACTION_ENRICHED_GENES, get_external_gene_name)] <- "#2CA02CFF"
  
  V(igr)$shape <- ifelse(names(V(igr)) == gene, "csquare", "circle")
  
  plot(igr,
       vertex.color=vcol,
       vertex.frame.color = vframecol,
       vertex.label.font = vfont,
       vertex.label.color = "black",
       vertex.label.family = "Helvetica",
       vertex.label.cex = 0.8,
       vertex.label.dist = 3,
       vertex.label.degree = -pi/2,
       # edge.curved = 0.1,
       edge.color=ecol,
       layout = layout.circle,
       # layout = layout.fruchterman.reingold
  )
  # title(gene)
  
}


# MB age STRING interactions ------------------------------------

MB_AGE_META <-
  readRDS("input/trap/MB_AGE_META.rds")

MB_AGE_DE_SPECIFIC_META <- MB_AGE_META %>%
  filter(padj < 0.01) %>%
  # filter(specific != "Shared" &
  #          status %in% c("Enriched", "Spliced")) %>%
  select(ensembl_gene_id,
         external_gene_name,
         description,
         baseMean,
         log2FoldChange,
         padj,
         status,
         pub_pd,
         specific,
         outcome_TRAP)


# I want to visualise the genes with DE to see how they cluster
# and whether clusters display similar changes in gene expression
# (e.g. a mitchondrial cluster goes down)
# retrieve the string IDs for MB age DEs

# USING STRING HOMOLOGS: More genes are retained in the end (around 75 %)

string_MB_AGE_HUMAN <- string_db$map(data.frame(ensembl_gene_id = MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
                                                external_gene_name = MB_AGE_DE_SPECIFIC_META$external_gene_name),
                                     "external_gene_name",
                                     removeUnmappedRows = FALSE)

# string_MB_AGE_MOUSE <- string_db_mouse$map(data.frame(ensembl_gene_id = MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
#                                                 external_gene_name = MB_AGE_DE_SPECIFIC_META$external_gene_name),
#                                "ensembl_gene_id",
#                                removeUnmappedRows = FALSE)

input_strings <- string_MB_AGE_HUMAN %>%
  drop_na %>%
  pull(STRING_id)

# % of AGE DE that are in STRING
length(input_strings) / nrow(MB_AGE_DE_SPECIFIC_META)

# retrieve the interactions
interactions <- string_db$get_interactions(input_strings) %>%
  mutate(from = string_MB_AGE_HUMAN$external_gene_name[match(from, string_MB_AGE_HUMAN$STRING_id)],
         to = string_MB_AGE_HUMAN$external_gene_name[match(to, string_MB_AGE_HUMAN$STRING_id)]) %>%
  drop_na() %>%
  distinct()

# % of STRING-mapped AGE DEs with known interactions
length(unique(c(interactions$from, interactions$to))) / length(input_strings)

# create the igraph
igr <- graph_from_data_frame(interactions,
                             directed = F)

# add ensembl_gene_id
V(igr)$ensembl_gene_id <- string_MB_AGE_HUMAN[match(V(igr)$name, string_MB_AGE_HUMAN$external_gene_name),]$ensembl_gene_id

# rename the main name to ensembl_gene_id
V(igr)$external_gene_name <- V(igr)$name
V(igr)$name <- V(igr)$ensembl_gene_id

# convert combined score to decimal
E(igr)$combined_score <- E(igr)$combined_score/1000
# plot the distribution of scores
plot(hist(E(igr)$combined_score))
# STRING says that a score of 0.5 means 1 in 2 are likely false positives
# I choose to filter so that 1 in 5 are false positives (score >= 0.8)
# filter edges for confidence and remove orphan nodes
igr <- delete.edges(igr, which(E(igr)$combined_score < 0.8 ))

# MB AGE DEG Degree (PLOT) ----

plot_data <- enframe(degree(igr)) %>%
  # left_join(anno_human,
  #           by = c("name" = "hsapiens_homolog_associated_gene_name")) %>%
  select(ensembl_gene_id = name,
         value) %>%
  filter(ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_META$ensembl_gene_id) %>%
  right_join(MB_AGE_DE_SPECIFIC_META) %>%
  # mutate(value = replace_na(value, 0)) %>%
  select(ensembl_gene_id:description,
         degree = value) %>%
  mutate(degree_bin =
           ifelse(is.na(degree),
                  "No\nSTRING\nIdentifier",
                  ifelse(degree == 0,
                         "0",
                         ifelse(degree < 5,
                                "1-4",
                                ifelse(degree < 10,
                                       "5-9",
                                       "10+")))),
         degree_bin = factor(degree_bin,
                             levels = c("No\nSTRING\nIdentifier",
                                        "0",
                                        "1-4",
                                        "5-9",
                                        "10+")))

MB_AGE_DEGREE <- plot_data %>%
  select(external_gene_name,
         degree)

p_mb_age_degree <- plot_data %>%
  ggplot(aes(x = degree_bin, fill = degree_bin)) +
  geom_bar(colour = "black") +
  geom_text(aes( label = ..count..,
                 y= ..count.. ), stat= "count", vjust = -.5) +
  scale_fill_d3() +
  labs(x = "Number of STRING\nInteractions",
       y = "Number of Genes") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 8))

# MB AGE TABLE OF HIGHEST DEGREE GENES ----

plot_data %>%
  filter(degree_bin %in% c("5-9", "10+"))

# MB AGE further STRING analysis ----

# add lfc values
V(igr)$log2FoldChange <- MB_AGE_DE_SPECIFIC_META[match(V(igr)$ensembl_gene_id, MB_AGE_DE_SPECIFIC_META$ensembl_gene_id),]$log2FoldChange
# check for NA values
sum(is.na(V(igr)$log2FoldChange))

# check degree
# degree(igr)

# remove 0 degree vertices
igr <- induced_subgraph(igr, degree(igr) > 0)

# remove quadruplets and smaller
igr <- induced_subgraph(igr, components(igr)$membership %in% which(components(igr)$csize > 4))

# define a layout
layout <- layout_nicely(igr)
# layout <- layout.fruchterman.reingold(igr)
# layout <- layout.kamada.kawai(igr)

plot(igr,
     layout = layout,
     edge.color = "grey70",
     vertex.label = NA,
     vertex.size = 5,
     vertex.color = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
)

# simplify igr
# igr <- simplify(igr)

# cluster with fast greedy
clusters <- cluster_fast_greedy(igr,
                                weights = E(igr)$combined_score)

plot(
  clusters,
  igr,
  layout = layout,
  edge.color = "grey70",
  vertex.label = NA,
  vertex.size = 4,
  # col = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
)

# percentage of edges that cross communities (I am deleting these)
crossing(clusters, igr) %>% sum / length(E(igr))

# delete crossing edges
igr <- delete.edges(igr, which(crossing(clusters, igr)))

# redo layout
set.seed(6)
layout <- layout_nicely(igr)

# MB AGE STRING Clusters ----
# plot
p_mb_age_deg_igraph <- as.grob(
  ~ plot(
    clusters,
    igr,
    layout = layout,
    edge.color = "grey70",
    edge.width = 0.5,
    # vertex.label = NA,
    vertex.size = 2,
    vertex.color = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
    vertex.frame.width = 0.5, 
    vertex.frame.color = "white"
  )
)

ggarrange(p_mb_age_deg_igraph)

ggsave(
  filename = "output/plots/figure_4/age_DE_DA_TRAP_clusters.png",
  width = 75,
  height = 75,
  dpi = 300, 
  units = "mm"
)

clusters %>% membership %>% table

membership(clusters)[order(membership(clusters))] %>% unique

cluster_order <- membership(clusters)[order(membership(clusters))]

cluster_order %>% table

anno <- readRDS("/active/paper/input/trap/anno.rds")

get_external_gene_name <- function(ensembl_gene_id) {
  anno[anno$ensembl_gene_id == ensembl_gene_id, ]$external_gene_name
}

sapply(split(names(cluster_order), cluster_order), function(x){sapply(x, get_external_gene_name)})

library(gprofiler2)

gp <- lapply(unique(cluster_order), 
             function(cl){
               gp <- gprofiler2::gost(query = names(cluster_order)[cluster_order == cl], 
                                      organism = "mmusculus", 
                                      exclude_iea = T)
               gp$result %>% 
                 arrange(p_value) %>%
                 return
               
             })

gp[[3]] %>%
  filter(p_value < 0.1,
         term_size < 500) %>% 
  slice_head(n = 20)

# MASC ----

source("scripts/masc.R")
library(lme4)
spatial_all_cells$age <- factor(spatial_all_cells$age, levels = c("YOUNG", "OLD"))
spatial_all_cells$genotype <- factor(spatial_all_cells$genotype, levels = c("WT", "OVX"))
spatial_all_cells$sample_name <- as.numeric(factor(spatial_all_cells$sample_name))
spatial_all_cells$cell_type_num <- as.numeric(factor(spatial_all_cells$cell_type))
masc_all <- MASC(data = spatial_all_cells, cluster = spatial_all_cells$cell_type_num, 
     contrast = "age", 
     random_effects = "sample_name", 
     fixed_effects = "genotype")

# export MASC data for shiny site
#
# bind_rows(lapply(masc_all, function(x) {
#   x[[1]]
# }), .id = 'name') %>%
#   filter(term == "ageOLD") %>%
#   mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
#   mutate(cl_num = as.numeric(str_extract(name, "(?<=cluster)[:digit:]{1,2}"))) %>%
#   left_join({
#     spatial_all_cells %>%
#       select("cl_num" = cell_type_publish_num,
#              cell_type_publish) %>%
#       distinct
#   }) %>%
#   select(cell_type_publish,
#          estimate,
#          fdr) %>%
#   inner_join(read_csv(
#     "/active/paper_shiny/input/setup/cell_type_names.csv",
#     col_names = c("cell_type_publish",
#                   "cell_type_full_publish")
#   )) %>%
#   select("cell_type" = cell_type_publish,
#          "cell_type_full" = cell_type_full_publish,
#          estimate,
#          fdr) %>%
#   saveRDS("/active/paper_shiny/input/ageing_cell_type_numbers/data.rds")

bind_rows(lapply(masc_all,function(x){x[[1]]}),.id='name') %>%
  filter(term == "ageOLD") %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  mutate(cl_num = as.numeric(str_extract(name, "(?<=cluster)[:digit:]{1,2}"))) %>%
  left_join({
    spatial_all_cells %>%
      select("cl_num" = cell_type_publish_num,
             cell_type_publish) %>%
      distinct
  }) %>%
  mutate(cell_type_publish = str_remove(cell_type_publish, "_NEU")) %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  ggplot(aes(x = estimate, 
             y = -log10(fdr), 
             label = ifelse(fdr < 0.005, 
                            cell_type_publish, 
                            ""), 
             color = ifelse(cell_type_publish == "DA SN", 
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

ggsave(
  "output/plots/figure_4/cell_types_ageing_OR.png",
  bg = "transparent",
  width = 80,
  height = 70,
  units = "mm",
  dpi = 300
)

{MB_AGE_META %>% filter(padj < 0.01) %>% pull(external_gene_name) %>% str_to_upper()}[{MB_AGE_META %>% filter(padj < 0.01) %>% pull(external_gene_name) %>% str_to_upper()} %in% pathways$KEGG_PARKINSONS_DISEASE]

