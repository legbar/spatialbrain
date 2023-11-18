library(tidyverse)
# Sys.setenv("MC_CORES"=18L)
library(parallel)
library(cowplot)
library(ggsci)
library(data.table)
library(ggrepel)
library(DRIMSeq)
library(readxl)
library(httr)
library(jsonlite)
library(homologene)
library(biomaRt)
library(DESeq2)
select <- dplyr::select
filter <- dplyr::filter

MB_FRACTION_META <- readRDS("input/trap/MB_FRACTION_META.rds")
anno <- readRDS("input/trap/anno.rds")

# number of TRAP enriched genes ----
MB_FRACTION_META %>%
  filter(sumz_adj < 0.01) %>%
  filter(log2FoldChange > 0) %>%
  filter(!conflict) %>%
  nrow

# DTU ----

# load the drimseq processed object: MB FRACTION (Cohort 1)
d_MB_FRACTION <- readRDS("input/trap/d_MB_FRACTION.rds")
# plot the number of transcripts per gene
counts_MB_FRACTION_DTU <- readRDS("input/trap/counts_MB_FRACTION_DTU.rds")

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
dds <- readRDS("/active/paper/input/trap/dds.rds")
dds_C1 <- dds[, colData(dds)$cohort == "C1"] %>% filter_zeros()
dds_C1_MB <- dds_C1[, colData(dds_C1)$compartment == "MB"] %>% filter_zeros()
dds_C1_MB_FILTER <- filter_genes(dds_C1_MB, grouping = c("fraction", "gene_id"))
single_isoform_genes <- single_isoform_genes[single_isoform_genes %in% dds_C1_MB_FILTER]

plot_data <- enframe(table(table(counts(d_MB_FRACTION)$gene_id))) %>%
  select(n_isoforms = name,
         n_genes = value) %>%
  mutate(n_genes = as.numeric(n_genes)) %>%
  add_row(.before = 1,
          n_isoforms = "1",
          n_genes = length(single_isoform_genes)) %>%
  filter(n_genes > 10)

p_mb_fraction_n_isoforms <- ggplot(plot_data,
                                   aes(x = n_isoforms,
                                       y = n_genes,
                                       fill = factor(n_isoforms,
                                                     levels = seq(max(n_isoforms), 1, -1)))) +
  geom_col(colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Number of Transcripts",
       y = "Number of \nGenes Analysed") +
  scale_fill_brewer() +
  theme(legend.position = "none")  +
  geom_text(aes(label = n_genes,
                x = n_isoforms,
                y = n_genes),
            vjust = -0.5)

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
res_MB_FRACTION$gene_id_simple <- substr(res_MB_FRACTION$gene_id, 1, 18)

MB_FRACTION_ENRICHED_GENES <- readRDS("input/trap/MB_FRACTION_ENRICHED_GENES.rds")

# Plot the number of DTU genes
plot_data <- res_MB_FRACTION %>%
  select(ensembl_gene_id = gene_id_simple,
         -gene_id,
         everything()) %>%
  left_join(anno) %>%
  arrange(adj_pvalue) %>%
  mutate(enriched = ifelse(ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
                           "Enriched", "Not enriched"),
         adj_pvalue = ifelse(pvalue == 1, 1, adj_pvalue), # correct pval 1 (NA) rows
         signif = adj_pvalue < 0.01)

plot_data %>%
  mutate(signif = ifelse(signif, "Differential Usage", "No Differential Usage"), 
         signif = factor(signif, levels = c("No Differential Usage", "Differential Usage"))) %>%
  ggplot(aes(x = enriched,
             fill = signif)) +
  geom_bar(colour = "black") +
  theme_cowplot(8) +
  scale_fill_d3() +
  labs(x = "Gene-level Enrichment",
       y = "Number of Genes",
       fill = "") +
  theme(legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ..count..),
             position=position_stack(vjust=0.5)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave("output/plots/figure_3/mb_dtu.png", 
       bg = "transparent", 
       width = 60, height = 80,
       units = "mm",
       dpi = 300)

# plot DTU examples ----

# load PD GWAS genes considered
homologs <- readRDS("input/trap/homologs_gwas.rds")

# load anno_tx
anno_tx <- readRDS("input/trap/anno_tx.rds")

res_MB_FRACTION %>%
  select(ensembl_gene_id = gene_id_simple,
         -gene_id,
         everything()) %>%
  left_join(anno) %>%
  arrange(adj_pvalue) %>%
  # filter(external_gene_name %in% homologs$mmusculus_homolog_associated_gene_name) %>% # to focus on GWAS candidate genes
  mutate(enriched = ifelse(ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
                           "Enriched", "Not enriched"),
         adj_pvalue = ifelse(pvalue == 1, 1, adj_pvalue), # correct pval 1 (NA) rows
         signif = adj_pvalue < 0.01)

  
# Snca
data <- counts(d_MB_FRACTION) %>%
  as_tibble %>%
  filter(str_detect(gene_id, "ENSMUSG00000025889")) %>%
  select(ensembl_gene_id = gene_id,
         ensembl_transcript_id = feature_id,
         everything()) %>%
  mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                       "[:alnum:]+(?=\\.[:digit:]+)"),
         ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                             "[:alnum:]+(?=\\.[:digit:]+)")) %>%
  left_join(anno_tx) %>%
  # mutate(feature_id = paste0("Transcript ", seq(1, nrow(.), 1))) %>%
  pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
               names_to = "sample_name",
               values_to = "count") %>%
  group_by(sample_name) %>%
  filter(median(count) > 0) %>%
  mutate(count = count / sum(count)) %>%
  inner_join(colData(dds), copy = TRUE) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"))

library(ggbeeswarm)
p_snca <- data %>%
  ggplot(aes_string(x = "fraction", 
          y = "count", 
          fill = "fraction")) +
  geom_boxplot() +
  geom_quasirandom(size = 0.5, shape = 21) +
  # stat_pvalue_manual(stat.test, label = "ast", y.position = c(1, 1)) +
  labs(x = "Fraction", y = "Proportion of \nTotal Expression", title = "Snca") +
  theme_cowplot(6) +
  theme(legend.position = "none") +
  facet_wrap(vars(external_transcript_name)) +
  panel_border() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_d3()

# Casr
data <- counts(d_MB_FRACTION) %>%
  as_tibble %>%
  filter(str_detect(gene_id, "ENSMUSG00000051980")) %>%
  select(ensembl_gene_id = gene_id,
         ensembl_transcript_id = feature_id,
         everything()) %>%
  mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                       "[:alnum:]+(?=\\.[:digit:]+)"),
         ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                             "[:alnum:]+(?=\\.[:digit:]+)")) %>%
  left_join(anno_tx) %>%
  # mutate(feature_id = paste0("Transcript ", seq(1, nrow(.), 1))) %>%
  pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
               names_to = "sample_name",
               values_to = "count") %>%
  group_by(sample_name) %>%
  filter(median(count) > 0) %>%
  mutate(count = count / sum(count)) %>%
  inner_join(colData(dds), copy = TRUE) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"))

p_casr <- data %>%
  ggplot(aes_string(x = "fraction", 
                    y = "count", 
                    fill = "fraction")) +
  geom_boxplot() +
  geom_quasirandom(size = 0.5, shape = 21) +
  # stat_pvalue_manual(stat.test, label = "ast", y.position = c(1, 1)) +
  labs(x = "Fraction", y = "Proportion of \nTotal Expression", title = "Casr") +
  theme_cowplot(6) +
  theme(legend.position = "none") +
  facet_wrap(vars(external_transcript_name), nrow = 1) +
  panel_border() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_d3()

# Stxbp1
data <- counts(d_MB_FRACTION) %>%
  as_tibble %>%
  filter(str_detect(gene_id, "ENSMUSG00000026797")) %>%
  select(ensembl_gene_id = gene_id,
         ensembl_transcript_id = feature_id,
         everything()) %>%
  mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                       "[:alnum:]+(?=\\.[:digit:]+)"),
         ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                             "[:alnum:]+(?=\\.[:digit:]+)")) %>%
  left_join(anno_tx) %>%
  # mutate(feature_id = paste0("Transcript ", seq(1, nrow(.), 1))) %>%
  pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
               names_to = "sample_name",
               values_to = "count") %>%
  group_by(sample_name) %>%
  filter(median(count) > 0) %>%
  mutate(count = count / sum(count)) %>%
  inner_join(colData(dds), copy = TRUE) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"))

p_stxbp1 <- data %>%
  ggplot(aes_string(x = "fraction", 
                    y = "count", 
                    fill = "fraction")) +
  geom_boxplot() +
  geom_quasirandom(size = 0.5, shape = 21) +
  # stat_pvalue_manual(stat.test, label = "ast", y.position = c(1, 1)) +
  labs(x = "Fraction", y = "Proportion of \nTotal Expression", title = "Stxbp1") +
  theme_cowplot(6) +
  theme(legend.position = "none") +
  facet_wrap(vars(external_transcript_name), nrow = 1) +
  panel_border() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_d3()

library(ggpubr)
p_splicing <- ggarrange(
  {ggarrange(p_snca, 
          p_stxbp1, 
          nrow = 1)}, 
  p_casr, 
  nrow = 2)

ggsave(
       filename = "output/figures/supplementary/splicing_examples.png", 
       plot = p_splicing, 
       # bg = "transparent", 
       width = 130, height = 80,
       units = "mm",
       dpi = 300)

# specificity index ----

# pSI_1.1.tar.gz copied to local cache
renv::install("pSI")
# Mouse to Human expression similarity: TRAP mouse vs sorted human cell equivalent
library(pSI)

# GABA (Ouwenga) ----

files <- list.files(
  "input/trap/gwas_priority/count_data/gabaergic/ouwenga",
  pattern = ".txt",
  recursive = T,
  full.names = T)

gabaergic_counts <- lapply(
  files, 
  function(x) {
    name <- str_extract(x, "GSM[:digit:]+")
    
    read_delim(x, 
               delim = "\t") %>%
      select(ensembl_gene_id = Feature, 
             !!quo_name(name) := Count)
  }
) %>%
  purrr::reduce(inner_join, by = "ensembl_gene_id")
gabaergic_counts0 <- as.matrix(gabaergic_counts[,-1])
rownames(gabaergic_counts0) <- gabaergic_counts$ensembl_gene_id
gabaergic_counts <- gabaergic_counts0
rm(gabaergic_counts0)

# GABA (Paul) ----
# Excitatory (Slc17a6 + (glutamatergic))

files <- list.files(
  "input/trap/gwas_priority/count_data/glutamatergic",
  pattern = ".txt",
  recursive = T,
  full.names = T)

glutamatergic_counts <- lapply(
  files, 
  function(x) {
    name <- str_extract(x, "GSM[:digit:]+")
    
    read_delim(x, 
               delim = "\t", 
               col_names = c("ensembl_gene_id", 
                             name))
  }
) %>%
  purrr::reduce(inner_join, by = "ensembl_gene_id")
glutamatergic_counts0 <- as.matrix(glutamatergic_counts[,-1])
rownames(glutamatergic_counts0) <- glutamatergic_counts$ensembl_gene_id
glutamatergic_counts <- glutamatergic_counts0
rm(glutamatergic_counts0)

# Oligodendrocytes (Voskuhl) ----

oligodendrocytes_counts <- read_delim("input/trap/gwas_priority/count_data/oligodendrocytic/GSE118451_Olig1.RiboTag_rawcounts.txt", 
                                      delim = "\t")
oligodendrocytes_counts0 <- as.matrix(oligodendrocytes_counts[,-1])
rownames(oligodendrocytes_counts0) <- oligodendrocytes_counts$...1
oligodendrocytes_counts <- oligodendrocytes_counts0
rm(oligodendrocytes_counts0)

# Microglia (Vasek) ----

files <- list.files(
  "input/trap/gwas_priority/count_data/microglial",
  pattern = ".txt",
  recursive = T,
  full.names = T)

microglia_counts <- lapply(
  files, 
  function(x) {
    name <- str_extract(x, "GSM[:digit:]+")
    
    read_delim(x, 
               delim = "\t", 
               col_names = c("ensembl_gene_id", 
                             name))
  }
) %>%
  purrr::reduce(inner_join, by = "ensembl_gene_id")
microglia_counts0 <- as.matrix(microglia_counts[,-1])
rownames(microglia_counts0) <- microglia_counts$ensembl_gene_id
microglia_counts <- microglia_counts0
rm(microglia_counts0)

# Astrocytes (Sakers) ----

files <- list.files(
  "input/trap/gwas_priority/count_data/astrocytic/",
  pattern = ".txt",
  recursive = T,
  full.names = T)

astrocytes_counts <- lapply(
  files, 
  function(x) {
    name <- str_extract(x, "GSM[:digit:]+")
    
    read_delim(x, 
               delim = "\t", 
               col_names = c("ensembl_gene_id", 
                             name))
  }
) %>%
  purrr::reduce(inner_join, by = "ensembl_gene_id")
astrocytes_counts0 <- as.matrix(astrocytes_counts[,-1])
rownames(astrocytes_counts0) <- astrocytes_counts$ensembl_gene_id
astrocytes_counts <- astrocytes_counts0
rm(astrocytes_counts0)

# Endothelia ----
# Pericytes ----
# Ependymal ----

# DA TRAP ----
dds_TRAP <- readRDS("input/trap/dds.rds")
dds_C1_MB_IP <- dds_TRAP[,colData(dds_TRAP)$region == "MB" &
                           colData(dds_TRAP)$fraction == "IP" &
                           colData(dds_TRAP)$cohort == "C1"]

# Combine rowMeans of each cell type ----

common_genes <- Reduce(intersect, list(rownames(dds_C1_MB_IP), 
                                       rownames(gabaergic_counts), 
                                       rownames(glutamatergic_counts), 
                                       rownames(oligodendrocytes_counts), 
                                       rownames(microglia_counts), 
                                       rownames(astrocytes_counts)))


# test <- glutamatergic_counts[which(common_genes %in% rownames(glutamatergic_counts)),]
# test <- test[sort(rownames(test)),]
# test
# glutamatergic_counts[match(common_genes %in% rownames(glutamatergic_counts)),]


cell_type_means <- tibble(ensembl_gene_id = rownames(dds_C1_MB_IP),
                          dopaminergic = round(rowMeans(counts(dds_C1_MB_IP)))) %>%
  left_join({
    tibble(ensembl_gene_id = rownames(gabaergic_counts),
           gabaergic = round(rowMeans(gabaergic_counts)))
  }) %>%
  left_join({
    tibble(ensembl_gene_id = rownames(glutamatergic_counts),
           glutamatergic = round(rowMeans(glutamatergic_counts)))
  }) %>%
  left_join({
    tibble(ensembl_gene_id = rownames(oligodendrocytes_counts),
           oligodendrocytes = round(rowMeans(oligodendrocytes_counts)))
  }) %>%
  left_join({
    tibble(ensembl_gene_id = rownames(microglia_counts),
           microglia = round(rowMeans(microglia_counts)))
  }) %>%
  left_join({
    tibble(ensembl_gene_id = rownames(astrocytes_counts),
           astrocytes = round(rowMeans(astrocytes_counts)))
  })

cell_type_means0 <- as.matrix(cell_type_means[, -1])
rownames(cell_type_means0) <- cell_type_means$ensembl_gene_id
cell_type_means <- cell_type_means0
rm(cell_type_means0)
cell_type_means[is.na(cell_type_means)] <- 0

cell_types_metadata <- tibble(sample_name = colnames(cell_type_means),
                              cell_type = colnames(cell_type_means))
cell_types_metadata <-
  data.frame(cell_type = cell_types_metadata$cell_type,
             row.names = cell_types_metadata$sample_name)

cell_types_dds <- DESeqDataSetFromMatrix(cell_type_means,
                                         cell_types_metadata,
                                         design = ~ 1)
cell_types_dds <- DESeq(cell_types_dds)

sizeFactors(cell_types_dds)

# Specificity Index ----

pSI.in <- DESeq2::counts(cell_types_dds,
                         normalized = T)
pSI.in <- pSI.in[rowMaxs(pSI.in) > 10,]
pSI.in <- as.data.frame(pSI.in, 
                        rownames = rownames(pSI.in))

pSI <- pSI::specificity.index(pSI.in = pSI.in)

pSI_dopaminergic <- pSI %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  # left_join(anno) %>% 
  filter(!is.na(dopaminergic)) %>% 
  arrange(dopaminergic)

pSI_gabaergic <- pSI %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  # left_join(anno) %>% 
  filter(!is.na(gabaergic)) %>% 
  arrange(gabaergic)

# Plot the number of specific genes identifed per cell type ----
enframe(colSums(!is.na(pSI))) %>%
  mutate(name = str_to_title(name)) %>%
  ggplot(aes(x = name, 
             y = value, 
             fill = name)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  theme_cowplot() +
  scale_y_continuous(expand = c(0, 0.05)) +
  labs(x = "Cell Type", 
       y = "Number of Specific Genes") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1))

# make a heatmap of specificity ----

specific_genes <- lapply(colnames(pSI), function(x){
  tib <- pSI %>%
    as_tibble(rownames = "ensembl_gene_id")
  # left_join(anno) %>% 
  tib <- filter(tib, !is.na(tib[[x]]))
  tib <- arrange(tib, tib[[x]])
  tib <- tib %>% select(c(ensembl_gene_id, x)) %>%
    slice_head(n = 500) %>%
    pull(ensembl_gene_id)
})
names(specific_genes) <- colnames(pSI)

specific_genes <- lapply(names(specific_genes), function(n){
  tibble(specificity = n, 
         ensembl_gene_id = specific_genes[[n]])
}) %>%
  bind_rows() %>%
  group_by(specificity) %>%
  mutate(group_id = 1000*cur_group_id()) %>%
  mutate(gene_order = group_id + row_number()) %>%
  ungroup %>%
  mutate(gene_order = rank(gene_order))

library(pheatmap)
library(RColorBrewer)
library(viridis)
pSI.in_heatmap <- pSI.in
colnames(pSI.in_heatmap) <- str_to_sentence(colnames(pSI.in_heatmap))
colnames(pSI.in_heatmap)[colnames(pSI.in_heatmap) == "Gabaergic"] <- "GABAergic"

pSI.in_heatmap[unique(Reduce(c, specific_genes$ensembl_gene_id)), ] %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id, 
               names_to = "cell_type", 
               values_to = "expr") %>%
  left_join(specific_genes) %>%
  # arrange(desc(expr), specificity) %>%
  # mutate(cell_type = factor(cell_type, levels = sort(unique(cell_type))), 
  #        ensembl_gene_id = factor(ensembl_gene_id)) %>%
  group_by(ensembl_gene_id) %>%
  mutate(expr = scale(expr, center = T)) %>% 
  ggplot(aes(x = gene_order, 
             y = cell_type, 
             fill = expr)) +
  geom_tile(colour = NA, 
            size = 1) +
  theme_cowplot(8) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0)), 
                     breaks = c(seq(0, 3000, 500))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_gradientn(values = scales::rescale(c(-2, -1, 0, 1, 2)), 
                       colours = c("#03071E", "#370617", "#6A040F", "#DC2F02", "#FFBA08")) +
  theme(axis.line = element_blank()) +
  labs(x = "Gene", 
       y = "", 
       fill = "Scaled\nExpression")
ggsave("output/plots/figure_3/psi_heatmap.png", 
       bg = "transparent", 
       width = 120, height = 45,
       units = "mm",
       dpi = 300)

# heatmap <- pSI.in_heatmap[unique(Reduce(c, specific_genes)), ] %>% t %>% pheatmap::pheatmap(scale = "column", 
#                                                                                             border_color = "grey", 
#                                                                                             cluster_rows = F, 
#                                                                                             cellheight = 15, 
#                                                                                             cellwidth = 0.1,
#                                                                                             show_colnames = F, 
#                                                                                             cluster_cols = F,
#                                                                                             color = inferno(10),
#                                                                                             gaps_row = c(1, 2, 3, 4, 5, 6), 
#                                                                                             gaps_col = seq(500, 3000, 500),
#                                                                                             cutree_cols = 6, 
#                                                                                             fontsize = 14, 
#                                                                                             legend_labels = "title")
# save_pheatmap_png <- function(x, filename, width=3600, height=3000, res = 300) {
#   png(filename, width = width, height = height, res = res)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# 
# save_pheatmap_png(heatmap, "/active/paper_md/input/figures/fig3/psi_heatmap.png")
# 
# dev.off()

# ----
# ----
# ----
# GWAS - Nalls 2019 ----

anno_human <- readRDS("/active/paper/input/trap/anno_human.rds")

# find all PD studies
# pd_gwas <- get_studies(efo_id = 'EFO_0002508')

# get variants
# pd_gwas_variants <- get_variants(study_id = "GCST009325") # Nalls, 2019
# pd_gwas_variants <- get_variants(study_id = "GCST004902") # Chang, 2017
# pd_gwas_variants <- pd_gwas_variants@variants
# pd_gwas_variants %>% View

# load Nalls summary statistics
nalls <- read_excel("input/trap/gwas_priority/Table S2. Detailed summary statistics on all nominated risk variants, known and novel_.xlsx")
nalls <- nalls %>% 
  filter(str_detect(SNP, "^rs") &
           !str_detect(CHR, "Signif|Hardy|Average"))

# simplify the nalls results ----
nalls <- nalls %>%
  select(SNP,
         nearest_select = "Nearest Gene",
         qtl_select = "QTL Nominated Gene (nearest QTL)")

# GET VARIANTS AND GENES WITHIN 500kb WINDOW WITH LD +0.5

# get_variant_genes_and_consequences <- function(lead_variant){
# 
#   # SET SERVER
#   server <- "https://rest.ensembl.org"
# 
#   # GET LD VARIANTS +0.5 LD
#    ext <- paste0("/ld/human/",
#                   lead_variant,
#                   "/1000GENOMES:phase_3:GBR?r2=0.5")
#     r <-
#       GET(paste(server, ext, sep = ""),
#           content_type("application/json"))
#     stop_for_status(r)
# 
#     linked_variants <- c(lead_variant,
#       unlist(fromJSON(toJSON(httr::content(
#         r
#       )))$variation2))
# 
#     print(paste("There are", length(linked_variants), "variants"))
# 
#     linked_variants <- split(linked_variants, ceiling(seq_along(linked_variants)/50))
# 
#     lapply(linked_variants, function(split_variants){
# 
#       # GET CONSEQUENCES
#       ext <- "/vep/human/id"
#       r <- POST(
#         paste(server, ext, sep = ""),
#         content_type("application/json"),
#         accept("application/json"),
#         body = paste0("{ \"ids\" : [\"",
#                       paste0(split_variants, sep = "", collapse = "\", \""),
#                       "\" ] }")
#       )
#       stop_for_status(r)
# 
#       variants_consequences <- fromJSON(toJSON(httr::content(r)))
# 
#       if("transcript_consequences" %in% colnames(variants_consequences)){
#         variants_consequences <- variants_consequences %>%
#           as_tibble %>%
#           select(id,
#                  most_severe_consequence,
#                  transcript_consequences) %>%
#           unnest(c(id, most_severe_consequence, transcript_consequences)) %>%
#           select(id,
#                  most_severe_consequence,
#                  gene_id,
#                  biotype) %>%
#           unnest(c(gene_id,
#                    biotype)) %>%
#           distinct() %>%
#           select(linked_variant = id,
#                  ensembl_gene_id = gene_id,
#                  biotype,
#                  consequence = most_severe_consequence
#           ) %>%
#           mutate(.before = everything(),
#                  lead_variant = lead_variant)
# 
#         print("Round Completed")
# 
#         return(variants_consequences)
# 
#       } else {
#         variants_consequences <- variants_consequences %>%
#           as_tibble %>%
#           select(id) %>%
#           unnest(id) %>%
#           distinct() %>%
#           select(linked_variant = id) %>%
#           mutate(ensembl_gene_id = "intergenic",
#                  biotype = "intergenic",
#                  consequence = "intergenic"
#           ) %>%
#           mutate(.before = everything(),
#                  lead_variant = lead_variant) %>%
#           distinct()
# 
#         print("Round Completed")
# 
#         return(variants_consequences)
# 
#       }
# 
#     }) %>%
#       bind_rows()
# 
# }
# 
# nalls_genes <- lapply(nalls$SNP,
#        function(x){
#          print(x)
#          get_variant_genes_and_consequences(x) %>%
#            left_join(anno_human,
#                      by = c("ensembl_gene_id" = "hsapiens_homolog_ensembl_gene"))
#        }) %>%
#   bind_rows()

# use genes from NALLS: +/- 1Mb, r2 > 0.5
nalls_genes <- readxl::read_excel("input/trap/gwas_priority/nalls_genes.xlsx") %>%
  select("external_gene_name" = "Nearest gene")

# add mouse homologs
homologs_ncbi <- nalls_genes %>%
  left_join({
    homologene(nalls_genes$external_gene_name, inTax = 9606, outTax = 10090) %>%
      select(external_gene_name = "9606",
             mmusculus_homolog_associated_gene_name = "10090")
  })

# set up marts for human genes and getting SNP coordinates
ensembl_human <- useMart("ensembl",
                         dataset="hsapiens_gene_ensembl")
ensembl_variation <- useMart("ENSEMBL_MART_SNP",
                             dataset = "hsapiens_snp")
homologs_ensembl <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "mmusculus_homolog_associated_gene_name"
  ),
  filters = "external_gene_name",
  values = nalls_genes$external_gene_name,
  mart = ensembl_human,
  uniqueRows = TRUE
)

# homologs_ensembl <- readRDS("R/objects/homologs_ensembl.rds")

all(nalls_genes$external_gene_name %in% homologs_ncbi$external_gene_name)
all(nalls_genes$external_gene_name %in% homologs_ensembl$external_gene_name)

# final homolog list
homologs <-left_join(homologs_ncbi,
                     homologs_ensembl,
                     by = "external_gene_name") %>%
  filter(!is.na(ensembl_gene_id)) %>%
  mutate(mmusculus_homolog_associated_gene_name.x = ifelse(
    is.na(mmusculus_homolog_associated_gene_name.x),
    mmusculus_homolog_associated_gene_name.y,
    mmusculus_homolog_associated_gene_name.x
  )) %>%
  select(ensembl_gene_id,
         external_gene_name,
         mmusculus_homolog_associated_gene_name = mmusculus_homolog_associated_gene_name.x) %>%
  filter(mmusculus_homolog_associated_gene_name != "") %>%
  distinct()

homologs %>% saveRDS("input/trap/homologs_gwas.rds")

# these humans genes are not in the database/data
nalls_genes$external_gene_name[!nalls_genes$external_gene_name %in% homologs$external_gene_name]

# add TRAP to nalls -----
# add MB and AXON FRACTION results
AXON_FRACTION_META <- readRDS("input/trap/AXON_FRACTION_META.rds")

nalls_trap <- homologs %>%
  left_join({
    MB_FRACTION_META %>%
      select(external_gene_name,
             log2FoldChange_MB = log2FoldChange,
             sumz_adj_MB = sumz_adj
      )},
    by = c("mmusculus_homolog_associated_gene_name" = "external_gene_name")) %>%
  left_join({
    AXON_FRACTION_META %>%
      select(external_gene_name,
             log2FoldChange_AXON = log2FoldChange,
             sumz_adj_AXON = sumz_adj
      )},
    by = c("mmusculus_homolog_associated_gene_name" = "external_gene_name")) %>%
  filter(mmusculus_homolog_associated_gene_name != "Ccdc62")

# remove rows with NA for both MB and AXON
nalls_trap <- nalls_trap[!c(is.na(nalls_trap$log2FoldChange_MB) &
                              is.na(nalls_trap$log2FoldChange_AXON)),] %>%
  distinct()

# get Nalls loci and genes within bounds -----
# GET NALLS SNP COORDINATES
# SET SERVER
# server <- "https://rest.ensembl.org"
# ext <- "/vep/human/id"
# r <- POST(
#   paste(server, ext, sep = ""),
#   content_type("application/json"),
#   accept("application/json"),
#   body = paste0("{ \"ids\" : [\"",
#                 paste0(nalls$SNP, sep = "", collapse = "\", \""),
#                 "\" ] }")
# )
# stop_for_status(r)

# nalls_consequences <- fromJSON(toJSON(httr::content(r)))
# nalls_consequences %>% saveRDS("input/trap/gwas_priority/nalls_consequences.rds")
nalls_consequences <- readRDS("input/trap/gwas_priority/nalls_consequences.rds")

nalls_loci <- nalls_consequences %>%
  as_tibble %>%
  unnest(c(id,
           seq_region_name,
           start)) %>%
  select(SNP = id,
         chr = seq_region_name,
         start = start) %>%
  filter(!str_detect(chr, "CHR")) %>%
  distinct() %>%
  rowwise() %>%
  mutate(end = start + 1e6,
         start = max(start - 1e6, 0))

# # GET GENES WITHIN +-1MB OF NALLS SNPS
# ext <- paste0("/overlap/region/human/",
#               paste0(nalls_loci$chr,
#                      ":",
#                      nalls_loci$start,
#                      "-",
#                      nalls_loci$end,
#                      "?feature=gene"))
# 
# nalls_all_genes <- lapply(ext, function(x){
# 
#   r <- GET(paste(server, x, sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   fromJSON(toJSON(httr::content(r))) %>%
#     pull(id) %>%
#     unlist()
# 
# })
# 
# names(nalls_all_genes) <- nalls_loci$SNP
# 
# nalls_all_genes <- nalls_all_genes %>%
#   enframe %>%
#   unnest(value) %>%
#   select(SNP = name,
#          hsapiens_homolog_ensembl_gene = value)

# nalls_all_genes %>% saveRDS("input/trap/gwas_priority/nalls_all_genes.rds")
nalls_all_genes <- readRDS("input/trap/gwas_priority/nalls_all_genes.rds")


# make some summary numbers ----
# total number of snps
n_nalls_snps <- length(unique(nalls$SNP))
# total number of qtls selected by nalls
n_nalls_qtl_selected <- length(unique(nalls$qtl_select[!is.na(nalls$qtl_select)]))
# total number of genes selected by nalls
n_nalls_genes_selected <- length(unique(nalls$nearest_select))
# all homologs measured in trap
n_nalls_homologs_values <- nalls_trap %>%
  pull(mmusculus_homolog_associated_gene_name) %>%
  unique
# total number of homologs measured in trap
n_nalls_homologs <- length(n_nalls_homologs_values)

# biotypes ----

# biotypes of all genes measured in nalls
# n_nalls_genes_biotypes <- getBM(
#   attributes = c("external_gene_name",
#                  "gene_biotype"),
#   filters = "external_gene_name",
#   values = unique(nalls_genes$external_gene_name),
#   mart = ensembl_human,
#   uniqueRows = TRUE
# ) %>%
# filter(!c(external_gene_name == "ZNRD1ASP" &
#            str_detect(gene_biotype, "nitary")))

# saveRDS(n_nalls_genes_biotypes, "input/trap/gwas_priority/n_nalls_genes_biotypes.rds")
n_nalls_genes_biotypes <- readRDS("input/trap/gwas_priority/n_nalls_genes_biotypes.rds")

n_nalls_genes <- length(unique(n_nalls_genes_biotypes$external_gene_name))

# make a gene summary table ----

# summary tibble
plot_data <- tibble(
  external_gene_name = nalls_genes$external_gene_name,
  "Nalls, 2019" = TRUE,
  "Mouse\nHomolog" = external_gene_name %in% homologs$external_gene_name,
  "Detected\nin TRAP" = external_gene_name %in% nalls_trap$external_gene_name,
  "Closest Gene" = external_gene_name %in% nalls$nearest_select,
  "QTL" = external_gene_name %in% nalls$qtl_select
) %>%
  left_join(n_nalls_genes_biotypes) %>%
  pivot_longer(-c(external_gene_name, gene_biotype),
               names_to = "category",
               values_to = "logical") %>%
  mutate(category = factor(category,
                           levels = c("Nalls, 2019", "Mouse\nHomolog",
                                      "Detected\nin TRAP", "Closest Gene",
                                      "QTL")))

# plot the number of homologs detected ----
p_gwas_homologs_trap_detected <- plot_data %>%
  filter(category %in% c("Nalls, 2019", "Mouse\nHomolog", "Detected\nin TRAP")) %>%
  filter(logical) %>%
  ggplot(aes(x = category,
             fill = category)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  labs(y = "Number of Genes",
       x = "Gene Conversion") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(y = ..count.., label = ..count..),
            stat = "count", vjust = -0.5) +
  theme(legend.position = "none")

# plot the missing outcomes ----
p_gwas_missing_outcomes <- plot_data %>%
  pivot_wider(names_from = category,
              values_from = logical) %>%
  filter(!`Detected\nin TRAP`) %>%
  pivot_longer(-c(external_gene_name, gene_biotype),
               names_to = "category",
               values_to = "logical") %>%
  filter(category %in% c("Nalls, 2019", "Closest Gene", "QTL")) %>%
  filter(logical) %>%
  mutate(category = factor(ifelse(category == "Nalls, 2019",
                                  "No Association",
                                  category),
                           levels = c("No Association",
                                      "Closest Gene",
                                      "QTL"))) %>%
  ggplot(aes(x = category,
             fill = category)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  geom_text(aes(y = ..count.., label = ..count..),
            stat = "count", vjust = -0.5) +
  theme(legend.position = "none") +
  labs(y = "Number of Genes",
       x = "Outcome") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


# plot the biotype of detected and undetected ----
# remove underscores from titles function for plotting
remove_underscores <- function(colnames) {
  gsub("_", " ", colnames)
}
p_gwas_biotype_breakdown <- plot_data %>%
  pivot_wider(names_from = "category",
              values_from = "logical") %>%
  mutate(gene_biotype = ifelse(str_detect(gene_biotype, "pseudogene"),
                               "pseudogene", gene_biotype),
         gene_biotype = factor(remove_underscores(replace_na(gene_biotype, "No annotation"))),
         gene_biotype = relevel(gene_biotype, ref = "protein coding"),
         "Detected\nin TRAP" = ifelse(`Detected\nin TRAP`, "Detected", "Not Detected")) %>%
  group_by(`Detected\nin TRAP`) %>%
  mutate(n = n()) %>%
  ggplot(aes(x = `Detected\nin TRAP`,
             fill = gene_biotype)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  theme(legend.position = c(0.8, 0.75)) +
  labs(y = "Number of Genes",
       x = "Detected in TRAP",
       fill = "Gene Biotype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(y = n, label = n),
            stat = "identity", vjust = -0.5, check_overlap = T)

# table the missing associated genes ----
missing_associated_genes <- plot_data %>%
  pivot_wider(names_from = "category",
              values_from = "logical") %>%
  filter(!`Detected\nin TRAP`) %>%
  pull(external_gene_name)

t_gwas_nalls_missing_assoc <- nalls %>%
  filter(nearest_select %in% missing_associated_genes |
           qtl_select %in% missing_associated_genes)

# combine the datasets ----
nalls_trap_combined <- nalls_trap %>%
  left_join(nalls_all_genes,
            by = c("ensembl_gene_id" = "hsapiens_homolog_ensembl_gene")) %>% # add snp info
  select(
    SNP,
    ensembl_gene_id,
    mmusculus_homolog_associated_gene_name,
    log2FoldChange_MB,
    sumz_adj_MB,
    log2FoldChange_AXON,
    sumz_adj_AXON
  ) %>%
  left_join(nalls) %>%
  select(SNP,
         nearest_select,
         qtl_select,
         everything())

# confirm no genes are lost by merging
length(unique(nalls_trap$ensembl_gene_id)) == length(unique(nalls_trap_combined$ensembl_gene_id))

# label the snps with no significant enrichment ----
aba_expression_axon <- readRDS("input/trap/aba_expression_axon.rds")

nalls_trap_axon_specific <- AXON_FRACTION_META %>%
  left_join(aba_expression_axon) %>%
  filter(external_gene_name %in% nalls_trap_combined$mmusculus_homolog_associated_gene_name) %>%
  select(external_gene_name,
         expression_energy) %>%
  filter(expression_energy < 0.5) %>% # not currently implemented
  pull(external_gene_name)

nalls_trap_combined <- nalls_trap_combined %>%
  rowwise() %>%
  mutate(
    # nonsig = min(sumz_adj_MB, sumz_adj_AXON) > 0.01,
    non_dopa_MB = sumz_adj_MB < 0.01 & log2FoldChange_MB < 0 |
      sumz_adj_MB > 0.01 | is.na(sumz_adj_MB),
    non_dopa_AXON = sumz_adj_AXON < 0.01 & log2FoldChange_AXON < 0 |
      sumz_adj_AXON > 0.01 | is.na(sumz_adj_AXON) | !mmusculus_homolog_associated_gene_name %in% nalls_trap_axon_specific,
    non_dopa = ifelse(non_dopa_MB & non_dopa_AXON,
                      TRUE, FALSE))

# Categorise Nalls selections by dopa status ----
MB_FRACTION_TRANSLATED_GENES <- readRDS("input/trap/MB_FRACTION_TRANSLATED_GENES.rds")

nalls_selection_dopa_status <-
  tibble(external_gene_name = unique(nalls$nearest_select, nalls$qtl_select)) %>%
  left_join(homologs) %>%
  select(
    external_gene_name = "mmusculus_homolog_associated_gene_name",
    hsapiens_homolog_ensembl_gene = ensembl_gene_id,
    hsapiens_homolog_associated_gene_name = external_gene_name
  ) %>%
  left_join(anno_human) %>%
  drop_na() %>%
  mutate(translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES)

nalls_selection_dopa_status <- nalls_selection_dopa_status %>%
  filter(translated)

# label the single gene snps ----
nalls_trap_combined <- nalls_trap_combined %>%
  group_by(SNP) %>%
  mutate(single = n() == 1) %>%
  ungroup()

# plot the single gene multi gene breakdown with non dopa information ----
plot_data <- nalls_trap_combined %>%
  select(
    SNP,
    single,
    non_dopa
  ) %>%
  group_by(SNP) %>%
  summarise(single = sum(single) == n(),
            non_dopa = sum(non_dopa) == n()) %>%
  filter(!is.na(SNP))

p_gwas_single_dopa_breakdown <- plot_data %>%
  mutate(multigene = !single) %>%
  pivot_longer(-c(SNP, non_dopa),
               names_to = "category",
               values_to = "logical") %>%
  filter(logical) %>%
  mutate(category = str_to_sentence(category),
         non_dopa = ifelse(non_dopa, "Not Dopaminergic", "Dopaminergic")) %>%
  ggplot(aes(category,
             fill = non_dopa)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  labs(y = "Number of loci",
       x = "Genes per SNP",
       fill = "Not Dopaminergic") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ..count..),
             position=position_stack(vjust=0.5))


# list of single gene snps ----
nalls_single_snps <- nalls_trap_combined %>%
  filter(single) %>%
  pull(SNP) %>%
  unique

# tabulate the final combined list ----

# remove single gene loci
gwas <- nalls_trap_combined %>%
  filter(!single) %>%
  select(SNP,
         mmusculus_homolog_associated_gene_name,
         log2FoldChange_MB,
         log2FoldChange_AXON,
         non_dopa,
         non_dopa_MB,
         non_dopa_AXON
  ) %>%
  distinct() %>%
  filter(!is.na(SNP))

n_gwas_snps <- length(unique(gwas$SNP))

n_gwas_nondopa <- gwas %>%
  group_by(SNP) %>%
  summarise(non_dopa = sum(non_dopa) == n()) %>%
  filter(non_dopa) %>%
  pull(SNP) %>%
  length

gwas <- gwas %>%
  filter(!non_dopa)

# gwas_single_dopa <- gwas %>%
#   select(SNP,
#          mmusculus_homolog_associated_gene_name) %>%
#   group_by(SNP) %>%
#   filter(n() == 1)
#
# gwas_mb_axon_agree <- gwas %>%
#   filter(!non_dopa) %>%
#   filter(!SNP %in% gwas_single_dopa$SNP) %>%
#   group_by(SNP) %>%
#   arrange(desc(log2FoldChange_MB)) %>%
#   mutate(rank_MB = row_number()) %>%
#   arrange(desc(log2FoldChange_AXON)) %>%
#   mutate(rank_AXON = row_number()) %>%
#   filter(rank_MB == 1 &
#            rank_AXON == 1)
#
# gwas %>%
#   filter(!non_dopa) %>%
#   filter(!SNP %in% gwas_single_dopa$SNP) %>%
#   filter(!SNP %in% gwas_mb_axon_agree$SNP) %>%
#   group_by(SNP) %>%
#   slice_max(order_by = log2FoldChange_MB,
#             n = 1) %>% View
#   arrange(desc())

# rank based on MB enrichment
gwas_mb <- gwas %>%
  filter(!non_dopa_MB) %>%
  group_by(SNP) %>%
  slice_max(order_by = log2FoldChange_MB,
            n = 1) %>%
  select(SNP,
         gene = mmusculus_homolog_associated_gene_name,
         lfc = log2FoldChange_MB)
# number of SNPS able to be ranked by MB data
gwas_mb$SNP %>% unique

# rank based on axonal enrichment
gwas_axon <- gwas %>%
  filter(!non_dopa_AXON) %>%
  group_by(SNP) %>%
  slice_max(order_by = log2FoldChange_AXON,
            n = 1) %>%
  select(SNP,
         gene = mmusculus_homolog_associated_gene_name,
         lfc = log2FoldChange_AXON)
# inner_join({
#   AXON_FRACTION_META %>%
#     select(external_gene_name,
#            expression_energy)
# },
# by = c("gene" = "external_gene_name")) %>%
# filter(expression_energy < 0.5) %>%
# select(-expression_energy)

# number of SNPS able to be ranked by MB data
gwas_axon$SNP %>% unique

t_gwas <- gwas_mb %>%
  full_join({
    gwas_axon
  },
  by = "SNP",
  suffix = c("_MB", "_AXON")) %>%
  mutate(lfc_MB = signif(lfc_MB, 3)) %>%
  mutate(lfc_AXON = signif(lfc_MB, 3)) %>%
  mutate(across(everything(), ~as.character(.))) %>%
  mutate(across(everything(), ~replace_na(., replace = ""))) %>%
  arrange(desc(lfc_MB), desc(lfc_AXON))

t_gwas_nalls_compare <- t_gwas %>%
  left_join(nalls) %>%
  select(-c(lfc_MB, lfc_AXON)) %>%
  mutate(across(everything(), ~replace_na(., "")))


# the number of gwas loci selected by TRAP ----
n_gwas_trap_loci <- nrow(t_gwas)

# plot soma axon driver ----
# Only multigene loci with at least one dopaminergic gene are included

plot_data <- t_gwas %>%
  mutate(agreement = gene_MB == gene_AXON,
         driver = ifelse(lfc_AXON == "",
                         "Soma",
                         ifelse(lfc_MB == "",
                                "Axon",
                                ifelse(as.numeric(lfc_MB) > as.numeric(lfc_AXON),
                                       "Soma", "Axon"))),
         trap_gene = ifelse(driver == "Soma", gene_MB, gene_AXON))

p_gwas_driving_compartment <- plot_data %>%
  ggplot(aes(x = driver,
             fill = driver)) +
  geom_bar(colour = "black") +
  scale_fill_d3()  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ..count..),
             position=position_stack(vjust=0.5)) +
  labs(x = "Driving Compartment",
       y = "Number of Loci") +
  theme(legend.position = "none")

# plot compartment agreement ----
p_gwas_compartment_agreement <- plot_data %>%
  filter(gene_MB != "" & gene_AXON != "") %>%
  mutate(agreement = ifelse(agreement, "Same Gene", "Different Gene")) %>%
  ggplot(aes(x = agreement,
             fill = agreement)) +
  geom_bar(colour = "black") +
  scale_fill_d3()  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ..count..),
             position=position_stack(vjust=0.5)) +
  labs(x = "Compartment Agreement",
       y = "Number of Loci") +
  theme(legend.position = "none")


# plot nalls trap agreement ----

gwas_trap_genes_human <- homologs %>%
  filter(mmusculus_homolog_associated_gene_name %in% plot_data$trap_gene) %>%
  pull(external_gene_name)

p_gwas_trap_nearest_qtl_compare <- tibble(trap_gene = gwas_trap_genes_human) %>%
  mutate(nearest_select = trap_gene %in% nalls$nearest_select,
         qtl_select = trap_gene %in% nalls$qtl_select) %>%
  pivot_longer(-trap_gene,
               names_to = "selection",
               values_to = "common") %>%
  mutate(common = factor(ifelse(common, "Common", "Different"),
                         levels = c("Different", "Common")),
         selection = ifelse(selection == "nearest_select",
                            "Nearest Gene",
                            "QTL")) %>%
  ggplot(aes(x = selection,
             fill = common)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ..count..),
             position=position_stack(vjust=0.5)) +
  labs(x = "Nalls Selection",
       y = "Number of Loci",
       fill = "TRAP Selection") +
  theme(legend.position = "top")

plot_data <- nalls %>%
  filter(SNP %in% t_gwas$SNP) %>%
  pivot_longer(-SNP,
               names_to = "selection",
               values_to = "external_gene_name") %>%
  left_join(plot_data) %>%
  mutate(external_gene_name = replace_na(external_gene_name, "No Selection")) %>%
  mutate(non_dopa = ifelse(external_gene_name == "No Selection",
                           "No Selection",
                           ifelse(external_gene_name %in% nalls_selection_dopa_status$hsapiens_homolog_associated_gene_name,
                                  "Enriched",
                                  "Not Enriched"))) %>%
  mutate(selection = ifelse(selection == "nearest_select",
                            "Nearest Gene", "QTL"),
         non_dopa = ifelse(non_dopa == "TRUE",
                           "Not Dopaminergic",
                           ifelse(non_dopa == "FALSE",
                                  "Dopaminergic",
                                  non_dopa)))

p_gwas_dopa_status_selection <- plot_data %>%
  ggplot(aes(x = selection,
             fill = non_dopa)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  # theme(legend.position = "top") +
  labs(y = "Number of loci",
       x = "Selection Type",
       fill = "Dopaminergic\nStatus") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ifelse(..count.. > 10, ..count.., "")),
             position=position_stack(vjust=0.5))

# find the percentage of nalls selections that were dopaminergic ----

nalls_selections <- tibble(external_gene_name = nalls %>% pull(nearest_select, qtl_select) %>% unique) %>%
  left_join(homologs)

n_nalls_selections_dopa_perc <- ((sapply(nalls_selections$mmusculus_homolog_associated_gene_name, get_ensembl_gene_id) %in% MB_FRACTION_ENRICHED_GENES)%>% sum )*100/nrow(nalls_selections)

# plot the distribution of significant enrichment/depletion among genes located ----
# within +- 1MB LD 0.5 of Nalls 2019 loci
plot_data <- nalls_trap %>%
  select(mmusculus_homolog_associated_gene_name,
         log2FoldChange = log2FoldChange_MB,
         sumz_MB = sumz_adj_MB
  ) %>%
  filter(sumz_MB < 0.01) %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rn = row_number(),
         n_enriched = sum(log2FoldChange > 0),
         n_depleted = sum(log2FoldChange < 0),
         enriched = ifelse(log2FoldChange > 0, "Enriched", "Depleted"),
         label = ifelse(rn <= 10, mmusculus_homolog_associated_gene_name, ""))

circle_data <- plot_data %>%
  filter(mmusculus_homolog_associated_gene_name %in% nalls_trap_axon_specific)

p_gwas_enrichment <- plot_data %>%
  ggplot(aes(x = rn,
             y = log2FoldChange)) +
  geom_point(aes(colour = enriched)) +
  geom_point(data = circle_data,
             size = 3,
             shape = 21) +
  geom_hline(yintercept = 0) +
  scale_color_d3() +
  geom_label_repel(aes(label = label),
                   box.padding = 0.5,
                   max.overlaps = 100,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20) +
  theme(legend.position = "none")  +
  geom_text(
    x = 100,
    y = 3,
    check_overlap = T,
    aes(label = paste(n_enriched, "genes\nenriched"),
        colour = enriched),
    size = 5,
    fontface = "bold"
  ) +
  geom_text(
    x = 100,
    y = -2,
    check_overlap = T,
    aes(label = paste(n_depleted, "genes\ndepleted")),
    colour = pal_d3()(2)[1],
    size = 5,
    fontface = "bold"
  ) +
  labs(x = "Enrichment/Depletion Rank",
       y = "Log2 Enrichment or Depletion")


# MAKE A LIST OF GWAS GENES + SINGLE ----

GWAS_GENES <- unique(
  c(
    t_gwas$gene_MB[t_gwas$gene_MB != ""],
    t_gwas$gene_AXON[t_gwas$gene_AXON != ""],
    nalls_trap_combined$mmusculus_homolog_associated_gene_name[nalls_trap_combined$single & !nalls_trap_combined$non_dopa]
  )
)

GWAS_GENES_BROAD <- nalls_trap_combined %>%
  filter(!non_dopa) %>%
  pull(mmusculus_homolog_associated_gene_name) %>%
  unique

# gene GWAS enrichment/depletion ----
# within +- 1MB LD 0.5 of Nalls 2019 loci
plot_data <- nalls_trap %>%
  select(mmusculus_homolog_associated_gene_name,
         log2FoldChange = log2FoldChange_MB,
         sumz_MB = sumz_adj_MB
  ) %>%
  filter(sumz_MB < 0.01) %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(rn = row_number(),
         n_enriched = sum(log2FoldChange > 0),
         n_depleted = sum(log2FoldChange < 0),
         enriched = ifelse(log2FoldChange > 0, "Enriched", "Depleted"),
         label = ifelse(rn <= 10, mmusculus_homolog_associated_gene_name, ""))

circle_data <- plot_data %>%
  filter(mmusculus_homolog_associated_gene_name %in% nalls_trap_axon_specific)

plot_data %>%
  ggplot(aes(x = rn,
             y = log2FoldChange)) +
  geom_point(aes(colour = enriched)) +
  geom_point(data = circle_data,
             size = 2,
             shape = 21) +
  geom_hline(yintercept = 0) +
  scale_color_d3() +
  theme_cowplot(8) +
  geom_label_repel(aes(label = label),
                   box.padding = 0.5,
                   max.overlaps = 100,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20, 
                   size = 2) +
  theme(legend.position = "none")  +
  geom_text(
    x = 100,
    y = 3,
    check_overlap = T,
    aes(label = paste(n_enriched, "genes\nenriched"),
        colour = enriched),
    # size = 5,
    fontface = "bold"
  ) +
  geom_text(
    x = 100,
    y = -2,
    check_overlap = T,
    aes(label = paste(n_depleted, "genes\ndepleted")),
    colour = pal_d3()(2)[1],
    # size = 5,
    fontface = "bold"
  ) +
  labs(x = "Enrichment/Depletion Rank",
       y = "Log2 Enrichment or Depletion")

ggsave("output/plots/figure_3/gwas_enrich_deplete.png", 
       bg = "transparent", 
       width = 90, height = 90,
       units = "mm",
       dpi = 300)

# plot R2 D' CASR exemplar ----

# GET 1000G POPULATION INFO

# library(httr)
# library(jsonlite)
# ext <- "/info/variation/populations/homo_sapiens?filter=LD"
# server <- "https://rest.ensembl.org"
# r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
# stop_for_status(r)
# populations_1000g <- fromJSON(toJSON(httr::content(r))) %>%
#   as_tibble %>%
#   unnest(everything())
# colnames(populations_1000g) <- str_to_title(colnames(populations_1000g)) 
# write_delim(populations_1000g,
#             "output/trap/populations_1000g.txt",
#             delim = "\t")
populations_1000g <- read_delim("output/trap/populations_1000g.txt")

linked_variants_and_lead_mapping <- lapply(nalls_loci$SNP[nalls_loci$SNP == "rs55961674"], function(lead_variant){
  
  # SET SERVER
  server <- "https://rest.ensembl.org"
  
  # SET LEAD VARIANT
  # lead_variant <- "rs55961674"
  
  # GET LEAD VARIANT INFO
  ext <- paste0("/variation/human/",
                lead_variant,
                "?")
  r <-
    GET(paste(server, ext, sep = ""),
        content_type("application/json"))
  stop_for_status(r)
  
  lead_variant_mapping <- tibble(lead_variant,
                                 lead_start = unlist(fromJSON(toJSON(httr::content(
                                   r
                                 )))$mappings$start))
  
  # GET LD VARIANTS +0.5 LD
  ext <- paste0("/ld/human/",
                lead_variant,
                "/1000GENOMES:phase_3:GBR")
  r <-
    GET(paste(server, ext, sep = ""),
        content_type("application/json"))
  stop_for_status(r)
  
  linked_variants <- fromJSON(toJSON(httr::content(r)))
  
  print(paste("There are", nrow(linked_variants), "variants"))
  
  return(list("linked_variants" = linked_variants, 
              "lead_variant_mapping" = lead_variant_mapping))
  
})

# for CASR
linked_variants <- linked_variants_and_lead_mapping[[1]]$linked_variants
lead_mapping <- linked_variants_and_lead_mapping[[1]]$lead_variant_mapping

loci_mapped <- nalls_loci$SNP[sapply(linked_variants, length) != 0]

# Split into groups of max 50
# if running multiple SNPs above
# linked_variants_list <-
#   split(unlist(linked_variants$variation2), ceiling(seq_along(unlist(
#     linked_variants$variation2
#   )) / 50))
# df <- linked_variants %>%
#   bind_rows()
# df <- unlist(df$variation2)
# df

# if running just a single (e.g. CASR)
linked_variants_list <-
split(unlist(linked_variants$variation2), ceiling(seq_along(unlist(
  linked_variants$variation2
)) / 50))
df <- linked_variants %>%
  bind_rows()
df <- unlist(df$variation2)
df

linked_variants_split <- split(df, (seq(length(df))-1) %/% 50) 

linked_variants_consequences <- lapply(linked_variants_split, function(x){
  # GET CONSEQUENCES (Gene info) of each variant within window
  ext <- "/vep/human/id"
  r <- POST(
    paste(server, ext, sep = ""),
    content_type("application/json"),
    accept("application/json"),
    body = paste0("{ \"ids\" : [\"",
                  paste0(x, sep = "", collapse = "\", \""),
                  "\" ] }")
  )
  
  stop_for_status(r)
  
  variants_consequences <- fromJSON(toJSON(httr::content(r))) %>%
    as_tibble() %>%
    unnest(c(seq_region_name, 
             start, 
             id, 
             allele_string, 
             assembly_name, 
             most_severe_consequence, 
             strand, 
             input, 
             end)) 
}) %>% 
  bind_rows() %>%
  mutate(transcript = transcript_consequences != "NULL", 
         regulatory = regulatory_feature_consequences != "NULL", 
         motif = motif_feature_consequences != "NULL",
         intergenic = intergenic_consequences != "NULL")

# Get the consequences for intragenic variants
linked_variants_consequences_transcripts <- linked_variants_consequences %>%
  filter(transcript) %>%
  select(id,
         start,
         most_severe_consequence,
         transcript_consequences) %>%
  unnest(c(id, most_severe_consequence, transcript_consequences)) %>%
  select(id,
         start,
         most_severe_consequence,
         gene_id,
         transcript_id,
         biotype) %>%
  unnest(c(gene_id,
           transcript_id,
           biotype)) %>%
  distinct() %>%
  select(linked_variant = id,
         start,
         ensembl_gene_id = gene_id,
         ensembl_transcript_id = transcript_id,
         biotype,
         consequence = most_severe_consequence
  )

# Get the consequences for intergenic variants
linked_variants_consequences_intergenic <- linked_variants_consequences %>%
  filter(!transcript) %>% 
  mutate(ensembl_gene_id = "intergenic",
         ensembl_transcript_id = "intergenic",
         biotype = "intergenic") %>%
  select(linked_variant = id, 
         start,
         ensembl_gene_id,
         ensembl_transcript_id,
         biotype, 
         consequence = most_severe_consequence)

# Combine linked variants consequences  
linked_variants <- linked_variants %>%
  select(lead_variant = variation1, linked_variant = variation2, r2, d_prime) %>%
  as_tibble %>%
  unnest(everything()) %>%
  left_join(lead_mapping) %>%
  left_join({
    bind_rows(linked_variants_consequences_transcripts, 
              linked_variants_consequences_intergenic)
  }) %>%
  mutate(diff = lead_start - start)

ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl",
                      version = 101)

# Get mouse homologs of linked variants
anno_homologs <- getBM(attributes = c("ensembl_gene_id",
                                      "external_gene_name",
                                      "description",
                                      "mmusculus_homolog_ensembl_gene",
                                      "mmusculus_homolog_associated_gene_name"),
                       filters = "ensembl_gene_id",
                       values = linked_variants$ensembl_gene_id,
                       mart = ensembl,
                       uniqueRows = TRUE)

plot_data <- linked_variants %>%
  left_join(anno_homologs) %>%
  select(-c(ensembl_transcript_id,
            biotype)) %>%
  distinct() %>%
  mutate(d_prime = as.double(d_prime),
         r2 = as.double(r2)) %>%
  left_join({
    MB_FRACTION_META %>%
      select(mmusculus_homolog_ensembl_gene = ensembl_gene_id,
             log2FoldChange,
             sumz_adj)
  }) %>%
  left_join({
    pSI %>%
      select(pSI = dopaminergic) %>%
      as_tibble(rownames = "ensembl_gene_id")
  }, 
  by = c("mmusculus_homolog_ensembl_gene" = "ensembl_gene_id")
  ) %>% 
  mutate(pSI = -log10(pSI)) %>%
  pivot_longer(c(d_prime, r2, log2FoldChange, pSI),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = ifelse(metric == "d_prime", 
                         "D*'\\''", 
                         ifelse(metric == "r2", 
                                "R^{2}", 
                                ifelse(metric == "log2FoldChange", 
                                       "Log[2]~Fold~Change", 
                                       "pSI"))), 
         metric = factor(metric, levels = c("R^{2}", 
                                            "D*'\\''", 
                                            "Log[2]~Fold~Change", 
                                            "pSI"))) %>%
  filter(!is.na(external_gene_name))

plot_data %>% 
  ggplot(aes(x = diff, 
             y = value, 
             colour = external_gene_name
  )) +
  geom_point(size = 0.5) +
  scale_color_d3(palette = "category20") +
  facet_wrap(vars(metric), 
             nrow = 2, 
             scales = "free_y", 
             labeller = "label_parsed") +
  geom_vline(xintercept = 0, 
             linetype = "dashed") +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  theme_cowplot(8) +
  panel_border() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Distance from Lead Variant (bp)", 
       colour = "Gene") +
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5), 
        strip.background = element_rect(colour = "black", 
                                        fill = "white"))

ggsave("output/plots/figure_3/gwas_casr.png", 
       bg = "transparent", 
       width = 90, height = 80,
       units = "mm",
       dpi = 300)


# plot stereoseq GWAS overlap ----

files <- list.files("output/markers/general", 
                    pattern = "*.csv", 
                    full.names = T)
names(files) <- str_extract(files, "(?<=general/)[:graph:]+(?=.csv)")

markers <- lapply(names(files), function(n){read_csv(files[[n]])})
names(markers) <- names(files)
markers <- bind_rows(markers, .id = "cell_type_publish")

# homologs not found in stereoseq
# homologs$mmusculus_homolog_associated_gene_name[!homologs$mmusculus_homolog_associated_gene_name %in% markers$names]

# filter for Nalls homologs
markers <- markers %>%
  filter(names %in% homologs$mmusculus_homolog_associated_gene_name)

plot_data <- markers %>% 
  filter(pvals_adj < 0.01) %>%
  group_by(cell_type_publish) %>%
  filter(n() > 5) %>%
  mutate(summary = median(logfoldchanges * -log10(pvals_adj))) %>%
  ungroup %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  mutate(cell_type_publish = factor(cell_type_publish), 
         cell_type_publish = fct_reorder(cell_type_publish, summary)) %>%
  mutate(neu_glia = ifelse(str_detect(cell_type_publish, "GLIA|MICRO|ASTRO|RBC|OLIG"), 
                                      "GLIAL", "NEURONAL"))

# intersection between DA enriched genes and glial enriched genes
DA_GWAS_SS_GENES <- plot_data %>%
  filter(logfoldchanges > 0) %>%
  filter(str_detect(cell_type_publish, "DA")) %>%
  pull(names) %>%
  unique %>%
  sort

GLIA_GWAS_SS_GENES <- plot_data %>%
  filter(logfoldchanges > 0) %>%
  filter(neu_glia == "GLIAL") %>%
  pull(names) %>%
  unique

intersect(DA_GWAS_SS_GENES, 
          GLIA_GWAS_SS_GENES)

# plot ss gwas gene enrichment
plot_data %>%
  ggplot(
    aes(
      x = cell_type_publish, 
      y = logfoldchanges
    )
  ) +
  # geom_boxplot(aes(fill = neu_glia), 
  #              outlier.size = 0.3, 
  #              lwd = 0.2) +
  # geom_violin(alpha = 0.2, 
  #             lwd = 0.2) +
  geom_jitter(
    aes(
      # size = -log10(pvals_adj),
      colour = neu_glia
    ),
    width = 0.1, 
    alpha = 0.7, 
    size = 0.3
  ) +
  # geom_point(data = plot_data %>% filter(str_detect(cell_type_publish,
  #                                                   "^DA")),
  #            aes(size = -log10(pvals_adj)),
  #            shape = 21, colour = "black") +
  # geom_label_repel(data = plot_data %>%
  #                    filter(logfoldchanges > 3),
  #   aes(label = names),
  #   size = 3
  # ) +
  # geom_text_repel(data = filter(plot_data, logfoldchanges > 3), 
  #                 aes(label = names), 
  #                 size = 3) +
  theme_cowplot(7) +
  geom_hline(yintercept = 0, 
             linetype = "dashed") +
  # coord_flip() +
  scale_color_d3() +
  scale_fill_d3() +
  labs(x = '', 
       y = expression(Log[2]~-Fold  ~ Enrichment), 
       colour = "Cell Class", 
       # title = "PD GWAS Candidate Gene Enrichment/Depletion", 
       # subtitle = "GWAS Loci from Nalls, 2019"
       ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top", 
        plot.margin = margin(7, 7, 7, 12, unit = "pt"))
  # scale_size(range = c(0.1, 0.75), 
  #            guide = "none")

ggsave("output/plots/figure_3/stereoseq_gwas.png", 
       bg = "transparent", 
       width = 110, height = 70,
       units = "mm",
       dpi = 300)

# plot specificity index GWAS genes ----

plot_data <- pSI %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  filter(external_gene_name %in% homologs$mmusculus_homolog_associated_gene_name) %>%
  pivot_longer(c(dopaminergic:astrocytes), 
               names_to = "cell_type", 
               values_to = "pSI") %>%
  mutate(pSI = replace_na(pSI, 1)) %>%
  filter(pSI != 1)
  
plot_data %>%
  filter(pSI < 0.01) %>% 
  # filter(cell_type %in% c("dopaminergic", "oligodendrocytes")) %>%
  arrange(external_gene_name) %>%
  mutate(cell_type = str_to_sentence(cell_type)) %>%
  mutate(cell_type = ifelse(cell_type == "Gabaergic", 
                            "GABAergic", cell_type)) %>%
  group_by(cell_type) %>%
  mutate(n = n()) %>%
  ggplot(aes(x = cell_type, 
             fill = cell_type)) +
  geom_bar(colour = "black") +
  theme_cowplot(7) +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), 
                     breaks = seq(0, 10, 2)) +
  # geom_text(aes(y = n, label = n),
  #           stat = "identity", vjust = -0.5, check_overlap = T, size = 3) +
  labs(x = "", 
       y = "Number of GWAS Candidate Genes\nSpecifically Expressed", 
       # title = "Specific Expression of GWAS Candidates", 
       # subtitle = "TRAP/RiboTag Datasets"
       ) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/plots/figure_3/trap_gwas.png", 
       bg = "transparent", 
       width = 60, height = 65,
       units = "mm",
       dpi = 300)

  ggplot(aes(x = cell_type,
             y = pSI)) +
  geom_boxplot() +
  geom_jitter(width = 0.2)

specific_genes <- lapply(colnames(pSI), function(x){
  tib <- pSI %>%
    as_tibble(rownames = "ensembl_gene_id")
  # left_join(anno) %>% 
  tib <- filter(tib, !is.na(tib[[x]]))
  tib <- arrange(tib, tib[[x]])
  tib <- tib %>% select(c(ensembl_gene_id, x)) %>%
    # slice_head(n = 500) %>%
    pull(ensembl_gene_id)
})
names(specific_genes) <- colnames(pSI)

specific_genes <- lapply(names(specific_genes), function(n){
  tibble(specificity = n, 
         ensembl_gene_id = specific_genes[[n]])
}) %>%
  bind_rows() %>%
  group_by(specificity) %>%
  mutate(group_id = 1000*cur_group_id()) %>%
  mutate(gene_order = group_id + row_number()) %>%
  ungroup %>%
  mutate(gene_order = rank(gene_order))

library(pheatmap)
library(RColorBrewer)
library(viridis)
pSI.in_heatmap <- pSI.in
colnames(pSI.in_heatmap) <- str_to_sentence(colnames(pSI.in_heatmap))
colnames(pSI.in_heatmap)[colnames(pSI.in_heatmap) == "Gabaergic"] <- "GABAergic"

pSI.in_heatmap[unique(Reduce(c, specific_genes$ensembl_gene_id)), ] %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id, 
               names_to = "cell_type", 
               values_to = "expr") %>%
  left_join(specific_genes) %>%
  group_by(ensembl_gene_id) %>%
  mutate(expr = scale(expr, center = T)) %>%
  left_join(anno) %>%
  filter(external_gene_name %in% homologs$mmusculus_homolog_associated_gene_name) %>% 
  select(-c(specificity, group_id, gene_order)) %>%
  distinct %>%
  # filter(expr > mean(expr)) %>%
  ggplot(aes(x = cell_type, 
             y = expr, 
             size = expr)) +
  # geom_violin() +
  geom_jitter(width = 0.2) +
  coord_flip()











# check how many of nalls nearest select genes are enriched in TRAP ----

nalls_trap %>%
  filter(external_gene_name %in% nalls$nearest_select) %>% 
  filter(sumz_adj_MB < 0.01 & log2FoldChange_MB > 0) %>% View
  
{homologs %>%
  filter(external_gene_name %in% nalls$nearest_select) %>%
  pull(ensembl_gene_id)} 

homologs_ensembl %>%
  filter(external_gene_name %in% nalls$nearest_select)

# examine Casr StereoSeq counts ----
counts_casr <- read_csv("output/counts/counts_casr.csv")

counts_casr %>%
  mutate(cell_type_publish = ifelse(str_detect(cell_type_publish, "DA_"), 
                                    "DA", cell_type_publish)) %>%
  filter(count_Casr > 0) %>%
  ggplot(aes(x = count_Casr)) +
  geom_histogram() +
  facet_wrap(vars(cell_type_publish))
