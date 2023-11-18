library(tidyverse)
library(ggsci)
library(cowplot)
library(parallel)
library(ggrepel)
setwd("/laune_zfs/scratch/peter/f_active/paper_23/")

pksapply <-
  purrr::partial(base::sapply, simplify = F, USE.NAMES = T)

# outlier iqr function
#Interquartile range outlier detection function (IROF)
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

# DA markers pval/logfc ----

# TRAP enrichment
MB_FRACTION_META <- readRDS("input/trap/MB_FRACTION_META.rds")

trap_enrichment <- select(MB_FRACTION_META, external_gene_name, pvalue_C1, log2FoldChange) %>%
  mutate(p_1s = pvalue_C1 / 2,
         p_1s = ifelse(log2FoldChange > 0, p_1s, 1-p_1s))


# da_sn_markers <- read_csv("output/markers/general/DA_SN.csv") %>%
#   filter(pvals_adj < 0.01)
# da_vta_markers <- read_csv("output/markers/general/DA_VTA.csv") %>%
#   filter(pvals_adj < 0.01)

da_markers <- read_csv('input/05_annotation/da_markers.csv') %>%
  filter(pvals_adj < 0.01)

stereoseq_da_enrichment <- select(da_markers, "external_gene_name" = names, 
                                  "log2FoldChange" = logfoldchanges, "p_2s" = pvals) %>%
  mutate(p_1s = p_2s / 2, 
         p_1s = ifelse(log2FoldChange > 0, p_1s, 1-p_1s))

da_enrichment_join <- inner_join(trap_enrichment, stereoseq_da_enrichment, by = "external_gene_name", 
           suffix = c("_trap", "_stereoseq"))

min_p_ss <- filter(da_enrichment_join, p_1s_stereoseq > 0) %>%
  slice_min(p_1s_stereoseq, n = 1) %>%
  pull(p_1s_stereoseq)

da_enrichment_join <- da_enrichment_join %>%
  group_by(external_gene_name) %>% 
  mutate(p_1s_stereoseq = ifelse(p_1s_stereoseq == 0, min_p_ss, p_1s_stereoseq)) %>%
  group_split() %>%
  pksapply(function(g){
    n <- g$external_gene_name[1]
    p <- c(g$p_1s_trap[1], g$p_1s_stereoseq[1])
    tibble(external_gene_name = g$external_gene_name[1], 
           p_fisher = metap::sumlog(p)$p)
  }) %>%
  bind_rows

da_enrichment_join <- da_enrichment_join %>%
  left_join(select(trap_enrichment, external_gene_name, log2FoldChange)) %>%
  left_join(select(stereoseq_da_enrichment, external_gene_name, log2FoldChange), 
            by = "external_gene_name", suffix = c("_trap", "_ss")) %>%
  rowwise() %>%
  mutate(lfc_mean = mean(c(log2FoldChange_trap, log2FoldChange_ss), na.rm = T))
  
filter(da_enrichment_join, p_fisher == 0) %>%
  ungroup %>%
  arrange(desc(lfc_mean)) %>%
  slice_head(n = 20) %>%
  pivot_longer(starts_with("log2FoldChange"), 
               names_to = "technology", 
               names_prefix = "log2FoldChange_", 
               values_to = "lfc") %>%
  mutate(external_gene_name = fct_reorder(factor(external_gene_name), lfc, min), 
         technology = ifelse(technology == "ss", "Stereo-seq", "TRAP")) %>%
  ggplot(aes(x = external_gene_name, 
             y = lfc, 
             fill = technology, 
             group = technology)) +
  geom_col(colour = "black", size = 0.15, position = position_dodge()) +
  theme_cowplot(9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  # scale_fill_gradient(low = 'lightgrey', 
  #                     # mid = 'tomato1',
  #                     high = 'royalblue4') +
  scale_fill_manual(values = pal_d3()(10)[c(3, 1)]) +
  coord_flip() +
  theme(axis.text = element_text(face = "bold", size = 10), 
        axis.title = element_text(face = "bold", size = 14), 
        legend.title = element_text(face = "bold", size = 12), 
        legend.text = element_text(face = "bold", size = 10), 
        legend.position = "top") +
  labs(x = "Gene", 
       y = expression(Log[2] ~ Fold ~ Enrichment), 
       fill = "")

ggsave2('output/plots/figure_2/da_markers_lfc_by_technology.png', width = 45 * 2, height = 65 * 2, 
        units = "mm")

# da_markers %>%
#   slice_head(n = 15) %>%
  mutate(names = fct_reorder(factor(names), logfoldchanges, min)) %>%
  ggplot(aes(x = names, 
             y = logfoldchanges, 
             fill = logfoldchanges)) +
  geom_col(colour = "black", size = 0.15) +
  theme_cowplot(9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_gradient(low = 'lightgrey', 
                      # mid = 'tomato1',
                      high = 'royalblue4') +
  coord_flip() +
  theme(legend.position = "none") +
  labs(x = "Gene", 
       y = expression(Log[2] ~ Fold ~ Enrichment))

ggsave2('output/plots/figure_2/da_markers_lfc.png', width = 65, height = 65, 
        units = "mm")

# DA Spatial ----

# da_spatial <- read_csv("/active/paper_backup/input/cell_type_lists/midbrain/da_spatial_context.csv")
# 
# p_da_spatial <- da_spatial %>%
#   mutate(da = !is.na(sn_vta)) %>%
#   group_by(sample_id) %>%
#   arrange(da) %>%
#   mutate(x = (x-mean(x))/sd(x), 
#          y = (y-mean(y))/sd(y)) %>%
#   mutate(alpha = ifelse(da, 0.65, 0.35)) %>%
#   ggplot(aes(x = x, y = y, colour = da, alpha = alpha)) +
#   geom_point() +
#   facet_wrap(vars(mouse_id)) +
#   scale_y_reverse() +
#   scale_color_manual(values = c("grey", "darkred")) +
#   theme_cowplot() +
#   panel_border() +
#   theme(axis.text = element_blank(), 
#         axis.ticks = element_blank(),
#         axis.line = element_blank(),
#         legend.position = "none", 
#         axis.title = element_blank(), 
#         strip.text = element_text(face = "bold", colour = "black", size = 8), 
#         strip.background = element_rect(fill = "white", colour = "black"))
# 
# saveRDS(p_da_spatial, "plots/p_da_spatial.rds")

spatial_all_cells <- read_csv("input/05_annotation/spatial_all_cells.csv")

palette_pk <- c("#1F77B4FF", 
                "#FF7F0EFF", 
                "#2CA02CFF", 
                "#D62728FF", 
                "#9467BDFF", 
                "#8C564BFF", 
                "#E377C2FF", 
                "#BCBD22FF",
                "#17BECFFF")

pksapply(unique(spatial_all_cells$sample_name), function(sn){
  spatial_all_cells %>%
    filter(sample_name == sn) %>%
    mutate(interest = ifelse(str_detect(cell_type, "DOPA"), 
                             cell_type, "OTHER")) %>%
    mutate(interest = factor(interest),
           interest = relevel(interest, ref = "OTHER")) %>%
    arrange(interest) %>%
    ggplot(aes(x = x, 
               y = -y, 
               colour = interest, 
               size = interest, 
               alpha = interest)) +
    geom_point() +
    labs(colour = "") +
    # facet_wrap(vars(sample_name), ncol = 6) +
    theme_cowplot() +
    # scale_color_manual(values = c("grey", sample(pal_d3(palette = "category20b")(20), n_ct-1))) +
    scale_color_manual(values = c("grey", palette_pk[c(1,4)])) +
    # scale_color_manual(values = c("grey", pal_d3(palette = "category10")(1))) +
    # scale_color_manual(values = c("grey", color_ct[[t]])) +
    scale_size_manual(values = c(0.0000002, 2, 2), guide = 'none') +
    scale_alpha_manual(values = c(0.8, 1, 1), guide = 'none') +
    theme(axis.text = element_blank(), 
          axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_blank(),
          # plot.title = element_text(hjust = 0.5, size = 24), 
          legend.key.size = unit(1, 'cm'), 
          legend.position = "none")
  
  ggsave(paste0("output/plots/figure_2/per_sample_name/da_spatial_exemplar_size1_", sn, ".png"),
         width = 8, height = 5, dpi = 300)
  
})


# DA Spatial markers ----

# markers <- c("Th", "Slc6a3", "Slc18a2", "Ret", "En1", "Slc10a4")
# 
# da_spatial_markers <- read_csv("input/markers/da_markers_spatial_counts.csv") %>%
#   pivot_longer(`6330403K07Rik`:Zwint, 
#                names_to = "gene", 
#                values_to = "count") %>%
#   filter(mouse_id == "OLD_WT_REL121.1b") %>%
#   filter(gene %in% markers) %>%
#   mutate(gene = factor(gene, levels = markers)) %>%
#   filter(count > 1) %>%
#   group_by(sample_id) %>%
#   mutate(alpha = count / quantile(count, 0.2)) %>%
#   arrange(da, count) %>%
#   mutate(x = (x-mean(x))/sd(x), 
#          y = (y-mean(y))/sd(y))
# 
# p_da_spatial_markers <- da_spatial_markers %>%
#   ggplot(aes(x = x, y = y, colour = count, alpha = count/max(count))) +
#   geom_point() +
#   facet_wrap(vars(gene)) +
#   scale_y_reverse() +
#   scale_color_continuous(type = "viridis") +
#   theme_cowplot() +
#   panel_border() +
#   theme(axis.text = element_blank(), 
#         axis.ticks = element_blank(),
#         axis.line = element_blank(),
#         legend.position = "none", 
#         axis.title = element_blank(), 
#         strip.text = element_text(face = "bold", colour = "black", size = 16), 
#         strip.background = element_rect(fill = "white", colour = "black"))
# 
# ggsave2('input/figures/fig2/da_markers_spatial.png')
# 
# saveRDS(p_da_spatial_markers, "plots/p_da_spatial_markers.rds")

# spatial DE genes ----

library(metap)

# DA SN VTA
files <- list.files("output/spatial_genes", 
                    pattern = "1[a-z].csv",
                    full.names = T)
names(files) <- str_extract(files, "(?<=spatial_genes/)[:graph:]+(?=.csv)")

results <- lapply(names(files), function(n){
  read_csv(files[n]) %>% 
    select(-...1) %>%
    mutate(mouse_id = n) %>%
    drop_na
}) %>%
  bind_rows() %>%
  group_by(gene) %>%
  group_split()

names(results) <- lapply(results, function(x){pull(x, gene) %>% unique})

results <- lapply(names(results), function(n){
  meta = results[[n]] %>% pull(pval) %>% sumlog
  p = meta$p
  tibble(gene = n, 
         p = p)
}) %>%
  bind_rows() %>%
  mutate(p = ifelse(is.na(p), 1, p)) %>%
  mutate(padj = p.adjust(p, method = "BH"))

sum(results$padj < 0.01)

results %>% write_csv('output/spatial_genes/results.csv')

da_spatial_counts <- read_csv("output/spatial_genes/da_spatial_counts_top_100.csv") %>%
  pivot_longer(-1, 
               names_to = "gene", 
               values_to = "count") %>%
  left_join(read_csv('output/spatial_genes/da_spatial_counts_top_100_meta.csv'))
colnames(da_spatial_counts)[1] <- 'sample_cell_id'

mclapply(c('Cplx1', 'Nrip3', 'Calb1', 'Hpcal1', 'Fxyd6', 'Sncg', 'Aldh1a1', 'Ahi1'), function(GENE){
  data <- da_spatial_counts %>%
    filter(gene == GENE) %>%
    filter(mouse_id == 'OLD_OVX_REL120.1b') %>%
    filter(y > 16750) %>%
    group_by(gene) %>%
    filter(!outlier_iqr(count))
  # mutate(count = ifelse(count > quantile(count, 0.9),
  #                       quantile(count, 0.9),
  #                       count)) %>%
  
  data %>%
    arrange(count) %>%
    ggplot(aes(x = x, y = y)) +
    # geom_point(shape = 21, size = 1, colour = 'black', alpha = 0.25) +
    geom_point(
      aes(colour = count, 
          size = count)) +
    scale_y_reverse() +
    theme_cowplot() +
    labs(title = GENE, 
         colour = "Count", size = "Count") +
    # panel_border() +
    # facet_wrap(vars(mouse_id), scales = 'free') +
    # scale_fill_manual(values = c('red', 'white')) +
    scale_colour_gradient(low = "lightgrey", high = "red", guide = "legend") +
    scale_size_continuous(range = c(0.25, 1)) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          # legend.position = "none",
          axis.title = element_blank(), 
          # strip.text = element_text(face = "bold", colour = "black", size = 16), 
          # strip.background = element_rect(fill = "white", colour = "black")
          plot.title = element_text(hjust = 0.5, size = 12)
    )
  
  ggsave2(paste0('output/spatial_genes/da_spatial_', GENE, '.png'),
          width = 50,
          height = 30,
          units = "mm")
}, mc.cores = 18)

results <- lapply(names(files), function(n){
  read_csv(files[n]) %>% 
    select(-...1) %>%
    mutate(mouse_id = n) %>%
    drop_na
}) %>%
  bind_rows() %>%
  group_by(gene) %>%
  group_split()

names(results) <- lapply(results, function(x){pull(x, gene) %>% unique})

results <- lapply(names(results), function(n){
  meta = results[[n]] %>% pull(pval) %>% sumlog
  p = meta$p
  tibble(gene = n, 
         p = p)
}) %>%
  bind_rows() %>%
  mutate(p = ifelse(is.na(p), 1, p)) %>%
  mutate(padj = p.adjust(p, method = "BH"))

sum(results$padj < 0.01)

results %>% write_csv('output/spatial_genes/results.csv')



# SN VTA Spatial ----

# sn_vta <- read_csv("input/clusters/sn_vta.csv")
read_csv('input/05_annotation/markers/da_spatial_meta.csv') %>%
  group_by(sample_id) %>%
  mutate(x = (x-mean(x))/sd(x), 
         y = (y-mean(y))/sd(y)) %>%
  ggplot(aes(x = x, y = y, fill = cell_type)) +
  geom_point(shape = 21, size = 2) +
  facet_wrap(vars(sample_name, cell_type)) +
  scale_y_reverse() +
  scale_fill_d3() +
  theme_cowplot() +
  panel_border() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none", 
        axis.title = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black", size = 10), 
        strip.background = element_rect(fill = "white", colour = "black"))

p_sn_vta_spatial <- dopaminergic_metadata %>%
  group_by(sample_id) %>%
  mutate(x = (x-mean(x))/sd(x), 
         y = (y-mean(y))/sd(y)) %>%
  ggplot(aes(x = x, y = y, fill = sn_vta)) +
  geom_point(shape = 21, size = 2) +
  facet_wrap(vars(mouse_id, sn_vta)) +
  scale_y_reverse() +
  scale_fill_d3() +
  theme_cowplot() +
  panel_border() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none", 
        axis.title = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black", size = 10), 
        strip.background = element_rect(fill = "white", colour = "black"))

# read_csv('input/05_annotation/markers/da_spatial_meta.csv') %>%
#   filter(sample_name == 'YOUNG_WT_4') %>%
#   # filter(y > 16750) %>%
#   mutate(cell_type = factor(cell_type, levels = c("DOPA SN", "DOPA VTA"))) %>%
#   arrange(desc(cell_type)) %>%
#   ggplot(aes(x = x, y = y, fill = cell_type)) +
#   geom_point(shape = 21, size = 1.5, stroke = 0.1) +
#   annotate(geom="text", 
#            x=15250, y=18500, 
#            label="VTA", 
#            size = 4, 
#            fontface = "bold") +
#   annotate(geom="text", 
#            x=11750, y=17750, 
#            label="SN", 
#            size = 4, 
#            fontface = "bold") +
#   annotate(geom="text", 
#            x=18600, y=17750, 
#            label="SN", 
#            size = 4, 
#            fontface = "bold") +
#   # facet_wrap(vars(sn_vta), nrow = 2) +
#   scale_y_reverse() +
#   scale_fill_d3() +
#   theme_cowplot() +
#   # panel_border() +
#   theme(axis.text = element_blank(), 
#         axis.ticks = element_blank(),
#         axis.line = element_blank(),
#         legend.position = "none",
#         axis.title = element_blank(), 
#         strip.text = element_text(face = "bold", colour = "black", size = 10), 
#         strip.background = element_rect(fill = "white", colour = "black"))
  # labs(fill = "Cluster")

ggsave2(paste0('output/plots/figure_2/da_spatial_clusters.png'), 
        width = 85, 
        height = 35, 
        units = "mm")

# SN VTA volcano ----

sn_markers <- read_csv('input/05_annotation/markers/sn_vta_markers.csv')

sn_markers %>%
  mutate(region = ifelse(logfoldchanges > 0 & pvals_adj < 0.001, 
                         "SN", 
                         ifelse(logfoldchanges < 0 & pvals_adj < 0.001, 
                                "VTA", "NS"))) %>%
  mutate(label = ifelse(abs(logfoldchanges) > 0.25 & -log10(pvals_adj) > 4, 
                        names, "")) %>%
  # filter(region != "NS") %>%
  ggplot(aes(x = logfoldchanges, 
             y = -log10(pvals_adj), 
             label = label)) +
  geom_point(shape = 21, 
             stroke = 0.1,
             aes(fill = region 
                 # size = -log10(pvals_adj)
                 )) +
  theme_cowplot(9) +
  scale_x_continuous(limits = c(-1.5, 1.5),
                     oob = scales::squish) +
  scale_y_continuous(limits = c(0, 25),
                     oob = scales::squish,
                     expand = expansion(mult = c(0, 0.1))) +
  scale_size_continuous(range(0.02, 1)) +
  scale_fill_manual(values = c(NA, pal_d3(palette = "category10")(2))) +
  labs(x = expression(Log[2] ~ Fold ~ Difference ~ In ~ Abundance), 
       y = expression(-Log[10] ~ Adjusted ~ italic(P))) +
  geom_text_repel(size = 2, aes(fontface = "bold")) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.position = "none", 
        # axis.title = element_text(size = 12)
        ) +
  theme(axis.text = element_text(face = "bold", size = 6), 
        axis.title = element_text(face = "bold", size = 6)
        # legend.title = element_text(face = "bold", size = 10), 
        # legend.text = element_text(face = "bold", size = 10)
  )

ggsave2(paste0('output/plots/figure_2/da_sn_vta_volcano.png'), 
        width = 80, 
        height = 70, 
        units = "mm")

# SN VTA spatial heatmap ----

da_spatial_counts <- read_csv("input/05_annotation/markers/da_spatial_counts.csv") %>%
  pivot_longer(-1, 
               names_to = "gene", 
               values_to = "count") %>%
  left_join(read_csv("input/05_annotation/markers/da_spatial_meta.csv")) %>%
  mutate(gene = factor(gene, levels = c("Calb1", "Cplx1", "Aldh1a1", "Rab3c")))
colnames(da_spatial_counts)[1] <- 'sample_cell_id'

# pksapply(unique(da_spatial_counts$sample_name), function(sn){
da_spatial_counts %>%
  filter(sample_name == "YOUNG_OVX_3") %>%
  group_by(gene) %>%
  # filter(!outlier_iqr(count)) %>% View
  group_by(gene) %>%
  mutate(
    # count = log2(count),
    # count = count > 1
    count_prop = count / max(count),
    Positive = count_prop > 0.25
    # count = count > 0.3
  ) %>%
  # filter(gene == GENE) %>%
  
  arrange(count) %>%
  ggplot(aes(x = x, y = y)) +
  # geom_point(shape = 21, size = 1, colour = 'black', alpha = 0.25) +
  geom_point(
    # aes(fill = Positive, 
    #     size = Positive, 
    #     alpha = Positive),
    aes(color = count_prop, 
        size = count_prop, 
        alpha = count_prop
        ),
    # colour = "black", 
    # shape = 21
    ) +
  scale_y_reverse() +
  theme_cowplot() +
  # labs(fill = "Count", size = "Count", alpha = "Count") +
  # labs(title = sn) +
  # panel_border() +
  facet_wrap(vars(gene), scales = 'free', ncol = 4) +
  # scale_fill_manual(values = c('red', 'white')) +
  # scale_colour_gradient2(low = "lightgrey", mid = "red", midpoint = 0.5, high = "red", guide = "legend") +
  scale_size_continuous(range = c(0.1, 5)) +
  # scale_size_manual(values = c(0.2, 3)) +
  scale_alpha_continuous(range = c(0.75, 0.9)) +
  # scale_fill_gradient2(low = "lightgrey", mid = "lightgrey", midpoint = 1, high = "red", guide = "legend") +
  scale_color_gradientn(colours = c("red", "red", "red", "lightgrey"),
                       values = c(1, 0)) +
  
  # scale_fill_manual(values = c("lightgrey", "red")) +
  # scale_size_manual(values = c(0.25, 1)) +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        axis.title = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black", size = 24),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)
  )
# })

ggsave2(paste0('output/plots/figure_2/da_spatial_marker_sn_vta.png'),
        width = 180*2,
        height = 30*2,
        units = "mm", 
        dpi = 300)

# sox6 otx2 ----

t <- pksapply(c("sox6", "otx2"), function(x){
  t <- read_csv(paste0("input/05_annotation/", x, "_positivity.csv"))[,-1]
  colnames(t)[ncol(t)] <- "count"
  
  t <- pivot_wider(t, names_from = GOI, values_from = count) %>%
    mutate(total = Positive + Negative, 
           prop = Positive / total, 
           age = str_extract(sample_name, "YOUNG|OLD"))
  
})

pksapply(t, function(x){
  glm(prop ~ cell_type, weights = total, data = x, family = "binomial") %>% summary
})

t <- bind_rows(t, .id = "Gene") %>%
  mutate(Gene = str_to_sentence(Gene))

t %>%
  ggplot(aes(x = cell_type, y = prop, 
             fill = cell_type)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_wrap(vars(Gene)) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.1))) +
  scale_fill_d3() +
  labs(x = "Cell Type", y = "Percentage of Cells %") +
  theme_cowplot() +
  panel_border(color = "black") +
  geom_bracket(data = filter(t, Gene == "Otx2"), 
               xmin = "DOPA SN", xmax = "DOPA VTA", y.position = 0.075, 
               label = "***", tip.length = c(0.3, 0.2), 
               size = 1, label.size = 6) +
  geom_bracket(data = filter(t, Gene == "Sox6"), 
               xmin = "DOPA SN", xmax = "DOPA VTA", y.position = 0.125, 
               label = "***", tip.length = c(0.1, 0.4), 
               size = 1, label.size = 6) +
  theme(legend.position = "none")

ggsave2(paste0('output/plots/figure_2/da_sox6_otx2_snv_vta.png'),
        width = 75 * 2,
        height = 65 * 2,
        units = "mm", 
        dpi = 300)


# Axon: Hobson ----
hobson_apex_genes <- readxl::read_excel("input/third_party_data/hobson_apex/elife-70921-fig2-data5-v2.xlsx", 
                                        sheet = 2) %>%
  mutate(Log2FC = as.numeric(Log2FC)) %>%
  filter(qvalue < 0.01) %>%
  filter(Log2FC > 0) %>%
  pull(Genes)

# protein-coding only
axon_trap_genes <- readRDS("input/trap/AXON_FRACTION_META.rds") %>%
  filter(!str_detect(external_gene_name, "^Gm")) %>%
  filter(!str_detect(external_gene_name, "Rik$")) %>%
  filter(sumz_adj < 0.01) %>%
  filter(log2FoldChange > 0) %>%
  filter(mb_translated) %>%
  filter(!conflict) %>%
  pull(external_gene_name)

total_genes <- readRDS("input/trap/AXON_FRACTION_META.rds") %>%
  filter(!str_detect(external_gene_name, "^Gm")) %>%
  filter(!str_detect(external_gene_name, "Rik$")) %>%
  pull(external_gene_name)

intersect(hobson_apex_genes, 
          axon_trap_genes)


overlap_significance <- function(genes_all, gene_sets, iterations) {
  observed <- length(reduce(gene_sets, intersect))
  simulated <- map_dbl(seq_len(iterations), function(x) {
    sim <- map(lengths(gene_sets), ~sample(genes_all, .x))
    sim <- length(reduce(sim, intersect))
    return(sim)
  })
  pval <- (sum(simulated >= observed) + 1) / (iterations + 1)
  return(list(pval=pval, simulated_values=simulated, observed=observed))
}

overlap_significance(total_genes, 
                     gene_sets = list(axon_trap_genes, 
                                      hobson_apex_genes), 
                     1e6)

observed <- length(intersect(hobson_apex_genes, 
                             axon_trap_genes))
sets = list(axon_trap_genes, 
            hobson_apex_genes)
all_genes <- total_genes

hyper_pval <- phyper(
  q=observed, m=length(sets[[1]]),
  n=length(all_genes) - length(sets[[1]]),
  k=length(sets[[2]]), lower.tail=FALSE
)

