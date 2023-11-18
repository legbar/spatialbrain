library(tidyverse)
library(parallel)
library(cowplot)
library(ggsci)
library(DESeq2)
library(ggpubr)
library(ggrepel)

pksapply <-
  purrr::partial(base::sapply, simplify = F, USE.NAMES = T)

# FIGURE A: MB fraction plot counts -------------------------------------------------

MB_FRACTION_META <- readRDS("input/trap/MB_FRACTION_META.rds")
anno <- readRDS("input/trap/anno.rds")

# MA-plot style
plot_data <- MB_FRACTION_META %>%
  # filter(abs(log2FoldChange) < 10) %>%
  mutate(enrichment = ifelse(sumz_adj > 0.01, "Unchanged",
                             ifelse(log2FoldChange > 0, "Enriched", "Depleted")))

# highlight markers
highlight_markers <- plot_data %>%
  left_join(anno) %>%
  filter(external_gene_name %in% c(
    c("Slc6a3", "Th", 
      # "Ddc", 
      "Slc18a2"),
    c("Gfap", 
      # "Gad2", 
      "S100b", 
      "Gad1"
      # "Aldh1l1"
    )
  ))

plot_data %>%
  ggplot(aes(x = baseMean_C1,
             y = log2FoldChange,
             colour = enrichment)) +
  theme_cowplot(10) +
  geom_point(size = 0.25,
             alpha = 0.5) +
  geom_point(data = highlight_markers,
             aes(fill = enrichment),
             size = 2,
             shape = 21,
             colour = "black") +
  geom_label_repel(data = highlight_markers,
                   aes(label = external_gene_name),
                   colour = "black",
                   arrow = arrow(length = unit(0.01, "npc")),
                   # box.padding = 1,
                   size = 5
                   # nudge_x = -1
  ) +
  scale_x_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6), 
                expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(-6, 6),
                     breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
                     oob = scales::squish) +
  scale_color_manual(values = c("#1F77B4FF",
                                "#2CA02CFF",
                                "#C1C1C1",
                                "#FF7F0EFF")) +
  scale_fill_manual(values = c("#1F77B4FF",
                               "#2CA02CFF",
                               "#C1C1C1",
                               "#FF7F0EFF")) +
  labs(x = "Mean Counts",
       y = "Log2 Fold Change",
       colour = "Enrichment") +
  theme(legend.position = "top",
        # legend.text = element_text(size = 14),
        # legend.title = element_text(size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ,
         fill = FALSE)

ggsave("output/plots/figure_3/ma_mb.png", 
       bg = "transparent", 
       width = 120, height = 120,
       units = "mm",
       dpi = 300)

# cell types ----

# load panglao DB
panglao <- read_delim("input/trap/PanglaoDB_markers_27_Mar_2020.tsv",
                      delim = "\t") %>%
  # filter(ifelse(str_detect(`official gene symbol`, "TH|SLC6A3|SLC18A2|DDC") & !str_detect(`cell type`, "drenergic|otonergic|nterneuron"),
  #               TRUE, ifelse(!str_detect(`official gene symbol`, "TH|SLC6A3|SLC18A2|DDC"), TRUE, FALSE ))) %>%
  group_by(`cell type`) %>%
  mutate(group_size = length(unique(`official gene symbol`)))

anno_human <- readRDS("input/trap/anno_human.rds")

# plot cell type enrichment
MB_FRACTION_META %>%
  left_join(anno_human) %>%
  left_join(panglao,
            by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
  mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
         score = -log10(sumz_adj) * log2FoldChange_C1) %>%
  group_by(`cell type`) %>%
  summarise(score = mean(score),
            group_proportion = n()/group_size,
            group_size = group_size,
            prop_score = score * group_proportion) %>%
  arrange(desc(prop_score)) %>%
  distinct() %>%
  filter(`cell type` %in% c("Astrocytes",
                            "Cholinergic neurons",
                            "Dopaminergic neurons",
                            "Endothelial cells",
                            "GABAergic neurons",
                            "Glutaminergic neurons",
                            "Interneurons",
                            "Microglia",
                            "Oligodendrocytes"
                            # "Oligodendrocyte progenitor cells"
  )) %>%
  mutate(`cell type` = factor(`cell type`)) %>%
  mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
  ggplot(aes(x = prop_score,
             y = `cell type`,
             fill = `cell type`)) +
  theme_cowplot(6) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  labs(x = "Enrichment Score",
       y = "Cell Type") +
  theme(legend.position = "none",
        # axis.text = element_text(size = 12),
        axis.title.y = element_blank())

ggsave("output/plots/figure_3/mb_enrichment.png", 
       bg = "transparent", 
       width = 60, height = 80,
       units = "mm",
       dpi = 300)

# pca function ----
make_pca <- function(vsd,
                     dds_object = dds,
                     # group,
                     # PC_x = PC1,
                     # PC_y = PC2,
                     identifier = "sample_name",
                     ntop = Inf,
                     add_meta = TRUE){
  
  
  # make the pca object
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = T)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select[-1], ])) #Not scaling as passing log transformed data: https://www.biostars.org/p/280615/#283126
  
  if(add_meta == TRUE){
    # calculate the variation per component
    var_explained0 <- pca$sdev^2/sum(pca$sdev^2) #Calculate PC variance
    names(var_explained0) <- colnames(pca$x)
    var_explained <- as_tibble(var_explained0,
                               rownames = "PC") %>%
      mutate(var = round(var_explained0*100)) %>%
      select(PC, var)
    
    pca$var_explained <- var_explained
    
    pca$x <- pca$x %>%
      as_tibble(rownames = identifier) %>%
      inner_join(colData(dds_object), copy = TRUE)
    
  }
  
  return(pca)
}
# FIG 1 PCA ----
dds <- readRDS("/active/paper/input/trap/dds.rds")
# C1
vsd <- vst(dds)
pca <- make_pca(vsd)
data <- plotPCA(vsd, intgroup = c("fraction", "cohort", "region"), ntop = Inf, returnData = TRUE)

data %>%
  mutate(fraction = ifelse(fraction == 'TOTAL', 'INPUT', 'IP'), 
         cohort = str_replace(cohort, '^C', 'Cohort ')) %>%
  ggplot(aes(PC1, PC2, colour = region, shape = fraction)) +
  geom_point(size = 2) +
  theme_cowplot() +
  scale_color_d3() +
  labs(colour = "Region",
       shape = "Fraction",
       x = paste("PC1:", pca$var_explained$var[1], "%"),
       y = paste("PC2:", pca$var_explained$var[2], "%")) +
  theme(legend.position = "top") +
  theme(legend.title = element_text(face = "bold")) +
  facet_wrap(vars(cohort), scales = "free")

ggsave("output/plots/figure_1/pca.png", 
       bg = "transparent", 
       width = 8, height = 3, 
       dpi = 300)

# FIG 1 UMAP ----

umap_all_cells <- read_csv("output/UMAP/all_cells.csv")

# GLIAL vs NEURONAL
umap_all_cells %>%
  group_by(cell_type_publish) %>%
  filter(n() > 500) %>%
  ungroup %>%
  mutate(NEU_GLIA = ifelse(str_detect(cell_type_publish, "GLIA|OLIG|RBC|MICRO|ASTRO"), 
                           "GLIAL", "NEURONAL")) %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             colour = NEU_GLIA)) +
  geom_point(alpha = 0.1, 
             size = 0.01) +
  # facet_wrap(vars(mouse_id)) +
  theme_cowplot() +
  scale_color_d3() +
  theme(legend.position = c(0.1, 0.85)) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  labs(colour = "")
ggsave("output/plots/figure_1/umap_all_cells.png", 
       bg = "transparent", 
       width = 4, height = 4, 
       dpi = 300)

# GLIAL UMAP
umap_all_cells %>%
  mutate(NEU_GLIA = ifelse(str_detect(cell_type_publish, "GLIA|OLIG|RBC|MICRO|ASTRO"), 
                           "GLIAL", "NEURONAL")) %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  filter(NEU_GLIA == "GLIAL") %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             colour = cell_type_publish)) +
  geom_point(alpha = 0.1, 
             size = 0.1) +
  # facet_wrap(vars(mouse_id)) +
  theme_cowplot() +
  scale_color_d3() +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  labs(colour = "") +
  scale_x_continuous(limits = c(5, 30))
ggsave("output/plots/figure_1/umap_glia.png", 
       bg = "transparent", 
       width = 4.5, height = 4, 
       dpi = 300)

# NEURONAL UMAP
set.seed(8)
umap_all_cells %>%
  mutate(NEU_GLIA = ifelse(str_detect(cell_type_publish, "GLIA|OLIG|RBC|ASTRO|MICRO"), 
                           "GLIAL", "NEURONAL")) %>%
  filter(NEU_GLIA == "NEURONAL") %>%
  group_by(cell_type_publish) %>%
  filter(n() > 500) %>%
  ungroup %>%
  mutate(cell_type_publish = ifelse(str_detect(cell_type_publish, "DA_"), "DOPAMINERGIC", cell_type_publish)) %>%
  mutate(cell_type_publish = str_remove(cell_type_publish, "NEU")) %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  mutate(cell_type_publish = ifelse(str_detect(cell_type_publish, "CA1|CA3|DG|SUBIC|LIMBIC"), 
                                    "HIPPOCAMPAL", cell_type_publish)) %>%
  group_by(cell_type_publish) %>%
  mutate(tally = n()) %>%
  # pull(cell_type_publish) %>% unique
  # mutate(cell_type_publish = ifelse(str_detect(cell_type_publish, 
  #                                              "^AGEING"), 
  #                                   "AGEING", 
  #                                   ifelse(str_detect(cell_type_publish,
  #                                                     "^NONAG"), 
  #                                          "NON-AGEING", 
  #                                          ifelse(str_detect(cell_type_publish, 
  #                                                            "MENING"),
  #                                                 "MENINGEAL", "RBC")))) %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             colour = cell_type_publish)) +
  geom_point(alpha = 0.1, 
             size = 0.5) +
  # facet_wrap(vars(mouse_id)) +
  scale_color_manual(values = sample(pal_d3(palette = "category20")(15), 
                                     replace = F)) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  labs(colour = "") +
  # facet_wrap(vars(cell_type_publish)) +
  # scale_x_continuous(limits = c(5, 30)) +
  theme_cowplot() +
  scale_y_continuous(limits = c(NA, 20))
ggsave("output/plots/figure_1/umap_neuronal.png", 
       bg = "transparent", 
       width = 8, height = 4, 
       dpi = 300)

# umap_all_cells %>%
#   filter(UMAP1 > 11.5 & UMAP1 < 14 & UMAP2 < -6 & UMAP2 > -8.5) %>% 
#   group_by(cell_type) %>%
#   tally
# 
# umap_all_cells %>%
#   filter(UMAP1 > 16 & UMAP1 < 18 & UMAP2 > 17 & UMAP2 < 18) %>% 
#   group_by(cell_type) %>%
#   tally

# FIG 1 SPATIAL ----

spatial_all_cells <- read_csv("input/05_annotation/spatial_all_cells.csv")

spatial_all_cells$cell_type %>% unique

palette_pk <- c("#1F77B4FF", 
                "#FF7F0EFF", 
                "#2CA02CFF", 
                "#D62728FF", 
                "#9467BDFF", 
                "#8C564BFF", 
                "#E377C2FF", 
                "#BCBD22FF",
                "#17BECFFF")

plot_examplar_cell_type_spatial <- function(sn = 'OLD_OVX_1', ct, t){
  df <- spatial_all_cells %>%
    filter(sample_name == sn) %>%
    mutate(interest = ifelse(str_detect(cell_type, ct), 
                             cell_type, "OTHER"), 
           size_interest = str_detect(cell_type, ct)) %>%
    mutate(interest = factor(interest),
           interest = relevel(interest, ref = "OTHER")) %>%
    arrange(interest)
  
  n_ct <- length(unique(df$interest))
  
  df %>%
    ggplot(aes(x = x, 
               y = -y, 
               colour = interest, 
               size = size_interest, 
               alpha = size_interest)) +
    geom_point() +
    labs(colour = "", 
         title = t) +
    theme_cowplot() +
    # scale_color_manual(values = c("grey", sample(pal_d3(palette = "category20b")(20), n_ct-1))) +
    scale_color_manual(values = c("grey", sample(palette_pk, n_ct-1))) +
    # scale_color_manual(values = c("grey", pal_d3(palette = "category10")(1))) +
    # scale_color_manual(values = c("grey", color_ct[[t]])) +
    scale_size_manual(values = c(0.02, 1), guide = 'none') +
    scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
    theme(axis.text = element_blank(), 
          axis.line = element_blank(), 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 36), 
          legend.key.size = unit(1, 'cm'), 
          legend.position = "none")
  
  ggsave(paste0("output/plots/figure_1/spatial_", t, ".png"), 
         width = 8, height = 5.5)
}

cts <- list("DOPAMINERGIC" = "DOPA SN|DOPA VTA", 
            "HIPPOCAMPUS" = "CA|SUBIC|DG",
            "THALAMUS" = "THALAMUS",
            "CORTEX EXCIT 1" = "L1|L4|L6|Tsh",
            "CORTEX EXCIT 2" = "L2/3|L5|AMYG",
            "CARTPT" = "Cartpt", 
            "OLIGO" = "OLIGO", 
            "ASTRO" = "ASTRO", 
            "CORTEX INHIB" = "IN Sst|IN Ad|IN Pva")

color_ct <- pal_d3(palette = "category20")(10)[1:length(cts)]
names(color_ct) <- names(cts)

set.seed(1)
pksapply(names(cts), function(n){
  plot_examplar_cell_type_spatial(sn = "YOUNG_WT_4", ct = cts[[n]], t = n)
})

spatial_all_cells %>%
  # filter(!cell_type_publish %in% c("RBC", 
  #                                  "GENERAL_MB_NEU"
  #                                  )) %>%
  # filter(!str_detect(
  #   cell_type_publish,
  #   "GLIA|OLIG|RBC")) %>%
  filter(mouse_id == "OLD_WT_REL121.1b") %>%
  mutate(cell_type_publish = ifelse(
    str_detect(cell_type_publish, "DA_"),
    "DOPAMINERGIC",
    cell_type_publish
  )) %>%
  mutate(cell_type_publish = str_remove(cell_type_publish, "NEU")) %>%
  mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  # mutate(cell_type_publish = ifelse(
  #   str_detect(cell_type_publish, "CA1|CA3|DG|SUBIC|LIMBIC"),
  #   "HIPPOCAMPAL",
  #   cell_type_publish
  # )) %>%
  # mutate(
  #   cell_type_publish = ifelse(str_detect(
  #     cell_type_publish,
  #     "GLIA|OLIG|RBC"),
  #     ifelse(
  #       str_detect(cell_type_publish,
#                  "^AGEING"),
#       "AGEING",
#       ifelse(
#         str_detect(cell_type_publish,
#                    "^NONAG"),
#         "NON-AGEING",
#         ifelse(str_detect(cell_type_publish,
#                           "MENING"),
#                "MENINGEAL", "RBC")
#       )
#     ), 
#     cell_type_publish)
#   ) %>% 
group_by(cell_type_publish, mouse_id) %>%
  filter(n() > 100) %>%
  ungroup %>%
  # group_by(cell_type_publish) %>% tally
  mutate(cell_type_publish = factor(cell_type_publish), 
         cell_type_publish = relevel(cell_type_publish, ref = "DOPAMINERGIC")) %>%
  arrange(cell_type_publish) %>%
  ggplot(aes(x = x, 
             y = -y, 
             colour = cell_type_publish)) +
  geom_point(size = 0.25) +
  # facet_wrap(vars(mouse_id)) +
  # scale_color_d3(palette = "category20") +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  labs(colour = "") +
  facet_wrap(vars(cell_type_publish)) +
  # scale_x_continuous(limits = c(5, 30)) +
  theme_cowplot() +
  facet_wrap(vars(cell_type_publish)) +
  panel_border() +
  # facet_wrap(vars(mouse_id), scales = "free") +
  theme(legend.position = "none",
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(size = 6))






































# 
# 
# 
# 
# 
# spatial_all_cells <- read_csv("output/spatial_coordinates/all_cells.csv")
# 
# spatial_all_cells$cell_type_publish %>% unique
# 
# spatial_all_cells %>%
#   # filter(!cell_type_publish %in% c("RBC", 
#   #                                  "GENERAL_MB_NEU"
#   #                                  )) %>%
#   # filter(!str_detect(
#   #   cell_type_publish,
#   #   "GLIA|OLIG|RBC")) %>%
#   filter(mouse_id == "OLD_WT_REL121.1b") %>%
#   mutate(cell_type_publish = ifelse(
#     str_detect(cell_type_publish, "DA_"),
#     "DOPAMINERGIC",
#     cell_type_publish
#   )) %>%
#   mutate(cell_type_publish = str_remove(cell_type_publish, "NEU")) %>%
#   mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
#   # mutate(cell_type_publish = ifelse(
#   #   str_detect(cell_type_publish, "CA1|CA3|DG|SUBIC|LIMBIC"),
#   #   "HIPPOCAMPAL",
#   #   cell_type_publish
#   # )) %>%
#   # mutate(
#   #   cell_type_publish = ifelse(str_detect(
#   #     cell_type_publish,
#   #     "GLIA|OLIG|RBC"),
#   #     ifelse(
#   #       str_detect(cell_type_publish,
#   #                  "^AGEING"),
#   #       "AGEING",
#   #       ifelse(
#   #         str_detect(cell_type_publish,
#   #                    "^NONAG"),
#   #         "NON-AGEING",
#   #         ifelse(str_detect(cell_type_publish,
#   #                           "MENING"),
#   #                "MENINGEAL", "RBC")
#   #       )
#   #     ), 
#   #     cell_type_publish)
#   #   ) %>% 
#   group_by(cell_type_publish, mouse_id) %>%
#   filter(n() > 100) %>%
#   ungroup %>%
#   # group_by(cell_type_publish) %>% tally
#   mutate(cell_type_publish = factor(cell_type_publish), 
#          cell_type_publish = relevel(cell_type_publish, ref = "DOPAMINERGIC")) %>%
#   arrange(cell_type_publish) %>%
#   ggplot(aes(x = x, 
#              y = -y, 
#              colour = cell_type_publish)) +
#   geom_point(size = 0.25) +
#   # facet_wrap(vars(mouse_id)) +
#   # scale_color_d3(palette = "category20") +
#   guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
#   labs(colour = "") +
#   facet_wrap(vars(cell_type_publish)) +
#   # scale_x_continuous(limits = c(5, 30)) +
#   theme_cowplot() +
#   facet_wrap(vars(cell_type_publish)) +
#   panel_border() +
#   # facet_wrap(vars(mouse_id), scales = "free") +
#   theme(legend.position = "none",
#         axis.text = element_blank(), 
#         axis.line = element_blank(), 
#         axis.title = element_blank(), 
#         axis.ticks = element_blank(), 
#         strip.text = element_text(size = 6))
# {
#   spatial_plot_aes <- tibble(cell_type_publish = spatial_all_cells$cell_type_publish %>% unique) %>%
#     mutate(size = 0.1, alpha = 0.5, colour = 'blue') %>%
#     mutate(cell_type_publish = ifelse(
#       str_detect(cell_type_publish, "DA_"),
#       "DOPAMINERGIC",
#       cell_type_publish
#     )) %>%
#     mutate(cell_type_publish = str_remove(cell_type_publish, "NEU")) %>%
#     mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
#     distinct
#   
#   pt_size = 0.5
#   
#   # DA
#   spatial_plot_aes[spatial_plot_aes$cell_type_publish == "DOPAMINERGIC",]$size <- 1
#   spatial_plot_aes[spatial_plot_aes$cell_type_publish == "DOPAMINERGIC",]$alpha <- 1
#   spatial_plot_aes[spatial_plot_aes$cell_type_publish == "DOPAMINERGIC",]$colour <- 'green'
#   
#   # THAL
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^THALAMIC D"),]$size <- 1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^THALAMIC D"),]$alpha <- 0.7
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^THALAMIC D"),]$colour <- 'red'
#   
#   # THAL
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^THALAMIC V"),]$size <- 1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^THALAMIC V"),]$alpha <- 0.7
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^THALAMIC V"),]$colour <- 'turquoise'
#   
#   # CA1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CA1"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CA1"),]$alpha <- 1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CA1"),]$colour <- 'green4'
#   
#   # DG
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^DG"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^DG"),]$alpha <- 1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^DG"),]$colour <- 'steelblue2'
#   
#   # CA3
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CA3"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CA3"),]$alpha <- 1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CA3"),]$colour <- 'orange'
#   
#   # CORTICAL_OUTER_EXCIT_NEU
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CORTICAL OUTER"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CORTICAL OUTER"),]$alpha <- 0.4
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CORTICAL OUTER"),]$colour <- 'slateblue1'
#   
#   # CORTICAL INNER
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CORTICAL INNER"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CORTICAL INNER"),]$alpha <- 0.75
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^CORTICAL INNER"),]$colour <- 'red4'
#   
#   # RN_PVALB_SNCG_NEU
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^NIGRAL"),]$size <- 1
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^NIGRAL"),]$alpha <- 0.5
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^NIGRAL"),]$colour <- 'lightcoral'
#   
#   # Dorsal subic
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC D"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC D"),]$alpha <- 0.5
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC D"),]$colour <- 'grey81'
#   
#   # Ventral subic
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC V"),]$size <- pt_size
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC V"),]$alpha <- 0.5
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC V"),]$colour <- 'darkolivegreen'
#   
#   # Meningeal glia
#   # spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^SUBIC V"),]$size <- 0.5
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^MENINGEAL"),]$alpha <- 0.5
#   spatial_plot_aes[str_detect(spatial_plot_aes$cell_type_publish, "^MENINGEAL"),]$colour <- 'cornsilk'
#   
#   spatial_colours <- spatial_plot_aes$colour
#   names(spatial_colours) <- spatial_plot_aes$cell_type_publish
#   spatial_sizes <- spatial_plot_aes$size
#   names(spatial_sizes) <- spatial_plot_aes$cell_type_publish
#   spatial_alphas <- spatial_plot_aes$alpha
#   names(spatial_alphas) <- spatial_plot_aes$cell_type_publish
#   
#   # set.seed(5)
#   spatial_all_cells %>%
#     # filter(!str_detect(
#     #   cell_type_publish,
#     #   "GLIA|OLIG|RBC|MICRO|ASTRO")) %>%
#     filter(mouse_id == "OLD_WT_REL121.1b") %>%
#     mutate(cell_type_publish = ifelse(
#       str_detect(cell_type_publish, "DA_"),
#       "DOPAMINERGIC",
#       cell_type_publish
#     )) %>%
#     mutate(cell_type_publish = str_remove(cell_type_publish, "NEU")) %>%
#     mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
#     # group_by(cell_type_publish, mouse_id) %>%
#     # filter(n() > 200) %>%
#     # tally
#     ungroup %>%
#     arrange(cell_type_publish) %>%
#     # group_by(cell_type_publish) %>% tally
#     mutate(cell_type_publish = factor(cell_type_publish), 
#            cell_type_publish = relevel(cell_type_publish, ref = "DOPAMINERGIC")) %>%
#     # arrange(size) %>%
#     ggplot(aes(x = x, 
#                y = -y, 
#                colour = cell_type_publish, 
#                size = cell_type_publish, 
#                alpha = cell_type_publish)) +
#     geom_point(position = position_jitter(height = 0.1, width = 0.1)) +
#     theme_cowplot() +
#     # facet_wrap(vars(mouse_id)) +
#     # scale_color_d3(palette = "category20") +
#     scale_color_manual(values = spatial_colours) +
#     scale_size_manual(values = spatial_sizes) +
#     # scale_size_manual(values = spatial_sizes*4) + # for A1 image
#     scale_alpha_manual(values = spatial_alphas) +
#     # guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
#     # labs(colour = "") +
#     # scale_color_manual(values = sample(pal_d3(palette = 'category20')(20), size = 20)) +
#     # facet_wrap(vars(cell_type_publish)) +
#     # scale_x_continuous(limits = c(5, 30)) +
#     # facet_wrap(vars(cell_type_publish)) +
#     # panel_border() +
#     # facet_wrap(vars(mouse_id), scales = "free") +
#     theme(legend.position = "none",
#           axis.text = element_blank(), 
#           axis.line = element_blank(), 
#           axis.title = element_blank(), 
#           axis.ticks = element_blank(), 
#           strip.text = element_text(size = 6), 
#           plot.background = element_rect(fill = 'black')) 
#   }
# 
# ggsave("output/plots/figure_1/spatial_full.png", 
#        bg = "transparent", 
#        width = 6.5, height = 4, 
#        dpi = 300)
# 
# ggsave("output/plots/figure_1/spatial_full_A1.png", 
#        bg = "transparent", 
#        width = 841, height = 594, 
#        units = "mm",
#        dpi = 300)
# 
# ggsave("output/plots/figure_1/spatial_full_A2.png", 
#        bg = "transparent", 
#        width = 594, height = 420, 
#        units = "mm",
#        dpi = 300)
# 
# # spatial_all_cells$cell_type_publish %>% unique %>% sort
# 
# # Hippocampus
# p_spatial_hipp <- spatial_all_cells %>%
#   filter(mouse_id == "OLD_WT_REL121.1b") %>%
#   # mutate(cell_type_publish = ifelse(
#   #   str_detect(cell_type_publish, "DA_"),
#   #   "DOPAMINERGIC",
#   #   cell_type_publish
#   # )) %>%
#   mutate(interest = ifelse(str_detect(cell_type_publish, "CA1|CA3|DG|SUBIC"), 
#                            cell_type_publish, "OTHER"), 
#          size_interest = str_detect(cell_type_publish, "CA1|CA3|DG|SUBIC")) %>%
#   # filter(str_detect(cell_type_publish, "CA1|CA3|DG|SUBIC")) %>%
#   mutate(interest = str_replace(interest, "CA1_NEU", "CA1"), 
#          interest = str_replace(interest, "CA3_NEU", "CA3"), 
#          interest = str_replace(interest, "DG_NEU", "DENTATE\nGYRUS"), 
#          interest = str_replace(interest, "SUBIC_DORSAL_NEU", "DORSAL\nSUBICULUM"), 
#          interest = str_replace(interest, "SUBIC_VENTRAL_CA1", "VENTRAL\nSUBICULUM")) %>%
#   mutate(interest = factor(interest), 
#          interest = relevel(interest, ref = "OTHER")) %>%
#   arrange(interest) %>%
#   ggplot(aes(x = x, 
#              y = -y, 
#              colour = interest, 
#              size = size_interest, 
#              alpha = size_interest)) +
#   geom_point() +
#   # facet_wrap(vars(mouse_id), scales = "free") +
#   # scale_color_d3(palette = "category20") +
#   # guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
#   labs(colour = "", 
#        title = "HIPPOCAMPUS") +
#   scale_x_continuous(limits = c(4250, 25000)) +
#   scale_y_continuous(limits = c(-15500, -2000)) +
#   theme_cowplot() +
#   scale_color_manual(values = c("grey", pal_d3()(5))) +
#   scale_size_manual(values = c(0.02, 1), guide = 'none') +
#   scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
#   theme(axis.text = element_blank(), 
#         axis.line = element_blank(), 
#         axis.title = element_blank(), 
#         axis.ticks = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         legend.key.size = unit(1, 'cm'))
#   
# ggsave("output/plots/figure_1/spatial_hipp.png", 
#        bg = "transparent", 
#        width = 6.5, height = 4, 
#        dpi = 300)
# 
# # Thalamus
# p_spatial_thal <- spatial_all_cells %>%
#   filter(mouse_id == "OLD_WT_REL121.1b") %>%
#   mutate(interest = ifelse(str_detect(cell_type_publish, "^THALA"), 
#                            cell_type_publish, "OTHER"), 
#          size_interest = str_detect(cell_type_publish, "^THALA")) %>%
#   mutate(interest = str_replace(interest, "THALAMIC_DORSOLAT_NEU", "DORSOLATERAL\nTHALAMUS"),
#          interest = str_replace(interest, "THALAMIC_VENTROMED_NEU", "VENTROMEDIAL\nTHALAMUS")) %>%
#   mutate(interest = factor(interest), 
#          interest = relevel(interest, ref = "OTHER")) %>%
#   arrange(interest) %>%
#   ggplot(aes(x = x, 
#              y = -y, 
#              colour = interest, 
#              size = size_interest, 
#              alpha = size_interest)) +
#   geom_point() +
#   # facet_wrap(vars(mouse_id), scales = "free") +
#   # scale_color_d3(palette = "category20") +
#   # guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
#   labs(colour = "", 
#        title = "THALAMUS") +
#   scale_x_continuous(limits = c(4250, 25000)) +
#   scale_y_continuous(limits = c(-15500, -2000)) +
#   theme_cowplot() +
#   scale_color_manual(values = c("grey", pal_d3()(5))) +
#   scale_size_manual(values = c(0.02, 1), guide = 'none') +
#   scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
#   theme(axis.text = element_blank(), 
#         axis.line = element_blank(), 
#         axis.title = element_blank(), 
#         axis.ticks = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         legend.key.size = unit(1, 'cm'))
# 
# ggsave("output/plots/figure_1/spatial_thal.png", 
#        bg = "transparent", 
#        width = 6.5, height = 4, 
#        dpi = 300)
# 
# # Reticulata/Red nucleus
# p_spatial_pvalb <- spatial_all_cells %>%
#   filter(mouse_id == "OLD_WT_REL121.1b") %>%
#   mutate(interest = ifelse(str_detect(cell_type_publish, "NIGRAL_RN|RN_PVALB"), 
#                            cell_type_publish, "OTHER"), 
#          size_interest = str_detect(cell_type_publish, "NIGRAL_RN|RN_PVALB")) %>%
#   mutate(interest = str_replace(interest, "NIGRAL_RN_PVALB_GABA_NEU", "PVALB+ GABAERGIC"),
#          interest = str_replace(interest, "RN_PVALB_SNCG_NEU", "PVALB+ SNCG+\nGABAERGIC")) %>%
#   mutate(interest = factor(interest), 
#          interest = relevel(interest, ref = "OTHER")) %>%
#   arrange(interest) %>%
#   ggplot(aes(x = x, 
#              y = -y, 
#              colour = interest, 
#              size = size_interest, 
#              alpha = size_interest)) +
#   geom_point() +
#   # facet_wrap(vars(mouse_id), scales = "free") +
#   # scale_color_d3(palette = "category20") +
#   # guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
#   labs(colour = "", 
#        title = "PVALB+ GABAERGIC") +
#   scale_x_continuous(limits = c(4250, 25000)) +
#   scale_y_continuous(limits = c(-15500, -2000)) +
#   theme_cowplot() +
#   scale_color_manual(values = c("grey", pal_d3()(5))) +
#   scale_size_manual(values = c(0.02, 1), guide = 'none') +
#   scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
#   theme(axis.text = element_blank(), 
#         axis.line = element_blank(), 
#         axis.title = element_blank(), 
#         axis.ticks = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         legend.key.size = unit(1, 'cm'))
# 
# ggarrange(p_spatial_hipp, 
#           p_spatial_thal, 
#           p_spatial_pvalb,
#           nrow = 3, 
#           align = 'v')
# 
# ggsave("output/plots/figure_1/spatial.png", 
#        bg = "transparent", 
#        width = 8, height = 12, 
#        dpi = 300)
# 
# # Dopaminergic neurons
# p_spatial_da <- spatial_all_cells %>%
#   filter(mouse_id == "OLD_WT_REL121.1b") %>%
#   # mutate(cell_type_publish = ifelse(
#   #   str_detect(cell_type_publish, "DA_"),
#   #   "DOPAMINERGIC",
#   #   cell_type_publish
#   # )) %>%
#   mutate(interest = ifelse(str_detect(cell_type_publish, "DA_SN|DA_VTA"), 
#                            cell_type_publish, "OTHER"), 
#          size_interest = str_detect(cell_type_publish, "DA_SN|DA_VTA")) %>%
#   mutate(interest = str_replace(interest, "DA_SN", "SN"), 
#          interest = str_replace(interest, "DA_VTA", "VTA")) %>%
#   mutate(interest = factor(interest), 
#          interest = relevel(interest, ref = "OTHER")) %>%
#   arrange(interest) %>%
#   ggplot(aes(x = x, 
#              y = -y, 
#              colour = interest, 
#              size = size_interest, 
#              alpha = size_interest)) +
#   geom_point() +
#   labs(colour = "", 
#        title = "DOPAMINERGIC NEURONS") +
#   scale_x_continuous(limits = c(4250, 25000)) +
#   scale_y_continuous(limits = c(-15500, -2000)) +
#   theme_cowplot() +
#   scale_color_manual(values = c("grey", pal_d3()(5))) +
#   scale_size_manual(values = c(0.02, 1), guide = 'none') +
#   scale_alpha_manual(values = c(0.4, 1), guide = 'none') +
#   theme(axis.text = element_blank(), 
#         axis.line = element_blank(), 
#         axis.title = element_blank(), 
#         axis.ticks = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         legend.key.size = unit(1, 'cm'))
# 
# # ggarrange(p_spatial_da, 
# #           p_spatial_hipp, 
# #           p_spatial_thal, 
# #           p_spatial_pvalb,
# #           nrow = 2, 
# #           ncol = 2,
# #           align = 'v')
# 
# p_spatial_da %>%
#   ggsave(filename = "output/plots/figure_1/spatial_da.png", 
#        bg = "transparent", 
#        width = 7, height = 4, 
#        dpi = 300)

# markers ----

files <- list.files("output/markers/general", 
                    full.names = T)
names(files) <- str_extract(files, "(?<=general/)[:graph:]+(?=.csv)")

markers <- lapply(names(files), function(n){
  read_csv(files[[n]]) %>%
    filter(pvals_adj < 0.01) %>%
    arrange(pvals_adj, desc(logfoldchanges)) %>%
    mutate(.before = everything(), 
           cell_type = n)
    # slice_head(n = 1) %>%
    # pull(names)
}) %>%
  bind_rows

markers %>%
  filter(cell_type == "OLIG")
