# Sys.setenv("MC_CORES"=18L)
# library(parallel)
# options(mc.cores = 18)
library(cowplot)
library(ggsci)
# library(Matrix)
# library(SingleCellExperiment)
# library(scater)
# library(MAST)
# library(Libra)
# library(Seurat)
library(ggh4x)
# library(ggrepel)
library(patchwork)
library(kableExtra)
library(tidyverse)

pksapply <-
  purrr::partial(base::sapply, simplify = F, USE.NAMES = T)

# load spatial obs ----
spatial_all_cells <-
  read_csv("input/05_annotation/spatial_all_cells.csv")
colnames(spatial_all_cells)[1] <- "cell_id"

# load DE results ----
genotype_pb_de <-
  read_csv("input/05_annotation/DE_pseudobulk_genotype.csv")
genotype_young_pb_de <-
  read_csv("input/05_annotation/DE_pseudobulk_genotype_by_age_YOUNG.csv")
genotype_old_pb_de <-
  read_csv("input/05_annotation/DE_pseudobulk_genotype_by_age_OLD.csv")
age_pb_de <- read_csv("input/05_annotation/DE_pseudobulk_age.csv")

meta <-
  tibble(cell_type = sort(unique(spatial_all_cells$cell_type)))

de_res <- list("GENOTYPE" = genotype_pb_de, 
     "GENOTYPE - YOUNG" = genotype_young_pb_de, 
     "GENOTYPE - OLD" = genotype_old_pb_de,
     "AGE" = age_pb_de)

de_res_n_per_cell_type <- pksapply(names(de_res), function(n){
  meta %>%
    left_join({
      de_res[[n]] %>%
        filter(gene != "Snca", 
               padj < 0.05) %>%
        group_by(cell_type) %>%
        summarise(n = n())
    }) %>%
    replace_na(list(n = 0))
})

pksapply(names(de_res_n_per_cell_type), function(n) {
  
  c <- list("AGE" = "red", 
            "GENOTYPE" = "darkgreen", 
            "GENOTYPE - YOUNG" = "darkgreen", 
            "GENOTYPE - OLD" = "darkblue")
  
  # a <- list("AGE" = 1, 
  #           "GENOTYPE" = 1, 
  #           "GENOTYPE - YOUNG" = 1, 
  #           "GENOTYPE - OLD" = 0.75)
  
  de_res_n_per_cell_type[[n]] %>%
    full_join(select(spatial_all_cells, sample_name, cell_type, x, y)) %>%
    filter(sample_name == "YOUNG_WT_4") %>%
    mutate(interest = n > 0,
           interest = factor(interest)) %>%
    arrange(desc(interest)) %>%
    ggplot(aes(
      x = x,
      y = -y,
      color = n,
      alpha = n
    )) +
    geom_point() +
    theme_cowplot(8) +
    scale_color_gradient(
      low = "lightgrey",
      high = c[[n]],
      guide = "legend",
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_alpha_continuous(
      range = c(0.2, 1),
      guide = "legend",
      breaks = scales::pretty_breaks(n = 4)
    ) +
    labs(title = n,
         color = "Number of\nDE genes",
         alpha = "Number of\nDE genes") +
    theme(
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 32, hjust = 0.5)
    ) +
    theme(
      legend.key.size = unit(2, 'cm'),
      #change legend key size
      legend.key.height = unit(2, 'cm'),
      #change legend key height
      legend.key.width = unit(2, 'cm'),
      #change legend key width
      legend.title = element_text(size = 24, hjust = 0.5),
      #change legend title font size
      legend.text = element_text(size = 24)
    ) + #change legend text font size
    guides(colour = guide_legend(override.aes = list(size = 10)))
  
  # ggsave(paste0("output/plots/figure_4/pb_de_", n, ".png"),
  #        width = 14.5, height = 9, dpi = 300)
})

# genotype supplementary volcanos ----


# young
ct <- de_res_n_per_cell_type$`GENOTYPE - YOUNG` %>%
  slice_max(order_by = n, n = 5) %>%
  pull(cell_type)

genotype_young_pb_de %>%
  filter(cell_type %in% ct) %>%
  filter(gene != "Snca") %>%
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05, 
                 alpha = padj < 0.05)) +
  geom_text_repel(data = filter(genotype_young_pb_de, 
                                cell_type %in% ct & padj < 0.05 & gene != "Snca"), 
                  aes(label = gene), 
                  size = 5) +
  facet_wrap(vars(cell_type)) +
  theme_cowplot() +
  scale_color_d3() +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_x_continuous(limits = c(-3, 3), oob = scales::squish) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Log2 Fold Change", 
       y = "-log10 Adjusted P Value", 
       color = "Adj P < 0.1",
       title = "Genotype: Young Brains") +
  theme(legend.position = "top", 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))

ggsave("output/plots/supp_axon_genotype/genotype_young_volcano.png", 
       width = 12, height = 12, dpi = 300)

# old
ct <- de_res_n_per_cell_type$`GENOTYPE - OLD` %>%
  slice_max(order_by = n, n = 5) %>%
  pull(cell_type)

genotype_old_pb_de %>%
  filter(cell_type %in% ct) %>%
  filter(gene != "Snca") %>%
  ggplot(aes(x = log2FoldChange, 
             y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05, 
                 alpha = padj < 0.05)) +
  geom_text_repel(data = filter(genotype_old_pb_de, 
                                cell_type %in% ct & padj < 0.05 & gene != "Snca"), 
                  aes(label = gene), 
                  size = 5) +
  facet_wrap(vars(cell_type)) +
  theme_cowplot() +
  scale_color_d3() +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_x_continuous(limits = c(-2, 2), oob = scales::squish) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Log2 Fold Change", 
       y = "-log10 Adjusted P Value", 
       color = "Adj P < 0.1",
       title = "Genotype: Old Brains") +
  theme(legend.position = "top", 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))

ggsave("output/plots/supp_axon_genotype/genotype_old_volcano.png", 
       width = 12, height = 12, dpi = 300)

pksapply(names(de_res_n_per_cell_type), function(n) {
  
  c <- list("AGE" = "red", 
            "GENOTYPE" = "darkgreen", 
            "GENOTYPE - YOUNG" = "darkgreen", 
            "GENOTYPE - OLD" = "darkblue")
  
  x <- de_res_n_per_cell_type[[n]]
  
  p1 <- x %>%
    full_join(select(spatial_all_cells, sample_name, cell_type, x, y)) %>%
    filter(sample_name == "YOUNG_WT_4") %>%
    mutate(interest = n > 0,
           interest = factor(interest)) %>%
    arrange(desc(interest)) %>%
    ggplot(aes(
      x = x,
      y = -y,
      color = n,
      alpha = n
    )) +
    geom_point() +
    theme_cowplot(8) +
    scale_color_gradient(
      low = "lightgrey",
      high = c[[n]],
      guide = "legend",
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_alpha_continuous(
      range = c(0.2, 1),
      guide = "legend",
      breaks = scales::pretty_breaks(n = 4)
    ) +
    labs(title = n,
         color = "Number of\nDE genes",
         alpha = "Number of\nDE genes") +
    theme(
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 42, hjust = 0.5)
    ) +
    theme(
      legend.position = "none"
      # legend.key.size = unit(2, 'cm'),
      # #change legend key size
      # legend.key.height = unit(2, 'cm'),
      # #change legend key height
      # legend.key.width = unit(2, 'cm'),
      # #change legend key width
      # legend.title = element_text(size = 24, hjust = 0.5),
      # #change legend title font size
      # legend.text = element_text(size = 24)
    ) #change legend text font size
    # guides(colour = guide_legend(override.aes = list(size = 10)))
  
  p2 <- slice_max(x, order_by = n, n = 5) %>%
    mutate(cell_type = factor(cell_type), 
           cell_type = fct_reorder(cell_type, n)) %>%
    ggplot(aes(x = n, 
               y = cell_type)) +
    geom_col(fill = c[[n]], color = "black") +
    geom_text(aes(label=n), vjust=0.5, hjust = -0.25, fontface = "bold", size = 12) +
    theme_cowplot(12) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "Number of DE genes", 
         y = "") +
    theme(axis.text = element_text(size = 36, face = "bold"), 
          axis.title = element_text(size = 42, face = "bold"), 
          axis.line = element_line(linewidth = 2))
  
  p1 / (plot_spacer() + p2 + plot_spacer() + plot_layout(widths = c(0.5, 4, 0.5)))
  
  ggsave(paste0("output/plots/figure_4/pb_de_", n, ".png"),
         width = 14.5, height = 18, dpi = 300)
  
})

# gprofiler ----

top_de_per_ct <- pksapply(names(de_res_n_per_cell_type), function(n) {
  
  x <- de_res_n_per_cell_type[[n]]
  
  x <- slice_max(x, order_by = n, n = 5) %>%
    mutate(cell_type = factor(cell_type), 
           cell_type = fct_reorder(cell_type, n))
  
})

# age_ct_most_de <- age_pb_de_n_per_cell_type %>% 
#   mutate(Comparison = "AGE") %>%
#   slice_max(order_by = n, n = 5) %>%
#   pull(cell_type)
# 
# age_pb_de %>%
#   filter(cell_type == "AMYGDALA")

temp <- de_res$AGE %>%
  filter(cell_type %in% top_de_per_ct$AGE$cell_type) %>%
  arrange(cell_type, pvalue) %>%
  group_by(cell_type) %>%
  group_split() %>%
  pksapply(function(x){
    bg = x %>%
      filter(!is.na(padj)) %>%
      pull(gene)
    g <- x %>%
      filter(padj < 0.05) %>%
      arrange(stat) %>%
      pull(gene) %>%
      gprofiler2::gost(organism = "mmusculus",
                       ordered_query = T,
                       significant = F,
                       correction_method = "fdr",
                       sources = "GO:BP",
                       exclude_iea = T,
                       # custom_bg = bg,
                       evcodes = T)
    g$result %>%
      filter(term_size <= 500) %>%
      mutate(.before = everything(),
             cell_type = x$cell_type[1])
  })

library(kableExtra)

order <- de_res_n_per_cell_type$AGE %>%
  slice_max(order_by = n, n = 5) %>%
  pull(cell_type)

t <- bind_rows(temp) %>%
  filter(!(cell_type == "EX L1" & str_detect(term_name, "inflammatory"))) %>%
  group_by(cell_type) %>%
  slice_head(n = 3) %>%
  mutate(term_name = str_to_sentence(term_name)) %>%
  select("Cell Type" = cell_type, 
         "Top Aging GO:BP Categories per Cell Type" = term_name) %>%
  ungroup %>%
    mutate(`Cell Type` = factor(`Cell Type`, levels = order)) %>%
  arrange(`Cell Type`)
  
t %>%
  select(-`Cell Type`) %>%
  kbl() %>%
  row_spec(0, bold = T) %>%
  pack_rows(index = table(t$`Cell Type`), indent = T) %>%
  # column_spec(1, bold = T, border_right = F) %>%
  kable_classic(full_width = F)

pksapply(names(de_res), function(n){
  de_res$AGE %>%
    filter(cell_type %in% top_de_per_ct$AGE$cell_type, 
           padj < 0.05)
})

go_age <- age_pb_de %>%
  filter(cell_type %in% age_ct_most_de) %>%
  group_by(cell_type) %>%
  group_split() %>%
  pksapply(function(x){
    bg = x %>%
      filter(!is.na(padj)) %>%
      pull(gene)
    g <- x %>%
      filter(padj < 0.05) %>%
      arrange(stat) %>%
      pull(gene) %>%
      gprofiler2::gost(organism = "mmusculus",
           ordered_query = T,
           significant = F,
           correction_method = "fdr",
           sources = "GO:BP",
           exclude_iea = T,
           # custom_bg = bg,
           evcodes = T)
    g$result %>%
      filter(term_size <= 500) %>%
      mutate(.before = everything(),
             cell_type = x$cell_type[1])
  })
names(go_age) <- pksapply(go_age, function(x){x$cell_type[1]})

# Microglia spatial ----
spatial_all_cells %>%
  filter(cell_type == "MICROGLIA") %>%
  group_by(sample_name) %>%
  summarise(n = n())

select(spatial_all_cells, sample_name, age, cell_type, x, y) %>%
  filter(sample_name %in% c("YOUNG_OVX_3", "OLD_OVX_2")) %>%
  mutate(interest = ifelse(cell_type == "MICROGLIA", "MICROGLIA", "Other"),
         interest = factor(interest, levels = c("Other", "MICROGLIA")), 
         age = factor(age, levels = c("YOUNG", "OLD"))) %>%
  arrange(interest) %>%
  ggplot(aes(
    x = x,
    y = -y,
    color = interest,
    size = interest,
    alpha = interest
  )) +
  geom_point() +
  theme_cowplot(8) +
  facet_grid(rows = vars(age)) +
  panel_border(color = "black") +
  scale_alpha_manual(
    values = c(0.2, 1),
    guide = "legend",
    # breaks = scales::pretty_breaks(n = 4)
  ) +
  scale_color_manual(values = c("lightgrey", pal_d3()(10)[4])) +
  scale_size_manual(values = c(1, 4)) +
  labs(title = "MICROGLIA") +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 46, hjust = 0.5),
    strip.text = element_text(size = 46)
  )

ggsave(paste0("output/plots/figure_4/microglia_spatial_age.png"),
       width = 14.5, height = 18, dpi = 300)

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

dopa_count1 <- spatial_all_cells %>%
  filter(cell_type %in% c("DOPA SN", "DOPA VTA")) %>%
  group_by(cell_type, age, sample_name) %>%
  summarise(n = n()) %>%
  summarise(n = mean(n))

dopa_count2 <- spatial_all_cells %>%
  filter(cell_type %in% c("DOPA SN", "DOPA VTA")) %>%
  group_by(cell_type, age, sample_name) %>%
  summarise(n = n())

library(ggpubr)

ggplot(dopa_count1, aes(x = age, y = n, fill = age)) +
  geom_col(position = "dodge", color = "black") +
  geom_jitter(data = dopa_count2, aes(group = age), position = position_dodge(width = 0.9)) +
  geom_bracket(
    data = filter(dopa_count1, cell_type == "DOPA SN"),
    xmin = "YOUNG", xmax = "OLD", y.position = 275,
    label = "*", tip.length = c(0.15, 0.2),
    size = 1, label.size = 12
  ) +
  scale_fill_manual(values = pal_d3()(10)[c(1, 2)]) +
  facet_wrap(vars(cell_type)) +
  theme_cowplot() +
  panel_border(color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "none") +
  labs(x = "Cell Type", y = "Number of cells per brain", 
       fill = "Age")

ggsave(
  "output/plots/figure_4/cell_types_ageing_DOPA.png",
  bg = "transparent",
  width = 130,
  height = 110,
  units = "mm",
  dpi = 300
)

bind_rows(lapply(masc_all,function(x){x[[1]]}),.id='name') %>%
  filter(term == "ageOLD") %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  mutate(cl_num = as.numeric(str_extract(name, "(?<=cluster)[:digit:]{1,2}"))) %>%
  left_join({
    spatial_all_cells %>%
      select("cl_num" = cell_type_num,
             cell_type) %>%
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
