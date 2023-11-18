library(cowplot)
library(ggsci)
library(ggh4x)
library(patchwork)
library(kableExtra)
library(tidyverse)

pksapply <-
  purrr::partial(base::sapply, simplify = F, USE.NAMES = T)

homologs <- readRDS("input/trap/homologs_gwas.rds")

markers <- read_csv("input/05_annotation/markers_per_cell_type.csv")
colnames(markers)[1] <- "cell_type"

markers <- markers %>%
  filter(names %in% homologs$mmusculus_homolog_associated_gene_name)

markers %>% 
  filter(!str_detect(cell_type, "ERYTH")) %>%
  filter(pvals_adj < 0.01) %>%
  group_by(cell_type) %>%
  # summarise(n = n()) %>%
  # arrange(n)
  # filter(n() > 5) %>%
  # summarise(summary = mean(logfoldchanges * -log10(pvals_adj))) %>%
  # summarise(summary = mean(scores)) %>%
  mutate(summary = mean(scores)) %>%
  ungroup %>%
  mutate(cell_type = factor(cell_type), 
         cell_type = fct_reorder(cell_type, summary)
         ) %>%
  mutate(neu_glia = ifelse(str_detect(cell_type, "ASTRO|OLIGO|OPC|MENINGES|GLIA"), 
                           "GLIAL", "NEURONAL"), 
         summary_alpha = rank(summary), 
         summary_alpha = summary_alpha / max(summary_alpha)) %>%
  ggplot(aes(x = cell_type, 
             y = logfoldchanges, 
             # alpha = -log10(pvals_adj),
             # alpha = summary_alpha,
             # fill = cell_type %in% c("DOPA SN", "DOPA VTA"),
             fill = neu_glia
             )) +
  # geom_col(color = "black") +
  geom_boxplot() +
  # geom_jitter() +
  scale_y_continuous(limits = c(-5, 5)) +
  scale_fill_manual(values = pal_d3()(10)[c(6, 10)]) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # geom_boxplot()
  # scale_fill_manual(values = c("grey", pal_d3()(10)[4])) +
  theme_cowplot(10) +
  labs(x = "", y = "Mean z-score") +
  theme(axis.text = element_text(angle = 45, hjust = 1, size = 14), 
        axis.title = element_text(size = 14),
        legend.position = "none", 
        plot.margin = margin(7, 7, 7, 12, unit = "pt"))

ggsave("output/plots/figure_2/stereoseq_gwas.png", 
       bg = "transparent", 
       width = 110*2, height = 70*2,
       units = "mm",
       dpi = 300)
