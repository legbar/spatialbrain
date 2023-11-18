library(tidyverse)
library(parallel)
library(cowplot)
library(ggsci)
library(ggpubr)
library(ggrepel)

spatial_all_cells <- read_csv("input/spatial_obs.csv") %>%
  group_by(cell_type, sample_name) %>%
  filter(n() > 100) %>%
  ungroup %>% 
  filter(!str_detect(cell_type, "small|lowerComplexity|iscard"))
  # mutate(cell_type_publish = factor(cell_type), 
  #        cell_type_publish = relevel(cell_type_publish, ref = "midbrainNeuronalAndCortexNonneuronal_other_dopaminergic_large")) %>%
  # arrange(cell_type_publish)

spatial_plot_aes <- tibble(cell_type = spatial_all_cells$cell_type %>% unique) %>%
  mutate(size = 0.1, alpha = 0.5, colour = 'blue') %>%
  # mutate(cell_type = ifelse(
  #   str_detect(cell_type, "opamine"),
  #   "DOPAMINERGIC",
  #   cell_type
  # )) %>%
  distinct

pt_size = 0.5

spatial_plot_aes$cell_type

# DA
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "opaminer"),]$size <- 1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "opaminer"),]$alpha <- 1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "opaminer"),]$colour <- 'green'

# THAL
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "halamus"),]$size <- 1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "halamus"),]$alpha <- 0.7
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "halamus"),]$colour <- 'red'

# # THAL
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "midbrainNeuronalAndCortexNonneuronal_other_other_GABAinterneurons_pvalb"),]$size <- 0.5
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "midbrainNeuronalAndCortexNonneuronal_other_other_GABAinterneurons_pvalb"),]$alpha <- 0.5
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "midbrainNeuronalAndCortexNonneuronal_other_other_GABAinterneurons_pvalb"),]$colour <- 'lightblue'

# CA1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "CA1"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "CA1"),]$alpha <- 1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "CA1"),]$colour <- 'green4'

# DG
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "entate"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "entate"),]$alpha <- 1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "entate"),]$colour <- 'steelblue2'

# CA3
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "CA3"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "CA3"),]$alpha <- 1
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "CA3"),]$colour <- 'orange'

# CORTICAL_OUTER_EXCIT_NEU
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "Tshz2"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "Tshz2"),]$alpha <- 0.4
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "Tshz2"),]$colour <- 'grey81'

# CORTICAL INNER
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_outer_retained_cortex_layer4Rorb_layer2slash3slash5"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_outer_retained_cortex_layer4Rorb_layer2slash3slash5"),]$alpha <- 0.75
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_outer_retained_cortex_layer4Rorb_layer2slash3slash5"),]$colour <- 'dodgerblue'

# RN_PVALB_SNCG_NEU
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_outer_retained_cortex_layer4Rorb_layer4"),]$size <- 0.5
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_outer_retained_cortex_layer4Rorb_layer4"),]$alpha <- 0.75
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_outer_retained_cortex_layer4Rorb_layer4"),]$colour <- 'lightcoral'

# Dorsal subic
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "optic"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "optic"),]$alpha <- 0.5
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "optic"),]$colour <- 'slateblue1'

# Ventral subic
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "mygdala"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "mygdala"),]$alpha <- 0.75
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "mygdala"),]$colour <- 'yellow'

spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_middle_dorsalSubicAndLayer5Etv1"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_middle_dorsalSubicAndLayer5Etv1"),]$alpha <- 0.75
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "cortexNeuronal_excitatory_middle_dorsalSubicAndLayer5Etv1"),]$colour <- 'red'

spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "l5slash6CtgfFoxp2Cdh18"),]$size <- pt_size
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "l5slash6CtgfFoxp2Cdh18"),]$alpha <- 0.75
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "l5slash6CtgfFoxp2Cdh18"),]$colour <- 'limegreen'

# Meningeal glia
# spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "^SUBIC V"),]$size <- 0.5
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "ening"),]$alpha <- 0.75
spatial_plot_aes[str_detect(spatial_plot_aes$cell_type, "ening"),]$colour <- 'cornsilk'

spatial_colours <- spatial_plot_aes$colour
names(spatial_colours) <- spatial_plot_aes$cell_type
spatial_sizes <- spatial_plot_aes$size
names(spatial_sizes) <- spatial_plot_aes$cell_type
spatial_alphas <- spatial_plot_aes$alpha
names(spatial_alphas) <- spatial_plot_aes$cell_type

spatial_all_cells %>% 
  # filter(sample_name == "OLD_OVX_1") %>%

  # filter(sample_name == "OLD_OVX_5") %>%
  # mutate(cell_type_publish = ifelse(
  #   str_detect(cell_type, "opaminergic"),
  #   "DOPAMINERGIC",
  #   cell_type
  # )) %>%
  # mutate(cell_type_publish = cell_type) %>%
  # mutate(cell_type_publish = str_remove(cell_type_publish, "NEU")) %>%
  # mutate(cell_type_publish = str_replace_all(cell_type_publish, "_", " ")) %>%
  ungroup %>%
  # group_by(cell_type_publish) %>% tally
  mutate(cell_type = factor(cell_type), 
         cell_type = relevel(cell_type, ref = "midbrainNeuronalAndCortexNonneuronal_other_dopaminergic_large")) %>%
  arrange(cell_type) %>%
  # arrange(size) %>%
  ggplot(aes(x = x, 
             y = -y, 
             colour = cell_type, 
             size = cell_type, 
             alpha = cell_type)) +
  geom_point(position = position_jitter(height = 0.1, width = 0.1)) +
  theme_cowplot() +
  facet_wrap(vars(sample_name)) +
  # scale_color_d3(palette = "category20") +
  scale_color_manual(values = spatial_colours) +
  scale_size_manual(values = spatial_sizes) +
  # scale_size_manual(values = spatial_sizes*4) + # for A1 image
  scale_alpha_manual(values = spatial_alphas) +
  # guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
  # labs(colour = "") +
  # scale_color_manual(values = sample(pal_d3(palette = 'category20')(20), size = 20)) +
  # facet_wrap(vars(cell_type_publish)) +
  # scale_x_continuous(limits = c(5, 30)) +
  # facet_wrap(vars(cell_type_publish)) +
  # panel_border() +
  # facet_wrap(vars(mouse_id), scales = "free") +
  theme(legend.position = "none",
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(size = 6), 
        plot.background = element_rect(fill = 'black')) 

ggsave("output/archive/plots/spatial_hetero.png", 
       bg = "black", dpi = 300, width = 40, height = 30)
