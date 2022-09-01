library(tidyverse)
library(cowplot)
library(ggrepel)

result <- read_csv("input/markers/malat1_da.csv")

result %>%
  ggplot(aes(x = Malat1_high_logfoldchanges, 
             y = -log10(Malat1_high_pvals_adj), 
             label = Malat1_high_names)) +
  geom_point(aes(colour = Malat1_high_logfoldchanges)) +
  geom_label_repel() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = expression(Log[2] ~ Fold ~ Change), 
       y = expression(-Log[10] ~ Adjusted ~ italic(P)), 
       title = "Malat1-positive DA neuron markers") +
  theme(legend.position = "none")
