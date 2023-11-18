library(tidyverse)
pksapply <-
  purrr::partial(base::sapply, simplify = F, USE.NAMES = T)
setwd("/laune_zfs/scratch/peter/f_active/paper_23/")

markers <- read_csv("input/05_annotation/markers_per_cell_type.csv")

markers %>%
  group_by(`Cell Type`) %>%
  slice_max(scores, n = 25) %>% 
  write_csv("input/top_markers_per_cell_type.csv")
