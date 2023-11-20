# RESULTS 3 ----
# ----
# ----
# ----


# SETUP ----
n <- dplyr::n
# load salmon data --------------------------------------------------------

# files <-
#   list.files(
#     "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/salmon",
#     pattern = "quant.sf",
#     recursive =  T,
#     full.names = T
#   )
#
# all(file.exists(files)) # everything exists
#
# names(files) <-
#   str_extract(files, "(?<=salmon\\/)[0-9|\\_]+(?=\\/quant.sf)") #name the files by the fastq_id from the enclosing folder
#
# files <- files[names(files) %in% metadata$fastq_id] # remove blanks
# files <-
#   files[order(match(names(files), metadata$fastq_id))] #reorder based on metadata
#
# files <- readRDS("R/objects/files.rds")
#
# metadata$files <- files
# metadata$names <- metadata$fastq_id
#
# se <- readRDS(file = "R/objects/se.rds")
# # if (!exists("se")) {
# #   se <- tximeta::tximeta(as.data.frame(metadata))
# # }
# # saveRDS(se, file = "R/objects/se.rds")
#
# gse <- readRDS("R/objects/gse.rds")
# # if (!exists("gse")) {
# #   gse <- tximeta::summarizeToGene(se)
# # }
# # saveRDS(gse, file = "R/objects/gse.rds")
#
#
# all(colnames(gse) == metadata$fastq_id)
#
# dds_salmon <- DESeqDataSet(gse, design = ~ 1)
#
# rownames(dds_salmon) <- gsub("\\..*$", "", rownames(dds_salmon))
# rowData(dds_salmon)$gene_id <-
#   gsub("\\..*$", "", rowData(dds_salmon)$gene_id)
#
#
# # anno_salmon <- get_anno(rownames(dds_salmon))
#
# dds_salmon <- collapseReplicates(dds_salmon,
#                                  groupby = dds_salmon$sample_name,
#                                  run = dds_salmon$lane_code)
#
# dds_salmon <- dds_salmon[rowSums(counts(dds_salmon)) > 0,]

dds_salmon <- readRDS("R/objects/dds_salmon_collapsed.rds")
# remove unused objects
rm(list = c("se", "gse"))


# ----
# ----
# ----
# END OF SETUP ----
# ----
# ----
# ----


# LIBRARY METRICS ----
# Make and enrichment tibble ----
# Make an enrichment status tibble
enrichment_status <- tibble(
  ensembl_gene_id = rownames(dds_C2_IP),
  enrichment = ifelse(
    ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
    ifelse(
      ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
      "Soma and Axon",
      "Soma only"
    ),
    ifelse(
      ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
      "Axon only",
      "No enrichment"
    )
  ),
) %>%
  mutate(
    enrichment = factor(
      enrichment,
      levels = c("Soma only",
                 "Axon only",
                 "Soma and Axon",
                 "No enrichment")
    )
  )

# Make a library proportion tibble - by compartment ----
# Make a library proportion tibble - split by compartment
library_proportion_info <- enrichment_status %>%
  left_join(enframe(rowMeans(counts(dds_C2_IP)[, colData(dds_C2_IP)$compartment == "MB"]),
                    name = "ensembl_gene_id",
                    value = "mb")) %>%
  left_join(enframe(rowMeans(counts(dds_C2_IP)[, colData(dds_C2_IP)$compartment == "AXON"]),
                    name = "ensembl_gene_id",
                    value = "axon"))


library_proportion_info_total <- enrichment_status %>%
  left_join(enframe(rowMeans(counts(dds_C2_TOTAL)[, colData(dds_C2_TOTAL)$compartment == "MB"]),
                    name = "ensembl_gene_id",
                    value = "mb")) %>%
  left_join(enframe(rowMeans(counts(dds_C2_TOTAL)[, colData(dds_C2_TOTAL)$compartment == "AXON"]),
                    name = "ensembl_gene_id",
                    value = "axon"))

# Plot just the proportion of library in each enrichment category by sample type (axon or mb) ----


p_library_proportion_compartment <- bind_rows({
  library_proportion_info %>%
    pivot_longer(c(mb, axon),
                 names_to = "compartment",
                 values_to = "count") %>%
    mutate(compartment = ifelse(compartment == "mb", "Soma", "Axon")) %>%
    group_by(compartment,
             enrichment) %>%
    summarise(count = sum(count)) %>%
    mutate(prop = count / sum(count),
           fraction = "TRAP")
}, {
  library_proportion_info_total %>%
    drop_na %>%
    pivot_longer(c(mb, axon),
                 names_to = "compartment",
                 values_to = "count") %>%
    mutate(compartment = ifelse(compartment == "mb", "Soma", "Axon")) %>%
    group_by(compartment,
             enrichment) %>%
    summarise(count = sum(count)) %>%
    mutate(prop = count / sum(count),
           fraction = "TOTAL")
}) %>%
  ggplot(aes(x = compartment,
             y = prop,
             fill = enrichment)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "TRAP Fraction and Compartment",
       y = "Percentage of Library",
       fill = "Enrichment") +
  geom_text(
    aes(
      x = compartment,
      y = prop,
      label = ifelse(prop > 0.05,
                     scales::percent(prop,
                                     accuracy = 1),
                     "")
    ),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = 2
  ) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(
    nrow = 1,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  ))  +
  facet_wrap(vars(fraction)) +
  panel_border()


# Plot the percentage of soma enrichment per cohort and age in soma TRAP samples ----------


p_mb_ip_age_library_composition <- bind_rows({
  enrichment_status %>%
    left_join(enframe(rowMeans(counts(dds_C1_MB_IP)[, colData(dds_C1_MB_IP)$age == "YOUNG"]),
                      name = "ensembl_gene_id",
                      value = "Young")) %>%
    left_join(enframe(rowMeans(counts(dds_C1_MB_IP)[, colData(dds_C1_MB_IP)$age == "OLD"]),
                      name = "ensembl_gene_id",
                      value = "Old")) %>%

    pivot_longer(c(Young, Old),
                 names_to = "age",
                 values_to = "count") %>%
    mutate(cohort = "C1")
}, {
  enrichment_status %>%
    left_join(enframe(rowMeans(counts(dds_C2_MB_IP)[, colData(dds_C2_MB_IP)$age == "YOUNG"]),
                      name = "ensembl_gene_id",
                      value = "Young")) %>%
    left_join(enframe(rowMeans(counts(dds_C2_MB_IP)[, colData(dds_C2_MB_IP)$age == "OLD"]),
                      name = "ensembl_gene_id",
                      value = "Old")) %>%

    pivot_longer(c(Young, Old),
                 names_to = "age",
                 values_to = "count") %>%
    mutate(cohort = "C2")
}, {
  enrichment_status %>%
    left_join(enframe(rowMeans(counts(dds_C3_MB_IP)[, colData(dds_C3_MB_IP)$age == "YOUNG"]),
                      name = "ensembl_gene_id",
                      value = "Young")) %>%
    left_join(enframe(rowMeans(counts(dds_C3_MB_IP)[, colData(dds_C3_MB_IP)$age == "OLD"]),
                      name = "ensembl_gene_id",
                      value = "Old")) %>%

    pivot_longer(c(Young, Old),
                 names_to = "age",
                 values_to = "count") %>%
    mutate(cohort = "C3")
}) %>%
  mutate(age = factor(age, levels = c("Young", "Old")),
         enrichment = ifelse(str_detect(enrichment, "oma"), "Enriched", "Not enriched"),
         enrichment = factor(enrichment, levels = c("Not enriched", "Enriched"))) %>%
  drop_na %>%
  group_by(age,
           cohort,
           enrichment) %>%
  summarise(count = sum(count)) %>%
  mutate(prop = count / sum(count)) %>%
  ggplot(aes(x = age,
             y = prop,
             fill = enrichment)) +
  geom_col(colour = "black") +
  scale_fill_manual(values = pal_d3()(2)[c(2, 1)]) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "TRAP Cohort and Age",
       y = "Percentage of Library",
       fill = "Enrichment") +
  geom_text(
    aes(
      x = age,
      y = prop,
      label = ifelse(prop > 0.05,
                     scales::percent(prop,
                                     accuracy = 1),
                     "")
    ),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = 2
  ) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(
    nrow = 1,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5
  )) +
  facet_wrap(vars(cohort))



# The number of genes in each category ----------
p_library_composition_enrichment_n_genes <- library_proportion_info %>%
  filter(enrichment != "No enrichment") %>%
  group_by(enrichment) %>%
  tally %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = enrichment,
             y = n,
             fill = enrichment)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(label = n),
            vjust = -0.5,
            fontface = 2) +
  geom_text(aes(label = scales::percent(prop),
                y = n/2),
            colour = "white",
            fontface = 2) +
  labs(x = "Enrichment",
       y = "Number of Genes") +
  theme(legend.position = "none")




# Plot the average expression count of each enrichment category ----
library_proportion_info %>%
  filter(enrichment != "No enrichment") %>%
  pivot_longer(c(mb, axon),
               names_to = "compartment",
               values_to = "count") %>%
  ggplot(aes(x = enrichment,
             y = count)) +
  scale_y_log10() +
  geom_violin(scale = "width",
              draw_quantiles = 0.5) +
  facet_wrap(vars(compartment), nrow = 1)

# ----
# ----
# ----
# AXON AGE - LIBRARY COMPOSITION CHECK ----


C2_AXON_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
      select(age, "Slc6a3" = count)
  }, {
    plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
      select("Th" = count)
  }, {
    #   plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
    #     select("Ddc" = count)
    # }, {
    plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-age,
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C2")


C3_AXON_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
      select(age, "Slc6a3" = count)
  }, {
    plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
      select("Th" = count)
  }, {
    #   plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
    #     select("Ddc" = count)
    # }, {
    plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-age,
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C3")



p_axon_ip_age_library_composition <- bind_rows(
  C2_AXON_IP_MARKERS,
  C3_AXON_IP_MARKERS) %>%
  mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Slc18a2"),
                         "Dopaminergic", "Other")) %>%
  ggplot(aes(x = age,
             y = log2(count + 1),
             colour = marker,
             group = gene)) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(cohort), nrow = 1) +
  panel_border() +
  scale_color_d3() +
  labs(x = "Cohort and Age",
       y = "Log2 Expression Count",
       colour = "Cell Type Marker") +
  theme(legend.position = "top")






C2_AXON_TOTAL_MARKERS <- bind_cols(
  {
    plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
      select(age, "Slc6a3" = count)
  }, {
    plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
      select("Th" = count)
  }, {
    #   plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
    #     select("Ddc" = count)
    # }, {
    plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-age,
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C2")


C3_AXON_TOTAL_MARKERS <- bind_cols(
  {
    plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
      select(age, "Slc6a3" = count)
  }, {
    plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
      select("Th" = count)
  }, {
    #   plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
    #     select("Ddc" = count)
    # }, {
    plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-age,
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C3")


p_axon_total_age_library_composition <- bind_rows(C2_AXON_TOTAL_MARKERS,
                                                  C3_AXON_TOTAL_MARKERS) %>%
  mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
                         "Dopaminergic", "Other")) %>%
  ggplot(aes(x = age,
             y = log2(count + 1),
             colour = marker,
             group = gene)) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(cohort), nrow = 1) +
  panel_border() +
  scale_color_d3() +
  labs(x = "Age",
       y = "Log2 Expression Count",
       colour = "Cell Type Marker") +
  theme(legend.position = "top")







# AXON AGE - ALTERNATIVE COMPOSITION CHECK ----
# This code is commented because I don't know which genes are reliably axonal in striatum:
# If I look at the proportion of library that is enriched, that will include
# other cell types, so the proportion doesnt change so much
# The first composition check (plotting DA markers) is informative.
# LIBRARY_PROPORTIONS <- bind_rows(
#   {
#     as_tibble(counts(dds_C2_AXON_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C2")
#   }, {
#     as_tibble(counts(dds_C3_AXON_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C3")
#   }
# ) %>%
#   pivot_longer(-c(ensembl_gene_id, cohort),
#                names_to = "sample_name",
#                values_to = "count") %>%
#   mutate(status = ifelse(ensembl_gene_id %in% AXON_FRACTION_ENRICHED_GENES,
#                          "Enriched",
#                          ifelse(ensembl_gene_id %in% AXON_FRACTION_ENRICHED_GENES,
#                                 "Depleted", "Unchanged")),
#          status = factor(status, levels = c("Depleted", "Unchanged", "Enriched"))) %>%
#   mutate(age = ifelse(str_detect(sample_name, "YOUNG"), "YOUNG", "OLD"),
#          age = factor(age, levels = c("YOUNG", "OLD"))) %>%
#   group_by(cohort, age, status) %>%
#   summarise(count = sum(count, na.rm = T)) %>%
#   mutate(prop = count/sum(count))
#
# # DON'T KNOW WHICH TO USE: Don't know which genes are definitely axon specific in striatum.
# # AXON_FRACTION_ENRICHED_GENES
# # AXON_FRACTION_ENRICHED_ABACONF_GENES
# # AXON_FRACTION_ENRICHED_AGNOSTIC_GENES
# # MB_FRACTION_ENRICHED_GENES
#
# p_age_mb_library_proportions <- LIBRARY_PROPORTIONS %>%
#   ggplot(aes(x = age,
#              y = prop,
#              fill = status)) +
#   geom_col() +
#   facet_wrap(vars(cohort), nrow = 1)  +
#   theme(legend.position = "top") +
#   geom_text(aes(y = prop, label = ifelse(prop > 0.1, scales::percent(prop, accuracy = 1), "")),
#             position = position_stack(vjust = 0.5),
#             show.legend = FALSE,
#             colour = "white",
#             fontface = 2) +
#   scale_fill_manual(values = pal_d3()(4)[c(4, 1, 3)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
#                      labels = scales::percent) +
#   panel_border() +
#   labs(x = "Age and Cohort",
#        y = "Percentage of library",
#        fill = "Enrichment")

# AXON AGE TOTAL DESeq2 -----------------------------------------------------

# Using C2 and C3 for total to maximise number of axonal samples.
# Will check whether combining worsens detection power like TRAP samples

# Filter for genes that are reliably detected in both cohorts
dds_C2_AXON_TOTAL_AGE_FILTER <- filter_genes(dds_C2_AXON_TOTAL,
                                             grouping = c("age", "gene_id"))
dds_C3_AXON_TOTAL_AGE_FILTER <- filter_genes(dds_C3_AXON_TOTAL,
                                             grouping = c("age", "gene_id"))

genes_AXON_TOTAL_AGE <- intersect(dds_C2_AXON_TOTAL_AGE_FILTER,
                                  dds_C3_AXON_TOTAL_AGE_FILTER)

dds_C2_AXON_TOTAL_AGE <- dds_C2_AXON_TOTAL[rownames(dds_C2_AXON_TOTAL) %in% genes_AXON_TOTAL_AGE,]
dds_C3_AXON_TOTAL_AGE <- dds_C3_AXON_TOTAL[rownames(dds_C3_AXON_TOTAL) %in% genes_AXON_TOTAL_AGE,]

dds_C2_AXON_TOTAL_AGE@design <- ~ age
dds_C3_AXON_TOTAL_AGE@design <- ~ age

colData(dds_C2_AXON_TOTAL_AGE) <- droplevels(colData(dds_C2_AXON_TOTAL_AGE))
colData(dds_C3_AXON_TOTAL_AGE) <- droplevels(colData(dds_C3_AXON_TOTAL_AGE))

dds_C2_AXON_TOTAL_AGE <- DESeq(
  dds_C2_AXON_TOTAL_AGE,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

dds_C3_AXON_TOTAL_AGE <- DESeq(
  dds_C3_AXON_TOTAL_AGE,
  fitType = "local",
  sfType = "poscounts",
  minReplicatesForReplace = Inf,
  parallel = TRUE,
)

# dds_C3_AXON_TOTAL_AGE <- dds_C3_AXON_TOTAL_AGE[rowData(dds_C3_AXON_TOTAL_AGE)$baseMean > 500 ,]

# results

res_C2_AXON_TOTAL_AGE <- DESeq2::results(
  dds_C2_AXON_TOTAL_AGE,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_AXON_TOTAL_AGE <- lfcShrink(
  dds_C2_AXON_TOTAL_AGE,
  res = res_C2_AXON_TOTAL_AGE,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)

res_C3_AXON_TOTAL_AGE <- DESeq2::results(
  dds_C3_AXON_TOTAL_AGE,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_AXON_TOTAL_AGE <- lfcShrink(
  dds_C3_AXON_TOTAL_AGE,
  res = res_C3_AXON_TOTAL_AGE,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)

AXON_AGE_TOTAL_META <-
  as_tibble(res_C2_AXON_TOTAL_AGE, rownames = "ensembl_gene_id") %>%
  left_join(
    as_tibble(res_C3_AXON_TOTAL_AGE, rownames = "ensembl_gene_id"),
    by = "ensembl_gene_id",
    suffix = c("_C2", "_C3")
  ) %>%
  left_join(anno) %>%
  mutate(
    across(starts_with("pvalue"), ~ replace_na(.x, 1 - 1e-10)),
    # replace NA with *almost* 1 (it can't be perfectly 1)
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)),
    # replace 0 values with minimum non-zero of the column
    across(starts_with("pvalue"), ~ .x / 2),
    # divide 2-sided pvalue by 2 into 1-sided
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)) # again replace 0 values with minimum non-zero of the column
  ) %>%
  select(sort(names(.))) %>%
  select(
    ensembl_gene_id,
    external_gene_name,
    description,
    everything(),
    -starts_with("padj")
  ) %>%
  # filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
  drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
  rowwise() %>%
  mutate(
    pvalue_C3 = ifelse(
      sign(log2FoldChange_C3) == sign(log2FoldChange_C2),
      pvalue_C3,
      1 - pvalue_C3)
  ) %>%
  mutate(across(starts_with("pvalue"),
                ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
  mutate(
    # sumlog = sumlog(c_across(pvalue_C2:pvalue_C3))$p,
    # calculate sumlog
    sumz = sumz(c_across(pvalue_C2:pvalue_C3))$p,
    # calculate stouffers
    # log2FoldChange = mean(c(
    #   log2FoldChange_C2,
    #   log2FoldChange_C3
    # ))
    log2FoldChange = log2FoldChange_C2
  ) %>%
  ungroup() %>%
  mutate(
    # sumlog_adj = p.adjust(sumlog, method = "fdr"),
    # correction for multiple comparisons
    sumz_adj = p.adjust(sumz, method = "fdr"),
    conflict = sign(log2FoldChange_C2) != sign(log2FoldChange_C3)) %>%
  mutate(outcome = ifelse(sumz_adj < 0.1,
                          ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
                          "Unchanged"))

# AXON AGE TOTAL OUTCOMES ----
AXON_AGE_TOTAL_OUTCOMES <- AXON_AGE_TOTAL_META %>%
  select(ensembl_gene_id,
         outcome)



# AXON AGE TRAP DESeq2 -----------------------------------------------------

# Just using Cohort 2: The dispersion in cohort 3 is just too large.
# Also, only 7000 genes are measurable in both
# versus 17000 in C2 alone
# Combining cohort 2 and 3 by meta analysis results in fewer detections,
# and the results of cohort 2 actually looking realistic (DA markers downregulated in old)

dds_C2_AXON_IP_AGE_FILTER <- filter_genes(dds_C2_AXON_IP,
                                          grouping = c("age", "gene_id"))

dds_C2_AXON_IP_AGE <- dds_C2_AXON_IP[dds_C2_AXON_IP_AGE_FILTER,]

dds_C2_AXON_IP_AGE@design <- ~ collection + age

colData(dds_C2_AXON_IP_AGE) <- droplevels(colData(dds_C2_AXON_IP_AGE))

dds_C2_AXON_IP_AGE <- DESeq(
  dds_C2_AXON_IP_AGE,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_C2_AXON_IP_AGE <- DESeq2::results(
  dds_C2_AXON_IP_AGE,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_AXON_IP_AGE <- lfcShrink(
  dds_C2_AXON_IP_AGE,
  res = res_C2_AXON_IP_AGE,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)

# AXON AGE: Use depleted genes to add to AXON_FRACTION_ENRICHED_GENES list ----
AXON_FRACTION_ENRICHED_GENES <- res_C2_AXON_IP_AGE %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  select(-c(lfcSE, pvalue)) %>%
  filter(padj < 0.01) %>%
  arrange(log2FoldChange) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
         axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
  filter(axon_translated &
           !mb_translated &
           log2FoldChange < 0) %>%
  pull(ensembl_gene_id)



# ----
# AXON AGE CONTINUED - MAKE THE METADATA TABLE ----


AXON_AGE_META <- res_C2_AXON_IP_AGE %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  # left_join((
  #   AXON_FRACTION_META %>%
  #     select("external_gene_name",
  #            expression_energy)
  # ),
  # by = c("external_gene_name")) %>%
  mutate(
    mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
    axon_enriched = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
    # aba_specific = expression_energy < 0.5,
    # aba_specific = replace_na(aba_specific, FALSE),
    # evidence = ifelse(
    #   axon_translated,
    #   ifelse(
    #     aba_specific & mb_translated,
    #     "Axon, Soma and ABA",
    #     ifelse(
    #       aba_specific,
    #       "Axon and ABA",
    #       ifelse(mb_translated,
    #              "Axon and Soma",
    #              "Axon only")
    #     )
    #   ),
    #   ifelse(aba_specific & mb_translated,
    #          "Soma and ABA",
    #          "Soma only")
    # ),
    # evidence = ifelse(!mb_translated & evidence == "Soma only",
    #                   "No support", evidence),
    # evidence = factor(
    #   evidence,
    #   levels = c(
    #     "Axon, Soma and ABA",
    #     "Axon and Soma",
    #     "Axon and ABA",
    #     "Axon only",
    #     "Soma and ABA",
    #     "Soma only",
    #     "No support"
    #   )
    # ),
    # mb_age = ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
    # mb_age_relationship = ifelse(
    #   ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id,
    #   ifelse(log2FoldChange > 0, "Common Upregulation",
    #          "Opposing"),
    #   ifelse(
    #     ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_DOWN$ensembl_gene_id,
    #     ifelse(log2FoldChange < 0, "Common Downregulation",
    #            "Opposing"),
    #     "Axon Specific"
    #   )
    # ),
    # mb_age_relationship = ifelse(padj < 0.01, mb_age_relationship, "NS in Axon"),
    outcome = ifelse(
      padj < 0.01,
      ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
      "Unchanged"
    )
  ) %>%
  left_join(AXON_AGE_TOTAL_OUTCOMES,
            by = "ensembl_gene_id",
            suffix = c("_TRAP", "_TOTAL")) %>%
  mutate(
    specific = ifelse(
      outcome_TRAP == outcome_TOTAL,
      "Shared",
      ifelse(
        outcome_TRAP == "Up",
        ifelse(outcome_TOTAL == "Unchanged",
               "Specific",
               "Specific Opposing"),
        ifelse(outcome_TOTAL == "Unchanged",
               "Specific",
               "Specific Opposing")
      )
    ),
    specific = replace_na(specific, "Specific")
  ) %>%
  left_join(publications_all_genes_axon) %>%
  left_join(publications_all_genes) %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name,
    description,
    # baseMean,
    log2FoldChange,
    padj,
    mb_translated,
    axon_enriched,
    # mb_age_relationship,
    pub_pd,
    pub_axon,
    # expression_energy,
    # aba_specific,
    # evidence,
    specific,
    outcome_TRAP,
    outcome_TOTAL
  ) %>%
  left_join(enrichment_status)

# AXON AGE - Number of genes in enrichment groups ----

p_axon_age_enrichment_groups <- AXON_AGE_META %>%
  filter(padj < 0.01 & log2FoldChange > 0) %>%
  mutate(enrichment = ifelse(enrichment == "Soma and Axon",
                             "Soma and\nAxon",
                             ifelse(enrichment == "Soma only",
                                    "Soma only",
                                    ifelse(enrichment == "Axon only",
                                           "Axon only",
                                           ifelse(enrichment == "No enrichment",
                                                  "No enrichment",
                                                  "No enrichment"))))) %>%
  ggplot(aes(x = enrichment,
             fill = enrichment %in% c("Soma only", "Soma and\nAxon"))) +
  geom_bar(colour = "black") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("grey", pal_d3()(1)[1])) +
  geom_text(aes(y = ..count..,
                label = ..count..),
            stat = "count",
            vjust = -0.5) +
  labs(x = "Enrichment Category",
       y = "Number of Genes")

# AXON AGE TRAP SPECIFICITY (PLOT) ----

p_axon_age_trap_specificity <- AXON_AGE_META %>%
  filter(padj < 0.01 & log2FoldChange > 0 & enrichment %in% c("Soma only", "Soma and Axon")) %>%
  ggplot(aes(x = specific,
             fill = specific)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(4, 1)]) +
  labs(x = "Enrichment Category",
       y = "Number of Genes",
       fill = "TRAP\nSpecificity") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(y = ..count..,
                label = ..count..),
            stat = "count",
            vjust = -0.5)

AXON_AGE_SHARED_DE_META <- AXON_AGE_META %>%
  filter(specific == "Shared") %>%
  nrow


# AXON_AGE_META %>% filter(padj < 0.01 &
#                            log2FoldChange > 0 &
#                            enrichment %in% c("Soma only", "Soma and Axon")) %>%
#   filter(ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id)

# AXON AGE EVIDENCE (PLOT) ----

p_axon_age_evidence <- AXON_AGE_META %>%
  left_join(enrichment_status) %>%
  filter(padj < 0.1) %>%
  filter(enrichment != "No enrichment") %>%
  mutate(outcome_TRAP = ifelse(log2FoldChange > 0, "Increased\nabundance", "Decreased\nabundance")) %>%
  # filter(outcome_TRAP != "Unchanged") %>%
  # mutate(evidence = ifelse(mb_translated & axon_translated,
  #                          "Soma and Axon",
  #                          ifelse(mb_translated,
  #                                 "Soma only",
  #                                 ifelse(axon_translated,
  #                                        "Axon only",
  #                                        "No enrichment"))),
  #        evidence = factor(evidence,
  #                          levels = c("Soma only",
  #                                     "Axon only",
#                                     "Soma and Axon",
#                                     "No enrichment"
#                                     ))) %>%
group_by(outcome_TRAP, enrichment) %>%
  tally %>%
  mutate(prop = n / sum(n),
         sum = sum(n)) %>%
  ggplot(aes(x = outcome_TRAP,
             y = prop,
             fill = enrichment)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  geom_text(aes(label = paste("Total:", sum),
                y = 1),
            check_overlap = T,
            vjust = -0.5,
            fontface = 2,
            size = 4) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Change in Age",
       y = "Percentage of Genes",
       fill = "Enrichment") +
  geom_text(
    aes(
      x = outcome_TRAP,
      y = prop,
      label = ifelse(prop > 0.1,
                     scales::percent(prop,
                                     accuracy = 1),
                     "")
    ),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = 2
  ) +
  theme(
    legend.position = "top"
  ) +
  guides(fill = guide_legend(nrow = 2,
                             byrow=TRUE,
                             title.position="top",
                             title.hjust = 0.5))

# Axon AGE: Plot log fold changes and putative dopaminergic selection ----


p_axon_age_putative_dopaminergic <- AXON_AGE_META %>%
  left_join(enrichment_status) %>%
  filter(padj < 0.01) %>%
  mutate(dopaminergic = ifelse(enrichment == "Soma only" |
                                 enrichment == "Soma and Axon" &
                                 log2FoldChange < 0 |
                                 enrichment == "Axon only" &
                                 log2FoldChange < 0 |
                                 enrichment == "No enrichment" &
                                 log2FoldChange < 0, "Dopaminergic", "Other Cell Types"),
         dopaminergic = ifelse(enrichment == "Soma and Axon" &
                                 log2FoldChange > 0,
                               "Potentially\nDopaminergic", dopaminergic)) %>%
  ggplot(aes(x = enrichment,
             y = log2FoldChange,
             colour = dopaminergic)) +
  geom_quasirandom() +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_color_d3() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  labs(x = "Enrichment Category",
       y = expression(Log[2] ~ Fold ~ Change)) +
  geom_hline(yintercept = 0,
             linetype = "dotted")



# AXON AGE: Assign putative dopaminergic based on aged axon TRAP downregulation ----


putative_dopaminergic <- AXON_AGE_META %>%
  left_join(enrichment_status) %>%
  filter(padj < 0.1) %>%
  mutate(dopaminergic = ifelse(enrichment == "Soma only" |
                                 enrichment == "Soma and Axon" &
                                 log2FoldChange < 0 |
                                 enrichment == "Axon only" &
                                 log2FoldChange < 0 |
                                 enrichment == "No enrichment" &
                                 log2FoldChange < 0, TRUE, FALSE)) %>%
  filter(dopaminergic) %>%
  pull(ensembl_gene_id)



putative_nondopaminergic <- AXON_AGE_META %>%
  left_join(enrichment_status) %>%
  filter(padj < 0.1) %>%
  mutate(dopaminergic = ifelse(enrichment == "Soma and Axon" &
                                 log2FoldChange > 0 |
                                 enrichment == "Axon only" &
                                 log2FoldChange > 0 |
                                 enrichment == "No enrichment" &
                                 log2FoldChange > 0, TRUE, FALSE)) %>%
  filter(dopaminergic) %>%
  pull(ensembl_gene_id)


# Using ABA specificity is no better than flipping a coin
# AXON_AGE_META %>%
#   left_join(enrichment_status) %>%
#   filter(padj < 0.1) %>%
#   # filter(enrichment != "No enrichment") %>%
#   mutate(outcome = ifelse(log2FoldChange > 0, "Increased\nabundance", "Decreased\nabundance")) %>%
#   group_by(outcome, enrichment, aba_specific) %>%
#   tally %>%
#   mutate(aba_specific = aba_specific,
#          prop = n / sum(n),
#          sum = sum(n)) %>%
#   ggplot(aes(x = enrichment,
#              y = prop,
#              fill = aba_specific)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   facet_wrap(vars(outcome),
#              nrow = 1)



# ----
# ----
# ----
# AXON VS MB - LIBRARY COMPOSITION CHECK ----

C2_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C2_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype", "region", "compartment"), returnData = T) %>%
      select(age, genotype, region, compartment, "Slc6a3" = count)
  }, {
    plotCounts(dds_C2_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
      select("Th" = count)
    # }, {
    #   plotCounts(dds_C2_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
    #     select("Ddc" = count)
  }, {
    plotCounts(dds_C2_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C2_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C2_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C2_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-c(age, genotype, region, compartment),
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C2")


C3_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C3_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype", "region", "compartment"), returnData = T) %>%
      select(age, genotype, region, compartment, "Slc6a3" = count)
  }, {
    plotCounts(dds_C3_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
      select("Th" = count)
    # }, {
    #   plotCounts(dds_C3_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
    #     select("Ddc" = count)
  }, {
    plotCounts(dds_C3_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C3_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C3_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C3_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-c(age, genotype, region, compartment),
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C3")


p_axon_vs_mb_ip_library_composition <- bind_rows(C2_IP_MARKERS,
                                                 C3_IP_MARKERS) %>%
  mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
                         "Dopaminergic", "Other"),
         compartment = ifelse(compartment == "MB", "Soma", "Axon")) %>%
  filter(count > 1) %>%
  # group_by(age, gene) %>%
  # summarise(count = mean(count)) %>%
  ggplot(aes(x = compartment,
             y = log2(count + 1),
             colour = marker,
             group = gene)) +
  geom_quasirandom() +
  # geom_smooth(method = "lm") +
  # facet_wrap(vars(cohort), nrow = 1) +
  panel_border() +
  scale_color_d3() +
  labs(x = "Compartment",
       y = expression(Log[2] ~ Expression ~ Count),
       colour = "Cell Type Marker") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5))

# AXON VS MB - ALTERNATIVE COMPOSITION CHECK ----

LIBRARY_PROPORTIONS <- bind_rows(
  {
    as_tibble(counts(dds_C2_IP), rownames = "ensembl_gene_id") %>%
      mutate(cohort = "C2")
  }, {
    as_tibble(counts(dds_C3_IP), rownames = "ensembl_gene_id") %>%
      mutate(cohort = "C3")
  }
) %>%
  pivot_longer(-c(ensembl_gene_id, cohort),
               names_to = "sample_name",
               values_to = "count") %>%
  left_join(colData(dds), copy = T) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
         enrichment = ifelse(mb_translated & axon_translated,
                             "Soma and Axon",
                             ifelse(mb_translated,
                                    "Soma only",
                                    ifelse(axon_translated,
                                           "Axon only",
                                           "No enrichment"))),
         enrichment = factor(enrichment,
                             levels = c("Soma only",
                                        "Axon only",
                                        "Soma and Axon",
                                        "No enrichment"
                             )),
         compartment = ifelse(compartment == "MB", "Soma", "Axon")) %>%
  group_by(compartment, enrichment) %>%
  summarise(count = sum(count, na.rm = T)) %>%
  mutate(prop = count/sum(count)) %>%
  drop_na()



# LIBRARY_PROPORTIONS %>%
#   ggplot(aes(x = compartment,
#              y = prop,
#              fill = enrichment)) +
#   geom_col(colour = "black") +
#   # facet_wrap(vars(cohort), nrow = 1)  +
#   theme(legend.position = "top") +
#   geom_text(aes(y = prop, label = ifelse(prop > 0.1, scales::percent(prop, accuracy = 1), "")),
#             position = position_stack(vjust = 0.5),
#             show.legend = FALSE,
#             colour = "white",
#             fontface = 2) +
#   scale_fill_d3() +
#   # scale_fill_manual(values = pal_d3()(4)[c(4, 1, 3)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
#                      labels = scales::percent) +
#   # panel_border() +
#   labs(x = "Compartment",
#        y = "Percentage of library",
#        fill = "Enrichment")  +
#   guides(fill = guide_legend(title.position="top",
#                              title.hjust = 0.5,
#                              nrow = 2,
#                              byrow = T))
























# AXON VS MB: THIS IS COMPLETE ----------------------
# AXON VS MB: TOTAL DESeq2 ----

dds_C2_TOTAL_COMPARTMENT_FILTER <- filter_genes(dds_C2_TOTAL,
                                                grouping = c("compartment", "gene_id"))
dds_C3_TOTAL_COMPARTMENT_FILTER <- filter_genes(dds_C3_TOTAL,
                                                grouping = c("compartment", "gene_id"))

genes_TOTAL_COMPARTMENT <- intersect(dds_C2_TOTAL_COMPARTMENT_FILTER,
                                     dds_C3_TOTAL_COMPARTMENT_FILTER)

dds_C2_TOTAL_COMPARTMENT <- dds_C2_TOTAL[rownames(dds_C2_TOTAL) %in% genes_TOTAL_COMPARTMENT,]
dds_C3_TOTAL_COMPARTMENT <- dds_C3_TOTAL[rownames(dds_C3_TOTAL) %in% genes_TOTAL_COMPARTMENT,]

dds_C2_TOTAL_COMPARTMENT@design <- ~ compartment
dds_C3_TOTAL_COMPARTMENT@design <- ~ compartment

colData(dds_C2_TOTAL_COMPARTMENT) <- droplevels(colData(dds_C2_TOTAL_COMPARTMENT))
colData(dds_C3_TOTAL_COMPARTMENT) <- droplevels(colData(dds_C3_TOTAL_COMPARTMENT))

dds_C2_TOTAL_COMPARTMENT <- DESeq(
  dds_C2_TOTAL_COMPARTMENT,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

dds_C3_TOTAL_COMPARTMENT <- DESeq(
  dds_C3_TOTAL_COMPARTMENT,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_C2_TOTAL_COMPARTMENT <- DESeq2::results(
  dds_C2_TOTAL_COMPARTMENT,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_TOTAL_COMPARTMENT <- lfcShrink(
  dds_C2_TOTAL_COMPARTMENT,
  res = res_C2_TOTAL_COMPARTMENT,
  contrast = c("compartment", "AXON", "MB"),
  type = "ashr",
  parallel = T
)

res_C3_TOTAL_COMPARTMENT <- DESeq2::results(
  dds_C3_TOTAL_COMPARTMENT,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_TOTAL_COMPARTMENT <- lfcShrink(
  dds_C3_TOTAL_COMPARTMENT,
  res = res_C3_TOTAL_COMPARTMENT,
  contrast = c("compartment", "AXON", "MB"),
  type = "ashr",
  parallel = T
)

summary(res_C2_TOTAL_COMPARTMENT)
summary(res_C3_TOTAL_COMPARTMENT)

AXON_VS_MB_TOTAL_META <-
  as_tibble(res_C2_TOTAL_COMPARTMENT, rownames = "ensembl_gene_id") %>%
  left_join(
    as_tibble(res_C3_TOTAL_COMPARTMENT, rownames = "ensembl_gene_id"),
    by = "ensembl_gene_id",
    suffix = c("_C2", "_C3")
  ) %>%
  left_join(anno) %>%
  mutate(
    across(starts_with("pvalue"), ~ replace_na(.x, 1 - 1e-10)),
    # replace NA with *almost* 1 (it can't be perfectly 1)
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)),
    # replace 0 values with minimum non-zero of the column
    across(starts_with("pvalue"), ~ .x / 2),
    # divide 2-sided pvalue by 2 into 1-sided
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)) # again replace 0 values with minimum non-zero of the column
  ) %>%
  select(sort(names(.))) %>%
  select(
    ensembl_gene_id,
    external_gene_name,
    description,
    everything(),
    -starts_with("padj")
  ) %>%
  filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
  drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
  rowwise() %>%
  mutate(across(starts_with("pvalue"),
                ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
  mutate(
    sumz = sumz(c_across(pvalue_C2:pvalue_C3))$p,
    log2FoldChange = log2FoldChange_C2
  ) %>%
  ungroup() %>%
  mutate(
    sumz_adj = p.adjust(sumz, method = "fdr"),
    conflict = sign(log2FoldChange_C2) != sign(log2FoldChange_C3)) # state whether there is a conflict in l2fc direction

# AXON VS MB: TOTAL outcomes plot ----

AXON_VS_MB_TOTAL_OUTCOMES <- AXON_VS_MB_TOTAL_META %>%
  mutate(outcome = ifelse(sumz_adj > 0.01,
                          "Unchanged",
                          ifelse(log2FoldChange > 0,
                                 "Axon", "Soma"))) %>%
  select(external_gene_name,
         outcome_total = outcome)

p_axon_vs_mb_total_outcome <- AXON_VS_MB_TOTAL_OUTCOMES %>%
  filter(outcome_total != "Unchanged") %>%
  ggplot(aes(x = outcome_total,
             fill = outcome_total)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2, 4)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "TOTAL Region",
       y = "Number of Genes") +
  geom_text(aes(label = paste("Total:", ..count..),
                y = ..count..),
            stat = "count",
            check_overlap = T,
            vjust = -0.5,
            fontface = 2,
            size = 4) +
  theme(legend.position = "none")

# AXON vs MB DESeq2 -------------------------------------------------------

dds_C2_IP@design <- ~ collection + age + compartment

colData(dds_C2_IP) <- droplevels(colData(dds_C2_IP))

dds_C2_IP_COMPARTMENT <- dds_C2_IP %>% filter_zeros()

dds_C2_IP_COMPARTMENT <- DESeq(
  dds_C2_IP_COMPARTMENT,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_C2_IP_COMPARTMENT <- DESeq2::results(
  dds_C2_IP_COMPARTMENT,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_IP_COMPARTMENT <- lfcShrink(
  dds_C2_IP_COMPARTMENT,
  res = res_C2_IP_COMPARTMENT,
  contrast = c("compartment", "AXON", "MB"),
  type = "ashr",
  parallel = T
)

# summary(res_C2_IP_COMPARTMENT)

AXON_VS_MB_META <-
  res_C2_IP_COMPARTMENT %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  select(-c(lfcSE, pvalue)) %>%
  arrange(log2FoldChange) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
  # left_join((AXON_FRACTION_META %>%
  #              select("external_gene_name",
  #                     expression_energy)),
  #           by = c("external_gene_name")) %>%
  mutate(compartment = ifelse(log2FoldChange > 0, "Axon", "Soma")

  ) %>%
  left_join(enrichment_status)

# evidence = ifelse(axon_translated,
#                   ifelse(aba_specific & mb_translated,
#                          "Axon, Soma and ABA",
#                          ifelse(aba_specific,
#                                 "Axon and ABA",
#                                 ifelse(mb_translated,
#                                        "Axon and Soma",
#                                        "Axon only"))),
#                   ifelse(aba_specific & mb_translated,
#                          "Soma and ABA",
#                          "Soma only")),
# evidence = ifelse(!mb_translated & evidence == "Soma only",
#                   "No support", evidence),
# evidence = factor(evidence, levels = c("Axon, Soma and ABA",
#                                        "Axon and Soma",
#                                        "Axon and ABA",
#                                        "Axon only",
#                                        "Soma and ABA",
#                                        "Soma only",
#                                        "No support")))


# AXON vs MB: Confidence stratified (PLOT) ----

plot_data <- AXON_VS_MB_META %>%
  filter(padj < 0.01)
# filter(mb_translated | axon_translated) # Commented because I am plotting all the evidence; even no support genes

p_axon_vs_mb_enrichment <- plot_data %>%
  left_join(enrichment_status) %>%
  group_by(compartment, enrichment) %>%
  tally %>%
  mutate(prop = n/sum(n),
         sum = sum(n)) %>%
  ggplot(aes(x = compartment,
             y = prop,
             fill = enrichment)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "TRAP Region",
       y = "Percentage of Genes",
       fill = "Axonal Enrichment") +
  geom_text(aes(x = compartment,
                y = 1,
                label = paste("Total:", sum)),
            check_overlap = T,
            vjust = -0.5) +
  geom_text(aes(x = compartment,
                y = prop,
                label = ifelse(prop > 0.05,
                               scales::percent(prop,
                                               accuracy = 0.1),
                               "")),
            position = position_stack(vjust = 0.5),
            colour = "white",
            fontface = 2) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 2,
                             byrow = T,
                             title.position = "top",
                             title.hjust = 0.5))



# AXON VS MB: Number of DEGS ----

AXON_VS_MB_META_SIGNIF_SUMMARY <- AXON_VS_MB_META %>%
  filter(padj < 0.01) %>%
  filter(mb_translated)

n_AXON_VS_MB_SIGNIF_GENES <- AXON_VS_MB_META %>%
  filter(padj < 0.01) %>%
  nrow

# require that genes be soma enriched.
n_AXON_VS_MB_SIGNIF_ENRICHED_GENES <- AXON_VS_MB_META %>%
  filter(padj < 0.01) %>%
  filter(mb_translated) %>%
  nrow


# AXON VS MB: AXON NUMBERS ----

n_AXON_VS_MB_AXON_ENRICHED <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
  filter(log2FoldChange > 0) %>%
  pull(external_gene_name) %>%
  length

n_AXON_VS_MB_AXON_ENRICHED_LOW_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
  filter(enrichment == "Axon only" &
           log2FoldChange > 0) %>%
  pull(external_gene_name) %>%
  length

n_AXON_VS_MB_AXON_ENRICHED_MEDIUM_EVIDENCE_SOMA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
  filter(enrichment == "Soma and Axon" &
           log2FoldChange > 0) %>%
  pull(external_gene_name) %>%
  length

# n_AXON_VS_MB_AXON_ENRICHED_MEDIUM_EVIDENCE_ABA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(enrichment == "Axon and ABA" &
#            log2FoldChange > 0) %>%
#   pull(external_gene_name) %>%
#   length

n_AXON_VS_MB_AXON_ENRICHED_SOMA_ONLY <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
  filter(enrichment == "Soma only" &
           log2FoldChange > 0) %>%
  pull(external_gene_name) %>%
  length

# n_AXON_VS_MB_AXON_ENRICHED_HIGH_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(enrichment == "Axon, Soma and ABA" &
#            log2FoldChange > 0) %>%
#   pull(external_gene_name) %>%
#   length

# AXON VS MB: SOMA NUMBERS ----
n_AXON_VS_MB_SOMA_ENRICHED <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
  filter(log2FoldChange < 0) %>%
  pull(external_gene_name) %>%
  length

# n_AXON_VS_MB_SOMA_ENRICHED_LOW_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon only" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length

# n_AXON_VS_MB_SOMA_ENRICHED_MEDIUM_EVIDENCE_SOMA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon and Soma" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length

# n_AXON_VS_MB_SOMA_ENRICHED_MEDIUM_EVIDENCE_ABA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon and ABA" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length

# n_AXON_VS_MB_SOMA_ENRICHED_HIGH_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon, Soma and ABA" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length


# AXON VS MB: Prioritising ----

axon_vs_mb_trap_total <- library_proportion_info_total %>%
  filter(enrichment != "No enrichment") %>%
  mutate(diff = log2(axon / mb)) %>%
  arrange(desc(diff)) %>%
  mutate(
    compartment = ifelse(diff > 1, "axon", "mb"),
    dopaminergic = ensembl_gene_id %in% putative_dopaminergic,
    non_dopaminergic = ensembl_gene_id %in% putative_nondopaminergic
  ) %>%
  left_join(anno) %>%
  left_join(publications_all_genes) %>%
  left_join(publications_all_genes_axon) %>%
  # filter(mb > 20 | axon > 20) %>%
  group_by(enrichment) %>%
  arrange(desc(diff)) %>%
  mutate(fraction = "TRAP") %>%
  ungroup() %>%
  select(ensembl_gene_id,
         diff_total = diff)


AXON_VS_MB_META_SIGNIF_PRIORITY <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
  left_join(enrichment_status) %>%
  filter(enrichment %in% c(
    "Soma and Axon",
    "Soma only") &
      log2FoldChange > 0) %>%
  left_join(publications_all_genes) %>%
  left_join(publications_all_genes_axon) %>%
  # filter(expression_energy < 1) %>%
  select(external_gene_name,
         description,
         log2FoldChange,
         padj,
         # expression_energy,
         enrichment,
         pub_pd,
         pub_axon,
         ensembl_gene_id) %>%
  distinct() %>%
  mutate(outcome_trap = ifelse(padj > 0.01,
                               "Unchanged",
                               ifelse(log2FoldChange > 0,
                                      "Axon", "Soma"))) %>%
  left_join(axon_vs_mb_trap_total) %>%
  filter(diff_total < 0.5) %>%
  left_join(AXON_VS_MB_TOTAL_OUTCOMES) %>%
  mutate(relationship = ifelse(outcome_trap == outcome_total,
                               "Shared",
                               ifelse(outcome_trap == "Axon",
                                      ifelse(outcome_total == "Unchanged",
                                             "Axon specific",
                                             "Axon opposing"),
                                      ifelse(outcome_total == "Unchanged",
                                             "Soma specific",
                                             "Soma opposing")))) %>%
  distinct


n_AXON_VS_MB_AXON_ENRICHED_PRIORITY <- nrow(AXON_VS_MB_META_SIGNIF_PRIORITY)

n_AXON_VS_MB_AXON_ENRICHED_HIGH_EVIDENCE_TOTAL_OPPOSING <- AXON_VS_MB_META_SIGNIF_PRIORITY %>%
  filter(relationship == "Axon opposing") %>%
  nrow

# AXON VS MB: TOTAL outcome comparison with HIGHEST CONFIDENCE TRAP GENES ----

p_axon_vs_mb_total_comparison <- AXON_VS_MB_META_SIGNIF_PRIORITY %>%
  filter(!is.na(relationship)) %>%
  mutate(relationship = ifelse(relationship == "Axon opposing",
                               "Axon\nopposing", relationship)) %>%
  ggplot(aes(x = relationship,
             fill = relationship)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(3, 1, 4)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(label = ..count..,
                y = ..count..),
            vjust = -0.5,
            stat = "count",
            fontface = 2) +
  theme(legend.position = "none") +
  labs(x = "Relationship to TOTAL Samples",
       y = "Number of Genes")

n_AXON_VS_MB_TOTAL_BREAKDOWN <- AXON_VS_MB_META_SIGNIF_PRIORITY %>%
  filter(!is.na(relationship)) %>%
  group_by(relationship) %>%
  tally() %>%
  mutate(prop = n / sum(n)) %>%
  pull(n)

# AXON VS MB: Example Atg7 ----
# p_axon_vs_mb_atg7 <- plotCounts(dds,
#                                 get_ensembl_gene_id("Atg7"),
#                                 c("cohort", "fraction", "compartment", "region", "age", "genotype"),
#                                 returnData = T) %>%
#   mutate(region = dplyr::recode(region,
#                                 MB = "Soma",
#                                 DS = "Dorsal\naxon",
#                                 VS = "Ventral\naxon")) %>%
#   group_by(compartment, fraction) %>%
#   filter(!outlier_iqr(count)) %>%
#   filter(count > 10) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   ggplot(aes(x = region,
#              y = log2(count + 1),
#              fill = region)) +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   facet_wrap(vars(fraction),
#              # scales = "free_y",
#              nrow = 1) +
#   labs(x = "TRAP Fraction and Compartment",
#        y = expression(Log[2] ~ Expression ~ Count),
#        title = "Atg7") +
#   panel_border() +
#   # scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none")

data <- plotCounts(dds,
                   get_ensembl_gene_id("Atg7"),
                   c("cohort", "fraction", "compartment", "region", "age", "genotype"),
                   returnData = T) %>%
  mutate(compartment = dplyr::recode(compartment,
                                     MB = "Soma",
                                     AXON = "Axon")) %>%
  group_by(compartment, fraction) %>%
  filter(!outlier_iqr(count)) %>%
  filter(count > 10) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"), 
         count = log2(count + 1))

p_axon_vs_mb_atg7 <- ggboxplot(data, 
                               x = "compartment", 
                               y = "count", 
                               fill = "compartment", 
                               facet.by = "fraction", 
                               outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  labs(x = "TRAP Fraction and Compartment", 
       y = expression(Log[2] ~ Expression ~ Count),
       title = "Atg7", 
       caption = "Enriched in Somal and Axonal TRAP Samples") +
  theme(legend.position = "none") +
  stat_pvalue_manual({
    bind_rows(
      {
        AXON_VS_MB_META_SIGNIF_PRIORITY %>%
          filter(external_gene_name == "Atg7") %>%
          mutate(fraction = "TRAP") %>%
          select(fraction, 
                 padj)
      },{
        AXON_VS_MB_TOTAL_META %>%
          filter(external_gene_name == "Atg7") %>%
          mutate(fraction = "TOTAL") %>%
          select(fraction, 
                 "padj" = sumz_adj)
      }
    ) %>%
      mutate(group1 = c("Soma", "Soma"), 
             group2 = c("Axon", "Axon"), 
             padj = signif(padj, 3))
  } %>%
    pval_asterisks(padj), 
  y.position = c(18, 11.5), 
  label = "ast",
  size = 3, 
  vjust = -0.5) +
  scale_y_continuous(limits = c(4, 20), 
                     n.breaks = 8)

# AXON VS MB: Example Snx1 ----

data <- plotCounts(dds,
                   get_ensembl_gene_id("Snx1"),
                   c("cohort", "fraction", "compartment", "region", "age", "genotype"),
                   returnData = T) %>%
  mutate(compartment = dplyr::recode(compartment,
                                MB = "Soma",
                                AXON = "Axon")) %>%
  group_by(compartment, fraction) %>%
  filter(!outlier_iqr(count)) %>%
  filter(count > 10) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"), 
         count = log2(count + 1))

p_axon_vs_mb_snx1 <- ggboxplot(data, 
                               x = "compartment", 
                               y = "count", 
                               fill = "compartment", 
                               facet.by = "fraction", 
                               outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  labs(x = "TRAP Fraction and Compartment", 
       y = expression(Log[2] ~ Expression ~ Count),
       title = "Snx1", 
       caption = "Enriched in Somal and Axonal TRAP Samples") +
  theme(legend.position = "none") +
  stat_pvalue_manual({
    bind_rows(
    {
      AXON_VS_MB_META_SIGNIF_PRIORITY %>%
        filter(external_gene_name == "Snx1") %>%
        mutate(fraction = "TRAP") %>%
        select(fraction, 
               padj)
    },{
      AXON_VS_MB_TOTAL_META %>%
        filter(external_gene_name == "Snx1") %>%
        mutate(fraction = "TOTAL") %>%
        select(fraction, 
               "padj" = sumz_adj)
    }
  ) %>%
    mutate(group1 = c("Soma", "Soma"), 
           group2 = c("Axon", "Axon"), 
           padj = signif(padj, 3))
    } %>%
    pval_asterisks(padj), 
  y.position = c(13.5, 13), 
  label = "ast",
  size = 3) +
  scale_y_continuous(limits = c(6, 15))

# p_axon_vs_mb_snx1 <- plotCounts(dds,
#                                 get_ensembl_gene_id("Snx1"),
#                                 c("cohort", "fraction", "compartment", "region", "age", "genotype"),
#                                 returnData = T) %>%
#   mutate(region = dplyr::recode(region,
#                                 MB = "Soma",
#                                 DS = "Dorsal\naxon",
#                                 VS = "Ventral\naxon")) %>%
#   group_by(compartment, fraction) %>%
#   filter(!outlier_iqr(count)) %>%
#   filter(count > 10) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   ggplot(aes(x = region,
#              y = log2(count + 1),
#              fill = region)) +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   facet_wrap(vars(fraction),
#              # scales = "free_y",
#              nrow = 1) +
#   labs(x = "TRAP Fraction and Compartment",
#        y = expression(Log[2] ~ Expression ~ Count),
#        title = "Snx1") +
#   panel_border() +
#   # scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none")

# AXON VS MB: Example Cul3 ----

data <- plotCounts(dds,
                   get_ensembl_gene_id("Cul3"),
                   c("cohort", "fraction", "compartment", "region", "age", "genotype"),
                   returnData = T) %>%
  mutate(compartment = dplyr::recode(compartment,
                                     MB = "Soma",
                                     AXON = "Axon")) %>%
  group_by(compartment, fraction) %>%
  filter(!outlier_iqr(count)) %>%
  filter(count > 10) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"), 
         count = log2(count + 1))

p_axon_vs_mb_cul3 <- ggboxplot(data, 
                               x = "compartment", 
                               y = "count", 
                               fill = "compartment", 
                               facet.by = "fraction", 
                               outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  labs(x = "TRAP Fraction and Compartment", 
       y = expression(Log[2] ~ Expression ~ Count),
       title = "Cul3", 
       caption = "Enriched in Somal and Axonal TRAP Samples") +
  theme(legend.position = "none") +
  stat_pvalue_manual({
    bind_rows(
      {
        AXON_VS_MB_META_SIGNIF_PRIORITY %>%
          filter(external_gene_name == "Cul3") %>%
          mutate(fraction = "TRAP") %>%
          select(fraction, 
                 padj)
      },{
        AXON_VS_MB_TOTAL_META %>%
          filter(external_gene_name == "Cul3") %>%
          mutate(fraction = "TOTAL") %>%
          select(fraction, 
                 "padj" = sumz_adj)
      }
    ) %>%
      mutate(group1 = c("Soma", "Soma"), 
             group2 = c("Axon", "Axon"), 
             padj = signif(padj, 3)) %>%
      pval_asterisks(padj)
  }, 
  y.position = c(13, 12), 
  label = "ast",
  size = 3, 
  vjust = -0.5) +
  scale_y_continuous(limits = c(6, 15))

# p_axon_vs_mb_cul3 <- plotCounts(dds,
#                                 get_ensembl_gene_id("Cul3"),
#                                 c("cohort", "fraction", "compartment", "region", "age", "genotype"),
#                                 returnData = T,
#                                 normalized = F) %>%
#   mutate(region = dplyr::recode(region,
#                                 MB = "Soma",
#                                 DS = "Dorsal\naxon",
#                                 VS = "Ventral\naxon")) %>%
#   group_by(compartment, fraction) %>%
#   filter(!outlier_iqr(count)) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   # filter(cohort != "C3") %>%
#   ggplot(aes(x = region,
#              y = log2(count + 1),
#              fill = region)) +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   facet_wrap(vars(fraction),
#              # scales = "free_y",
#              nrow = 1) +
#   labs(x = "TRAP Fraction and Compartment",
#        y = expression(Log[2] ~ Expression ~ Count),
#        title = "Cul3") +
#   panel_border() +
#   scale_y_continuous(limits = c(5, NA)) +
#   theme(legend.position = "none")

# AXON VS MB: Enrichment distribution ----

p_axon_vs_mb_distribution <- AXON_VS_MB_META %>%
  filter(mb_translated | axon_translated) %>%
  select(ensembl_gene_id,
         external_gene_name,
         log2FoldChange,
         padj,
         enrichment) %>%
  mutate(compartment = ifelse(
    log2FoldChange > 0 & padj < 0.01,
    "Axonal",
    ifelse(log2FoldChange < 0 &
             padj < 0.01, "Somal",
           "Non-significant")
  )) %>%
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    colour = compartment
  )) +
  geom_point(size = 0.2) +
  scale_x_continuous(limits = c(-5, 5),
                     oob = scales::squish) +
  scale_y_continuous(limits = c(0, 100),
                     oob = scales::squish) +
  geom_hline(yintercept = -log10(0.01),
             linetype = "dotted") +
  scale_color_d3() +
  facet_wrap(vars(enrichment),
             nrow = 1) +
  labs(
    x = expression(Log[2] ~ Expression ~ Count),
    y = expression(-Log[10] ~ Adjusted ~ italic(P)),
    colour = "Compartment of Greater Abundance", 
    caption = "Facet: TRAP Enrichment Category"
  ) +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  panel_border()

# p_axon_vs_mb_distribution <- bind_rows({
#   library_proportion_info %>%
#     filter(enrichment != "No enrichment") %>%
#     mutate(diff = log2(axon/mb)) %>%
#     arrange(desc(diff)) %>%
#     mutate(compartment = ifelse(diff > 1, "axon", "mb"),
#            dopaminergic = ensembl_gene_id %in% putative_dopaminergic) %>%
#     left_join(anno) %>%
#     # filter(mb > 20 | axon > 20) %>%
#     group_by(enrichment) %>%
#     arrange(desc(diff)) %>%
#     mutate(fraction = "TRAP")
# }, {
#   library_proportion_info_total %>%
#     filter(enrichment != "No enrichment") %>%
#     mutate(diff = log2(axon/mb)) %>%
#     arrange(desc(diff)) %>%
#     mutate(compartment = ifelse(diff > 1, "axon", "mb"),
#            dopaminergic = ensembl_gene_id %in% putative_dopaminergic) %>%
#     left_join(anno) %>%
#     # filter(mb > 20 | axon > 20) %>%
#     group_by(enrichment) %>%
#     arrange(desc(diff)) %>%
#     mutate(fraction = "TOTAL")
# }) %>%
#   group_by(fraction, enrichment) %>%
#   mutate(rank = row_number()/dplyr::n()) %>%
#   ggplot(aes(x = rank,
#              y = diff,
#              colour = enrichment,
#              group = enrichment)) +
#   geom_line() +
#   geom_hline(yintercept = 0.5, linetype = "dotted") +
#   geom_vline(xintercept = 0.5, linetype = "dotted") +
#   facet_wrap(vars(fraction), nrow = 1) +
#   labs(x = "Number of Genes",
#        y = "Compartmental Enrichment",
#        colour = "TRAP\nEnrichment") +
#   scale_x_continuous(breaks = c(0, 1),
#                      labels = c("1", "8000")) +
#   scale_y_continuous(breaks = c(-10, 10),
#                      labels = c("Soma", "Axon")) +
#   theme(legend.position = "top") +
#   guides(color = guide_legend(nrow = 2,
#                               byrow = T,
#                               title.position = "top",
#                               title.hjust = 0.5))

# AXON VS MB: FGSEA ----
# # get entrez IDs
# anno_fgsea <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "entrezgene_id"),
#   filters = "ensembl_gene_id",
#   values = AXON_VS_MB_META$ensembl_gene_id,
#   mart = ensembl,
#   uniqueRows = TRUE
# )
#
#
# # pathways <- fgsea::gmtPathways("R/objects/msigdb.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c2.cp.kegg.v7.4.symbols.gmt.txt")
# pathways <- fgsea::gmtPathways("R/objects/c5.go.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.bp.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.mf.v7.4.symbols.gmt.txt")
#
# ranks <- AXON_VS_MB_META %>%
#   filter(mb_translated) %>%
#   mutate(rank = -log10(padj) * log2FoldChange) %>%
#   inner_join(anno_human) %>%
#   select(hsapiens_homolog_associated_gene_name,
#          rank) %>%
#   filter(hsapiens_homolog_associated_gene_name != "") %>%
#   distinct(across(-rank), .keep_all = T) %>%
#   deframe
# names(ranks)
# AXON_VS_MB_FGSEA <- fgsea(pathways,
#                           ranks,
#                           eps = 0,
#                           minSize = 10,
#                           maxSize = 500,
#                           nPermSimple = 10000,
#                           BPPARAM = MulticoreParam(),
#                           nproc = 22)
#
# AXON_VS_MB_FGSEA %>% View
#
# AXON_VS_MB_FGSEA$leadingEdge[[3272]]

# ----
# ----
# ----
# AXON DS VS VS Subset common genes ----
dds_C2_AXON_IP_REGION_FILTER <- filter_genes(dds_C2_AXON_IP,
                                             grouping = c("region", "gene_id"))

# AXON DS VS VS TOTAL DESEQ2 ----

dds_C2_AXON_TOTAL_REGION <- dds_C2_AXON_TOTAL[rownames(dds_C2_AXON_TOTAL) %in% dds_C2_AXON_IP_REGION_FILTER,]

dds_C2_AXON_TOTAL_REGION@design <- ~ region

colData(dds_C2_AXON_TOTAL_REGION) <- droplevels(colData(dds_C2_AXON_TOTAL_REGION))

colData(dds_C2_AXON_TOTAL_REGION)$region <- relevel(colData(dds_C2_AXON_TOTAL_REGION)$region, ref = "VS")

dds_C2_AXON_TOTAL_REGION <- DESeq(
  dds_C2_AXON_TOTAL_REGION,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_C2_AXON_TOTAL_REGION <- DESeq2::results(
  dds_C2_AXON_TOTAL_REGION,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_AXON_TOTAL_REGION <- lfcShrink(
  dds_C2_AXON_TOTAL_REGION,
  res = res_C2_AXON_TOTAL_REGION,
  contrast = c("region", "DS", "VS"),
  type = "ashr",
  parallel = T
)

summary(res_C2_AXON_TOTAL_REGION)

DS_VS_VS_TOTAL_META <- res_C2_AXON_TOTAL_REGION %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  mutate(outcome = ifelse(padj < 0.1,
                          ifelse(log2FoldChange > 0, "Dorsal", "Ventral"),
                          "Unchanged"))

DS_VS_VS_TOTAL_OUTCOMES <- DS_VS_VS_TOTAL_META %>%
  select(ensembl_gene_id,
         outcome)

# AXON DS VS VS TRAP DESEQ2 ----

dds_C2_AXON_IP_REGION <- dds_C2_AXON_IP[dds_C2_AXON_IP_REGION_FILTER,]

dds_C2_AXON_IP_REGION@design <- ~ collection + region

colData(dds_C2_AXON_IP_REGION) <- droplevels(colData(dds_C2_AXON_IP_REGION))

colData(dds_C2_AXON_IP_REGION)$region <- relevel(colData(dds_C2_AXON_IP_REGION)$region, ref = "VS")

dds_C2_AXON_IP_REGION <- DESeq(
  dds_C2_AXON_IP_REGION,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_C2_AXON_IP_REGION <- DESeq2::results(
  dds_C2_AXON_IP_REGION,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)

res_C2_AXON_IP_REGION <- lfcShrink(
  dds_C2_AXON_IP_REGION,
  res = res_C2_AXON_IP_REGION,
  contrast = c("region", "DS", "VS"),
  type = "ashr",
  parallel = T
)

DS_VS_VS_META <- res_C2_AXON_IP_REGION %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
  left_join(enrichment_status) %>%
  mutate(outcome = ifelse(padj < 0.01,
                          ifelse(log2FoldChange > 0, "Dorsal", "Ventral"),
                          "Unchanged")) %>%
  left_join(DS_VS_VS_TOTAL_OUTCOMES,
            by = "ensembl_gene_id",
            suffix = c("_TRAP", "_TOTAL")) %>%
  mutate(specific = ifelse(outcome_TRAP == outcome_TOTAL,
                           "Shared",
                           ifelse(outcome_TRAP == "Dorsal",
                                  ifelse(outcome_TOTAL == "Unchanged",
                                         "Specific",
                                         "Specific Opposing"),
                                  ifelse(outcome_TOTAL == "Unchanged",
                                         "Specific",
                                         "Specific Opposing"))),
         specific = replace_na(specific, "Specific")) %>%
  left_join(publications_all_genes_axon) %>%
  left_join(publications_all_genes) %>%
  dplyr::select(ensembl_gene_id,
                external_gene_name,
                description,
                baseMean,
                log2FoldChange,
                padj,
                mb_translated,
                axon_translated,
                # mb_age_relationship,
                pub_pd,
                pub_axon,
                # expression_energy,
                # aba_specific,
                # evidence,
                specific,
                outcome_TRAP,
                outcome_TOTAL)

# left_join((AXON_FRACTION_META %>%
#              select("external_gene_name",
#                     expression_energy)),
#           by = c("external_gene_name")) %>%
# mutate(aba_specific = expression_energy < 0.5,
#        aba_specific = replace_na(aba_specific, FALSE),
#        # compartment = ifelse(log2FoldChange > 0, "DS", "VS"),
#        evidence = ifelse(axon_translated,
#                          ifelse(aba_specific & mb_translated,
#                                 "Axon, Soma and ABA",
#                                 ifelse(aba_specific,
#                                        "Axon and ABA",
#                                        ifelse(mb_translated,
#                                               "Axon and Soma",
#                                               "Axon only"))),
#                          ifelse(aba_specific & mb_translated,
#                                 "Soma and ABA",
#                                 "Soma only")),
#        evidence = ifelse(!mb_translated & evidence == "Soma only",
#                          "No support", evidence),
#        evidence = factor(evidence, levels = c("Axon, Soma and ABA",
#                                               "Axon and Soma",
#                                               "Axon and ABA",
#                                               "Axon only",
#                                               "Soma and ABA",
#                                               "Soma only",
#                                               "No support")),



# AXON DS VS VS: EVIDENCE STRATIFIED (2 PLOTS) ----

# FILTER FOR MB TRANSLATED GENES

DS_VS_VS_META_SIGNIF_SUMMARY <- DS_VS_VS_META %>%
  filter(padj < 0.01) %>%
  mutate(region = outcome_TRAP)

n_DS_VS_VS_DE_NOFILTER <- nrow(DS_VS_VS_META_SIGNIF_SUMMARY)

# JUST THE NUMBER OF GENES FIRST
p_ds_vs_vs_deg <- DS_VS_VS_META_SIGNIF_SUMMARY %>%
  group_by(region) %>%
  mutate(n = dplyr::n()) %>%
  ggplot(aes(x = region,
             fill = region
  )) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2, 4)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "TRAP Region",
       y = "Number of Genes") +
  geom_text(aes(label = paste("Total:", n),
                y = n),
            check_overlap = T,
            vjust = -0.5,
            fontface = 2,
            size = 4) +
  theme(legend.position = "none") +
  guides(fill = guide_legend(nrow = 2,
                             byrow=TRUE))

# THEN EVIDENCE BREAKDOWN
p_ds_vs_vs_deg_evidence <- DS_VS_VS_META_SIGNIF_SUMMARY %>%
  left_join(enrichment_status) %>%
  select(everything(),
         evidence = enrichment) %>%
  group_by(region, evidence) %>%
  tally %>%
  mutate(prop = n/sum(n),
         sum = sum(n)) %>%
  ggplot(aes(x = region,
             y = prop,
             fill = evidence)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "TRAP Region",
       y = "Percentage of Genes",
       fill = "Axonal Evidence") +
  geom_text(aes(x = region,
                y = 1.05,
                label = paste("Total:", sum)),
            check_overlap = T) +
  geom_text(aes(x = region,
                y = prop,
                label = ifelse(prop > 0.05,
                               scales::percent(prop,
                                               accuracy = 0.1),
                               "")),
            position = position_stack(vjust = 0.5),
            colour = "white",
            fontface = 2)

# AXON DS VS VS: Prioritising ----

# Selecting genes that are enriched in soma, and up in dorsal striatum.
# There are still interesting genes that are soma enriched, not axon enriched

DS_VS_VS_META_SIGNIF_PRIORITY <- DS_VS_VS_META_SIGNIF_SUMMARY %>%
  filter(mb_translated) %>%
  left_join(publications_all_genes_axon) %>%
  left_join(publications_all_genes) %>%
  # filter(expression_energy < 1) %>%
  select(ensembl_gene_id,
         external_gene_name,
         description,
         log2FoldChange,
         padj,
         axon_translated,
         # expression_energy,
         # evidence,
         pub_pd,
         pub_axon) %>%
  distinct() %>%
  filter(log2FoldChange > 0)

n_DS_VS_VS_DE_RETAINED <- nrow(DS_VS_VS_META_SIGNIF_PRIORITY)

# AXON DS VS VS: TOTAL outcomes plot ----

DS_VS_VS_TOTAL_OUTCOMES <- DS_VS_VS_TOTAL_META %>%
  mutate(outcome = ifelse(pvalue > 0.1,
                          "Unchanged",
                          ifelse(log2FoldChange > 0,
                                 "Dorsal", "Ventral"))) %>%
  select(external_gene_name,
         outcome_total = outcome,
         pvalue_TOTAL = pvalue)

p_ds_vs_vs_total_outcome <- DS_VS_VS_TOTAL_OUTCOMES %>%
  filter(outcome_total != "Unchanged") %>%
  ggplot(aes(x = outcome_total,
             fill = outcome_total)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2, 4)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "TOTAL Region",
       y = "Number of Genes") +
  geom_text(aes(label = paste("Total:", ..count..),
                y = ..count..),
            stat = "count",
            check_overlap = T,
            vjust = -0.5,
            fontface = 2,
            size = 4) +
  theme(legend.position = "none")

# AXON DS VS VS: TOTAL outcome comparison with HIGHEST CONFIDENCE TRAP GENES ----

plot_data <- DS_VS_VS_META_SIGNIF_PRIORITY %>%
  mutate(outcome_trap = ifelse(padj > 0.01,
                               "Unchanged",
                               ifelse(log2FoldChange > 0,
                                      "Dorsal", "Ventral"))) %>%
  left_join(DS_VS_VS_TOTAL_OUTCOMES) %>%
  mutate(relationship = ifelse(outcome_trap == outcome_total,
                               "Shared",
                               ifelse(outcome_trap == "Dorsal",
                                      ifelse(outcome_total == "Unchanged",
                                             "Dorsal specific",
                                             "Dorsal opposing"),
                                      ifelse(outcome_total == "Unchanged",
                                             "Ventral specific",
                                             "Ventral opposing"))))

DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON <- plot_data

p_ds_vs_vs_total_comparison <- plot_data %>%
  filter(!is.na(relationship)) %>%
  ggplot(aes(x = relationship,
             fill = relationship)) +
  geom_bar(colour = "black") +
  # scale_fill_manual(values = pal_d3()(4)[c(3, 1, 4)]) +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(label = ..count..,
                y = ..count..),
            vjust = -0.5,
            stat = "count",
            fontface = 2) +
  theme(legend.position = "none") +
  labs(x = "Relationship to TOTAL Samples",
       y = "Number of Genes")

n_DS_VS_VS_DE_PRIORITY_SPECIFIC <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
  filter(relationship == "Dorsal specific") %>%
  nrow

# n_DS_VS_VS_DE_PRIORITY_SPECIFIC_ABACONF <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific" &
#            expression_energy < 0.5) %>%
#   nrow
#
# n_DS_VS_VS_DE_PRIORITY_SPECIFIC_AXON_ENRICHED <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific" &
#            axon_translated) %>%
#   nrow
#
# n_DS_VS_VS_DE_PRIORITY_SPECIFIC_ABACONF_AXON_ENRICHED <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific" &
#            expression_energy < 0.5 &
#            axon_translated) %>%
#   nrow


# DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific" &
#            expression_energy < 0.5) %>% View

DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
  filter(relationship == "Dorsal specific")

# No GWAS targets
# DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific") %>%
#   mutate(gwas = external_gene_name %in% GWAS_GENES_BROAD) %>%
#   filter(gwas)

# plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Chrna4"), "region")

# AXON DS VS VS: Examples ----

# Oxr1
data <- view_counts(dds,
            "Oxr1") %>%
  left_join(colData(dds), copy = T) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
  group_by(cohort, fraction, region) %>%
  filter(!outlier_iqr(count)) %>%
  filter(count != 0 & region != "MB") %>%
  mutate(count = log2(count + 1))

p_ds_vs_vs_oxr1 <- ggboxplot(data, 
          x = "region", 
          y = "count", 
          fill = "region",
          outlier.shape = NA, 
          facet.by = "fraction") +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  theme(legend.position = "none") +
  scale_fill_d3() +
  labs(x = "TRAP Fraction and Striatal Region", 
       y = expression(Log[2] ~ Expression ~ Count),
       title = "Oxr1") +
  stat_pvalue_manual({
    bind_rows(
      {
        DS_VS_VS_META %>%
          filter(external_gene_name == "Oxr1") %>%
          mutate(fraction = "TRAP") %>%
          select(fraction, 
                 padj)
      },{
        DS_VS_VS_TOTAL_META %>%
          filter(external_gene_name == "Oxr1") %>%
          mutate(fraction = "TOTAL") %>%
          select(fraction, 
                 padj)
      }
    ) %>%
      mutate(group1 = c("DS", "DS"), 
             group2 = c("VS", "VS"), 
             padj = signif(padj, 3)) %>%
      pval_asterisks(padj)
  }, 
  y.position = c(14, 13), 
  label = "ast",
  vjust = -0.25) +
  scale_y_continuous(limits = c(9, 15))

# p_ds_vs_vs_oxr1 <- view_counts(dds,
#                                "Oxr1") %>%
#   left_join(colData(dds), copy = T) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   filter(count != 0) %>%
#   # filter(fraction == "IP") %>%
#   ggplot(aes(x = region,
#              y = log2(count+1),
#              fill = region)) +
#   geom_violin() +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Fraction and Region",
#        y = "Log2 Expression Count",
#        title = "Oxr1") +
#   panel_border() +
#   facet_wrap(vars(fraction),
#              nrow = 1) +
#   # scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none")


# p_ds_vs_vs_ptgs2 <- view_counts(dds_C2,
#                                 "Ptgs2") %>%
#   left_join(colData(dds_C2), copy = T) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   filter(count != 0) %>%
#   ggplot(aes(x = region,
#              y = log2(count+1),
#              fill = region)) +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Fraction and Region",
#        y = "Log2 Expression Count",
#        title = "Ptgs2 (COX-2)") +
#   panel_border() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none")


# Micu1
data <- view_counts(dds,
                    "Micu1") %>%
  left_join(colData(dds), copy = T) %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
  group_by(cohort, fraction, region) %>%
  filter(!outlier_iqr(count)) %>%
  filter(count != 0 & region != "MB") %>%
  mutate(count = log2(count + 1))

p_ds_vs_vs_micu1 <- ggboxplot(data, 
                             x = "region", 
                             y = "count", 
                             fill = "region",
                             outlier.shape = NA, 
                             facet.by = "fraction") +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  theme(legend.position = "none") +
  scale_fill_d3() +
  labs(x = "TRAP Fraction and Striatal Region", 
       y = expression(Log[2] ~ Expression ~ Count),
       title = "Micu1") +
  stat_pvalue_manual({
    bind_rows(
      {
        DS_VS_VS_META %>%
          filter(external_gene_name == "Micu1") %>%
          mutate(fraction = "TRAP") %>%
          select(fraction, 
                 padj)
      },{
        DS_VS_VS_TOTAL_META %>%
          filter(external_gene_name == "Micu1") %>%
          mutate(fraction = "TOTAL") %>%
          select(fraction, 
                 padj)
      }
    ) %>%
      mutate(group1 = c("DS", "DS"), 
             group2 = c("VS", "VS"), 
             padj = signif(padj, 3)) %>%
      pval_asterisks(padj)
  }, 
  y.position = c(14, 12.5), 
  label = "ast",
  vjust = -0.25) +
  scale_y_continuous(limits = c(9, 15))

# p_ds_vs_vs_micu1 <- view_counts(dds,
#                                 "Micu1") %>%
#   inner_join(colData(dds), copy = T) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   ggplot(aes(x = region,
#              y = log2(count+1),
#              fill = region)) +
#   geom_violin() +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Fraction and Region",
#        y = "Log2 Expression Count",
#        title = "Micu1") +
#   panel_border() +
#   facet_wrap(vars(fraction),
#              nrow = 1) +
#   scale_y_continuous(limits = c(8, NA)) +
#   theme(legend.position = "none")


# view_counts(dds,
#             "Cul3") %>%
#   inner_join(colData(dds_C2), copy = T) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   ggplot(aes(x = region,
#              y = log2(count+1),
#              fill = region)) +
#   geom_violin() +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Fraction and Region",
#        y = "Log2 Expression Count",
#        title = "Micu1") +
#   panel_border() +
#   facet_wrap(vars(fraction, age),
#              nrow = 1) +
#   scale_y_continuous(limits = c(8, NA)) +
#   theme(legend.position = "none")


# ----
# ----
# ----
# AXON vs MB DTU: Counts for IP and TOTAL SETUP ------------

# Keep this code, but start from the readRDS function (the commented code has already been run)
# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/salmon"
# files <-
#   list.files(
#     "input/data/salmon",
#     pattern = "quant.sf",
#     recursive =  T,
#     full.names = T
#   )
# 
# all(file.exists(files)) # everything exists
# 
# names(files) <-
#   str_extract(files, "(?<=salmon\\/)[0-9|\\_]+(?=\\/quant.sf)")
# 
# files <- files[names(files) %in% metadata$fastq_id] # remove blanks
# files <-
#   files[order(match(names(files), metadata$fastq_id))] #reorder based on metadata
# 
# metadata$files_salmon <- files
# metadata$names <- metadata$fastq_id
# 
# # formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/input/index/references/annotation.gtf"
# txdb <-
#   makeTxDbFromGFF(file = "input/index/references/annotation.gtf",
#                   dataSource = "GENCODE",
#                   organism = "Mus musculus")
# 
# tx2gene <-
#   AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
# 
# # import counts
# txi <- tximport(
#   files,
#   txOut = TRUE,
#   type = "salmon",
#   tx2gene = tx2gene,
#   ignoreTxVersion = F,
#   countsFromAbundance = "dtuScaledTPM"
# )
# 
# # remove zero count isoforms
# cts <- txi$counts
# 
# # collapse lanes
# all(colnames(cts) == metadata$fastq_id)
# sample_name <- as.factor(metadata$sample_name)
# length(sample_name) == ncol(cts)
# sp <- split(seq(along = sample_name), sample_name)
# countdata <- sapply(sp, function(i) rowSums(cts[,i, drop = FALSE]))
# dim(countdata)
# mode(countdata) <- "integer"
# # colsToKeep <- sapply(sp, `[`, 1)
# # collapsed <- cts[, colsToKeep]
# # dimnames(countdata) <- dimnames(collapsed)
# collapsed <- countdata
# all(colnames(collapsed) == levels(sample_name))
# cts <- collapsed
# rm(list = c("collapsed", "countdata"))
# cts <- cts[rowSums(cts) > 0, ]
# # check all rownames of cts are in tx2gene
# all(rownames(cts) %in% tx2gene$TXNAME)
# # match the rownames of tx2gene to cts
# tx2gene <- tx2gene[match(rownames(cts), tx2gene$TXNAME), ]
# all(rownames(cts) == tx2gene$TXNAME)
# 
# # collapse lane replicates
# # cts0 <- cts %>%
# #   as_tibble(rownames = "tx_id") %>%
# #   pivot_longer(
# #     -tx_id,
# #     names_to = c("lane", "sample_code"),
# #     names_sep = "_",
# #     values_to = "count"
# #   ) %>%
# #   group_by(tx_id,
# #            sample_code) %>%
# #   summarise(count = sum(count))
# 
# # remove duplicated features
# length(rownames(cts)) == length(unique(rownames(cts)))
# # create a data.frame for DRIMSeq
# counts <- data.frame(gene_id = tx2gene$GENEID,
#                      feature_id = tx2gene$TXNAME,
#                      cts,
#                      check.names = F) # it changes hyphens to full stops! C2 TOTAL names get changed
# 
# # filter counts for genes enriched in AXON and MB (matching the substring name - version number removed temporarily)
# # doing this for all IP vs IP comparisons because we need confidence that the gene
# # is dopaminergic in origin. For fraction comparisons, not filtering on enrichment
# # bexause there could be DTU without enrichment
# # Update: 21st May: Including axon enriched AGNOSTIC GENES: Can filter on ABA striatal intensity afterwards
# # Update: 24th May: Back to axon enriched + soma enriched list : The only genes I
# # prioritise from AXON VS MB gene-level are all soma + axon enriched. so stick with this level of evidence
# # Update 25th: Now including all detectable genes and filtering afterwards:
# # There is an unbalanced shift of differential expression towards the soma in
# # soma/axon enriched genes. In AXON VS MB DESeq2, all detectable genes were
# # compared and then filtered by enrichment evidence. So now doing the same approach here.
# 
# # counts <- counts[substr(counts$gene_id, 1, 18) %in% AXON_FRACTION_ENRICHED_GENES,]
# 
# # filter samples
# # MB_YOUNG_WT_MIXED_C2.TOTAL1.5
# # MB_YOUNG_WT_MIXED_C2.TOTAL1-5
# 
# # JUST USING COHORT 2 SAMPLES!
# # EDIT: 24th MAY: Adding C1 MB samples so that there are more TOTAL samples for comparison
# counts_IP <-
#   counts[, colnames(counts) %in% c("gene_id",
#                                    "feature_id",
#                                    (
#                                      metadata %>% filter(
#                                        fraction == "IP" &
#                                          cohort %in% c("C1", "C2") &
#                                          collection != "C1.IPPOOL" &
#                                          !sample_name %in% c("DS_OLD_WT_MIXED_C2.3",
#                                                              "VS_YOUNG_WT_MIXED_C2.3",
#                                                              "MB_OLD_WT_MIXED_C2.6",
#                                                              "MB_YOUNG_WT_MIXED_C2.7",
#                                                              "DS_OLD_WT_MIXED_C2.TOTAL1-3")
#                                      ) %>% pull(sample_name) %>% unique
#                                    ))]
# 
# counts_TOTAL <-
#   counts[, colnames(counts) %in% c("gene_id",
#                                    "feature_id",
#                                    (
#                                      metadata %>% filter(
#                                        fraction == "TOTAL" &
#                                          cohort %in% c("C1", "C2") &
#                                          !sample_name %in% c("DS_OLD_WT_MIXED_C2.3",
#                                                              "VS_YOUNG_WT_MIXED_C2.3",
#                                                              "MB_OLD_WT_MIXED_C2.6",
#                                                              "MB_YOUNG_WT_MIXED_C2.7",
#                                                              "DS_OLD_WT_MIXED_C2.TOTAL1-3")
#                                      ) %>% pull(sample_name) %>% unique
#                                    ))]
# 
# 
# saveRDS(counts_IP, "R/objects/counts_IP.rds")
# saveRDS(counts_TOTAL, "R/objects/counts_TOTAL.rds")

counts_IP <- readRDS("R/objects/counts_IP.rds")
counts_TOTAL <- readRDS("R/objects/counts_TOTAL.rds")

# AXON VS MB DTU: TOTAL DRIMSeq -----

counts <- counts_TOTAL

# create sample metadata for DRIMSeq
samps <- metadata %>%
  filter(sample_name %in% colnames(counts)) %>%
  select("sample_id" = sample_name,
         cohort,
         compartment) %>%
  distinct
samps <- as.data.frame(samps)
samps$cohort <- relevel(samps$cohort, ref = "C2")
samps$compartment <- relevel(samps$compartment, ref = "MB")

# create drimseq object
d <- dmDSdata(counts = counts,
              samples = samps)

# DRIMSeq filter
n <- nrow(samps)

table(samples(d)$cohort, samples(d)$compartment)
# n.small <- sum(samps$compartment == "MB" & samps$cohort == "C3") # all the cohort 2 MB samples multiplied by 2
# n.small <- 44
n.small <- sum(samps$compartment == "AXON")
d <- dmFilter(
  d,
  min_samps_gene_expr = n * 0.75,
  min_gene_expr = 10,
  min_samps_feature_expr = n.small,
  min_feature_expr = 10,
  min_samps_feature_prop = n.small,
  min_feature_prop = 0.1,
)

saveRDS(d, "R/objects/d_TOTAL_COMPARTMENT.rds")

# code to run in isolate env
# d <- readRDS("R/objects/d_TOTAL_COMPARTMENT.rds")
# 
# design_full <-
#   model.matrix( ~ cohort + compartment, data = DRIMSeq::samples(d))
# colnames(design_full)
# 
# set.seed(821196)
# system.time({
#   d <-
#     dmPrecision(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
#   d <-
#     dmFit(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
#   d <-
#     dmTest(d, coef = "compartmentAXON", BPPARAM = BiocParallel::MulticoreParam())
# })
# 
# saveRDS(d,
#         "R/objects/d_TOTAL_COMPARTMENT_PROCESSED.rds")

d_TOTAL_COMPARTMENT <- readRDS("R/objects/d_TOTAL_COMPARTMENT_PROCESSED.rds")

# get results
res_TOTAL_COMPARTMENT <- DRIMSeq::results(d_TOTAL_COMPARTMENT)
res.txp_TOTAL_COMPARTMENT <- DRIMSeq::results(d_TOTAL_COMPARTMENT, level = "feature")

# replace NA values with 1
no.na <- function(x){
  ifelse(is.na(x), 1, x)
}

res_TOTAL_COMPARTMENT$pvalue <- no.na(res_TOTAL_COMPARTMENT$pvalue)
res.txp_TOTAL_COMPARTMENT$pvalue <- no.na(res.txp_TOTAL_COMPARTMENT$pvalue)

# remove gene id version number
res_TOTAL_COMPARTMENT$gene_id_simple <- substr(res_TOTAL_COMPARTMENT$gene_id, 1, 18)


# AXON VS MB DTU: IP DRIMSeq -----

counts <- counts_IP

# create sample metadata for DRIMSeq
samps <- metadata %>%
  filter(sample_name %in% colnames(counts)) %>%
  select("sample_id" = sample_name,
         cohort,
         compartment) %>%
  distinct
samps <- as.data.frame(samps)
samps$cohort <- relevel(samps$cohort, ref = "C2")
samps$compartment <- relevel(samps$compartment, ref = "MB")

# create drimseq object
d <- dmDSdata(counts = counts,
              samples = samps)

methods(class = class(d))
#
substr(counts(d)$gene_id, 1, 18) %>% unique

# DRIMSeq filter
n <- nrow(samps)

table(samples(d)$cohort, samples(d)$compartment)
# n.small <- sum(samps$compartment == "MB" & samps$cohort == "C3") # all the cohort 2 MB samples multiplied by 2
# n.small <- 44
n.small <- sum(samps$compartment == "AXON")
d <- dmFilter(
  d,
  min_samps_gene_expr = n * 0.75,
  min_gene_expr = 10,
  min_samps_feature_expr = n.small,
  min_feature_expr = 10,
  min_samps_feature_prop = n.small,
  min_feature_prop = 0.1,
)

saveRDS(d, "R/objects/d_IP_COMPARTMENT.rds")

# code to run in isolate env
# d <- readRDS("R/objects/d_IP_COMPARTMENT.rds")
# 
# design_full <-
#   model.matrix( ~ cohort + compartment, data = DRIMSeq::samples(d))
# colnames(design_full)
# 
# set.seed(821196)
# system.time({
#   d <-
#     dmPrecision(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
#   d <-
#     dmFit(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
#   d <-
#     dmTest(d, coef = "compartmentAXON", BPPARAM = BiocParallel::MulticoreParam())
# })
# 
# saveRDS(d,
#         "R/objects/d_IP_COMPARTMENT_PROCESSED.rds")

d_IP_COMPARTMENT <- readRDS("R/objects/d_IP_COMPARTMENT_PROCESSED.rds")

# get results
res_IP_COMPARTMENT <- DRIMSeq::results(d_IP_COMPARTMENT)
res.txp_IP_COMPARTMENT <- DRIMSeq::results(d_IP_COMPARTMENT, level = "feature")

# replace NA values with 1
no.na <- function(x){
  ifelse(is.na(x), 1, x)
}

res_IP_COMPARTMENT$pvalue <- no.na(res_IP_COMPARTMENT$pvalue)
res.txp_IP_COMPARTMENT$pvalue <- no.na(res.txp_IP_COMPARTMENT$pvalue)

# remove gene id version number
res_IP_COMPARTMENT$gene_id_simple <- substr(res_IP_COMPARTMENT$gene_id, 1, 18)

# NOT DOING STAGE-R FOR THESIS - EDIT (2021-09-29): Now including!
# IP #############################
# get annotation
anno_drimseq <- get_anno(res_IP_COMPARTMENT$gene_id_simple, get_human = TRUE)

# stageR
# create a vector of gene pvalues, stripped of version number
pScreen_IP_COMPARTMENT <- res_IP_COMPARTMENT$pvalue
strp <- function(x) substr(x,1,18)
names(pScreen_IP_COMPARTMENT) <- strp(res_IP_COMPARTMENT$gene_id)

# create a vector of transcript pvalues
pConfirmation_IP_COMPARTMENT <- matrix(res.txp_IP_COMPARTMENT$pvalue, ncol=1)
rownames(pConfirmation_IP_COMPARTMENT) <- strp(res.txp_IP_COMPARTMENT$feature_id)

# create a df with transcript and gene ids (tx2gene)
tx2gene <- res.txp_IP_COMPARTMENT[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# run stageR at alpha of 0.05
stageRObj_IP_COMPARTMENT <- stageRTx(pScreen=pScreen_IP_COMPARTMENT, pConfirmation=pConfirmation_IP_COMPARTMENT,
                                  pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj_IP_COMPARTMENT <- stageWiseAdjustment(stageRObj_IP_COMPARTMENT, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj_IP_COMPARTMENT <- getAdjustedPValues(stageRObj_IP_COMPARTMENT, order=FALSE,
                                              onlySignificantGenes=TRUE)
})
drim.padj_IP_COMPARTMENT <- drim.padj_IP_COMPARTMENT %>%
  as_tibble %>%
  left_join(anno,
            by = c("geneID" = "ensembl_gene_id"))

# TOTAL #############################
# get annotation
anno_drimseq <- get_anno(res_TOTAL_COMPARTMENT$gene_id_simple, get_human = TRUE)

# stageR
# create a vector of gene pvalues, stripped of version number
pScreen_TOTAL_COMPARTMENT <- res_TOTAL_COMPARTMENT$pvalue
strp <- function(x) substr(x,1,18)
names(pScreen_TOTAL_COMPARTMENT) <- strp(res_TOTAL_COMPARTMENT$gene_id)

# create a vector of transcript pvalues
pConfirmation_TOTAL_COMPARTMENT <- matrix(res.txp_TOTAL_COMPARTMENT$pvalue, ncol=1)
rownames(pConfirmation_TOTAL_COMPARTMENT) <- strp(res.txp_TOTAL_COMPARTMENT$feature_id)

# create a df with transcript and gene ids (tx2gene)
tx2gene <- res.txp_TOTAL_COMPARTMENT[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# run stageR at alpha of 0.05
stageRObj_TOTAL_COMPARTMENT <- stageRTx(pScreen=pScreen_TOTAL_COMPARTMENT, pConfirmation=pConfirmation_TOTAL_COMPARTMENT,
                                     pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj_TOTAL_COMPARTMENT <- stageWiseAdjustment(stageRObj_TOTAL_COMPARTMENT, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj_TOTAL_COMPARTMENT <- getAdjustedPValues(stageRObj_TOTAL_COMPARTMENT, order=FALSE,
                                                 onlySignificantGenes=TRUE)
})
drim.padj_TOTAL_COMPARTMENT <- drim.padj_TOTAL_COMPARTMENT %>%
  as_tibble %>%
  left_join(anno,
            by = c("geneID" = "ensembl_gene_id"))

# post-hoc filtering on the standard deviation in proportions
# res.txp.filt <- DRIMSeq::results(d, level="feature")
# smallProportionSD <- function(d, filter=0.1) {
#   cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
#   gene.cts <- rowsum(cts, counts(d)$gene_id)
#   total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
#   props <- cts/total.cts
#   propSD <- sqrt(rowVars(props))
#   propSD < filter
# }
# filt <- smallProportionSD(d)
# res.txp.filt$pvalue[filt] <- 1
# res.txp.filt$adj_pvalue[filt] <- 1

# res_COMPARTMENT %>%
#   as_tibble %>%
#   arrange(adj_pvalue) %>%
#   filter(adj_pvalue < 0.01) %>%
#   inner_join(anno, by = c("gene_id_simple" = "ensembl_gene_id"))

# mutate(gwas = hsapiens_homolog_ensembl_gene %in% pd_gwas_genes_human) %>%
# filter(external_gene_name %in% t_gwas_nalls$`TRAP Top Candidate`) %>%

# AXON VS MB DTU: PLOTTING FUNCTION ----

# rabgap1l
data <- bind_rows(
  {
    counts(d_IP_COMPARTMENT) %>%
      as_tibble %>%
      filter(str_detect(gene_id, get_ensembl_gene_id("Rabgap1l"))) %>%
      select(ensembl_gene_id = gene_id,
             ensembl_transcript_id = feature_id,
             everything()) %>%
      mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                           "[:alnum:]+(?=\\.[:digit:]+)"),
             ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                                 "[:alnum:]+(?=\\.[:digit:]+)")) %>%
      left_join(anno_tx) %>%
      pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
                   names_to = "sample_name",
                   values_to = "count") %>%
      group_by(sample_name) %>%
      filter(median(count) > 0) %>%
      mutate(count = count / sum(count)) %>%
      inner_join(colData(dds), copy = TRUE) %>%
      mutate(fraction = "TRAP")
  }, {
    counts(d_TOTAL_COMPARTMENT) %>%
      as_tibble %>%
      filter(str_detect(gene_id, get_ensembl_gene_id("Rabgap1l"))) %>%
      select(ensembl_gene_id = gene_id,
             ensembl_transcript_id = feature_id,
             everything()) %>%
      mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                           "[:alnum:]+(?=\\.[:digit:]+)"),
             ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                                 "[:alnum:]+(?=\\.[:digit:]+)")) %>%
      left_join(anno_tx) %>%
      pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
                   names_to = "sample_name",
                   values_to = "count") %>%
      group_by(sample_name) %>%
      filter(median(count) > 0) %>%
      mutate(count = count / sum(count)) %>%
      mutate(compartment = ifelse(str_detect(sample_name, "MB"), "MB", "AXON")) %>%
      mutate(fraction = "TOTAL")
  }
) %>%
  mutate(fraction = factor(fraction, levels = c("TRAP", "TOTAL")),
         compartment = ifelse(compartment == "MB", "Soma", "Axon")) 



pval_IP <- res_IP_COMPARTMENT %>%
  left_join(anno, by = c("gene_id_simple" = "ensembl_gene_id")) %>%
  filter(external_gene_name == "Rabgap1l") %>%
  mutate(p = signif(adj_pvalue, 3)) %>%
  pull(p)

pval_TOTAL <- res_TOTAL_COMPARTMENT %>%
  left_join(anno, by = c("gene_id_simple" = "ensembl_gene_id")) %>%
  filter(external_gene_name == "Rabgap1l") %>%
  mutate(p = signif(adj_pvalue, 3)) %>%
  pull(p)

stat.test <- drim.padj_IP_COMPARTMENT %>%
  filter(external_gene_name == "Rabgap1l") %>%
  left_join(anno_tx, by = c("txID" = "ensembl_transcript_id")) %>%
  select(external_transcript_name, gene, transcript)

stat.test <- stat.test %>%
  mutate(group1 = "Soma", group2 = "Axon") %>%
  arrange(external_transcript_name) %>%
  mutate(fraction = paste0("TRAP: *P* = ", pval_IP)) %>%
  pval_asterisks(transcript)

data <- data %>%
  mutate(fraction = ifelse(fraction == "TRAP", 
                           paste0("TRAP: *P* = ", pval_IP), 
                           paste0("TOTAL: *P* = ", pval_TOTAL)))

p_axon_vs_mb_dtu_rabgap1l <- ggboxplot(data, 
                                   x = "compartment", 
                                   y = "count", 
                                   fill = "compartment",
                                   facet.by = c("fraction", "external_transcript_name"), 
                                   outlier.shape = NA) +
  theme_PK(font_size = 11) +
  scale_fill_d3() +
  theme(legend.position = "none", 
        strip.text.y = ggtext::element_markdown(size = 8)) +
  scale_y_continuous(limits = c(0, 1), 
                     sec.axis = sec_axis(~ . , name = "TRAP Fraction", breaks = NULL, labels = NULL)) +
  labs(y = "Proportion of Total Expression", 
       x = "Compartment", 
       title = "Rabgap1l") +
  stat_pvalue_manual(stat.test, 
                     y.position = c(0.8, 0.9, 0.45), 
                     label = "ast", 
                     size = 3, 
                     vjust = -0.25) +
  panel_border()

# p_axon_vs_mb_dtu_rabgap1l <- bind_rows(
#   {
#     counts(d_IP_COMPARTMENT) %>%
#       as_tibble %>%
#       filter(str_detect(gene_id, "ENSMUSG00000026721")) %>%
#       select(ensembl_gene_id = gene_id,
#              ensembl_transcript_id = feature_id,
#              everything()) %>%
#       mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
#                                            "[:alnum:]+(?=\\.[:digit:]+)"),
#              ensembl_transcript_id = str_extract(ensembl_transcript_id,
#                                                  "[:alnum:]+(?=\\.[:digit:]+)")) %>%
#       left_join(anno_tx) %>%
#       pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
#                    names_to = "sample_name",
#                    values_to = "count") %>%
#       group_by(sample_name) %>%
#       filter(median(count) > 0) %>%
#       mutate(count = count / sum(count)) %>%
#       inner_join(colData(dds), copy = TRUE) %>%
#       mutate(fraction = "TRAP")
#   }, {
#     counts(d_TOTAL_COMPARTMENT) %>%
#       as_tibble %>%
#       filter(str_detect(gene_id, "ENSMUSG00000026721")) %>%
#       select(ensembl_gene_id = gene_id,
#              ensembl_transcript_id = feature_id,
#              everything()) %>%
#       mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
#                                            "[:alnum:]+(?=\\.[:digit:]+)"),
#              ensembl_transcript_id = str_extract(ensembl_transcript_id,
#                                                  "[:alnum:]+(?=\\.[:digit:]+)")) %>%
#       left_join(anno_tx) %>%
#       pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
#                    names_to = "sample_name",
#                    values_to = "count") %>%
#       group_by(sample_name) %>%
#       filter(median(count) > 0) %>%
#       mutate(count = count / sum(count)) %>%
#       mutate(compartment = ifelse(str_detect(sample_name, "MB"), "MB", "AXON")) %>%
#       mutate(fraction = "TOTAL")
#   }
# ) %>%
#   mutate(fraction = factor(fraction, levels = c("TRAP", "TOTAL")),
#          compartment = ifelse(compartment == "MB", "Soma", "Axon")) %>%
#   ggplot(aes(x = compartment,
#              y = count,
#              fill = compartment)) +
#   geom_violin() +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Compartment",
#        y = "Proportion of Total Expression",
#        title = "Rabgap1l") +
#   theme(legend.position = "none") +
#   facet_grid(rows = vars(fraction),
#              cols = vars(external_transcript_name)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   panel_border()


# oxr1
data <- bind_rows(
  {
    counts(d_IP_COMPARTMENT) %>%
      as_tibble %>%
      filter(str_detect(gene_id, get_ensembl_gene_id("Oxr1"))) %>%
      select(ensembl_gene_id = gene_id,
             ensembl_transcript_id = feature_id,
             everything()) %>%
      mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                           "[:alnum:]+(?=\\.[:digit:]+)"),
             ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                                 "[:alnum:]+(?=\\.[:digit:]+)")) %>%
      left_join(anno_tx) %>%
      pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
                   names_to = "sample_name",
                   values_to = "count") %>%
      group_by(sample_name) %>%
      filter(median(count) > 0) %>%
      mutate(count = count / sum(count)) %>%
      inner_join(colData(dds), copy = TRUE) %>%
      mutate(fraction = "TRAP")
  }, {
    counts(d_TOTAL_COMPARTMENT) %>%
      as_tibble %>%
      filter(str_detect(gene_id, get_ensembl_gene_id("Oxr1"))) %>%
      select(ensembl_gene_id = gene_id,
             ensembl_transcript_id = feature_id,
             everything()) %>%
      mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
                                           "[:alnum:]+(?=\\.[:digit:]+)"),
             ensembl_transcript_id = str_extract(ensembl_transcript_id,
                                                 "[:alnum:]+(?=\\.[:digit:]+)")) %>%
      left_join(anno_tx) %>%
      pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
                   names_to = "sample_name",
                   values_to = "count") %>%
      group_by(sample_name) %>%
      filter(median(count) > 0) %>%
      mutate(count = count / sum(count)) %>%
      mutate(compartment = ifelse(str_detect(sample_name, "MB"), "MB", "AXON")) %>%
      mutate(fraction = "TOTAL")
  }
) %>%
  mutate(fraction = factor(fraction, levels = c("TRAP", "TOTAL")),
         compartment = ifelse(compartment == "MB", "Soma", "Axon")) 

pval_IP <- res_IP_COMPARTMENT %>%
  left_join(anno, by = c("gene_id_simple" = "ensembl_gene_id")) %>%
  filter(external_gene_name == "Oxr1") %>%
  mutate(p = signif(adj_pvalue, 3)) %>%
  pull(p)

pval_TOTAL <- res_TOTAL_COMPARTMENT %>%
  left_join(anno, by = c("gene_id_simple" = "ensembl_gene_id")) %>%
  filter(external_gene_name == "Oxr1") %>%
  mutate(p = signif(adj_pvalue, 3)) %>%
  pull(p)

stat.test <- drim.padj_IP_COMPARTMENT %>%
  filter(external_gene_name == "Oxr1") %>%
  left_join(anno_tx, by = c("txID" = "ensembl_transcript_id")) %>%
  select(external_transcript_name, gene, transcript)

stat.test <- stat.test %>%
  mutate(group1 = "Soma", group2 = "Axon") %>%
  arrange(external_transcript_name) %>%
  mutate(fraction = paste0("TRAP: *P* = ", pval_IP)) %>%
  pval_asterisks(transcript)

data <- data %>%
  mutate(fraction = ifelse(fraction == "TRAP", 
                           paste0("TRAP: *P* = ", pval_IP), 
                           paste0("TOTAL: *P* = ", pval_TOTAL)))

p_axon_vs_mb_dtu_oxr1 <- ggboxplot(data, 
                                       x = "compartment", 
                                       y = "count", 
                                       fill = "compartment",
                                       facet.by = c("fraction", "external_transcript_name"), 
                                       outlier.shape = NA) +
  theme_PK(font_size = 11) +
  scale_fill_d3() +
  theme(legend.position = "none", 
        strip.text.y = ggtext::element_markdown(size = 8)) +
  scale_y_continuous(limits = c(0, 1), 
                     sec.axis = sec_axis(~ . , name = "TRAP Fraction", breaks = NULL, labels = NULL)) +
  labs(y = "Proportion of Total Expression", 
       x = "Compartment", 
       title = "Oxr1") +
  stat_pvalue_manual(stat.test, 
                     y.position = c(0.4, 0.4, 0.9, 0.7, 0.45), 
                     label = "ast", 
                     size = 3, 
                     vjust = -0.25) +
  panel_border()

# p_axon_vs_mb_dtu_oxr1 <- bind_rows(
#   {
#     counts(d_IP_COMPARTMENT) %>%
#       as_tibble %>%
#       filter(str_detect(gene_id, get_ensembl_gene_id("Oxr1"))) %>%
#       select(ensembl_gene_id = gene_id,
#              ensembl_transcript_id = feature_id,
#              everything()) %>%
#       mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
#                                            "[:alnum:]+(?=\\.[:digit:]+)"),
#              ensembl_transcript_id = str_extract(ensembl_transcript_id,
#                                                  "[:alnum:]+(?=\\.[:digit:]+)")) %>%
#       left_join(anno_tx) %>%
#       pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
#                    names_to = "sample_name",
#                    values_to = "count") %>%
#       group_by(sample_name) %>%
#       filter(median(count) > 0) %>%
#       mutate(count = count / sum(count)) %>%
#       inner_join(colData(dds), copy = TRUE) %>%
#       mutate(fraction = "TRAP")
#   }, {
#     counts(d_TOTAL_COMPARTMENT) %>%
#       as_tibble %>%
#       filter(str_detect(gene_id, get_ensembl_gene_id("Oxr1"))) %>%
#       select(ensembl_gene_id = gene_id,
#              ensembl_transcript_id = feature_id,
#              everything()) %>%
#       mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
#                                            "[:alnum:]+(?=\\.[:digit:]+)"),
#              ensembl_transcript_id = str_extract(ensembl_transcript_id,
#                                                  "[:alnum:]+(?=\\.[:digit:]+)")) %>%
#       left_join(anno_tx) %>%
#       pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
#                    names_to = "sample_name",
#                    values_to = "count") %>%
#       group_by(sample_name) %>%
#       filter(median(count) > 0) %>%
#       mutate(count = count / sum(count)) %>%
#       mutate(compartment = ifelse(str_detect(sample_name, "MB"), "MB", "AXON")) %>%
#       mutate(fraction = "TOTAL")
#   }
# ) %>%
#   mutate(fraction = factor(fraction, levels = c("TRAP", "TOTAL")),
#          compartment = ifelse(compartment == "MB", "Soma", "Axon")) %>%
#   ggplot(aes(x = compartment,
#              y = count)) +
#   # geom_violin() +
#   geom_boxplot(aes(fill = compartment),
#                outlier.shape = NA) +
#   # geom_quasirandom(size = 0.5, shape = 21) +
#   # geom_quasirandom(shape = 21,
#   #                  size = 2,
#   #                  colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Compartment",
#        y = "Proportion of Total Expression",
#        title = "Oxr1") +
#   theme(legend.position = "none") +
#   facet_grid(rows = vars(fraction),
#              cols = vars(external_transcript_name)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   panel_border()





# AXON VS MB DTU: TOTAL SUMMARY ----

AXON_VS_MB_TOTAL_DTU <- res_TOTAL_COMPARTMENT %>%
  as_tibble %>%
  arrange(adj_pvalue) %>%
  select(ensembl_gene_id = gene_id_simple, adj_pvalue) %>%
  left_join(anno) %>%
  left_join(aba_expression_axon) %>%
  left_join(publications_all_genes) %>%
  mutate(signif = adj_pvalue < 0.01,
         mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
         # axon_abundant = ensembl_gene_id %in% AXON_VS_MB_META_SIGNIF_PRIORITY$ensembl_gene_id,
         aba_specific = expression_energy < 0.5,
         aba_specific = replace_na(aba_specific, FALSE),
         evidence = ifelse(axon_translated,
                           ifelse(aba_specific & mb_translated,
                                  "Axon, Soma and ABA",
                                  ifelse(aba_specific,
                                         "Axon and ABA",
                                         ifelse(mb_translated,
                                                "Axon and Soma",
                                                "Axon only"))),
                           ifelse(aba_specific & mb_translated,
                                  "Soma and ABA",
                                  "Soma only")),
         evidence = ifelse(!mb_translated & evidence == "Soma only",
                           "No support", evidence),
         evidence = factor(evidence, levels = c("Axon, Soma and ABA",
                                                "Axon and Soma",
                                                "Axon and ABA",
                                                "Axon only",
                                                "Soma and ABA",
                                                "Soma only",
                                                "No support"))) %>%
  filter(signif)

AXON_VS_MB_TOTAL_DTU %>%
  ggplot(aes(x = evidence)) +
  geom_bar()

# AXON VS MB DTU: IP SUMMARY ----

AXON_VS_MB_TRAP_DTU <- res_IP_COMPARTMENT %>%
  as_tibble %>%
  arrange(adj_pvalue) %>%
  select(ensembl_gene_id = gene_id_simple, adj_pvalue) %>%
  left_join(anno) %>%
  left_join(aba_expression_axon) %>%
  left_join(publications_all_genes) %>%
  left_join(publications_all_genes_axon) %>%
  mutate(signif = adj_pvalue < 0.01,
         mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
         # axon_abundant = ensembl_gene_id %in% AXON_VS_MB_META_SIGNIF_PRIORITY$ensembl_gene_id,
         aba_specific = expression_energy < 0.5,
         aba_specific = replace_na(aba_specific, FALSE),
         evidence = ifelse(axon_translated,
                           ifelse(aba_specific & mb_translated,
                                  "Axon, Soma and ABA",
                                  ifelse(aba_specific,
                                         "Axon and ABA",
                                         ifelse(mb_translated,
                                                "Axon and Soma",
                                                "Axon only"))),
                           ifelse(aba_specific & mb_translated,
                                  "Soma and ABA",
                                  "Soma only")),
         evidence = ifelse(!mb_translated & evidence == "Soma only",
                           "No support", evidence),
         evidence = factor(evidence, levels = c("Axon, Soma and ABA",
                                                "Axon and Soma",
                                                "Axon and ABA",
                                                "Axon only",
                                                "Soma and ABA",
                                                "Soma only",
                                                "No support")),
         specific = !ensembl_gene_id %in% AXON_VS_MB_TOTAL_DTU$ensembl_gene_id)

# AXON VS MB DTU: NUMBERS ----

# sort this out: how many genes did we really consider to begin with?
n_AXON_VS_MB_DTU_CONSIDERED <- length(unique(counts_IP$gene_id))
n_AXON_VS_MB_DTU_CONSIDERED <- nrow(AXON_VS_MB_META)

n_AXON_VS_MB_DTU_RETAINED <- nrow(AXON_VS_MB_TRAP_DTU)

n_AXON_VS_MB_DTU_FILTER_PERCENT <- 100*signif(n_AXON_VS_MB_DTU_RETAINED/n_AXON_VS_MB_DTU_CONSIDERED, 2)

n_AXON_VS_MB_DTU_ANYCONFIDENCE <- AXON_VS_MB_TRAP_DTU %>%
  filter(adj_pvalue < 0.01) %>%
  nrow

# AXON VS MB DTU: IP vs TOTAL PLOT ----

# Filtered for genes enriched in either axon or soma
p_axon_vs_mb_dtu_specificity <- AXON_VS_MB_TRAP_DTU %>%
  filter(signif) %>%
  filter(mb_translated) %>%
  left_join(enrichment_status) %>%
  select(-evidence) %>%
  rename(evidence = enrichment) %>%
  # mutate(evidence = ifelse(evidence == "Axon, Soma and ABA",
  #                          "Axon, Soma\nand ABA",
  #                          ifelse(evidence == "Axon and Soma",
  #                                 "Axon and\nSoma",
  #                                 ifelse(evidence == "Soma and ABA",
  #                                        "Soma and\nABA",
  #                                        "Soma only")))) %>%
  group_by(evidence) %>%
  mutate(n = dplyr::n()) %>%
  mutate(specific = ifelse(ensembl_gene_id %in% {AXON_VS_MB_TOTAL_DTU %>%
      filter(signif) %>%
      pull(ensembl_gene_id)},
      "Shared", "Specific")) %>%
  ggplot(aes(x = evidence,
             fill = specific)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_text(aes(label = n,
                y = n),
            check_overlap = T,
            vjust = -0.5,
            fontface = 2,
            size = 4) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ifelse(..count.. > 50, ..count.., "")),
             position=position_stack(vjust=0.5)) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 9)) +
  # guides(fill = guide_legend(nrow = 4,
  #                            byrow=TRUE)) +
  labs(x = "Enrichment Category",
       y = "Number of Genes with DTU",
       fill = "Specificity")


# AXON VS MB DTU: Isoforms BAR CHART -----------------------

# counts <- counts_IP
#
# # plot number of genes DTU by number of isoforms
# n_tx <- tibble(gene_id = unique(counts$gene_id),
#                n_tx = sapply(unique(counts$gene_id), function(x) sum(counts(d_IP_COMPARTMENT)$gene_id == x))) %>%
#   mutate(signif = ifelse(gene_id %in% res_IP_COMPARTMENT$gene_id[res_IP_COMPARTMENT$adj_pvalue < 0.01], "DTU", "NS"))
# n_tx$n_tx[table(counts$gene_id)[match(unique(counts$gene_id), names(table(counts$gene_id)))] == 1] <- 1
# n_tx$n_tx[n_tx$n_tx == 0] <- "Filtered"
# n_tx$n_tx <- factor(n_tx$n_tx, levels = c("Filtered", seq(1, 6, 1)))

n_tx <- readRDS("R/objects/n_tx.rds")

# Filtered for genes enriched in either axon or soma
p_axon_vx_mb_dtu_bar <- n_tx %>%
  filter(substr(gene_id, 1, 18) %in% unique(MB_FRACTION_TRANSLATED_GENES,
                                            AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)) %>%
  filter(!is.na(n_tx) &
           n_tx != "6") %>%
  ggplot(aes(x = n_tx,
             fill = signif)) +
  geom_bar(colour = "black") +
  stat_count(geom = "text", colour = "white", size = 3.5, fontface = 2,
             aes(label = ifelse(..count.. < 100, "", ..count..)),
             position=position_stack(vjust=0.5)) +
  scale_fill_manual(values = pal_d3()(2)[c(2, 1)]) +
  labs(x = "Number of Isoforms",
       y = "Number of Genes",
       fill = "Outcome") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = c(0.75, 0.75))


# ----
# ----
# ----
# MB AGE - Plan ----
# C3 ONT is now being used in place of C3 short read, based on correlation findings.
# MB AGE - LIBRARY COMPOSITION CHECK ----

C1_MB_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
      select(age, genotype, "Slc6a3" = count)
  }, {
    plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
      select("Th" = count)
    # }, {
    #   plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
    #     select("Ddc" = count)
  }, {
    plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-c(age, genotype),
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C1")


C2_MB_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
      select(age, genotype, "Slc6a3" = count)
  }, {
    plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
      select("Th" = count)
    # }, {
    #   plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
    #     select("Ddc" = count)
  }, {
    plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-c(age, genotype),
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C2")


C3_MB_IP_MARKERS <- bind_cols(
  {
    plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
      select(age, genotype, "Slc6a3" = count)
  }, {
    plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
      select("Th" = count)
    # }, {
    #   plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
    #     select("Ddc" = count)
  }, {
    plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-c(age, genotype),
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C3")





C3_MB_IP_ONT_MARKERS <- bind_cols(
  {
    plotCounts(dds_ONT, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
      select(age, genotype, "Slc6a3" = count)
  }, {
    plotCounts(dds_ONT, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
      select("Th" = count)
    # }, {
    #   plotCounts(dds_ONT, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
    #     select("Ddc" = count)
  }, {
    plotCounts(dds_ONT, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
      select("Slc18a2" = count)
  }, {
    plotCounts(dds_ONT, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
      select("Gfap" = count)
  }, {
    plotCounts(dds_ONT, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
      select("Gad2" = count)
  }, {
    plotCounts(dds_ONT, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
      select("S100b" = count)
  }
) %>%
  pivot_longer(-c(age, genotype),
               names_to = "gene",
               values_to = "count") %>%
  mutate(cohort = "C3 ONT")



# Take home points from this plot:
#
# Enrichment can be seen (and depletion of non-DA markers)
#
# In some cohorts (C1 and a little bit C2/3), DA markers go up
# while in C3 ONT, markers go down
# this creates a problem in each cohort individually, as other DA genes
# may be classed as DE
#
# The different collection strategy of C1 shows: DA markers are definitely higher there
# and non DA are lower (this is the only cohort where this is the case)
#
# By using meta-analysis, differences in DA markers and associated genes that follow
# are screened, preserving only non-conflicting changes across cohorts
#
# Other point: non-DA markers are really going up in all cohorts apart from C1
# This emphasises the importance of holding imited value in results from
# DEPLETED genes.
p_mb_ip_age_markers <- bind_rows(C1_MB_IP_MARKERS,
                                 C2_MB_IP_MARKERS,
                                 C3_MB_IP_MARKERS,
                                 C3_MB_IP_ONT_MARKERS) %>%
  mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
                         "Dopaminergic", "Other")) %>%
  # group_by(age, gene) %>%
  # summarise(count = mean(count)) %>%
  ggplot(aes(x = age,
             y = log2(count + 1),
             colour = marker,
             group = gene)) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(cohort), nrow = 1) +
  panel_border() +
  scale_color_d3() +
  labs(x = "Age and Cohort",
       y = "Log2 Expression Count",
       colour = "Cell Type Marker") +
  theme(legend.position = "top")






# C1_MB_TOTAL_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#   }, {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#       select("Ddc" = count)
#   }, {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
#
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C1")
#
#
# C2_MB_TOTAL_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#   }, {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#       select("Ddc" = count)
#   }, {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C2")
#
#
# C3_MB_TOTAL_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#   }, {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#       select("Ddc" = count)
#   }, {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C3")
#
#
# # Not very revealing: Or perhaps it is: DA and non DA markers with overlapping abundance,
# # some up in age, some down, so not a clear picture.
# # Probably not evidence of a difference in library composition here.
# # The number of detected genes is no different, either (First results chapter).
# # This is good!
# p_mb_total_age_library_composition <- bind_rows(C1_MB_TOTAL_MARKERS,
#                                                 C2_MB_TOTAL_MARKERS,
#                                                 C3_MB_TOTAL_MARKERS) %>%
#   mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
#                          "Dopaminergic", "Other")) %>%
#   # group_by(age, gene) %>%
#   # summarise(count = mean(count)) %>%
#   ggplot(aes(x = age,
#              y = log2(count + 1),
#              colour = marker,
#              group = gene)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(vars(cohort), nrow = 1) +
#   panel_border() +
#   scale_color_d3() +
#   labs(x = "Age",
#        y = "Log2 Expression Count",
#        colour = "Cell Type Marker") +
#   theme(legend.position = "top")







# MB AGE - ALTERNATIVE COMPOSITION CHECK ----

LIBRARY_PROPORTIONS <- bind_rows(
  {
    as_tibble(counts(dds_C1_MB_IP), rownames = "ensembl_gene_id") %>%
      mutate(cohort = "C1")
  }, {
    as_tibble(counts(dds_C2_MB_IP), rownames = "ensembl_gene_id") %>%
      mutate(cohort = "C2")
  }, {
    as_tibble(counts(dds_C3_MB_IP), rownames = "ensembl_gene_id") %>%
      mutate(cohort = "C3")
  }, {
    as_tibble(counts(dds_ONT), rownames = "ensembl_gene_id") %>%
      mutate(cohort = "C3 ONT")
  }
) %>%
  pivot_longer(-c(ensembl_gene_id, cohort),
               names_to = "sample_name",
               values_to = "count") %>%
  mutate(status = ifelse(ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
                         "Enriched",
                         ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
                                "Depleted", "Unchanged")),
         status = factor(status, levels = c("Depleted", "Unchanged", "Enriched"))) %>%
  mutate(age = ifelse(str_detect(sample_name, "YOUNG"), "YOUNG", "OLD"),
         age = factor(age, levels = c("YOUNG", "OLD"))) %>%
  group_by(cohort, age, status) %>%
  summarise(count = sum(count, na.rm = T)) %>%
  mutate(prop = count/sum(count))

p_age_mb_library_proportions <- LIBRARY_PROPORTIONS %>%
  ggplot(aes(x = age,
             y = prop,
             fill = status)) +
  geom_col(colour = "black") +
  facet_wrap(vars(cohort), nrow = 1)  +
  theme(legend.position = "top") +
  geom_text(aes(y = prop, label = ifelse(prop > 0.1, scales::percent(prop, accuracy = 1), "")),
            position = position_stack(vjust = 0.5),
            show.legend = FALSE,
            colour = "white",
            fontface = 2) +
  scale_fill_manual(values = pal_d3()(4)[c(4, 1, 3)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::percent) +
  panel_border() +
  labs(x = "Age and Cohort",
       y = "Percentage of library",
       fill = "Enrichment")

# MB AGE TOTAL DESeq2 -----------------------------------------------------
# TRAP will not have been performed on these samples, so there isnt a batch difference to worry about

# C1 TOTAL MB
dds_MB_TOTAL_AGE <- dds_C1_MB_TOTAL %>% filter_zeros()

dds_MB_TOTAL_AGE@design <- ~ age

colData(dds_MB_TOTAL_AGE) <- droplevels(colData(dds_MB_TOTAL_AGE))

dds_MB_TOTAL_AGE <- DESeq(
  dds_MB_TOTAL_AGE,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_MB_TOTAL_AGE <- DESeq2::results(
  dds_MB_TOTAL_AGE,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_MB_TOTAL_AGE <- lfcShrink(
  dds_MB_TOTAL_AGE,
  res = res_MB_TOTAL_AGE,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)

summary(res_MB_TOTAL_AGE)

# dds_C3_MB_TOTAL_AGE <- dds_C3_MB_TOTAL[rownames(dds_C3_MB_TOTAL) %in% genes_MB_IP_AGE,]
#
# dds_C3_MB_TOTAL_AGE@design <- ~ age
#
# colData(dds_C3_MB_TOTAL_AGE) <- droplevels(colData(dds_C3_MB_TOTAL_AGE))
#
# dds_C3_MB_TOTAL_AGE <- DESeq(
#   dds_C3_MB_TOTAL_AGE,
#   sfType = "poscounts",
#   fitType = "local",
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_C3_MB_TOTAL_AGE <- DESeq2::results(
#   dds_C3_MB_TOTAL_AGE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C3_MB_TOTAL_AGE <- lfcShrink(
#   dds_C3_MB_TOTAL_AGE,
#   res = res_C3_MB_TOTAL_AGE,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_C3_MB_TOTAL_AGE)


MB_AGE_TOTAL_META <-
  as_tibble(res_MB_TOTAL_AGE, rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  select(
    ensembl_gene_id,
    external_gene_name,
    description,
    everything()
  ) %>%
  mutate(outcome = ifelse(log2FoldChange > 0,
                          "Up", "Down"),
         outcome = ifelse(padj < 0.1,
                          outcome, "Unchanged"),
         mb_translation = ifelse(ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
                                 "Translated",
                                 ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
                                        "Depleted", "Unchanged")))



MB_AGE_TOTAL_OUTCOMES <- MB_AGE_TOTAL_META %>%
  select(ensembl_gene_id,
         outcome)


# IF ALSO USING C3 DATA
# MB_AGE_TOTAL_META <-
#   as_tibble(res_C1_MB_TOTAL_AGE, rownames = "ensembl_gene_id") %>%
#   left_join(
#     as_tibble(res_C3_MB_TOTAL_AGE, rownames = "ensembl_gene_id"),
#     by = "ensembl_gene_id",
#     suffix = c("_C1", "_C3")
#   ) %>%
#   left_join(anno) %>%
#   select(
#     ensembl_gene_id,
#     external_gene_name,
#     description,
#     everything()
#   ) %>%
#   mutate(
#     across(starts_with("pvalue"), ~ replace_na(.x, 1 - 1e-10)),
#     # replace NA with *almost* 1 (it can't be perfectly 1)
#     across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)),
#     # replace 0 values with minimum non-zero of the column
#     across(starts_with("pvalue"), ~ .x / 2),
#     # divide 2-sided pvalue by 2 into 1-sided
#     across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)) # again replace 0 values with minimum non-zero of the column
#   ) %>%
#   select(sort(names(.))) %>%
#   select(
#     ensembl_gene_id,
#     external_gene_name,
#     description,
#     everything(),
#     -starts_with("padj")
#   ) %>%
#   filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
#   drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
#   rowwise() %>%
#   mutate(
#     pvalue_C3 = ifelse(
#       # correct pvalues where there is a conflict in log2foldchange between groups
#       sign(log2FoldChange_C1) == sign(log2FoldChange_C3),
#       pvalue_C3,
#       1 - pvalue_C3
#     )) %>%
#   mutate(across(starts_with("pvalue"),
#                 ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
#   mutate(
#     sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
#     log2FoldChange = log2FoldChange_C1
#   ) %>%
#   ungroup() %>%
#   mutate(
#     sumz_adj = p.adjust(sumz, method = "fdr"),
#     conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
#   arrange(desc(log2FoldChange)) %>%
#   mutate(outcome = ifelse(log2FoldChange > 0,
#                           "Up", "Down"),
#          outcome = ifelse(sumz_adj < 0.1,
#                           outcome, "Unchanged"),
#          mb_translation = ifelse(ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#                                  "Translated",
#                                  ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                         "Depleted", "Unchanged")))


# MB AGE - C1 C2 C3_ONT MERGED DATASET ----

dds_C1_C2_C3ONT_MB_IP <-
  dds[rownames(dds) %in% rownames(dds_ONT),
      colData(dds)$cohort != "C3" &
        colData(dds)$fraction == "IP" &
        colData(dds)$region == "MB"]
dds_C1_C2_C3ONT_MB_IP_AGE_FILTER <-
  filter_genes(dds_C1_C2_C3ONT_MB_IP,
               grouping = c("cohort", "age", "gene_id"))
dds_ONT_AGE_FILTER <- filter_genes(dds_ONT,
                                   grouping = c("age", "gene_id"))

genes_MB_IP_AGE <- intersect(
  dds_C1_C2_C3ONT_MB_IP_AGE_FILTER,
  dds_ONT_AGE_FILTER
)

counts_C1_C2_C3ONT_MB_IP <- counts(dds_C1_C2_C3ONT_MB_IP)[rownames(dds_C1_C2_C3ONT_MB_IP) %in% genes_MB_IP_AGE,]
counts_ONT <- counts(dds_ONT)[rownames(dds_ONT) %in% genes_MB_IP_AGE,]

all(rownames(counts_ONT) == rownames(counts_C1_C2_C3ONT_MB_IP))

dim(counts_C1_C2_C3ONT_MB_IP)
dim(counts_ONT)

counts_C1_C2_C3ONT_MB_IP <- cbind(counts_C1_C2_C3ONT_MB_IP,
                                  counts_ONT)

coldata_C1_C2_C3ONT_MB_IP <- colData(dds_C1_C2_C3ONT_MB_IP)[,c("sample_name", "cohort", "collection", "age", "genotype")]
coldata_ONT <- colData(dds_ONT)[,c("sample_name", "cohort", "collection", "age", "genotype")]
coldata_C1_C2_C3ONT_MB_IP <- rbind(coldata_C1_C2_C3ONT_MB_IP, coldata_ONT)

dds_C1_C2_C3ONT_MB_IP <- DESeqDataSetFromMatrix(counts_C1_C2_C3ONT_MB_IP,
                                                coldata_C1_C2_C3ONT_MB_IP,
                                                design = ~ collection + age)

# dds_C1_C2_C3ONT_MB_IP <- dds_C1_C2_C3ONT_MB_IP[rownames(dds_C1_C2_C3ONT_MB_IP) %in% MB_FRACTION_TRANSLATED_GENES, ]

dds_C1_C2_C3ONT_MB_IP <- DESeq(dds_C1_C2_C3ONT_MB_IP,
                               minReplicatesForReplace = Inf,
                               parallel = TRUE)

res_C1_C2_C3ONT_MB_IP <- DESeq2::results(
  dds_C1_C2_C3ONT_MB_IP,
  cooksCutoff = Inf,
  filterFun = ihw,
  alpha = 0.01,
  lfcThreshold = log2(1.05),
  parallel = T
)
res_C1_C2_C3ONT_MB_IP <- lfcShrink(
  dds_C1_C2_C3ONT_MB_IP,
  res = res_C1_C2_C3ONT_MB_IP,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  lfcThreshold = log2(1.05),
  parallel = T
)

summary(res_C1_C2_C3ONT_MB_IP)

MB_AGE_META <- res_C1_C2_C3ONT_MB_IP %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  left_join(anno) %>%
  # filter(padj < 0.01) %>%
  mutate(outcome = ifelse(log2FoldChange > 0,
                          "Up", "Down"),
         outcome = ifelse(padj < 0.01, outcome, "Unchanged"),
         status = ifelse(ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
                         "Enriched",
                         ifelse(ensembl_gene_id %in% MB_FRACTION_DTU_NOT_ENRICHED_GENES &
                                  !ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
                                "Spliced",
                                ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
                                       "Depleted", "Unchanged")))) %>%
  left_join(publications_all_genes) %>%
  left_join(MB_AGE_TOTAL_OUTCOMES, by = "ensembl_gene_id", suffix = c("_TRAP", "_TOTAL")) %>%
  mutate(specific = ifelse(outcome_TRAP == outcome_TOTAL,
                           "Shared",
                           ifelse(outcome_TRAP == "Up",
                                  ifelse(outcome_TOTAL == "Unchanged",
                                         "Specific",
                                         "Specific Opposing"),
                                  ifelse(outcome_TOTAL == "Unchanged",
                                         "Specific",
                                         "Specific Opposing"))),
         specific = replace_na(specific, "Specific"))

# MB AGE ENRICHMENT CATEGORY ----

p_age_translation_status <- MB_AGE_META %>%
  filter(padj < 0.01) %>%
  left_join(enrichment_status) %>%
  ggplot(aes(x = enrichment,
             fill = enrichment)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  # scale_fill_manual(values = pal_d3()(5)[c(4, 3, 1, 5)]) +
  labs(x = "Enrichment Category",
       y = "Number of Genes") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7)) +
  geom_text(aes(label = ..count..,
                y = ..count..),
            stat = "count",
            fontface = 2,
            vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# MB_AGE_META %>%
#  filter(padj < 0.01) %>%
#  ggplot(aes(x = status,
#             fill = status)) +
#  geom_bar(colour = "black") +
#  scale_fill_manual(values = pal_d3()(5)[c(4, 3, 1, 5)]) +
#  labs(x = "Translation Category",
#       y = "Number of Genes") +
#  theme(legend.position = "none",
#        axis.text.x = element_text(size = 7)) +
#  geom_text(aes(label = ..count..,
#                y = ..count..),
#            stat = "count",
#            fontface = 2,
#            vjust = -0.5) +
#  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# MB AGE TRAP ENRICHMENT OUTCOME ----

plot_data <- MB_AGE_META %>%
  filter(padj < 0.01) %>%
  left_join(enrichment_status)

p_mb_age_enrichment_outcome <- plot_data %>%
  ggplot(aes(x = outcome_TRAP)) +
  geom_bar(position = "fill",
           colour = "black",
           aes(fill = enrichment)) +
  scale_fill_d3() +
  # scale_fill_manual(values = pal_d3()(5)[c(4, 3, 1, 5)]) +
  labs(x = "Change in Age",
       y = "Percentage of Genes",
       fill = "TRAP\nSpecificity") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  geom_text(data = . %>%
              group_by(outcome_TRAP, enrichment) %>%
              tally() %>%
              mutate(p = n / sum(n)) %>%
              ungroup(),
            aes(y = p,
                group = enrichment,
                label = ifelse(p > 0.05, scales::percent(signif(p, 2)), "")),
            position = position_stack(vjust = 0.5),
            show.legend = FALSE,
            colour = "white",
            fontface = 2) +
  geom_text(data = plot_data %>%
              group_by(outcome_TRAP) %>%
              tally(),
            aes(y = 1,
                label = paste("Total:", n)),
            vjust = -0.5) +
  guides(fill = guide_legend(nrow = 2,
                             byrow=TRUE))

# MB AGE TRAP SPECIFICITY ----

p_mb_age_specificity <- MB_AGE_META %>%
  filter(padj < 0.01) %>%
  filter(status %in% c("Enriched", "Spliced")) %>%
  group_by(outcome_TRAP) %>%
  mutate(n = dplyr::n()) %>%
  ggplot(aes(x = outcome_TRAP,
             fill = specific)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(4, 3, 1)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  geom_text(aes(label = n,
                y = n),
            stat = "identity",
            vjust = -0.5,
            check_overlap = T,
            fontface = 2) +
  labs(x = "Change in Age",
       y = "Number of Genes",
       fill = "TRAP\nSpecificity") +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 3,
                             byrow=TRUE))



# MB AGE - MAKE A LIST OF ALL CONFIDENT DEGS -------

MB_AGE_DE_SPECIFIC_META <- MB_AGE_META %>%
  filter(padj < 0.01) %>%
  filter(specific != "Shared" &
           status %in% c("Enriched", "Spliced")) %>%
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

MB_AGE_DE_SPECIFIC_UP <- MB_AGE_DE_SPECIFIC_META %>%
  filter(log2FoldChange > 0)

MB_AGE_DE_SPECIFIC_DOWN <- MB_AGE_DE_SPECIFIC_META %>%
  filter(log2FoldChange < 0)

# MB AGE Specificity Examples (PLOT) ----
gene_ids <- bind_rows({
  #   MB_AGE_DE_SPECIFIC_META %>%
  #     filter(baseMean > 50 &
  #              specific == "Shared") %>%
  #     filter(external_gene_name == "Ass1") %>% # Or Tyrobp
  #     select(specific,
  #            ensembl_gene_id,
  #            external_gene_name)
  # }, {
  MB_AGE_DE_SPECIFIC_META %>%
    filter(baseMean > 50 &
             specific == "Specific") %>%
    filter(external_gene_name == "Mtfr1") %>%
    select(specific,
           ensembl_gene_id,
           external_gene_name)
}, {
  MB_AGE_DE_SPECIFIC_META %>%
    filter(baseMean > 50 &
             specific == "Specific Opposing") %>%
    filter(external_gene_name == "Scg2") %>%
    select(specific,
           ensembl_gene_id,
           external_gene_name)
}
)

data <- assay(vst(dds_C1_MB)) %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  filter(ensembl_gene_id %in% gene_ids$ensembl_gene_id) %>%
  pivot_longer(-ensembl_gene_id,
               names_to = "sample_name",
               values_to = "count") %>%
  left_join(gene_ids) %>%
  left_join(colData(dds_C1_MB),
            copy = T) %>%
  mutate(specific = factor(specific,
                           levels = c("Shared",
                                      "Specific",
                                      "Specific Opposing"))) %>%
  arrange(specific) %>%
  group_by(specific) %>%
  mutate(number = cur_group_id()) %>%
  ungroup() %>%
  unite(col = group,
        fraction,
        age,
        sep = "\n",
        remove = F) %>%
  unite(col = facet,
        specific,
        external_gene_name,
        sep = ": ",
        remove = F) %>%
  mutate(group = factor(group, levels = c("TOTAL\nYOUNG",
                                          "TOTAL\nOLD",
                                          "IP\nYOUNG",
                                          "IP\nOLD")),
         facet = fct_reorder(facet, number)) 

p_mb_age_deg_specificity_examples <- ggboxplot(data, 
          x = "group", 
          y = "count", 
          fill = "facet", 
          facet.by = "facet", 
          scales = "free", 
          outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 0.5) +
  theme_PK(font_size = 11) +
  panel_border() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 8),
        strip.text = element_text(face = 2)) +
  scale_fill_d3() +
  labs(x = "Fraction & Age",
       y = "Log2 Expression Count") +
  stat_pvalue_manual(bind_rows({
    MB_AGE_META %>%
      filter(external_gene_name %in% c("Mtfr1", "Scg2")) %>%
      select(external_gene_name, padj) %>%
      mutate(group1 = "IP\nOLD", group2 = "IP\nYOUNG", 
             padj = signif(padj, 3))
  }, {
    MB_AGE_TOTAL_META %>%
      filter(external_gene_name %in% c("Mtfr1", "Scg2")) %>%
      select(external_gene_name, padj) %>%
      mutate(group1 = "TOTAL\nOLD", group2 = "TOTAL\nYOUNG", 
             padj = signif(padj, 3))
  }) %>%
    mutate(facet = ifelse(external_gene_name == "Mtfr1", 
                          "Specific: Mtfr1", 
                          "Specific Opposing: Scg2")) %>%
    arrange(facet, group1) %>%
    pval_asterisks(padj), 
  y.position = c(14, 13.6, 9.55, 9), 
  vjust = -0.25, 
  label = "ast",
  size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# p_mb_age_deg_specificity_examples <- assay(vst(dds_C1_MB)) %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(ensembl_gene_id %in% gene_ids$ensembl_gene_id) %>%
#   pivot_longer(-ensembl_gene_id,
#                names_to = "sample_name",
#                values_to = "count") %>%
#   left_join(gene_ids) %>%
#   left_join(colData(dds_C1_MB),
#             copy = T) %>%
#   mutate(specific = factor(specific,
#                            levels = c("Shared",
#                                       "Specific",
#                                       "Specific Opposing"))) %>%
#   arrange(specific) %>%
#   group_by(specific) %>%
#   mutate(number = cur_group_id()) %>%
#   ungroup() %>%
#   unite(col = group,
#         fraction,
#         age,
#         sep = "\n",
#         remove = F) %>%
#   unite(col = facet,
#         specific,
#         external_gene_name,
#         sep = ": ",
#         remove = F) %>%
#   mutate(group = factor(group, levels = c("TOTAL\nYOUNG",
#                                           "TOTAL\nOLD",
#                                           "IP\nYOUNG",
#                                           "IP\nOLD")),
#          facet = fct_reorder(facet, number)) %>%
#   ggplot(aes(x = group,
#              y = count)) +
#   geom_violin(colour = "black",
#               fill = NA,
#               scale = "width") +
#   geom_quasirandom(shape = 21,
#                    size = 2,
#                    colour = "black",
#                    aes(fill = specific)) +
#   facet_wrap(vars(facet),
#              nrow = 1,
#              scales = "free_y") +
#   scale_fill_d3() +
#   labs(x = "Fraction & Age",
#        y = "Log2 Expression Count") +
#   theme(panel.grid.major = element_blank(),
#         legend.position = "none",
#         # axis.ticks.y = element_blank(),
#         axis.text.x = element_text(size = 8),
#         strip.text = element_text(face = 2)) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
#   panel_border()




# MB AGE NUMBERS ----

n_MB_AGE_DEG <- MB_AGE_META %>%
  filter(padj < 0.01) %>%
  nrow

n_MB_AGE_DEG_SPECIFIC_OPPOSING_TRANSLATED <- MB_AGE_DE_SPECIFIC_META %>%
  filter(specific == "Specific Opposing") %>%
  nrow

# MB age STRING interactions ------------------------------------

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
igr <- induced_subgraph(igr, components(igr)$membership %in% which(components(igr)$csize > 2))

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
layout <- layout_nicely(igr)

# MB AGE STRING Clusters ----
# plot
p_mb_age_deg_igraph <- as.grob(
  ~ plot(
    clusters,
    igr,
    layout = layout,
    edge.color = "grey70",
    vertex.label = NA,
    vertex.size = 3,
    vertex.color = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
  )
)

cluster_order <- membership(clusters)[order(membership(clusters))]

sapply(split(names(cluster_order), cluster_order), function(x){sapply(x, get_external_gene_name)})

# Trying to automate GO search and labelling is a nightmare
# Opting to label group with the gene of highest degree per group instead
# temp <- gost(query = split(names(cluster_order), cluster_order),
#      organism = "mmusculus",
#      multi_query = T,
#      custom_bg = NULL)
# temp$result %>%
#   mutate(significant = str_remove_all(significant, "c\\(|\\)")) %>%
#   separate(col = significant,
#            into = paste0("signif_", unique(cluster_order)),
#            sep = ", ") %>%
#   pivot_longer(starts_with("signif"),
#                names_to = "group",
#                values_to = "signif",
#                names_prefix = "signif_") %>%
#   mutate(signif = as.logical(signif)) %>%
#   filter(signif) %>%
#   arrange(group) %>%
#   filter(source == "REAC") %>% View

# MB AGE Clusters Heatmap ----

# get genes with highest degree
gene_labels <- sapply(split(degree(igr), membership(clusters)), function(x) names(sort(x, decreasing = T)[1]) )
# Convert to readable names
gene_labels <- sapply(gene_labels, get_external_gene_name)

# expand to the number of genes
gene_labels <- gene_labels[membership(clusters)]

# plot an ageing heatmap
plot_data <- t(scale(t(assay(vst(dds_C1_MB_IP))[match(names(membership(clusters)), rownames(dds_C1_MB_IP)),order(colData(dds_C1_MB_IP)$age)])))
plot_data_quantiles <- quantile(plot_data, probs = seq(0, 1, length.out = 10))

# make an annotation column df
annot_cols <- data.frame(Age = str_extract(colnames(plot_data), pattern = "(?<=MB_)[:alpha:]+"))
rownames(annot_cols) <- colnames(plot_data)
# make an annotation row df
annot_rows <- data.frame(Group = gene_labels)
# correct non-unique rownames
rownames(plot_data) <- make.unique(rownames(plot_data))
rownames(annot_rows) <- rownames(plot_data)

mat_colors <- list(Age = pal_d3()(10)[c(3, 2)],
                   Group = pal_d3(palette = "category20")(20)[c(1:length(unique(gene_labels)))])
names(mat_colors$Age) <- unique(annot_cols$Age)
names(mat_colors$Group) <- unique(gene_labels)

p_mb_age_deg_heatmap <- as.grob(
  ~ pheatmap(
    mat = plot_data[order(names(gene_labels)),],
    scale = "none",
    treeheight_row = 0,
    show_rownames = F,
    show_colnames = F,
    annotation_col = annot_cols,
    annotation_row = annot_rows,
    annotation_colors = mat_colors,
    # color = brewer.pal(n = length(plot_data_quantiles) - 1,
    #                    name = "PiYG"),
    color = colorspace::diverging_hcl(11,
                                      h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7)),
    # color = colorspace::diverging_hcl(11,
    #                                   h = c(180, 50), c = 80, l = c(20, 95), power = c(0.7, 1.3)),
    # color = brewer.pal(n = 8,
    #                    name = "RdBu"),
    breaks = plot_data_quantiles,
    border_color = NA,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = which(unlist(
      lapply(seq(2:length(names(gene_labels)[order(names(gene_labels))])), function(n)
        names(gene_labels)[order(names(gene_labels))][n] > names(gene_labels)[order(names(gene_labels))][n - 1])
    )),
    gaps_col = which(unlist(
      lapply(seq(2:ncol(plot_data)), function(n)
        colData(dds_C1_MB_IP)$age[n] != colData(dds_C1_MB_IP)$age[n-1])
    ))
    # cellwidth = 6,
    # cellheight = 2
  )
)
# ----
# ----
# ----
# MB AGE DE SPECIFIC SIMPLE VIEWING TABLE ----
MB_AGE_DE_SPECIFIC_META

# MB AGE CLUSTER TABLE ----

MB_AGE_TRAP_OUTCOMES <- MB_AGE_DE_SPECIFIC_META %>%
  select(ensembl_gene_id,
         status,
         specific,
         outcome_TRAP,
         pub_pd)

t_mb_age_cluster_membership <- enframe(membership(clusters), name = "ensembl_gene_id", value = "cluster") %>%
  mutate(external_gene_name = sapply(ensembl_gene_id, get_external_gene_name),
         cluster = as.character(cluster)) %>%
  left_join(
    enframe(gene_labels, name = "cluster", value = "lead") %>% distinct()) %>%
  left_join(MB_AGE_DEGREE) %>%
  left_join(MB_AGE_TRAP_OUTCOMES) %>%
  left_join(anno) %>%
  mutate(description = str_extract(description, ".+(?=\\[)"),
         description = ifelse(str_detect(description,
                                         "^[:upper:]"),
                              description,
                              str_to_sentence(description)),
         cluster = factor(paste0("Cluster ", cluster), levels = paste0("Cluster ", seq(1, length(unique(cluster)))))) %>%
  arrange(cluster, desc(degree))

# t_mb_age_cluster_membership %>%
#   mutate(GWAS = external_gene_name %in% GWAS_GENES_BROAD) %>%
#   filter(GWAS) %>% View


# ----
# ----
# ----
# AXON AGE DA MARKER DOWNREG EXAMPLE (PLOT) ----

p_axon_age_ds_vs <- view_counts(dds_C2_IP,
                                gene = "Slc6a3") %>%
  left_join(colData(dds_C2_IP),
            copy = T) %>%
  # filter(region != "MB") %>%
  filter(fraction == "IP") %>%
  group_by(region) %>%
  filter(!outlier_iqr(count, "high")) %>%
  ggplot(aes(x = age,
             y = count)) +
  geom_violin(colour = "black") +
  geom_quasirandom(shape = 21,
                   size = 3,
                   colour = "black",
                   aes(fill = region)) +
  facet_wrap(vars(region),
             scales = "free_y") +
  scale_fill_d3() +
  labs(x = "Region and Age",
       y = "Expression Count",
       title = "Slc6a3") +
  theme(legend.position = "none")


# AXON AGE ATG7 EXAMPLE (PLOT) ----

p_axon_age_ds_vs_atg7 <- view_counts(dds_C2_IP,
                                     gene = "Atg7") %>%
  left_join(colData(dds_C2_IP),
            copy = T) %>%
  # filter(region != "MB") %>%
  filter(fraction == "IP") %>%
  group_by(region) %>%
  filter(!outlier_iqr(count, "high")) %>%
  ggplot(aes(x = age,
             y = count)) +
  geom_violin(colour = "black") +
  geom_quasirandom(shape = 21,
                   size = 3,
                   colour = "black",
                   aes(fill = region)) +
  facet_wrap(vars(region),
             scales = "free_y") +
  scale_fill_d3() +
  labs(x = "Region and Age",
       y = "Expression Count",
       title = "Atg7") +
  theme(legend.position = "none")




# AXON AGE PANGLAO CELL TYPES ----

p_axon_age_panglao <- AXON_AGE_META %>%
  left_join(anno_human) %>%
  left_join(panglao,
            by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
  mutate(padj = ifelse(padj == 0, min(padj[padj > 0]), padj),
         score = -log10(padj) * log2FoldChange) %>%
  group_by(`cell type`) %>%
  summarise(score = mean(score),
            group_proportion = dplyr::n()/group_size,
            group_size = group_size,
            prop_score = score * group_proportion) %>%
  arrange(desc(prop_score)) %>%
  distinct() %>%
  filter(`cell type` %in% c("Astrocytes",
                            "Dopaminergic neurons",
                            "Microglia",
                            "Oligodendrocytes",
                            "Interneurons"
  )) %>%
  mutate(`cell type` = factor(`cell type`)) %>%
  mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
  ggplot(aes(x = prop_score,
             y = `cell type`,
             fill = `cell type`)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  labs(x = "Differential Expression\nScore",
       y = "Cell Type") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = 0,
             linetype = "dotted")

# AXON AGE GENES THAT GO UP ----



n_AXON_AGE_DE <- AXON_AGE_META %>% filter(padj < 0.01) %>% nrow

n_AXON_AGE_DE_UP_SOMA_ONLY <- AXON_AGE_META %>% filter(padj < 0.01 &
                                                         mb_translated &
                                                         !axon_enriched &
                                                         log2FoldChange > 0) %>%
  nrow()

n_AXON_AGE_DE_UP_SOMA_AND_AXON <- AXON_AGE_META %>% filter(padj < 0.01 &
                                                             mb_translated &
                                                             axon_enriched &
                                                             log2FoldChange > 0) %>%
  nrow()

n_AXON_AGE_DE_UP_SOMA_AND_AXON_MB_AGREE <- AXON_AGE_META %>% filter(padj < 0.01 &
                                                                      log2FoldChange > 0 &
                                                                      enrichment %in% c("Soma and Axon")) %>%
  filter(ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id) %>%
  nrow()


# n_AXON_AGE_DE_UP_DOPAMINERGIC <- AXON_AGE_META %>%
#   left_join(enrichment_status) %>%
#  filter(padj < 0.01 &
#           log2FoldChange > 0 &
#           enrichment != "Axon only") %>%
#  filter(mb_translated |
#           axon_translated) %>%
#  nrow
#
# n_AXON_AGE_DE_UP_DOPAMINERGIC_COMMON <- AXON_AGE_META %>%
#  filter(padj < 0.01 &
#           log2FoldChange > 0 &
#           confidence != "Axon TRAP and ABA") %>%
#  filter(mb_translated |
#           axon_enriched & aba_specific) %>%
#  filter(mb_age_relationship == "Common Upregulation") %>%
#  nrow


# AXON AGE SARM1 DS/VS ----

# DS
dds_C2_DS_IP_AGE <- dds_C2_AXON_IP[,colData(dds_C2_AXON_IP)$region == "DS"] %>%
  filter_zeros()
dds_C2_DS_IP_AGE <- dds_C2_DS_IP_AGE[dds_C2_AXON_IP_AGE_FILTER,]
dds_C2_DS_IP_AGE@design <- ~ collection + age
colData(dds_C2_DS_IP_AGE) <- droplevels(colData(dds_C2_DS_IP_AGE))
dds_C2_DS_IP_AGE <- DESeq(
  dds_C2_DS_IP_AGE,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)
res_C2_DS_IP_AGE <- DESeq2::results(
  dds_C2_DS_IP_AGE,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_DS_IP_AGE <- lfcShrink(
  dds_C2_DS_IP_AGE,
  res = res_C2_DS_IP_AGE,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)
DS_AGE_META <- as_tibble(res_C2_DS_IP_AGE, rownames = "ensembl_gene_id") %>%
  mutate(outcome_DS_AGE = ifelse(padj < 0.01,
                                 ifelse(log2FoldChange < 0,
                                        "Downregulated", "Upregulated"),
                                 "NS in DS")) %>%
  left_join(anno) %>%
  arrange(padj)
# select(ensembl_gene_id, outcome_DS_AGE)

# VS
dds_C2_VS_IP_AGE <- dds_C2_AXON_IP[,colData(dds_C2_AXON_IP)$region == "VS"] %>%
  filter_zeros()
dds_C2_VS_IP_AGE <- dds_C2_VS_IP_AGE[dds_C2_AXON_IP_AGE_FILTER,]
dds_C2_VS_IP_AGE@design <- ~ collection + age
colData(dds_C2_VS_IP_AGE) <- droplevels(colData(dds_C2_VS_IP_AGE))
dds_C2_VS_IP_AGE <- DESeq(
  dds_C2_VS_IP_AGE,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)
res_C2_VS_IP_AGE <- DESeq2::results(
  dds_C2_VS_IP_AGE,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_VS_IP_AGE <- lfcShrink(
  dds_C2_VS_IP_AGE,
  res = res_C2_VS_IP_AGE,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)
VS_AGE_META <- as_tibble(res_C2_VS_IP_AGE, rownames = "ensembl_gene_id") %>%
  mutate(outcome_VS_AGE = ifelse(padj < 0.01,
                                 ifelse(log2FoldChange < 0,
                                        "Downregulated", "Upregulated"),
                                 "NS in VS")) %>%
  left_join(anno) %>%
  arrange(padj)

# PLOT GENE TESTERS -----


# AGE IP TESTER
plotCounts(
  dds,
  get_ensembl_gene_id("Sarm1"),
  c("cohort", "compartment", "age", "region", "fraction"),
  returnData = T
) %>%
  filter(
    fraction == "IP"
    # cohort != "C3"
  ) %>%
  ggplot(aes(x = age, y = count)) +
  geom_violin() +
  geom_quasirandom(aes(colour = age)) +
  facet_grid(rows = vars(cohort), cols = vars(region), scales = "free_y") +
  panel_border() +
  labs(x = "Age and Region",
       y = "Expression Count",
       title = "micu1") +
  scale_color_d3() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none")

# ENRICHMENT TESTER
plotCounts(
  dds,
  get_ensembl_gene_id("Sarm1"),
  c("cohort", "compartment", "age", "region", "fraction"),
  returnData = T
) %>%
  # filter(fraction == "IP" &
  #          cohort != "C3") %>%
  ggplot(aes(x = fraction, y = count)) +
  geom_violin() +
  geom_quasirandom(aes(colour = fraction)) +
  facet_grid(rows = vars(cohort), cols = vars(region), scales = "free_y") +
  panel_border() +
  labs(x = "Age and Region",
       y = "Expression Count",
       title = "Sarm1") +
  scale_color_d3() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none")

# GENOTYPE IP TESTER
plotCounts(
  dds_C1_C2_C3ONT_MB_IP,
  get_ensembl_gene_id("Caln1"),
  c("cohort", "age", "genotype"),
  returnData = T
) %>%
  filter(
    # fraction == "IP" &
    cohort != "C2"
  ) %>%
  ggplot(aes(x = genotype, y = count)) +
  geom_violin() +
  geom_quasirandom(aes(colour = genotype)) +
  facet_grid(rows = vars(cohort), cols = vars(age), scales = "free_y") +
  panel_border() +
  labs(x = "Genotype and Region",
       y = "Expression Count",
       title = "Cpne6") +
  scale_color_d3() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "none")

# AXON AGE COPINE 6 ----

data <- plotCounts(
  dds,
  get_ensembl_gene_id("Cpne6"),
  c("cohort", "compartment", "age", "region", "fraction"),
  returnData = T
) %>%
  filter(fraction == "IP" &
           cohort != "C3") %>%
  mutate(count = log2(count + 1), 
         compartment = factor(ifelse(compartment == "MB", "Soma", "Axon"), 
                              levels = c("Soma", "Axon")))

p_axon_age_cpne6 <- ggboxplot(data, 
          x = "age", y = "count", fill = "age", 
          facet.by = "compartment", 
          outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  theme(legend.position = "none") +
  labs(title = "Cpne6", 
       x = "Age and Region",
       y = expression(Log[2] ~ Expression ~ Count)
       ) +
  stat_pvalue_manual(bind_rows({
    MB_AGE_META %>%
      filter(external_gene_name == "Cpne6") %>%
      select(external_gene_name, padj) %>%
      mutate(compartment = "Soma", 
             padj = signif(padj, 3))
  }, {
    AXON_AGE_META %>%
      filter(external_gene_name == "Cpne6") %>%
      select(external_gene_name, padj) %>%
      mutate(compartment = "Axon", 
             padj = signif(padj, 3))
  }) %>%
    mutate(group1 = "YOUNG", group2 = "OLD", 
           compartment = factor(compartment, levels = c("Soma", "Axon"))) %>%
    pval_asterisks(padj), 
  y.position = c(13.25, 13.7), 
  vjust = -0.25, 
  label = "ast",
  size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# p_axon_age_cpne6 <- plotCounts(
#   dds,
#   get_ensembl_gene_id("Cpne6"),
#   c("cohort", "compartment", "age", "region", "fraction"),
#   returnData = T
# ) %>%
#   filter(fraction == "IP" &
#            cohort != "C3") %>%
#   ggplot(aes(x = age, y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = age),
#                    shape = 21,
#                    size = 2) +
#   facet_grid(cols = vars(region), scales = "free_y") +
#   panel_border() +
#   labs(x = "Age and Region",
#        y = "Expression Count",
#        title = "Cpne6") +
#   scale_fill_d3() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 8))

# AXON AGE NMNAT2 ----

data <- plotCounts(
  dds,
  get_ensembl_gene_id("Nmnat2"),
  c("cohort", "compartment", "age", "region", "fraction"),
  returnData = T
) %>%
  filter(fraction == "IP" &
           cohort != "C3") %>%
  mutate(count = log2(count + 1), 
         compartment = factor(ifelse(compartment == "MB", "Soma", "Axon"), 
                              levels = c("Soma", "Axon")))

p_axon_age_nmnat2 <- ggboxplot(data, 
                              x = "age", y = "count", fill = "age", 
                              facet.by = "compartment", 
                              outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  theme(legend.position = "none") +
  labs(title = "Nmnat2", 
       x = "Age and Region",
       y = expression(Log[2] ~ Expression ~ Count)
  ) +
  stat_pvalue_manual(bind_rows({
    MB_AGE_META %>%
      filter(external_gene_name == "Nmnat2") %>%
      select(external_gene_name, padj) %>%
      mutate(compartment = "Soma", 
             padj = signif(padj, 3))
  }, {
    AXON_AGE_META %>%
      filter(external_gene_name == "Nmnat2") %>%
      select(external_gene_name, padj) %>%
      mutate(compartment = "Axon", 
             padj = signif(padj, 3))
  }) %>%
    mutate(group1 = "YOUNG", group2 = "OLD", 
           compartment = factor(compartment, levels = c("Soma", "Axon"))) %>%
    pval_asterisks(padj), 
  y.position = c(13.75, 12.9), 
  label = "ast",
  vjust = -0.25, 
  size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# p_axon_age_nmnat2 <- plotCounts(
#   dds,
#   get_ensembl_gene_id("Nmnat2"),
#   c("cohort", "compartment", "age", "region", "fraction"),
#   returnData = T
# ) %>%
#   filter(fraction == "IP" &
#            cohort != "C3") %>%
#   ggplot(aes(x = age, y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = age),
#                    shape = 21,
#                    size = 2) +
#   facet_grid(cols = vars(region), scales = "free_y") +
#   panel_border() +
#   labs(x = "Age and Region",
#        y = "Expression Count",
#        title = "Nmnat2") +
#   scale_fill_d3() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 8))


# p_axon_fraction_nmnat2 <- plotCounts(
#   dds,
#   get_ensembl_gene_id("Nmnat2"),
#   c("cohort", "compartment", "age", "region", "fraction"),
#   returnData = T
# ) %>%
#   filter(count > 1) %>%
#   # filter(fraction == "IP" &
#   #          cohort != "C3") %>%
#   ggplot(aes(x = fraction, y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = age),
#                    shape = 21,
#                    size = 2) +
#   facet_grid(cols = vars(region), scales = "free_y") +
#   panel_border() +
#   labs(x = "Fraction and Region",
#        y = "Expression Count",
#        title = "Sarm1") +
#   scale_fill_d3() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none")


# AXON AGE SARM1 ----

data <- plotCounts(
  dds,
  get_ensembl_gene_id("Sarm1"),
  c("cohort", "compartment", "age", "region", "fraction"),
  returnData = T
) %>%
  filter(fraction == "IP" &
           cohort != "C3") %>%
  mutate(count = log2(count + 1), 
         region = factor(ifelse(region == "MB", "Soma", 
                                ifelse(region == "DS", "Dorsal Axon", "Ventral Axon")), 
                              levels = c("Soma", "Dorsal Axon", "Ventral Axon")))

p_axon_age_sarm1 <- ggboxplot(data, 
                               x = "age", y = "count", fill = "age", 
                               facet.by = "region", 
                               outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  theme(legend.position = "none") +
  labs(title = "Sarm1", 
       x = "Age and Region",
       y = expression(Log[2] ~ Expression ~ Count)
  ) +
  scale_y_continuous(limits = c(8.5, 11.5), 
                     oob = scales::squish, 
                     expand = expansion(mult = c(0.05, 0.15))) +
  stat_pvalue_manual({
    bind_rows({
      MB_AGE_META %>%
        filter(external_gene_name == "Sarm1") %>%
        select(external_gene_name, padj) %>%
        mutate(region = "Soma", 
               padj = signif(padj, 3))
    }, {
      DS_AGE_META %>%
        filter(external_gene_name == "Sarm1") %>%
        select(external_gene_name, padj) %>%
        mutate(region = "Dorsal Axon", 
               padj = signif(padj, 3))
    }, {
      VS_AGE_META %>%
        filter(external_gene_name == "Sarm1") %>%
        select(external_gene_name, padj) %>%
        mutate(region = "Ventral Axon", 
               padj = signif(padj, 3))
    }) %>%
      mutate(group1 = "YOUNG", group2 = "OLD", 
             region = factor(region, levels = c("Soma", "Dorsal Axon", 
                                                "Ventral Axon")))
  } %>%
    pval_asterisks(padj), 
  y.position = c(11, 11, 11.25), 
  label = "ast",
  vjust = -0.25, 
  size = 3)

# p_axon_age_sarm1 <- plotCounts(
#   dds,
#   get_ensembl_gene_id("Sarm1"),
#   c("cohort", "compartment", "age", "region", "fraction"),
#   returnData = T
# ) %>%
#   filter(count > 1) %>%
#   filter(fraction == "IP" &
#            cohort != "C3") %>%
#   ggplot(aes(x = age, y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = age),
#                    shape = 21,
#                    size = 2) +
#   facet_grid(cols = vars(region), scales = "free_y") +
#   panel_border() +
#   labs(x = "Age and Region",
#        y = "Expression Count",
#        title = "Sarm1") +
#   scale_fill_d3() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 8))

# p_axon_fraction_sarm1 <- plotCounts(
#   dds,
#   get_ensembl_gene_id("Sarm1"),
#   c("cohort", "compartment", "age", "region", "fraction"),
#   returnData = T
# ) %>%
#   filter(count > 1) %>%
#   # filter(fraction == "IP" &
#   #          cohort != "C3") %>%
#   ggplot(aes(x = fraction, y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = age),
#                    shape = 21,
#                    size = 2) +
#   facet_grid(cols = vars(region), scales = "free_y") +
#   panel_border() +
#   labs(x = "Fraction and Region",
#        y = "Expression Count",
#        title = "Sarm1") +
#   scale_fill_d3() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 8))

# AXON AGE PTGES2 ----

data <- plotCounts(
  dds,
  get_ensembl_gene_id("Ptges2"),
  c("cohort", "compartment", "age", "region", "fraction"),
  returnData = T
) %>%
  filter(fraction == "IP" &
           cohort != "C3") %>%
  mutate(count = log2(count + 1), 
         compartment = factor(ifelse(compartment == "MB", "Soma", "Axon"), 
                              levels = c("Soma", "Axon")))

p_axon_age_ptges2 <- ggboxplot(data, 
                               x = "age", y = "count", fill = "age", 
                               facet.by = "compartment", 
                               outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  theme_PK(font_size = 11) +
  panel_border() +
  scale_fill_d3() +
  theme(legend.position = "none") +
  labs(title = "Ptges2", 
       x = "Age and Region",
       y = expression(Log[2] ~ Expression ~ Count)
  ) +
  stat_pvalue_manual(bind_rows({
    MB_AGE_META %>%
      filter(external_gene_name == "Ptges2") %>%
      select(external_gene_name, padj) %>%
      mutate(compartment = "Soma", 
             padj = signif(padj, 3))
  }, {
    AXON_AGE_META %>%
      filter(external_gene_name == "Ptges2") %>%
      select(external_gene_name, padj) %>%
      mutate(compartment = "Axon", 
             padj = signif(padj, 3))
  }) %>%
    mutate(group1 = "YOUNG", group2 = "OLD", 
           compartment = factor(compartment, levels = c("Soma", "Axon"))) %>%
    pval_asterisks(padj), 
  y.position = c(11.4, 11.75), 
  label = "ast",
  vjust = -0.25, 
  size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# p_axon_age_ptges2 <- plotCounts(
#   dds,
#   get_ensembl_gene_id("Ptges2"),
#   c("cohort", "compartment", "age", "region", "fraction"),
#   returnData = T
# ) %>%
#   filter(fraction == "IP" &
#            cohort != "C3") %>%
#   ggplot(aes(x = age, y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = age),
#                    shape = 21,
#                    size = 2) +
#   facet_grid(cols = vars(region), scales = "free_y") +
#   panel_border() +
#   labs(x = "Age and Region",
#        y = "Expression Count",
#        title = "Ptges2") +
#   scale_color_d3() +
#   scale_fill_d3() +
#   scale_y_continuous(limits = c(0, NA))+
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 8))




# ----
# ----
# ----
# GENOTYPE START HERE --------
# SNCA OVX: Transcripts ----

# get transcript external gene names
# anno_SNCA <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "ensembl_transcript_id",
#                  "external_transcript_name"),
#   filters = "external_gene_name",
#   values = "SNCA",
#   mart = ensembl_human,
#   uniqueRows = TRUE
# )

# anno_Snca <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "ensembl_transcript_id",
#                  "external_transcript_name"),
#   filters = "external_gene_name",
#   values = "Snca",
#   mart = ensembl,
#   uniqueRows = TRUE
# )

# anno_SNCA_joined <- bind_rows(anno_SNCA, anno_Snca)

anno_SNCA_joined <- read_delim("R/objects/snca_transcripts.txt", delim = "\t") %>%
  select(ensembl_gene_id = "Gene stable ID",
         ensembl_transcript_id = "Transcript stable ID",
         external_transcript_name = "Transcript name")

# TRAP transcripts
# formerly "/zfs/analysis/trap/active/testing_env/snca/input_lr"
files <-
  list.files(
    "input/data/snca",
    pattern = "quant.sf",
    recursive =  T,
    full.names = T
  )

all(file.exists(files)) # everything exists

names(files) <-
  str_extract(files, "(?<=input_lr\\/)[:alnum:]+(?=\\/quant.sf)")

snca_trap <- lapply(files, function(x){
  read_delim(x,
             delim = "\t") %>%
    mutate(ensembl_transcript_id = str_extract(Name, "[:alnum:]+[\\.][:digit:]+(?=\\|)")) %>%
    select(ensembl_transcript_id,
           everything(),
           -Name)
}) %>%
  bind_rows(.id = "barcode") %>%
  mutate(barcode = paste0("barcode", str_pad(barcode, 2, "0", side = "left"))) %>% 
  left_join(colData(dds_ONT),
            copy = T) %>%
  mutate(ensembl_transcript_id = substr(ensembl_transcript_id, 1, 15)) %>%
  filter(ensembl_transcript_id %in% anno_SNCA_joined$ensembl_transcript_id) %>%
  mutate(species = ifelse(str_detect(ensembl_transcript_id, "ENST"),
                          "SNCA", "Snca"),
         genotype = factor(ifelse(genotype == "WT", "WT", "SNCA-OVX"),
                           levels = c("WT", "SNCA-OVX")))



# set the transcript order
snca_order <- snca_trap %>%
  filter(barcode == "barcode01") %>%
  pull(ensembl_transcript_id)

# plot the OVX effect gene-level
p_ovx_gene_level <- snca_trap %>%
  group_by(sample_name, genotype, species) %>%
  summarise(reads = sum(TPM)) %>%
  ggplot(aes(x = genotype,
             y = reads,
             fill = species)) +
  geom_quasirandom(shape = 21,
                   size = 2) +
  facet_wrap(vars(species)) +
  panel_border() +
  scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
  labs(x = "Genotype",
       y = "Counts") +
  theme(legend.position = "none")

# p_ovx_transcript_level <- snca_trap %>%
#   group_by(species, sample_name) %>%
#   mutate(ensembl_transcript_id = factor(paste0("T", row_number()),
#                                         levels = paste0("T", row_number()))) %>%
# ggplot(aes(x = ensembl_transcript_id,
#              y = TPM,
#            fill = species)) +
#   geom_quasirandom(shape = 21,
#                    size = 2) +
#   facet_grid(rows = vars(genotype),
#              cols = vars(species),
#              scales = "free_x",
#              space = "free_x") +
#   panel_border() +
#   labs(x = "Transcript",
#        y = "Counts") +
#   scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
#   theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 90,
#                                    vjust = 0.5))

# compare with human counts from ips neurons

# SNCA_snca_transcripts <- c(anno_SNCA$ensembl_transcript_id,
#                            anno_Snca$ensembl_transcript_id)

SNCA_snca_transcripts <- anno_SNCA_joined$ensembl_transcript_id

# formerly "/zfs/analysis/rwm_transcriptomics/data/combined_kallisto/final_counts_20200207"
files <-
  list.files(
    "input/data/snca/ips",
    pattern = "abundance.tsv",
    recursive =  T,
    full.names = T
  )

names(files) <-
  str_extract(files, "(?<=ips\\/).+(?=\\/abundance.tsv)")

snca_ips <- lapply(files, function(x){
  read_delim(x,
             delim = "\t") %>%
    mutate(ensembl_transcript_id = substr(target_id, 1, 15)) %>%
    filter(ensembl_transcript_id %in% SNCA_snca_transcripts) %>%
    select(ensembl_transcript_id,
           everything(),
           -target_id)
}) %>%
  bind_rows(.id = "sample_name") %>%
  mutate(species = "SNCA")

# compare with human counts from ips neurons - nicheterwitz

# formerly  "/zfs/analysis/bm_rna/2021/nichterwitz"
files <-
  list.files(
    "input/data/snca/nichterwitz",
    pattern = "quant.sf",
    recursive =  T,
    full.names = T
  )

names(files) <-
  str_extract(files, "(?<=nichterwitz\\/).+(?=\\/quant.sf)")

snca_lcm <- lapply(files, function(x){
  read_delim(x,
             delim = "\t") %>%
    mutate(ensembl_transcript_id = str_extract(Name, "[:alnum:]+(?=[\\.][:digit:]+)")) %>%
    select(ensembl_transcript_id,
           everything(),
           -Name) %>%
    filter(ensembl_transcript_id %in% SNCA_snca_transcripts)
})  %>%
  bind_rows(.id = "sample_name") %>%
  mutate(species = "SNCA")

# compare with human counts from ips neurons - parkinnen
# formerly "/zfs/analysis/bm_rna/2021/parkinnen/"
files <-
  list.files(
    "input/data/snca/parkinnen/",
    pattern = "abundance.tsv",
    recursive =  T,
    full.names = T
  )

names(files) <-
  str_extract(files, "(?<=parkinnen\\/).+(?=\\/abundance.tsv)")

snca_lcm <- lapply(files, function(x){
  read_delim(x,
             delim = "\t") %>%
    mutate(ensembl_transcript_id = str_extract(target_id, "[:alnum:]+(?=[\\.][:digit:]+)")) %>%
    select(ensembl_transcript_id,
           everything(),
           -target_id) %>%
    filter(ensembl_transcript_id %in% SNCA_snca_transcripts)
})  %>%
  bind_rows(.id = "sample_name") %>%
  mutate(species = "SNCA")



# plot tpm per transcript
p_snca_wt_ovx_ips_lcm <- bind_rows(
  {
    snca_trap %>%
      select(ensembl_transcript_id,
             species,
             sample_name,
             genotype,
             counts = "TPM")
  }, {
    snca_ips[order(factor(snca_ips$ensembl_transcript_id, levels = snca_order)),] %>%
      mutate(species = "SNCA",
             genotype = "Human iPSC\nDopamine Neurons") %>%
      select(ensembl_transcript_id,
             species,
             sample_name,
             genotype,
             counts = "tpm")
  }, {
    snca_lcm[order(factor(snca_lcm$ensembl_transcript_id, levels = snca_order)),] %>%
      mutate(species = "SNCA",
             genotype = "Human LCM\nDopamine Neurons") %>%
      select(ensembl_transcript_id,
             species,
             sample_name,
             genotype,
             counts = "tpm")
  }
) %>%
  arrange(sample_name) %>%
  group_by(species, sample_name) %>%
  left_join(anno_SNCA_joined) %>%
  mutate(ensembl_transcript_id = factor(paste0("T", row_number()),
                                        levels = paste0("T", row_number()))) %>%
  mutate(genotype = factor(genotype,
                           levels = c("WT",
                                      "SNCA-OVX",
                                      "Human iPSC\nDopamine Neurons",
                                      "Human LCM\nDopamine Neurons"))) %>%
  # filter(counts > 20) %>%
  ggplot(aes(x = external_transcript_name,
             y = log2(counts+1),
             # y = counts,
             fill = species)) +
  geom_quasirandom(
    # shape = 21,
    size = 1
  ) +
  facet_grid(rows = vars(genotype),
             cols = vars(species),
             scales = "free_x",
             space = "free_x") +
  panel_border() +
  labs(x = "Transcript",
       y = expression(Log[2] ~ Expression ~ Count)) +
  scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5))
  # scale_y_continuous(limits = c(0, 500),
  #                    oob = scales::squish)




# Young OVX vs Young WT ----

# Cohort 1
dds_C1_MB_IP_YOUNG_GENOTYPE <- dds_C1_MB_IP[,colData(dds_C1_MB_IP)$age == "YOUNG"]
colData(dds_C1_MB_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_C1_MB_IP_YOUNG_GENOTYPE))
dds_C1_MB_IP_YOUNG_GENOTYPE@design <- ~ genotype
dds_C1_MB_IP_YOUNG_GENOTYPE <- DESeq(
  dds_C1_MB_IP_YOUNG_GENOTYPE,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_C1_MB_IP_YOUNG_GENOTYPE <- DESeq2::results(
  dds_C1_MB_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C1_MB_IP_YOUNG_GENOTYPE <- lfcShrink(
  dds_C1_MB_IP_YOUNG_GENOTYPE,
  res = res_C1_MB_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)

# Cohort 3 ONT
dds_ONT_MB_IP_YOUNG_GENOTYPE <- dds_ONT[,colData(dds_ONT)$age == "YOUNG"]
colData(dds_ONT_MB_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_YOUNG_GENOTYPE))
dds_ONT_MB_IP_YOUNG_GENOTYPE@design <- ~ genotype
dds_ONT_MB_IP_YOUNG_GENOTYPE <- DESeq(
  dds_ONT_MB_IP_YOUNG_GENOTYPE,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_ONT_MB_IP_YOUNG_GENOTYPE <- DESeq2::results(
  dds_ONT_MB_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_ONT_MB_IP_YOUNG_GENOTYPE <- lfcShrink(
  dds_ONT_MB_IP_YOUNG_GENOTYPE,
  res = res_ONT_MB_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)

# Meta-analyis
MB_GENOTYPE_YOUNG_META <-
  as_tibble(res_C1_MB_IP_YOUNG_GENOTYPE, rownames = "ensembl_gene_id") %>%
  left_join(
    as_tibble(res_ONT_MB_IP_YOUNG_GENOTYPE, rownames = "ensembl_gene_id"),
    by = "ensembl_gene_id",
    suffix = c("_C1", "_C3")
  ) %>%
  left_join(anno) %>%
  mutate(
    across(starts_with("pvalue"), ~ replace_na(.x, 1 - 1e-10)),
    # replace NA with *almost* 1 (it can't be perfectly 1)
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)),
    # replace 0 values with minimum non-zero of the column
    across(starts_with("pvalue"), ~ .x / 2),
    # divide 2-sided pvalue by 2 into 1-sided
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)) # again replace 0 values with minimum non-zero of the column
  ) %>%
  select(sort(names(.))) %>%
  select(
    ensembl_gene_id,
    external_gene_name,
    description,
    everything(),
    -starts_with("padj")
  ) %>%
  # filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
  drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
  rowwise() %>%
  mutate(
    pvalue_C3 = ifelse(
      sign(log2FoldChange_C3) == sign(log2FoldChange_C1),
      pvalue_C3,
      1 - pvalue_C3)
  ) %>%
  mutate(across(starts_with("pvalue"),
                ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
  mutate(
    # sumlog = sumlog(c_across(pvalue_C1:pvalue_C3))$p,
    # calculate sumlog
    sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
    # calculate stouffers
    log2FoldChange = mean(c(
      log2FoldChange_C1,
      log2FoldChange_C3
    ))
  ) %>%
  ungroup() %>%
  mutate(
    # correction for multiple comparisons
    sumz_adj = p.adjust(sumz, method = "fdr"),
    conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
  # mutate(score = -log10(sumz_adj) * abs(log2FoldChange)) %>%
  # arrange(desc(score), desc(log2FoldChange)) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_enriched = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
  group_by(external_gene_name) %>%
  filter(row_number() == 1) # Remove Snca, keep SNCA

# Old OVX vs Old WT ----
# Cohort 1
dds_C1_MB_IP_OLD_GENOTYPE <- dds_C1_MB_IP[,colData(dds_C1_MB_IP)$age == "OLD"]
colData(dds_C1_MB_IP_OLD_GENOTYPE) <- droplevels(colData(dds_C1_MB_IP_OLD_GENOTYPE))
dds_C1_MB_IP_OLD_GENOTYPE@design <- ~ genotype
dds_C1_MB_IP_OLD_GENOTYPE <- DESeq(
  dds_C1_MB_IP_OLD_GENOTYPE,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_C1_MB_IP_OLD_GENOTYPE <- DESeq2::results(
  dds_C1_MB_IP_OLD_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C1_MB_IP_OLD_GENOTYPE <- lfcShrink(
  dds_C1_MB_IP_OLD_GENOTYPE,
  res = res_C1_MB_IP_OLD_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)

# Cohort 3 ONT
dds_ONT_MB_IP_OLD_GENOTYPE <- dds_ONT[,colData(dds_ONT)$age == "OLD"]
colData(dds_ONT_MB_IP_OLD_GENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_OLD_GENOTYPE))
dds_ONT_MB_IP_OLD_GENOTYPE@design <- ~ genotype
dds_ONT_MB_IP_OLD_GENOTYPE <- DESeq(
  dds_ONT_MB_IP_OLD_GENOTYPE,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_ONT_MB_IP_OLD_GENOTYPE <- DESeq2::results(
  dds_ONT_MB_IP_OLD_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_ONT_MB_IP_OLD_GENOTYPE <- lfcShrink(
  dds_ONT_MB_IP_OLD_GENOTYPE,
  res = res_ONT_MB_IP_OLD_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)

# Meta-analyis
MB_GENOTYPE_OLD_META <-
  as_tibble(res_C1_MB_IP_OLD_GENOTYPE, rownames = "ensembl_gene_id") %>%
  left_join(
    as_tibble(res_ONT_MB_IP_OLD_GENOTYPE, rownames = "ensembl_gene_id"),
    by = "ensembl_gene_id",
    suffix = c("_C1", "_C3")
  ) %>%
  left_join(anno) %>%
  mutate(
    across(starts_with("pvalue"), ~ replace_na(.x, 1 - 1e-10)),
    # replace NA with *almost* 1 (it can't be perfectly 1)
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)),
    # replace 0 values with minimum non-zero of the column
    across(starts_with("pvalue"), ~ .x / 2),
    # divide 2-sided pvalue by 2 into 1-sided
    across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)) # again replace 0 values with minimum non-zero of the column
  ) %>%
  select(sort(names(.))) %>%
  select(
    ensembl_gene_id,
    external_gene_name,
    description,
    everything(),
    -starts_with("padj")
  ) %>%
  # filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
  drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
  rowwise() %>%
  mutate(
    pvalue_C3 = ifelse(
      sign(log2FoldChange_C3) == sign(log2FoldChange_C1),
      pvalue_C3,
      1 - pvalue_C3)
  ) %>%
  mutate(across(starts_with("pvalue"),
                ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
  mutate(
    # sumlog = sumlog(c_across(pvalue_C1:pvalue_C3))$p,
    # calculate sumlog
    sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
    # calculate stouffers
    log2FoldChange = mean(c(
      log2FoldChange_C1,
      log2FoldChange_C3
    ))
  ) %>%
  ungroup() %>%
  mutate(
    # correction for multiple comparisons
    sumz_adj = p.adjust(sumz, method = "fdr"),
    conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
  # mutate(score = -log10(sumz_adj) * abs(log2FoldChange)) %>%
  # arrange(desc(score), desc(log2FoldChange)) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
         axon_enriched = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
  group_by(external_gene_name) %>%
  filter(row_number() == 1) # Remove Snca, keep SNCA

# Numbers
n_MB_GENOTYPE_OLD_NGENES <- MB_GENOTYPE_OLD_META %>% filter(sumz_adj < 0.01) %>% nrow
n_MB_GENOTYPE_OLD_NGENES_FDR0.1 <- MB_GENOTYPE_OLD_META %>% filter(sumz_adj < 0.1) %>% nrow


# MA plot of genotype, by age ----
plot_data <- MB_GENOTYPE_OLD_META %>%
  mutate(signif = sumz_adj < 0.01) %>%
  arrange(signif) %>%
  mutate(age = "OLD") %>%
  bind_rows({
    MB_GENOTYPE_YOUNG_META %>%
      mutate(signif = sumz_adj < 0.1) %>%
      arrange(signif)  %>%
      mutate(age = "YOUNG")
  })

labels <- plot_data %>%
  filter(signif) %>%
  select(external_gene_name,
         baseMean_C1,
         log2FoldChange,
         mb_translated,
         age) %>%
  mutate(fontface = ifelse(mb_translated, "bold", "plain"),
         label_colour = mb_translated,
         external_gene_name = ifelse(external_gene_name == "Snca", "SNCA", external_gene_name)) %>%
  group_by(external_gene_name) %>%
  filter(row_number() == 1)

p_mb_genotype_ma <- plot_data %>%
  ggplot(aes(x = baseMean_C1,
             y = log2FoldChange)) +
  geom_point(aes(fill = signif),
             shape = 21,
             colour = "black") +
  geom_point(data = {plot_data %>% filter(signif)},
             size = 2,
             colour = pal_d3()(3)[3]) +
  geom_label_repel(data = labels,
                   aes(label = external_gene_name,
                       fontface = fontface,
                       colour = label_colour)) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-1, 1),
                     oob = scales::squish) +
  scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
  scale_color_manual(values = c("black", pal_d3()(3)[3])) +
  theme(legend.position = "none") +
  labs(x = "Expression Counts",
       y = "Log2 Fold Change"
  ) +
  facet_wrap(vars(factor(age, levels = c("YOUNG", "OLD"))),
             nrow = 1) +
  panel_border()

# Direction of change: OVX vs WT, by age ----
p_mb_genotype_old_enrichment <- plot_data %>%
  filter(!external_gene_name %in% c("Snca", "SNCA")) %>%
  left_join(enrichment_status) %>%
  filter(sumz_adj < 0.1) %>%
  mutate(outcome = ifelse(log2FoldChange > 0, "Increase", "Decrease")) %>%
  ggplot(aes(x = outcome)) +
  geom_bar(aes(fill = enrichment),
           colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "top") +
  geom_text(data = . %>%
              group_by(outcome) %>%
              tally(),
            aes(label = paste("Total:", n),
                y = n),
            vjust = -0.5) +
  labs(x = "Change in SNCA-OVX",
       y = "Number of Genes",
       fill = "Enrichment")
  # guides(fill = guide_legend(nrow = 2,
  #                     byrow = T,
  #                     title.position = "top",
  #                     title.hjust = 0.5))

# Old WT vs Young ALL ----

# subset samples
dds_ONT_MB_IP_AGE_WT <- dds_ONT[,colData(dds_ONT)$age == "YOUNG" |
                                         colData(dds_ONT)$age == "OLD" &
                                         colData(dds_ONT)$genotype == "WT"]
# process results
dds_ONT_MB_IP_AGE_WT@design <- ~ collection + age
colData(dds_ONT_MB_IP_AGE_WT) <- droplevels(colData(dds_ONT_MB_IP_AGE_WT))
dds_ONT_MB_IP_AGE_WT <- DESeq(
  dds_ONT_MB_IP_AGE_WT,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_ONT_MB_IP_AGE_WT <- DESeq2::results(
  dds_ONT_MB_IP_AGE_WT,
  contrast = c("age", "OLD", "YOUNG"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_ONT_MB_IP_AGE_WT <- lfcShrink(
  dds_ONT_MB_IP_AGE_WT,
  res = res_ONT_MB_IP_AGE_WT,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)

# Old OVX vs Young ALL ----

# subset samples
dds_ONT_MB_IP_AGE_OVX <- dds_ONT[,colData(dds_ONT)$age == "YOUNG" |
                                  colData(dds_ONT)$age == "OLD" &
                                  colData(dds_ONT)$genotype == "OVX"]
# process results
dds_ONT_MB_IP_AGE_OVX@design <- ~ collection + age
colData(dds_ONT_MB_IP_AGE_OVX) <- droplevels(colData(dds_ONT_MB_IP_AGE_OVX))
dds_ONT_MB_IP_AGE_OVX <- DESeq(
  dds_ONT_MB_IP_AGE_OVX,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_ONT_MB_IP_AGE_OVX <- DESeq2::results(
  dds_ONT_MB_IP_AGE_OVX,
  contrast = c("age", "OLD", "YOUNG"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_ONT_MB_IP_AGE_OVX <- lfcShrink(
  dds_ONT_MB_IP_AGE_OVX,
  res = res_ONT_MB_IP_AGE_OVX,
  contrast = c("age", "OLD", "YOUNG"),
  type = "ashr",
  parallel = T
)

# Direction of change: Old vs Young, by genotype ----

# combine the results of ageing from both genotypes
res_ONT_AGE_GENOTYPE_SPLIT <- bind_rows({
  res_ONT_MB_IP_AGE_WT %>%
    as_tibble(rownames = "ensembl_gene_id") %>%
    left_join(anno) %>%
    left_join(enrichment_status) %>%
    mutate(outcome = ifelse(log2FoldChange < 0, "Decrease", "Increase")) %>%
    mutate(genotype = "WT")
}, {
  res_ONT_MB_IP_AGE_OVX %>%
    as_tibble(rownames = "ensembl_gene_id") %>%
    left_join(anno) %>%
    left_join(enrichment_status) %>%
    mutate(outcome = ifelse(log2FoldChange < 0, "Decrease", "Increase")) %>%
    mutate(genotype = "OVX")
}) %>% drop_na

p_age_genotype_outcome <- res_ONT_AGE_GENOTYPE_SPLIT %>%
  filter(padj < 0.01 & external_gene_name != "Snca") %>%
  ggplot(aes(x = outcome)) +
  geom_bar(aes(fill = enrichment),
           colour = "black") +
  scale_fill_d3() +
  facet_wrap(vars(genotype)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme(legend.position = c(0.8, 0.8)) +
  geom_text(data = . %>%
              group_by(outcome, genotype) %>%
              tally(),
            aes(label = paste("Total:", n),
                y = n),
            vjust = -0.5) +
  labs(x = "Change in Age",
       y = "Number of Genes",
       fill = "Enrichment") +
  panel_border()

# Magnitude of change: Old vs Young, by genotype ----

p_age_genotype_lfc <- res_ONT_AGE_GENOTYPE_SPLIT %>%
  filter(ensembl_gene_id %in% {res_ONT_MB_IP_AGE_WT %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      left_join(anno) %>%
      filter(padj < 0.1) %>%
      pull(ensembl_gene_id)}) %>%
  group_by(ensembl_gene_id) %>%
  filter(dplyr::n() == 2) %>%
  mutate(score = -log10(padj) * log2FoldChange,
         genotype = ifelse(genotype == "OVX", "SNCA-OVX", "WT")) %>%
  ggplot(aes(x = genotype,
             y = log2FoldChange)) +
  geom_quasirandom(alpha = 0.5) +
  geom_hline(yintercept = 0,
             linetype = "dotted") +
  labs(x = "Genotype",
       y = "Log2 Fold Change in Age") +
  facet_wrap(vars(enrichment),
             nrow = 1) +
  panel_border()


res_ONT_AGE_GENOTYPE_SPLIT %>%
  filter(ensembl_gene_id %in% {res_ONT_MB_IP_AGE_WT %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      left_join(anno) %>%
      filter(padj < 0.1) %>%
      pull(ensembl_gene_id)}) %>%
  group_by(ensembl_gene_id) %>%
  filter(dplyr::n() == 2) %>%
  mutate(score = -log10(padj) * log2FoldChange,
         genotype = ifelse(genotype == "OVX", "SNCA-OVX", "WT")) %>%
  ggplot(aes(x = genotype,
             y = log2FoldChange)) +
  geom_quasirandom(alpha = 0.5) +
  geom_hline(yintercept = 0,
             linetype = "dotted") +
  labs(x = "Genotype",
       y = "Log2 Fold Change in Age") +
  facet_wrap(vars(enrichment),
             nrow = 1) +
  panel_border()

# Zfp672 ----

p_genotype_zfp672 <- bind_rows({
  plotCounts(dds_C1_MB_IP,
             get_ensembl_gene_id("Zfp672"),
             c("age", "genotype"),
             returnData = T) %>% mutate(cohort = "C1")
},
{
  plotCounts(dds_ONT,
             get_ensembl_gene_id("Zfp672"),
             c("age", "genotype"),
             returnData = T) %>% mutate(cohort = "C3 ONT")
}) %>%
  unite(c(age, genotype), col = age_genotype, sep = "\n") %>%
  mutate(age_genotype = factor(age_genotype,
                               levels = c("YOUNG\nWT",
                                          "YOUNG\nOVX",
                                          "OLD\nWT",
                                          "OLD\nOVX"))) %>%
  ggplot(aes(x = age_genotype,
             y = count,
             colour = age_genotype)) +
  geom_quasirandom() +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(vars(cohort),
             scales = "free",
             nrow = 1) +
  scale_color_d3() +
  panel_border() +
  labs(x = "Age & Genotype",
       y = "Expression Count",
       title = "Zfp672") +
  theme(legend.position = "none")

# OVX vs WT: Age independent ----
colData(dds_ONT) <- droplevels(colData(dds_ONT))
dds_ONT@design <- ~ collection + age + genotype
dds_ONT <- DESeq(
  dds_ONT,
  minReplicatesForReplace = Inf,
  parallel = T
)
res_ONT_MB_IP_GENOTYPE <- DESeq2::results(
  dds_ONT,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_ONT_MB_IP_GENOTYPE <- lfcShrink(
  dds_ONT,
  res = res_ONT_MB_IP_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)

# OVX vs WT: Axon age independent ----

# dds_C3_AXON_IP_GENOTYPE <- dds_C3_AXON_IP
# dds_C3_AXON_IP_GENOTYPE <- dds_C3_AXON_IP_GENOTYPE[filter_genes(dds_C3_AXON_IP_GENOTYPE, grouping = c("genotype", "age", "gene_id")), ]
# colData(dds_C3_AXON_IP_GENOTYPE) <- droplevels(colData(dds_C3_AXON_IP_GENOTYPE))
# dds_C3_AXON_IP_GENOTYPE@design <- ~ collection + age + region + genotype
# dds_C3_AXON_IP_GENOTYPE <- DESeq(
#   dds_C3_AXON_IP_GENOTYPE,
#   sfType = "poscounts",
#   fitType = "local",
#   minReplicatesForReplace = Inf,
#   parallel = T
# )
# res_C3_AXON_IP_GENOTYPE <- DESeq2::results(
#   dds_C3_AXON_IP_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C3_AXON_IP_GENOTYPE <- lfcShrink(
#   dds_C3_AXON_IP_GENOTYPE,
#   res = res_C3_AXON_IP_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   type = "ashr",
#   parallel = T
# )

# Old OVX vs Old WT: Axon ----

# dds_C3_AXON_IP_OLD_GENOTYPE <- dds_C3_AXON_IP[,colData(dds_C3_AXON_IP)$age == "OLD"]
# dds_C3_AXON_IP_OLD_GENOTYPE <- dds_C3_AXON_IP_OLD_GENOTYPE[filter_genes(dds_C3_AXON_IP_OLD_GENOTYPE, grouping = c("genotype", "gene_id")), ]
# colData(dds_C3_AXON_IP_OLD_GENOTYPE) <- droplevels(colData(dds_C3_AXON_IP_OLD_GENOTYPE))
# dds_C3_AXON_IP_OLD_GENOTYPE@design <- ~ collection + genotype
# dds_C3_AXON_IP_OLD_GENOTYPE <- DESeq(
#   dds_C3_AXON_IP_OLD_GENOTYPE,
#   sfType = "poscounts",
#   fitType = "local",
#   minReplicatesForReplace = Inf,
#   parallel = T
# )
# res_C3_AXON_IP_OLD_GENOTYPE <- DESeq2::results(
#   dds_C3_AXON_IP_OLD_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C3_AXON_IP_OLD_GENOTYPE <- lfcShrink(
#   dds_C3_AXON_IP_OLD_GENOTYPE,
#   res = res_C3_AXON_IP_OLD_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   type = "ashr",
#   parallel = T
# )

# Young OVX vs Young WT: Axon ----

dds_C3_AXON_IP_YOUNG_GENOTYPE <- dds_C3_AXON_IP[,colData(dds_C3_AXON_IP)$age == "YOUNG"]
dds_C3_AXON_IP_YOUNG_GENOTYPE <- dds_C3_AXON_IP_YOUNG_GENOTYPE[filter_genes(dds_C3_AXON_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id")), ]
colData(dds_C3_AXON_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_C3_AXON_IP_YOUNG_GENOTYPE))
dds_C3_AXON_IP_YOUNG_GENOTYPE@design <- ~ collection + genotype
dds_C3_AXON_IP_YOUNG_GENOTYPE <- DESeq(
  dds_C3_AXON_IP_YOUNG_GENOTYPE,
  sfType = "poscounts",
  fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = T
)
res_C3_AXON_IP_YOUNG_GENOTYPE <- DESeq2::results(
  dds_C3_AXON_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_AXON_IP_YOUNG_GENOTYPE <- lfcShrink(
  dds_C3_AXON_IP_YOUNG_GENOTYPE,
  res = res_C3_AXON_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)

# Young OVX vs Young WT: DS ----

dds_C3_DS_IP_YOUNG_GENOTYPE <- dds_C3_AXON_IP[,colData(dds_C3_AXON_IP)$age == "YOUNG" & colData(dds_C3_AXON_IP)$region == "DS"]
filter_C3_DS_IP_YOUNG_GENOTYPE <- filter_genes(dds_C3_DS_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id"), prop = 5/5)
dds_C3_DS_IP_YOUNG_GENOTYPE <- dds_C3_DS_IP_YOUNG_GENOTYPE[filter_C3_DS_IP_YOUNG_GENOTYPE,]
dds_C3_DS_IP_YOUNG_GENOTYPE <- dds_C3_DS_IP_YOUNG_GENOTYPE[filter_genes(dds_C3_DS_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id")), ]
colData(dds_C3_DS_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_C3_DS_IP_YOUNG_GENOTYPE))
dds_C3_DS_IP_YOUNG_GENOTYPE@design <- ~ collection + genotype
dds_C3_DS_IP_YOUNG_GENOTYPE <- DESeq(
  dds_C3_DS_IP_YOUNG_GENOTYPE,
  sfType = "poscounts",
  fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = T
)
res_C3_DS_IP_YOUNG_GENOTYPE <- DESeq2::results(
  dds_C3_DS_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_DS_IP_YOUNG_GENOTYPE <- lfcShrink(
  dds_C3_DS_IP_YOUNG_GENOTYPE,
  res = res_C3_DS_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)
summary(res_C3_DS_IP_YOUNG_GENOTYPE, alpha = 0.01)



# res_C3_DS_IP_YOUNG_GENOTYPE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(padj < 0.01) %>%
#   left_join(anno) %>%
#   left_join(enrichment_status) %>%
#   ggplot(aes(x = log2FoldChange > 0,
#              fill = enrichment)) +
#   geom_bar()



# plotCounts(dds_C3_DS_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Cav2"), "genotype")



# Young OVX vs Young WT: VS ----

dds_C3_VS_IP_YOUNG_GENOTYPE <- dds_C3_AXON_IP[,colData(dds_C3_AXON_IP)$age == "YOUNG" & colData(dds_C3_AXON_IP)$region == "VS"]
filter_C3_VS_IP_YOUNG_GENOTYPE <- filter_genes(dds_C3_VS_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id"), prop = 5/5)
dds_C3_VS_IP_YOUNG_GENOTYPE <- dds_C3_VS_IP_YOUNG_GENOTYPE[filter_C3_VS_IP_YOUNG_GENOTYPE,]
dds_C3_VS_IP_YOUNG_GENOTYPE <- dds_C3_VS_IP_YOUNG_GENOTYPE[filter_genes(dds_C3_VS_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id")), ]
colData(dds_C3_VS_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_C3_VS_IP_YOUNG_GENOTYPE))
dds_C3_VS_IP_YOUNG_GENOTYPE@design <- ~ collection + genotype
dds_C3_VS_IP_YOUNG_GENOTYPE <- DESeq(
  dds_C3_VS_IP_YOUNG_GENOTYPE,
  sfType = "poscounts",
  fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = T
)
res_C3_VS_IP_YOUNG_GENOTYPE <- DESeq2::results(
  dds_C3_VS_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_VS_IP_YOUNG_GENOTYPE <- lfcShrink(
  dds_C3_VS_IP_YOUNG_GENOTYPE,
  res = res_C3_VS_IP_YOUNG_GENOTYPE,
  contrast = c("genotype", "OVX", "WT"),
  type = "ashr",
  parallel = T
)
summary(res_C3_VS_IP_YOUNG_GENOTYPE, alpha = 0.01)



res_C3_VS_IP_YOUNG_GENOTYPE %>%
  as_tibble(rownames = "ensembl_gene_id") %>%
  filter(padj < 0.01) %>%
  left_join(anno) %>%
  left_join(enrichment_status) %>%
  ggplot(aes(x = log2FoldChange > 0,
             fill = enrichment)) +
  geom_bar()

# res_C3_VS_IP_YOUNG_GENOTYPE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   left_join(enrichment_status) %>%
#   left_join(publications_all_genes) %>%
#   left_join(publications_all_genes_axon) %>%
#   filter(padj < 0.01) %>%
#   filter(str_detect(enrichment, "oma")) %>%
#   View

# Young OVX DS vs Young WT ALL ----

dds_C3_DSOVX_IP_YOUNG_GENOTYPE <- dds_C3_AXON_IP[,colData(dds_C3_AXON_IP)$age == "YOUNG"]
filter_C3_DSOVX_IP_YOUNG_GENOTYPE <- filter_genes(dds_C3_DSOVX_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id"), prop = 5/5)
dds_C3_DSOVX_IP_YOUNG_GENOTYPE <- dds_C3_DSOVX_IP_YOUNG_GENOTYPE[filter_C3_DSOVX_IP_YOUNG_GENOTYPE,]
dds_C3_DSOVX_IP_YOUNG_GENOTYPE <- dds_C3_DSOVX_IP_YOUNG_GENOTYPE[filter_genes(dds_C3_DSOVX_IP_YOUNG_GENOTYPE, grouping = c("genotype", "gene_id")), ]
colData(dds_C3_DSOVX_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_C3_DSOVX_IP_YOUNG_GENOTYPE))
colData(dds_C3_DSOVX_IP_YOUNG_GENOTYPE)$group <- factor(ifelse(colData(dds_C3_DSOVX_IP_YOUNG_GENOTYPE)$region == "DS" &
                                                          colData(dds_C3_DSOVX_IP_YOUNG_GENOTYPE)$genotype == "OVX",
                                                        "DS_OVX", "Other"), levels = c("Other", "DS_OVX"))
dds_C3_DSOVX_IP_YOUNG_GENOTYPE@design <- ~ collection + region + group
dds_C3_DSOVX_IP_YOUNG_GENOTYPE <- DESeq(
  dds_C3_DSOVX_IP_YOUNG_GENOTYPE,
  sfType = "poscounts",
  fitType = "local",
  minReplicatesForReplace = Inf,
  # parallel = T
)
resultsNames(dds_C3_DSOVX_IP_YOUNG_GENOTYPE)
res_C3_DSOVX_IP_YOUNG_GENOTYPE <- DESeq2::results(
  dds_C3_DSOVX_IP_YOUNG_GENOTYPE,
  contrast = c("group", "DS_OVX", "Other"),
  cooksCutoff = Inf,
  filterFun = ihw,
  # parallel = T
)
res_C3_DSOVX_IP_YOUNG_GENOTYPE <- lfcShrink(
  dds_C3_DSOVX_IP_YOUNG_GENOTYPE,
  res = res_C3_DSOVX_IP_YOUNG_GENOTYPE,
  contrast = c("group", "DS_OVX", "Other"),
  type = "ashr",
  # parallel = T
)
summary(res_C3_DSOVX_IP_YOUNG_GENOTYPE, alpha = 0.01)


# res_C3_DSOVX_IP_YOUNG_GENOTYPE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   left_join(enrichment_status) %>%
#   left_join(publications_all_genes) %>%
#   left_join(publications_all_genes_axon) %>%
#   filter(padj < 0.01) %>%
#   filter(str_detect(enrichment, "oma")) %>%
#   filter(abs(log2FoldChange) < 2)

  plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Rbl2"), c("region", "genotype"))

# Plot Young OVX vs Young WT: DS DEGs ----

# res_C3_DS_IP_YOUNG_GENOTYPE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   left_join(enrichment_status) %>%
#   left_join(publications_all_genes) %>%
#   left_join(publications_all_genes_axon) %>%
#   filter(padj < 0.01) %>%
#   filter(str_detect(enrichment, "oma")) %>%
#   filter(abs(log2FoldChange) < 2)

  
data <- bind_rows(
  {
    plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Mapre2"), c("region", "genotype"),
               returnData = T) %>%
      mutate(external_gene_name = "Mapre2")
  }, {
    plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Csnk1a1"), c("region", "genotype"),
               returnData = T) %>%
      mutate(external_gene_name = "Csnk1a1")
  }, {
    plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Dusp8"), c("region", "genotype"),
               returnData = T) %>%
      mutate(external_gene_name = "Dusp8")
  }
) %>%
  unite(c(region, genotype), col = region_genotype, sep = "\n", remove = F) %>%
  mutate(region_genotype = factor(region_genotype,
                                  levels = c("VS\nWT",
                                             "VS\nOVX",
                                             "DS\nWT",
                                             "DS\nOVX")), 
         count = log2(count + 1))

stat.test <- bind_rows(
  {
    res_C3_VS_IP_YOUNG_GENOTYPE %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      mutate(region = "VS")
  }, {
    res_C3_DS_IP_YOUNG_GENOTYPE %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      mutate(region = "DS")
  }
) %>% 
  left_join(anno) %>%
  filter(external_gene_name %in% c("Csnk1a1", "Dusp8", "Mapre2")) %>%
  mutate(group1 = "WT", group2 = "OVX", 
         padj = signif(padj, 3), 
         external_gene_name = factor(external_gene_name, 
                                     levels = c("Csnk1a1", 
                                                "Dusp8", 
                                                "Mapre2")), 
         region = factor(region, levels = c("VS", "DS"))) %>%
  pval_asterisks(padj)
  
p_genotype_ds <- ggboxplot(data, 
            x = "genotype", 
            y = "count", 
            fill = "genotype", 
            outlier.shape = NA, 
            facet.by = c("external_gene_name", "region"), 
            scales = "free", 
            nrow = 3) +
    geom_quasirandom(shape = 21, size = 1) +
    scale_fill_d3() +
    theme_PK(font_size = 11) +
    panel_border() +
    theme(legend.position = "none") +
    labs(x = "Genotype & Region", 
         y = expression(Log[2] ~ Expression ~ Count)) +
    stat_pvalue_manual(stat.test, 
                       y.position = c(13.5,13.5,13.65,13.5,13.75,14.25), 
                       label = "ast", 
                       vjust = -0.25, 
                       size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))

# p_genotype_ds <- bind_rows(
#   {
#     plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Mapre2"), c("region", "genotype"),
#                returnData = T) %>%
#       mutate(external_gene_name = "Mapre2")
#   }, {
#     plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Csnk1a1"), c("region", "genotype"),
#                returnData = T) %>%
#       mutate(external_gene_name = "Csnk1a1")
#   }, {
#     plotCounts(dds_C3_AXON_IP_YOUNG_GENOTYPE, get_ensembl_gene_id("Dusp8"), c("region", "genotype"),
#                returnData = T) %>%
#       mutate(external_gene_name = "Dusp8")
#   }
# ) %>%
#   unite(c(region, genotype), col = region_genotype, sep = "\n") %>%
#   mutate(region_genotype = factor(region_genotype,
#                                levels = c("VS\nWT",
#                                           "VS\nOVX",
#                                           "DS\nWT",
#                                           "DS\nOVX"))) %>%
#   ggplot(aes(x = region_genotype,
#              y = count,
#              colour = region_genotype)) +
#   geom_quasirandom() +
#   scale_y_continuous(limits = c(0, NA),
#                      expand = expansion(mult = c(0.05, 0.1))) +
#   facet_wrap(vars(external_gene_name),
#              scales = "free",
#              nrow = 3) +
#   scale_color_d3() +
#   panel_border() +
#   labs(x = "Region & Genotype",
#        y = "Expression Count") +
#   theme(legend.position = "none")






c("Mapre2",
  "Csnk1a1",
  "Rpn2",
  "Dusp8")

# ZFP672

data <- bind_rows({
  plotCounts(dds_C1_MB_IP,
             get_ensembl_gene_id("Zfp672"),
             c("age", "genotype"),
             returnData = T) %>% mutate(cohort = "C1")
},
{
  plotCounts(dds_ONT,
             get_ensembl_gene_id("Zfp672"),
             c("age", "genotype"),
             returnData = T) %>% mutate(cohort = "C3 ONT")
}) %>%
  unite(c(age, genotype), col = age_genotype, sep = "\n") %>%
  mutate(age_genotype = factor(age_genotype,
                               levels = c("YOUNG\nWT",
                                          "YOUNG\nOVX",
                                          "OLD\nWT",
                                          "OLD\nOVX")), 
         count = log2(count + 1))
  # filter(cohort )

stat.test <- bind_rows(
  {
    MB_GENOTYPE_YOUNG_META %>%
      filter(external_gene_name == "Zfp672") %>%
      select(external_gene_name, pvalue_C1, pvalue_C3, sumz_adj) %>%
      pivot_longer(starts_with("pvalue"), 
                   names_to = "cohort", 
                   values_to = "pvalue", 
                   names_prefix = "pvalue_") %>%
      mutate(group1 = "YOUNG\nWT", group2 = "YOUNG\nOVX")
  }, {
    MB_GENOTYPE_OLD_META %>%
      filter(external_gene_name == "Zfp672") %>%
      select(external_gene_name, pvalue_C1, pvalue_C3, sumz_adj) %>%
      pivot_longer(starts_with("pvalue"), 
                   names_to = "cohort", 
                   values_to = "pvalue", 
                   names_prefix = "pvalue_") %>%
      mutate(group1 = "OLD\nWT", group2 = "OLD\nOVX")
  }
) %>%
  mutate(cohort = ifelse(cohort == "C3", "C3 ONT", "C1"), 
         pvalue = signif(pvalue, 3)) %>%
  pval_asterisks(pvalue)

sumz_values <- stat.test %>%
  pull(sumz_adj) %>% 
  unique %>%
  signif(3)

p_genotype_zfp672 <- ggboxplot(data, 
          x = "age_genotype", y = "count", fill = "age_genotype", 
          facet.by = "cohort", scales = "free", outlier.shape = NA) +
  geom_quasirandom(shape = 21, size = 1) +
  scale_fill_d3() +
  labs(x = "Age & Genotype", 
       y = expression(Log[2] ~ Expression ~ Count), 
       title = "Zfp672", 
       caption = paste0("Stouffer's *P*: Young Genotype = ", sumz_values[1], ", Old Genotype = ", sumz_values[2])) +
  stat_pvalue_manual(stat.test, 
                     label = "ast", 
                     y.position = c(10, 5.3, 10, 5.3), 
                     size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_PK(font_size = 11) +
  panel_border() +
  theme(legend.position = "none", 
        plot.caption = ggtext::element_markdown())


#   bind_rows({
#   plotCounts(dds_C1_MB_IP,
#              get_ensembl_gene_id("Zfp672"),
#              c("age", "genotype"),
#              returnData = T) %>% mutate(cohort = "C1")
# },
# {
#   plotCounts(dds_ONT,
#              get_ensembl_gene_id("Zfp672"),
#              c("age", "genotype"),
#              returnData = T) %>% mutate(cohort = "C3 ONT")
# }) %>%
#   unite(c(age, genotype), col = age_genotype, sep = "\n") %>%
#   mutate(age_genotype = factor(age_genotype,
#                                levels = c("YOUNG\nWT",
#                                           "YOUNG\nOVX",
#                                           "OLD\nWT",
#                                           "OLD\nOVX"))) %>%
#   ggplot(aes(x = age_genotype,
#              y = count,
#              colour = age_genotype)) +
#   geom_quasirandom() +
#   scale_y_continuous(limits = c(0, NA),
#                      expand = expansion(mult = c(0.05, 0.1))) +
#   facet_wrap(vars(cohort),
#              scales = "free",
#              nrow = 1) +
#   scale_color_d3() +
#   panel_border() +
#   labs(x = "Age & Genotype",
#        y = "Expression Count",
#        title = "Zfp672") +
#   theme(legend.position = "none")
# Testing ----





# bind_rows(
#   {
#     MB_AGE_META %>%
#       left_join(enrichment_status) %>%
#       filter(str_detect(enrichment, "oma"),
#              padj < 0.1) %>%
#       mutate(genotype = "WT")
#   }, {
#     res_ONT_MB_IP_OLD_GENOTYPE %>%
#       as_tibble(rownames = "ensembl_gene_id") %>%
#       left_join(anno) %>%
#       left_join(enrichment_status) %>%
#       filter(padj < 0.1 & log2FoldChange > 0 & str_detect(enrichment, "oma")) %>%
#       mutate(genotype = "OVX")
#   }
# ) %>%
#   select(external_gene_name,
#          log2FoldChange,
#          padj,
#          genotype) %>%
#   arrange(external_gene_name) %>%
#   group_by(external_gene_name) %>%
#   filter(sum(genotype == "OVX") == 1) %>% View
#   filter(dplyr::n() == 1) %>% View
# pivot_wider(values_from = score,
#               names_from = genotype) %>%
#   unnest(c(WT, OVX)) %>%
#   # filter(WT != 0 & OVX != 0) %>%
#   filter(sign(WT) != sign(OVX)) %>%
#   arrange(desc(OVX)) %>%
#   ggplot(aes(x = WT, y = OVX)) +
#   geom_point() +
#   coord_cartesian(xlim = c(-1, 1),
#                   ylim = c(-1, 1))








# ----
# ----
# ----

# Remove large objects ----
# sort( sapply(ls(),function(x){object.size(get(x))}))

rm(talon_qc)
rm(talon_annot)
rm(txi)

# Save output ----
save(list = c(ls(pattern = "^p_"),
              ls(pattern = "^t_"),
              ls(pattern = "^n_")),
     file = "R/objects/04-05-results-objects.RData")

save.image("R/objects/RESULTS_2_and_3.RData")

# load("R/objects/03-results-code.Rdata")
