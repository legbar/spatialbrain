getwd()

# SETUP ----
# set seed for reproducibility
set.seed(821196)
# load packages and save package versions
source("R/renv.R")
# load sample metadata
source("R/metadata.R")
# load outlier tracking functions
source("R/outlier_tracking.R")
# load annotation function
source("R/anno.R")
# PCA function
source("R/pca.R")
# Pvalue asterisk function
source("R/pval_asterisks.R")
# remove underscores from titles function for plotting
remove_underscores <- function(colnames) {
  gsub("_", " ", colnames)
}
# fix function associations
rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
n <- dplyr::n
# set theme
source("R/theme_set.R")
# ----
# ----
# ----
# RESULTS 1 ----
# ----
# ----
# ----
# eGFP ELISA --------------------------------------------------------------
egfp_elisa <-
  read_delim("input/experiments/TRAP_GFP_RUN1_RUN2.csv", delim = ",") %>%
  pivot_longer(everything(),
               names_to = "region",
               values_to = "concentration") %>%
  filter(!is.na(concentration)) %>%
  mutate(genotype = ifelse(region == "MB fs/fs", "RELfs/fs", "RELfs/-::DATcre/-"),
         region = ifelse(region == "MB fs/fs", "MB", region)) %>%
  group_by(region) %>%
  filter(!outlier_iqr(concentration)) %>% 
  # filter(region != "CTX") %>% 
  mutate(region = factor(region, levels = c("MB", "DS", "VS", "CTX", "MB fs/fs")))

egfp_elisa_comparisons <- egfp_elisa %>% 
  filter(genotype != "RELfs/fs")

# perform one-way ANOVA
egfp_elisa_comparisons_aov <- aov(log10(concentration) ~ region, data = egfp_elisa_comparisons)
# summarise the results of one-way ANOVA
summary(egfp_elisa_comparisons_aov)
# plot residuals vs fitted line to check for homogeneity of variances
plot(egfp_elisa_comparisons_aov, 1)
# test Levene's test as an alternative: http://www.sthda.com/english/wiki/one-way-anova-test-in-r?utm_content=buffer3431c&utm_medium=social&utm_source=facebook.com&utm_campaign=buffer#check-anova-assumptions-test-validity
library(car)
# the p-value is larger than 0.05, so we can assume homogeneity of variances
leveneTest(log10(concentration) ~ region, data = egfp_elisa_comparisons)  
# if this assumption could not be met, we could use Welch's ANOVA
oneway.test(log10(concentration) ~ region, data = egfp_elisa_comparisons)
# or pairwise t-tests with corrections for mult comparisons
pairwise.t.test(log10(egfp_elisa_comparisons$concentration), 
                egfp_elisa_comparisons$region,
                p.adjust.method = "BH", pool.sd = FALSE)

# also need to check normality of residuals
plot(egfp_elisa_comparisons_aov, 2)
# most of the points fall along the line representing quantiles of the normal distribution
# we can also use the Shapiro-Wilk test
egfp_elisa_comparisons_aov_residuals <- residuals(egfp_elisa_comparisons_aov)
shapiro.test(egfp_elisa_comparisons_aov_residuals)
# if we failed the assumptions of normality, we could use the Kruskal-Wallis rank sum test
kruskal.test(log10(concentration) ~ region, data = egfp_elisa_comparisons)

# going back to the one-way ANOVA, we need to test for post-hoc differences
# we can do this with using Tukey multiple pairwise-comparisons
TukeyHSD(egfp_elisa_comparisons_aov)
# you could also use pairwise t-tests using the exact code as above, but
# this approach would not use the variance of all the data, only the variance
# of the two groups being compared

# if we had used the Kruskal-wallis test, we could use Dunn's test or 
# Wilcoxon's test for multiple pairwise non-parametric comparisons, 
# correcting for multiple testing: https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/#multiple-pairwise-comparisons

# to plot the data with post hoc test results
library(rstatix)
stat.test <- egfp_elisa_comparisons_aov %>% 
  tukey_hsd() %>%
  pval_asterisks(p.adj)
# mutate(group1 = paste0(group1, "_REL"), 
#        group2 = paste0(group2, "_REL"))

p_egfp_elisa <- egfp_elisa %>%
  ungroup %>%
  mutate(concentration = log10(concentration + 1)) %>%
  mutate(region_genotype = factor(ifelse(genotype == "RELfs/fs", 
                                         paste0(region, "\nRELfs/fs"), 
                                         as.character(region)), 
                                  levels = c("MB", "DS", "VS", "CTX", 
                                             "MB\nRELfs/fs"))) %>%
  ggboxplot(x = "region_genotype", y = "concentration", 
            fill = "region_genotype", 
            outlier.shape = NA) +
  geom_quasirandom(shape = 21,
                   colour = "black",
                   size = 1.5) +
  stat_pvalue_manual(stat.test,
                     label = "ast",
                     y.position = c(3.4, 3.6, 3.8, 1.9, 2.1, 2.45)) +
  scale_fill_d3() +
  labs(x = "Region",
       y = expression(Log[10] ~ eGFP ~ Content ~ (pg))) +
  theme(legend.position = "none")

rm(list = ls(pattern = "^egfp"))

# qPCR marker -------------------------------------------------------------

label <-
  tibble(
    region = c("MB", "MB"),
    gene = c("Slc6a3", "Slc6a3"),
    fold_change = c(2, 0.5),
    label = c("Enrichment", "Depletion")
  ) %>%
  mutate(region = factor(region, levels = c("MB", "DS", "VS")))

marker_qpcr <-
  read_delim("input/experiments/trap_qpcr_markers_20210223.csv",
             delim = ",") %>%
  mutate(region = factor(region, levels = c("MB", "DS", "VS")),
         gene = factor(gene, levels = c("Th", "Slc6a3", "Gad2", "Gfap")), 
         fold_change = log2(fold_change + 1e-10))

stat.test <- marker_qpcr %>%
  group_by(region, gene) %>%
  t_test(fold_change ~ 1, mu = 0) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  pval_asterisks(padj) %>%
  left_join(marker_qpcr) %>%
  group_by(gene, region) %>%
  summarise(fold_change = max(fold_change), 
            ast = ast)

p_marker_qpcr <- marker_qpcr %>%
  ggboxplot(x = "gene", y = "fold_change", fill = "gene", 
            outlier.shape = NA) +
  geom_quasirandom(shape = 21,
                   size = 0.5,
                   colour = "black") +
  facet_wrap(vars(region)) +
  panel_border() +
  scale_y_continuous(limits = c(-4, 12), n.breaks = 10) +
  scale_fill_d3() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             size = 1.2) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 8, 1, 1), "lines"),
    strip.text = element_text(face = 2), 
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  labs(
    x = "Gene",
    fill = "Gene",
    y = expression(Log[2] ~ Fold ~ Change ~ "in" ~ TRAP)
  ) +
  geom_text(
    x = 15.5,
    aes(label = label),
    data = label,
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  geom_text(data = stat.test, 
            aes(label = ast, 
                y = fold_change + 1), 
            check_overlap = TRUE, 
            size = 3)

rm(list = c("label",
            "marker_qpcr"))


# cohort design -----------------------------------------------------------

t_ip_conditions <- metadata %>%
  select(sample_name, cohort, fraction, region, age, genotype) %>%
  distinct() %>%
  filter(fraction == "IP") %>%
  group_by(region, age, genotype, cohort) %>%
  summarise(n = n()) %>%
  mutate(animals = ifelse(cohort == "C2" & region != "MB", 
                    n*4, n)) %>%
  group_by(region, age, genotype) %>%
  summarise(n = sum(n), 
            animals = sum(animals)) %>%
  rename_with(str_to_title)

t_total_conditions <- metadata %>%
  select(sample_name, cohort, fraction, region, age, genotype) %>%
  distinct() %>%
  filter(fraction == "TOTAL") %>%
  group_by(region, age, genotype) %>%
  summarise(n = n()) %>%
  rename_with(str_to_title)

t_conditions <- t_ip_conditions %>%
  left_join(
    t_total_conditions,
    by = c("Region", "Age", "Genotype"),
    suffix = c(" TRAP Samples", " TOTAL Samples")
  ) %>%
  select(Region, Age, Genotype, `N TRAP Samples`, `N TOTAL Samples`, "N Animals" = Animals)


# rna yield ----------------------------------

rna_yield_IP <- metadata %>%
  select(
    cohort,
    collection,
    fraction,
    compartment,
    region,
    age,
    genotype,
    sample_name,
    rna_yield
  ) %>%
  distinct() %>%
  filter(!is.na(rna_yield) &
           fraction == "IP")

rna_yield_outliers_high <- rna_yield_IP %>%
  group_by(cohort, compartment) %>%
  filter(outlier_iqr(rna_yield, low_high = "high")) %>%
  group_by(cohort)

rna_yield_outliers_low <- rna_yield_IP %>%
  group_by(cohort, compartment) %>%
  filter(outlier_iqr(rna_yield, low_high = "low")) %>%
  group_by(cohort)

outlier_tracking <- update_outlier_tracking(rna_yield_outliers_high,
                                            reason = "High RNA Yield")
outlier_tracking <- update_outlier_tracking(rna_yield_outliers_low,
                                            reason = "Low RNA Yield")

rm(list = c("rna_yield_outliers_high",
            "rna_yield_outliers_low"))

n_rna_yield_outliers <- outlier_tracking %>%
  filter(str_detect(reason, "Yield")) %>%
  nrow

rna_yield_IP <- rna_yield_IP %>%
  mutate(outlier = sample_name %in% (outlier_tracking %>% filter(str_detect(reason, "Yield")))$sample_name) %>%
  rename_with(remove_underscores) %>%
  rename_with(str_to_title) %>%
  rename("RNA Yield (pg/uL)" = "Rna Yield")

rna_yield_IP_noOutliers <- rna_yield_IP %>%
  filter(!Outlier)

t_outlier_rna_yield <- outlier_tracking %>%
  filter(str_detect(reason, "Yield")) %>%
  rename_with(remove_underscores) %>%
  rename_with(str_to_title)

# cohort

# variances are not homogenous
rna_yield_IP_noOutliers %>%
  group_by(Region) %>%
  rstatix::levene_test(`RNA Yield (pg/uL)` ~ Cohort)
# the data are not normal
rna_yield_IP_noOutliers %>%
  group_by(Region) %>%
  select(rna_yield = "RNA Yield (pg/uL)") %>%
  rstatix::shapiro_test(rna_yield)

# so we use the Kruskal-Wallis test
stat.test <- rna_yield_IP_noOutliers %>%
  group_by(Region) %>%
  rstatix::kruskal_test(`RNA Yield (pg/uL)` ~ Cohort)

# find the signif comparisons
stat.test <- rna_yield_IP_noOutliers %>%
  group_by(Region) %>% 
  rstatix::pairwise_wilcox_test(`RNA Yield (pg/uL)` ~ Cohort, p.adjust.method = "BH") %>%
  mutate(p.adj = signif(p.adj, 3)) %>%
  pval_asterisks(p.adj)

p_rna_yield_IP_cohort <- rna_yield_IP %>% 
  mutate(Region = factor(Region, levels = c("MB", "DS", "VS"))) %>%
  ggboxplot(x = "Cohort", y = "RNA Yield (pg/uL)", fill = "Cohort", 
            outlier.shape = NA, facet.by = "Region", scales = "free_x") +
  geom_quasirandom(aes(group = Cohort, 
                       shape = Outlier), 
                   dodge.width = 0.75, 
                   width = 0.1, 
                   size = 1) +
  scale_shape_manual(values = c(21, 4)) +
  scale_color_d3() +
  scale_fill_d3() +
  # facet_grid(vars(Region), scales = "free_x") +
  stat_pvalue_manual(stat.test, 
                     label = "ast",
                     y.position = c(8e3, 1.275e4, 1.4e4, 8e3, 9.1e3)) +
  panel_border() +
  scale_y_continuous(n.breaks = 8, 
                     limits = c(0, 15000))

# rna_yield_IP %>%
#   filter(!sample_name %in% (outlier_tracking %>%
#                               filter(str_detect(reason, "Yield")))$sample_name) %>%
#   do(tidy(lm(
#     rna_yield ~ cohort + age + region + genotype, data = .
#   )))

# region difference

# variances are not homogenous
rna_yield_IP_noOutliers %>%
  filter(Cohort != "C1") %>% # cohort 1 only has MB
  group_by(Cohort) %>%
  rstatix::levene_test(`RNA Yield (pg/uL)` ~ Region)

# the data are not normal
rna_yield_IP_noOutliers %>%
  filter(Cohort != "C1") %>% # cohort 1 only has MB
  group_by(Cohort) %>%
  select(rna_yield = "RNA Yield (pg/uL)") %>%
  rstatix::shapiro_test(rna_yield)

# so we use the Kruskal-Wallis test
stat.test <- rna_yield_IP_noOutliers %>%
  filter(Cohort != "C1") %>% # cohort 1 only has MB
  group_by(Cohort) %>%
  rstatix::kruskal_test(`RNA Yield (pg/uL)` ~ Region)

# find the signif comparisons
stat.test <- rna_yield_IP_noOutliers %>%
  filter(Cohort != "C1") %>% # cohort 1 only has MB
  group_by(Cohort) %>% 
  rstatix::pairwise_wilcox_test(`RNA Yield (pg/uL)` ~ Region, p.adjust.method = "BH") %>%
  pval_asterisks(p.adj)

p_rna_yield_IP_region <- rna_yield_IP %>% 
  mutate(Region = factor(Region, levels = c("MB", "DS", "VS"))) %>%
  ggboxplot(x = "Region", y = "RNA Yield (pg/uL)", fill = "Region", 
            outlier.shape = NA, facet.by = "Cohort", scales = "free_x") +
  geom_quasirandom(aes(group = Region, 
                       shape = Outlier), 
                   dodge.width = 0.75, 
                   width = 0.1, 
                   size = 0.5) +
  scale_shape_manual(values = c(21, 4)) +
  scale_color_d3() +
  scale_fill_d3() +
  stat_pvalue_manual(stat.test, 
                     y.position = c(8.5e3, 1.05e4, 1.25e4, 1.25e4, 1.35e4, 1.5e4), 
                     label = "ast", 
                     size = 2) +
  panel_border() +
  theme(legend.position = "top", 
        legend.box = "vertical", 
        axis.text = element_text(size = 8)) +
  scale_y_continuous(n.breaks = 8, 
                     expand = expansion(mult = c(0.05, 0.1))) +
  guides(shape = "none")

# age difference
stat.test <- rna_yield_IP_noOutliers %>%
  group_by(Cohort, Region) %>%
  rstatix::kruskal_test(`RNA Yield (pg/uL)` ~ Age)

# find the signif comparisons
stat.test <- rna_yield_IP_noOutliers %>%
  group_by(Cohort, Region) %>%
  rstatix::pairwise_wilcox_test(`RNA Yield (pg/uL)` ~ Age, p.adjust.method = "BH") %>%
  pval_asterisks(p.adj)

p_rna_yield_IP_age <- rna_yield_IP %>%
  ggboxplot(x = "Age", y = "RNA Yield (pg/uL)", fill = "Age", 
            facet.by = c("Cohort", "Region"), scales = "free_x", 
            outlier.shape = NA) +
  geom_quasirandom(aes(group = Age, 
                       shape = Outlier), 
                   dodge.width = 0.75, 
                   width = 0.1, 
                   shape = 21, 
                   size = 0.1) +
  scale_shape_manual(values = c(21, 4)) +
  scale_color_d3() +
  scale_fill_d3() +
  stat_pvalue_manual(stat.test, 
                     y.position = c(9e3, 9e3, 
                                    1.1e4, 9.5e3, 
                                    1.3e4, 5e3, 9e3), 
                     label = "ast", 
                     size = 3, 
                     vjust = -0.5) +
  scale_y_continuous(n.breaks = 6, expand = expansion(mult = c(0.1, 0.15)), 
                     limits = c(0, NA)) +
  panel_border() +
  theme(legend.position = "top", 
        axis.text = element_text(size = 8))

# genotype difference
stat.test <- rna_yield_IP_noOutliers %>%
  filter(Cohort != "C2") %>%
  group_by(Cohort, Region) %>%
  rstatix::kruskal_test(`RNA Yield (pg/uL)` ~ Genotype)

# find the signif comparisons
stat.test <- rna_yield_IP_noOutliers %>%
  group_by(Cohort, Region) %>%
  rstatix::pairwise_wilcox_test(`RNA Yield (pg/uL)` ~ Age, p.adjust.method = "BH") %>%
  pval_asterisks(p.adj)


# collection difference
# rna_yield_IP %>%
#   filter(!sample_name %in% (outlier_tracking %>%
#                               filter(str_detect(reason, "Yield")))$sample_name) %>%
#   group_by(cohort) %>%
#   do(tidy(lm(rna_yield ~ collection, data = .))) %>%
#   filter(p.value < 0.05)





# Once I start talking about sequencing QC, the first thing to talk about is
# adapter content: this justifies read trimming.
# There will be a plot of adapter content by base position for cohort 3,
# a table of the adapter sequences trimmed
#
# Then I can go on to talk about alignment percentage
# alignment category
# fragment length
#
# Then reasons why reads failed to map after trimming (from FASTQ QC report of unmapped reads)

# Although adapter content comes before trimming in results1,
# the list of trimming outliers is required for p_adapter_content
# so trimming code is run first




# trimming percentage base ------------------------------------------------

# cutadapt-multiqc does not report read percentage trimmed separately.
# it reports the number of bases trimmed for r1 and r2 however

#import multiqc's summary of fastqc data - need to fix snakemake to output report to file and not stdout
# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/testing_environment/cutadapt_multiqc/multiqc_data/multiqc_cutadapt.txt"
multiqc_cutadapt <-
  read_delim(
    "input/qc/multiqc_cutadapt.txt",
    delim = "\t"
  ) %>%
  mutate(fastq_id = str_extract(Sample, "(?<=WTCHG_)[0-9|\\_]+(?=\\_1|\\_2|\\_val)")) %>%
  inner_join(metadata) %>%
  filter(!str_detect(sample_name, "BLANK"))

# there is no difference in lane, so averaging across lanes
multiqc_cutadapt %>%
  group_by(cohort) %>%
  do(tidy(lm(percent_trimmed ~ lane_code, data = .)))

# average across lanes
multiqc_cutadapt <- multiqc_cutadapt %>%
  group_by(cohort, sample_name) %>%
  summarise(percent_trimmed = mean(percent_trimmed)) %>%
  group_by(cohort) %>%
  mutate(outlier = outlier_iqr(percent_trimmed, low_high = "high")) #We are only interested in samples with excessive trimming. So abs() is not used here


# trimming percentage cohort ----------------------------------------------

p_trim_percent_cohort <- multiqc_cutadapt %>%
  ggplot(aes(x = cohort,
             y = percent_trimmed)) +
  geom_quasirandom(aes(fill = outlier),
                   shape = 21,
                   colour = "black",
                   size = 2) +
  scale_fill_d3() +
  # ggtitle("Read Trimming: Cohort and Region Differences") +
  labs(x = "Cohort",
       y = "Percentage of bp trimmed (%)",
       fill = "Outlier") +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 3)))


# trimming percentage outliers --------------------------------------------

t_trim_percent_outliers <- multiqc_cutadapt %>%
  group_by(cohort) %>%
  filter(outlier_iqr(percent_trimmed, low_high = "high")) %>%
  select(sample_name, percent_trimmed) %>%
  filter(percent_trimmed > 1)

outlier_tracking <- update_outlier_tracking(t_trim_percent_outliers,
                                            reason = "High bp trimming percentage")

t_trim_percent_outliers <- outlier_tracking %>%
  filter(str_detect(reason, "trimming percentage")) %>%
  rename_with(remove_underscores) %>%
  rename_with(str_to_title)



# trimming percentage condition (region) ----------------------------------------------

# region the only condition with an effect
multiqc_cutadapt %>% inner_join(metadata) %>%
  filter(cohort == "C3" &
           !sample_name %in% t_trim_percent_outliers$`Sample Name`) %>%
  do(tidy(lm(
    percent_trimmed ~ region + age + genotype, data = .
  )))

data <- multiqc_cutadapt %>% inner_join(metadata) %>%
  filter(cohort == "C3" &
           !sample_name %in% t_trim_percent_outliers$`Sample Name`) 

# shapiro: fail
data %>%
  group_by(region, fraction) %>%
  shapiro_test(percent_trimmed)

# levene: fail
data %>%
  group_by(fraction) %>%
  levene_test(percent_trimmed ~ region)

# kruskal
data %>%
  group_by(fraction) %>%
  kruskal_test(percent_trimmed ~ region)

# wilcoxon
stat.test <- data %>% group_by(fraction) %>%
  pairwise_wilcox_test(percent_trimmed ~ region, p.adjust.method = "BH") %>%
  pval_asterisks(p.adj)

p_trim_percent_region <- data %>%
  ggboxplot(x = "region", 
            y = "percent_trimmed", 
            fill = "region", 
            outlier.shape = NA, 
            facet.by = "fraction") +
  geom_quasirandom(aes(group = region), 
                   dodge.width = 0.75, 
                   width = 0.1, 
                   shape = 21, 
                   size = 0.5) +
  scale_color_d3() +
  scale_fill_d3() +
  stat_pvalue_manual(stat.test, 
                     y.position = c(35.5, 36.5, 33, 31.5, 34, 35.5), 
                     label = "ast", 
                     vjust = -0.5,
                     size = 3) +
  panel_border() +
  labs(x = "Region", fill = "Region") +
  theme(axis.title.y = element_blank(), 
        legend.position = "top", 
        axis.text = element_text(size = 8)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
  
# 
# p_trim_percent_region <-
#   multiqc_cutadapt %>% inner_join(metadata) %>%
#   filter(cohort == "C3" &
#            !sample_name %in% t_trim_percent_outliers$`Sample Name`) %>%
#   ggplot(aes(region, percent_trimmed,
#              fill = region)) +
#   geom_quasirandom(shape = 21,
#                    colour = "black",
#                    size = 2) +
#   scale_fill_d3() +
#   stat_compare_means(comparisons = list(c("MB", "DS"), c("MB", "VS")),
#                      label = "p.signif") +
#   labs(x = "Region",
#        # y = "Percentage of bp trimmed (%)",
#        fill = "Region") +
#   guides(fill = guide_legend(override.aes = list(size = 3))) +
#   theme(legend.position = "top",
#         axis.title.y = element_blank(),
#         panel.grid.major.x = element_blank()) +
#   scale_y_continuous(limits = c(NA, 38))

# cohort 3 trim percent range ---------------------------------------------

n_trim_percent_c3_range <-
  multiqc_cutadapt %>% filter(cohort == "C3") %>% filter(!outlier) %>% pull(percent_trimmed) %>% range %>% round(digits = 0)








# cohort 3 trimmed read length --------------------------------------------
# this is not being used because rnaseqc reports the MAXIMUM read length aligned per sample, not the mean read length
# multiqc_rnaseqc <- read_delim("/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/qc/multiqc/multiqc_data/multiqc_rna_seqc.txt",
#                               delim = "\t") %>%
#   mutate(fastq_id = str_extract(Sample, "[0-9|\\_]+(?=\\_csorted)")) %>%
#   inner_join(metadata) %>%
#   group_by(sample_name, cohort, plate, collection, fraction, region, age, genotype) %>% # average lanes
#   summarise(across(everything(), mean))
#
# multiqc_rnaseqc %>%
#   group_by(cohort) %>%
#   # filter(cohort == "C3") %>%
#   do(tidy(lm(`Read Length`~fraction, data = .)))
#
# p_read_length_c3 <- multiqc_rnaseqc %>%
#   filter(cohort == "C3" ) %>%
#   ggplot(aes(x = `Read Length`)) +
#   stat_ecdf(geom = "step", pad = FALSE) +
#   scale_color_d3() +
#   labs(x = "Read Length",
#        y = "Cumulative \nProportion")
#
# p_read_length_c3_fraction <- multiqc_rnaseqc %>%
# filter(cohort == "C3") %>%
#   ggplot(aes(x = `Read Length`,
#              colour = fraction)) +
#   stat_ecdf(geom = "step", pad = FALSE) +
#   scale_color_d3() +
#   labs(x = "Read Length",
#        y = "Cumulative \nProportion",
#        colour = "Fraction") +
#   theme(legend.position = "top")
#
# p_read_length_c3_region_ip <- multiqc_rnaseqc %>%
#   filter(cohort == "C3" &
#            fraction == "IP") %>%
#   ggplot(aes(x = `Read Length`,
#              colour = region)) +
#   stat_ecdf(geom = "step", pad = FALSE) +
#   scale_color_d3() +
#   labs(x = "Read Length",
#        # y = "Cumulative Proportion",
#        colour = "Region") +
#   theme(legend.position = "top",
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank())
#
# p_read_length_c3_age_ip <- multiqc_rnaseqc %>%
#   filter(cohort == "C3" &
#            fraction == "IP") %>%
#   ggplot(aes(x = `Read Length`,
#              colour = age)) +
#   stat_ecdf(geom = "step", pad = FALSE) +
#   scale_color_d3() +
#   labs(x = "Read Length",
#        # y = "Cumulative Proportion",
#        colour = "Age") +
#   theme(legend.position = "top",
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank())













# load ngsReports ---------------------------------------------------------
#Using ngsReports to obtain fastqc data for plotting. https://bioconductor.org/packages/release/bioc/vignettes/ngsReports/inst/doc/ngsReportsIntroduction.html#mean-sequence-quality-per-read
#I don't want to load ngsReports fully, because it masks so many packages

# fileDir <- "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/qc"
# files <- list.files(fileDir, pattern = "fastqc.zip$", full.names = TRUE,
#                     recursive = TRUE)
# fdl <- ngsReports::FastqcDataList(files)
# saveRDS(fdl, file = "R/objects/fdl.rds")
fdl <- readRDS("R/objects/fdl.rds")
# add_ngsReports function -------------------------------------------------
add_ngsReports_metadata <- function(report) {
  report %>%
    mutate(
      read_mate = ifelse(str_detect(Filename, "_1"),
                         "Read 1", "Read 2"),
      fastq_id = str_extract(
        Filename,
        "(?<=WTCHG_)[0-9|\\_]+(?=\\_1|\\_2|\\_val|\\_mapped|\\_unmapped)"
      ),
      mapped = ifelse(
        str_detect(Filename, "unmapped"),
        "unmapped",
        ifelse(str_detect(Filename, "mapped"), "mapped", "input")
      ),
      trimmed = ifelse(
        str_detect(Filename, "val"),
        "Trimmed",
        ifelse(mapped != "input", "Trimmed", "Untrimmed")
      )
    ) %>%
    inner_join(metadata)
}

# Adapter content base ---------------------------------------------------------
adapter_levels <- ngsReports::getModule(fdl, "Adapter_Content") %>%
  add_ngsReports_metadata() %>%
  filter(cohort == "C3" &
           trimmed == "Untrimmed") %>%
  pull(Position) %>%
  unique()

adapter_content <- ngsReports::getModule(fdl, "Adapter_Content") %>%
  add_ngsReports_metadata() %>%
  filter(cohort == "C3") %>%
  filter(read_mate == "Read 1") %>% # read 1 and 2 are virtually the same
  mutate(Position = factor(Position, levels = adapter_levels)) %>%
  pivot_longer(
    cols = c(
      "Illumina_Universal_Adapter",
      "Illumina_Small_RNA_3'_Adapter",
      "Illumina_Small_RNA_5'_Adapter",
      "Nextera_Transposase_Sequence",
      "SOLID_Small_RNA_Adapter"
    ),
    names_to = "adapter",
    values_to = "percentage"
  ) %>%
  filter(adapter == "Illumina_Universal_Adapter")

# adapter content by lane - no difference --------------------------------
# p_adapter_content_lane <- adapter_content %>%
#   filter(cohort == "C3") %>%
#   ggplot(aes(x = Position,
#              y = percentage)) +
#   geom_line(aes(group = Filename),
#             alpha = 0.5) +
#   facet_wrap(vars(lane_code),
#              nrow = 2)



# adapter content by condition: no difference -----------------------------
# t_adapter_content_condition <- adapter_content %>%
#   filter(cohort == "C3" &
#            Position == "140" &
#            trimmed == "Untrimmed") %>%
#   group_by(fraction, region, age, genotype, sample_name) %>%
#   summarise(percentage = mean(percentage)) %>%
#   ungroup() %>%
#   do(tidy(lm(percentage~fraction+region+age+genotype, data = .)))


# adapter content cohort 3 ------------------------------------------------
adapter_content <- adapter_content %>%
  group_by(Position, trimmed, sample_name) %>%
  summarise(percentage = mean(percentage)) %>%
  ungroup()

p_adapter_content <- adapter_content %>%
  filter(trimmed == "Untrimmed") %>%
  mutate(Outlier = sample_name %in% t_trim_percent_outliers$`Sample Name`) %>%
  ggplot(aes(x = Position,
             y = percentage)) +
  geom_line(aes(group = sample_name,
                colour = Outlier),
            alpha = 0.5) +
  scale_color_d3() +
  # ggtitle("Adapter content per base: Cohort 3") +
  labs(y = "Percentage") +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ),
  legend.position = "top",
  panel.grid.major.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))









# clontech sequence detection ---------------------------------------------

n_percent_clontech <-
  suppressMessages(ngsReports::getModule(fdl, "Overrepresented_sequences")) %>%
  add_ngsReports_metadata() %>%
  filter(read_mate == "Read 1" &
           trimmed == "Untrimmed" &
           mapped == "input") %>%
  filter(str_detect(Possible_Source, "Clontech")) %>%
  pull(sample_name) %>%
  unique() %>%
  length() * 100 / 120




# trimming per base quality --------------------------------------------------------

#Per base quality - cohort
base_levels <-
  ngsReports::getModule(fdl, "Per_base_sequence_quality") %>%
  add_ngsReports_metadata() %>%
  filter(cohort == "C3" &
           mapped == "input") %>%
  pull(Base) %>%
  unique()

p_per_base_sequence_quality <-
  ngsReports::getModule(fdl, "Per_base_sequence_quality") %>%
  add_ngsReports_metadata() %>%
  filter(cohort == "C3" &
           mapped == "input") %>%
  mutate(Base = factor(Base, levels = base_levels)) %>%
  group_by(sample_name, trimmed, Base) %>% # To average across lanes
  mutate(Mean = mean(Mean)) %>%
  select(sample_name,
         trimmed,
         Base,
         Mean,
         cohort,
         region,
         age,
         genotype,
         sex,
         fraction) %>%
  distinct() %>%
  unite(
    col = sample_name_trimmed,
    c(sample_name, trimmed),
    sep = "_",
    remove = F
  ) %>%
  ggplot(aes(x = Base,
             y = Mean)) +
  geom_line(alpha = 0.2,
            aes(group = sample_name_trimmed,
                colour = trimmed)) +
  labs(y = "Q-Score",
       colour = "Trimming",
       x = "Position") +
  scale_color_d3() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ),
  legend.position = "top",
  panel.grid.major.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))








# alignment - import data --------------------------------------------------------------------

# STAR
# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/testing_environment/star_multiqc/multiqc_data/multiqc_star.txt"
multiqc_star <-
  read_delim(
    "input/qc/multiqc_star.txt",
    delim = "\t"
  ) %>%
  mutate(
    trimmed = ifelse(!str_detect(Sample, "untrimmed"), "Trimmed", "Untrimmed"),
    fastq_id = str_extract(Sample, "^[0-9]+_[0-9]+")
  ) %>%
  inner_join(metadata) %>%
  select(fastq_id, trimmed, everything())

# Salmon
# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/multiqc_temp/multiqc_data/multiqc_salmon.txt"
multiqc_salmon <-
  read_delim(
    "input/qc/multiqc_salmon.txt",
    delim = "\t"
  ) %>%
  mutate(
    trimmed = ifelse(!str_detect(Sample, "untrimmed"), "Trimmed", "Untrimmed"),
    fastq_id = str_extract(Sample, "^[0-9]+_[0-9]+")
  ) %>%
  select(fastq_id, trimmed, everything()) %>%
  inner_join(metadata)

# Kallisto
# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/multiqc_temp/multiqc_data/multiqc_kallisto.txt"
multiqc_kallisto <-
  read_delim(
    "input/qc/multiqc_kallisto.txt",
    delim = "\t"
  ) %>%
  mutate(
    trimmed = ifelse(str_detect(Sample, "val"), "Trimmed", "Untrimmed"),
    fastq_id = str_extract(Sample, "[0-9]+_[0-9]+")
  ) %>%
  select(fastq_id, trimmed, everything()) %>%
  inner_join(metadata)


# alignment percentage - merged -----------------------------------------

alignment_by_aligners <- multiqc_star %>%
  rename(alignment_percentage = uniquely_mapped_percent,
         aligned_reads = uniquely_mapped) %>%
  group_by(sample_name, trimmed, cohort) %>%
  summarise(
    alignment_percentage = mean(alignment_percentage),
    aligned_reads = sum(aligned_reads)
  ) %>%
  mutate(aligner = "STAR") %>%
  bind_rows((
    multiqc_salmon %>%
      rename(alignment_percentage = percent_mapped,
             aligned_reads = num_mapped) %>%
      group_by(sample_name, trimmed, cohort) %>%
      summarise(
        alignment_percentage = mean(alignment_percentage),
        aligned_reads = sum(aligned_reads)
      ) %>%
      mutate(aligner = "Salmon")
  ),
  (
    multiqc_kallisto %>%
      rename(alignment_percentage = percent_aligned,
             aligned_reads = pseudoaligned_reads) %>%
      group_by(sample_name, trimmed, cohort) %>%
      summarise(
        alignment_percentage = mean(alignment_percentage),
        aligned_reads = sum(aligned_reads)
      ) %>%
      mutate(aligner = "Kallisto")
  )
  ) %>%
  inner_join((
    metadata %>%
      select(sample_name,
             collection,
             fraction,
             region,
             age,
             genotype,
             sex) %>%
      distinct()
  )) %>%
  mutate(trimmed = factor(trimmed, levels = c("Untrimmed", "Trimmed")))

# alignment percentage - outliers -----------------------------------------

t_alignment_percentage_outliers <- alignment_by_aligners %>%
  filter(trimmed == "Trimmed") %>%
  group_by(cohort, sample_name) %>%
  summarise(alignment_percentage = mean(alignment_percentage)) %>%
  filter(outlier_iqr(alignment_percentage, low_high = "low")) %>%
  filter(alignment_percentage < 80)

outlier_tracking <-
  update_outlier_tracking(t_alignment_percentage_outliers, reason = "Low alignment percentage")

t_alignment_percentage_outliers <- outlier_tracking %>%
  filter(str_detect(reason, "alignment")) %>%
  rename_with(remove_underscores) %>%
  rename_with(str_to_title)

# remove outliers from plots
alignment_by_aligners <- alignment_by_aligners %>%
  filter(!sample_name %in% t_alignment_percentage_outliers$sample_name)


# alignment percentage - trimming -------------------------------------------

# two groups in each comparison
# unequal variance
alignment_by_aligners %>%
  group_by(cohort, aligner) %>%
  levene_test(formula = alignment_percentage ~ trimmed)

# data are not normal
alignment_by_aligners %>%
  group_by(cohort, aligner) %>%
  rstatix::shapiro_test(alignment_percentage)

# wilcoxon rank sum test
stat.test <- alignment_by_aligners %>%
  group_by(cohort, aligner) %>%
  rstatix::wilcox_test(formula = alignment_percentage ~ trimmed) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  pval_asterisks(padj)

p_alignment_percentage_aligner_trimming <- ggboxplot(alignment_by_aligners, 
          x = "trimmed", y = "alignment_percentage", 
          fill = "cohort", facet.by = c("cohort", "aligner"), 
          scales = "free_y", outlier.shape = NA) +
  geom_quasirandom(shape = 21, 
                   size = 0.5) +
  scale_y_continuous(n.breaks = 6, 
                     expand = expansion(mult = c(0.1, 0.15)), 
                     limits = c(NA, 100)) +
  theme(axis.text = element_text(size = 8), 
        legend.position = "none",
        axis.title.x = element_blank()) +
  stat_pvalue_manual(stat.test, 
                     y.position = c(97, 91, 91, 99, 95, 92, 100, 100, 100), 
                     label = "ast", 
                     size = 3, 
                     vjust = -0.5) +
  labs(x = "Trimmed", y = "Alignment Percentage") +
  panel_border()

# p_alignment_percentage_aligner_trimming <- alignment_by_aligners %>%
#   ggplot(aes(x = trimmed,
#              y = alignment_percentage,
#              fill = cohort)) +
#   geom_quasirandom(shape = 21,
#                    colour = "black",
#                    size = 2) +
#   scale_fill_d3() +
#   facet_wrap(vars(cohort, aligner)) +
#   theme(legend.position = "none",
#         axis.title.x = element_blank()) +
#   labs(x = "Trimmed", y = "Alignment Percentage")
# # ggtitle("Alignment by Cohort and Aligner: Effect of Trimming")

n_alignment_percentage_cohort <- alignment_by_aligners %>%
  filter(trimmed == "Trimmed") %>%
  group_by(cohort) %>%
  summarise(alignment_percentage = mean(alignment_percentage)) %>%
  pull(alignment_percentage) %>%
  round(1)

n_alignment_percentage_cohort_mean <- mean(n_alignment_percentage_cohort)

# remove untrimmed -------------------------------------------

alignment_by_aligners <- alignment_by_aligners %>%
  filter(trimmed == "Trimmed")

# alignment percentage - aligner differences ------------------------------

symnums <- list(cutpoints = c(0, 0.001, 0.01, 0.1, 1), symbols = c("***", "**", "*", "NS"))

p_alignment_percentage_aligners <- alignment_by_aligners %>%
  ggboxplot(x = "aligner", y = "alignment_percentage", 
            fill = "aligner", outlier.shape = NA) +
  # ggplot(aes(x = aligner,
  #            y = alignment_percentage,
  #            fill = aligner)) +
  geom_quasirandom(shape = 21,
                   size = 1,
                   colour = "black") +
  scale_fill_d3()  +
  theme(legend.position = "none") +
  labs(x = "Aligner", y = "Alignment Percentage") +
  stat_compare_means(comparisons = list(c("Salmon", "STAR"), 
                                        c("Kallisto", "Salmon"),
                                        c("Kallisto", "STAR")),
                     label = "p.signif",
                     method = "wilcox.test", 
                     hide.ns = F,
                     symnum.args = symnums) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
# ggtitle("Alignment Percentage: By Aligner, Fraction and Region")


# alignment percentage - lane differences ---------------------------------
# no major difference across lanes in alignment percentage
multiqc_star %>%
  rename(alignment_percentage = uniquely_mapped_percent) %>%
  group_by(sample_name, lane_code, trimmed, cohort) %>%
  summarise(alignment_percentage = mean(alignment_percentage)) %>%
  mutate(aligner = "STAR") %>%
  bind_rows((
    multiqc_salmon %>%
      rename(alignment_percentage = percent_mapped) %>%
      group_by(sample_name, lane_code, trimmed, cohort) %>%
      summarise(alignment_percentage = mean(alignment_percentage)) %>%
      mutate(aligner = "Salmon")
  ),
  (
    multiqc_kallisto %>%
      rename(alignment_percentage = percent_aligned) %>%
      group_by(sample_name, lane_code, trimmed, cohort) %>%
      summarise(alignment_percentage = mean(alignment_percentage)) %>%
      mutate(aligner = "Kallisto")
  )
  ) %>%
  group_by(cohort, trimmed) %>%
  do(tidy(lm(alignment_percentage ~ lane_code, data = .))) %>%
  filter(p.value < 0.05)

# alignment percentage - conditions ---------------------------------------

# check for region
alignment_by_aligners %>%
  filter(cohort != "C1") %>%
  group_by(aligner, cohort) %>%
  do(tidy(lm(
    alignment_percentage ~ fraction + age + region, data = .
  ))) %>%
  filter(p.value < 0.05)
# check for collection, and age (collection overlaps with fraction)
alignment_by_aligners %>%
  filter(fraction == "IP") %>%
  group_by(aligner, cohort) %>%
  do(tidy(lm(
    alignment_percentage ~ collection + age, data = .
  ))) %>%
  filter(p.value < 0.05)
# check for genotype
alignment_by_aligners %>%
  filter(cohort != "C2") %>%
  group_by(aligner, cohort) %>%
  do(tidy(lm(
    alignment_percentage ~ fraction + age + genotype, data = .
  ))) %>%
  filter(p.value < 0.05)
# check for fraction
alignment_by_aligners %>%
  group_by(aligner, cohort) %>%
  do(tidy(lm(alignment_percentage ~ fraction, data = .))) %>%
  filter(p.value < 0.05)


# fraction by aligner
p_alignment_percentage_fraction_aligner <- alignment_by_aligners %>%
  ggboxplot(x = "fraction", y = "alignment_percentage", 
            fill = "fraction", outlier.shape = NA) +
  geom_quasirandom(shape = 21,
                   size = 1,
                   colour = "black") +
  scale_fill_d3() +
  facet_wrap(vars(aligner)) +
  # stat_compare_means(label.x = 1.5,
  #                    label = "p.signif",
  #                    label.y = 98, 
  #                    symnum.args = symnums) +
  stat_compare_means(label = "p.signif",
                     method = "wilcox.test", 
                     label.x = 1.5,
                     hide.ns = F,
                     symnum.args = symnums) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), 
                     limits = c(NA, 100)) +
  theme(legend.position = "none") +
  labs(x = "Fraction", y = "Alignment Percentage") +
  panel_border()

# region in IP by aligner


alignment_by_aligners %>%
  filter(fraction == "IP") %>%
  group_by(aligner) %>%
  levene_test(formula = alignment_percentage ~ region)

stat.test <- alignment_by_aligners %>%
  filter(fraction == "IP") %>%
  group_by(aligner) %>%
  rstatix::wilcox_test(formula = alignment_percentage ~ region) %>%
  mutate(padj = p.adjust(p, "BH")) %>%
  pval_asterisks(padj)

p_alignment_percentage_ip_region_aligner <-
  alignment_by_aligners %>%
  filter(fraction == "IP") %>%
  ggboxplot(x = "region", y = "alignment_percentage", 
            fill = "region", outlier.shape = NA) +
  geom_quasirandom(shape = 21,
                   size = 1,
                   colour = "black") +
  scale_fill_d3() +
  facet_wrap(vars(aligner)) +
  stat_pvalue_manual(stat.test, 
                     label = "ast", 
                     size = 3, 
                     y.position = c(97, 99, 98, 95.5, 97.5, 96.5, 92.5, 94.5, 93.5)) +
  theme(legend.position = "none") +
  labs(x = "Region", y = "Alignment Percentage") +
  scale_y_continuous(limits = c(NA, 100)) +
  panel_border()












# alignment category - trimming - to help diagnose mapping failure cause -------------------------------------------------
p_star_mapping_categories_trimming <- multiqc_star %>%
  filter(cohort == "C3") %>%
  pivot_longer(
    cols = c(
      "uniquely_mapped_percent",
      "multimapped_percent",
      "multimapped_toomany_percent",
      "unmapped_mismatches_percent",
      "unmapped_tooshort_percent",
      "unmapped_other_percent"
    ),
    names_to = "Mapping Type",
    values_to = "Proportion"
  ) %>%
  mutate(`Mapping Type` = str_to_title(remove_underscores(gsub(
    "_percent", "", `Mapping Type`
  ))),
  trimmed = factor(trimmed, levels = c("Untrimmed", "Trimmed"))) %>%
  group_by(sample_name, cohort, trimmed, fraction, `Mapping Type`) %>%
  summarise(Proportion = mean(Proportion)) %>%
  ggplot(aes(x = sample_name,
             y = Proportion)) +
  geom_col(aes(fill = `Mapping Type`),
           colour = "black",
           size = 0.2) +
  scale_fill_d3() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    axis.line.x = element_blank()
  ) +
  facet_wrap(vars(trimmed),
             scales = "free",
             nrow = 1) +
  scale_y_continuous(expand = c(0, 2)) +
  labs(x = "Sample",
       y = "Percentage of reads (%)",
       fill = "Alignment Type")


# alignment category - cohort ----------------------------------------

p_star_mapping_categories_cohort <- multiqc_star %>%
  filter(trimmed == "Trimmed") %>%
  pivot_longer(
    cols = c(
      "uniquely_mapped_percent",
      "multimapped_percent",
      "multimapped_toomany_percent",
      "unmapped_mismatches_percent",
      "unmapped_tooshort_percent",
      "unmapped_other_percent"
    ),
    names_to = "Mapping Type",
    values_to = "Proportion"
  ) %>%
  mutate(`Mapping Type` = str_to_title(remove_underscores(gsub(
    "_percent", "", `Mapping Type`
  )))) %>%
  group_by(sample_name, cohort, trimmed, fraction, `Mapping Type`) %>%
  summarise(Proportion = mean(Proportion)) %>%
  ggplot(aes(x = sample_name,
             y = Proportion)) +
  geom_col(aes(fill = `Mapping Type`),
           colour = "black",
           size = 0.2) +
  scale_fill_d3() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10),
    axis.line.x = element_blank()
  ) +
  facet_wrap(vars(cohort),
             scales = "free",
             nrow = 1) +
  scale_y_continuous(expand = c(0, 2)) +
  labs(x = "Sample",
       y = "Percentage (%)")

# avg_input_read_length
# this is plotted against uniquely mapped percentage, because short read length is not a justifcation to exclude samples on its own
# it can provide some explanation why a sample might have a low mapping percentage, however
p_avg_input_read_length_mapping_rate <- multiqc_star %>%
  ggplot(aes(x = avg_input_read_length,
             y = uniquely_mapped_percent)) +
  geom_point()

# number of splice junctions reliably detected (in SJ.out.tab)
p_num_splices_fraction <- multiqc_star %>%
  ggplot(aes(x = fraction,
             y = num_splices,
             colour = fraction)) +
  geom_boxplot() +
  stat_compare_means(label.x = 1.5) +
  scale_color_d3()

#mismatch_rate
p_mismatch_rate <- multiqc_star %>%
  ggplot(aes(x = fraction,
             y = mismatch_rate,
             colour = fraction)) +
  geom_boxplot() +
  stat_compare_means(label.x = 1.5) +
  scale_color_d3()

outlier_tracking <- multiqc_star %>%
  filter(outlier_iqr(mismatch_rate)) %>%
  update_outlier_tracking(reason = "mismatch_rate")




#deletion_rate
#deletion_length
#insertion_rate
#insertion_length















# fragment length by aligner and cohort ---------------------------------------------------------

# p_fragment_length_aligner_cohort <- multiqc_salmon %>%
#   select(fastq_id, trimmed, lane_code, sample_name, "fragment_length_mean" = frag_length_mean,
#          cohort, collection, fraction, region, age, genotype, sex) %>%
#   mutate(aligner = "salmon") %>%
#   bind_rows(
#     (
#       multiqc_kallisto %>%
#         select(fastq_id, trimmed, lane_code, sample_name, "fragment_length_mean" = fragment_length,
#                cohort, collection, fraction, region, age, genotype, sex) %>%
#         mutate(aligner = "kallisto")
#     )
#   ) %>%
#   group_by(sample_name, trimmed, aligner) %>%
#   mutate(fragment_length_mean = mean(fragment_length_mean)) %>%
#   select(-lane_code) %>%
#   pivot_wider(names_from = "aligner",
#               values_from = "fragment_length_mean") %>%
#   filter(trimmed == "Trimmed") %>%
#   ggplot(aes(x = salmon,
#              y = kallisto,
#              fill = cohort)) +
#   geom_point(shape = 21, colour = "black", size = 2) +
#   geom_abline(intercept = 0,
#               slope = 1,
#               linetype = "dotted") +
#   scale_x_continuous(limits = c(130, 350)) +
#   scale_y_continuous(limits = c(130, 350)) +
#   labs(x = "Fragment Length: Salmon (bp)",
#        y = "Fragment Length: Kallisto (bp)",
#        fill = "Cohort") +
#   scale_fill_d3() +
#   theme(legend.position = "top") +
#   guides(fill = guide_legend(override.aes = list(size = 3) ) ) +
#   annotate("text", x = 340, y = 350, label = "paste(bold(y), \" = \", bold(x))",
#            parse = TRUE, size = 5)


# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/qc/multiqc/multiqc_data/multiqc_rna_seqc.txt"
multiqc_rnaseqc <-
  read_delim(
    "input/qc/multiqc_rna_seqc.txt",
    delim = "\t"
  ) %>%
  mutate(fastq_id = str_extract(Sample, "[0-9|\\_]+(?=\\_csorted)")) %>%
  inner_join(metadata)


p_fragment_length_aligner_cohort <- multiqc_salmon %>%
  select(
    fastq_id,
    trimmed,
    lane_code,
    sample_name,
    "fragment_length_mean" = frag_length_mean,
    cohort,
    collection,
    fraction,
    region,
    age,
    genotype,
    sex
  ) %>%
  mutate(aligner = "salmon") %>%
  bind_rows((
    multiqc_kallisto %>%
      select(
        fastq_id,
        trimmed,
        lane_code,
        sample_name,
        "fragment_length_mean" = fragment_length,
        cohort,
        collection,
        fraction,
        region,
        age,
        genotype,
        sex
      ) %>%
      mutate(aligner = "kallisto")
  )) %>%
  filter(trimmed == "Trimmed") %>%
  bind_rows((
    multiqc_rnaseqc %>%
      select(
        fastq_id,
        lane_code,
        sample_name,
        "fragment_length_mean" = `Average Fragment Length`,
        cohort,
        collection,
        fraction,
        region,
        age,
        genotype,
        sex
      ) %>%
      mutate(aligner = "STAR")
  )) %>%
  group_by(sample_name, aligner) %>%
  mutate(
    fragment_length_mean = mean(fragment_length_mean),
    aligner = ifelse(
      aligner == "STAR",
      "STAR",
      ifelse(aligner == "salmon", "Salmon", "Kallisto")
    )
  ) %>%
  select(-lane_code) %>%
  # pivot_wider(names_from = "aligner",
  #             values_from = "fragment_length_mean") %>%
  ggplot(aes(x = fragment_length_mean,
             fill = cohort)) +
  geom_density(alpha = 0.5) +
  facet_wrap(vars(aligner),
             nrow = 3) +
  scale_fill_d3() +
  labs(x = "Mean Fragment Length",
       y = "Density",
       fill = "Cohort") +
  theme(legend.position = "top")







# successful vs failed reads - FASTQC ---------------------------------------------------


# FASTQC summary
fastqc_summary <- ngsReports::getModule(fdl, "Summary") %>%
  add_ngsReports_metadata()

fastqc_summary %>%
  filter(cohort == "C3" &
           mapped != "input") %>%
  group_by(Category, Status, cohort, mapped) %>%
  summarise(status = n()) %>%
  filter(Status == "FAIL")


# failed reads - per base sequence content base -----------------------------------------------

per_base_sequence_content <-
  ngsReports::getModule(fdl, "Per_base_sequence_content") %>%
  add_ngsReports_metadata()


# per base sequence content plot  -----------------------------------------

# cohort 3
# cohort 3 base levels
base_levels <-
  ngsReports::getModule(fdl, "Per_base_sequence_content") %>%
  add_ngsReports_metadata() %>%
  filter(cohort == "C3" &
           mapped != "input") %>%
  pull(Base) %>%
  unique()

p_per_base_sequence_content_c3 <- per_base_sequence_content %>%
  filter(cohort == "C3" &
           mapped != "input",
         read_mate == "Read 1") %>%
  pivot_longer(cols = c(A, G, C, T),
               names_to = "Nucleotide",
               values_to = "Percentage") %>%
  mutate(
    Base = factor(Base, levels = base_levels),
    mapped = ifelse(mapped == "mapped", "Aligned", "Failed")
  ) %>%
  group_by(sample_name, mapped, Base, Nucleotide) %>% # To average across lanes
  mutate(Percentage = mean(Percentage)) %>%
  ggplot(aes(x = Base,
             y = Percentage)) +
  geom_line(alpha = 0.2,
            aes(group = interaction(sample_name, Nucleotide),
                colour = Nucleotide)) +
  scale_color_d3() +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  alpha = 1))) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    title = element_text(size = 12),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Position",
       y = "Percentage (%)") +
       # title = "Per Base Sequence Content: Aligned vs Failed Reads") +
       facet_wrap(vars(mapped), nrow = 2)





# failed reads - quality --------------------------------------------------

base_levels <-
 ngsReports::getModule(fdl, "Per_base_sequence_quality") %>%
 add_ngsReports_metadata() %>%
 filter(cohort == "C3" &
          mapped != "input") %>%
 pull(Base) %>%
 unique()

p_per_base_sequence_quality_failed <-
 ngsReports::getModule(fdl, "Per_base_sequence_quality") %>%
 add_ngsReports_metadata() %>%
 filter(cohort == "C3" &
          mapped != "input") %>%
 mutate(
   Base = factor(Base, levels = base_levels),
   mapped = ifelse(mapped == "mapped", "Aligned", "Failed")
 ) %>%
 group_by(sample_name, mapped, Base) %>%
 mutate(Mean = mean(Mean)) %>%
 select(sample_name,
        mapped,
        Base,
        Mean,
        cohort,
        region,
        age,
        genotype,
        sex,
        fraction) %>%
 distinct() %>%
 unite(
   col = sample_name_mapped,
   c(sample_name, mapped),
   sep = "_",
   remove = F
 ) %>%
 ggplot(aes(x = Base,
            y = Mean)) +
 geom_line(alpha = 0.2,
           aes(group = sample_name_mapped,
               colour = mapped)) +
 labs(y = "Q-Score",
      colour = "Alignment",
      x = "Position") +
 scale_color_d3() +
 theme(
   axis.text.x = element_text(
     angle = 90,
     vjust = 0.5,
     hjust = 1,
     size = 10
   ),
   legend.position = "top",
   panel.grid.major.x = element_blank()
 ) +
 guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))



# remove fdl ----
rm(fdl)

# library size - cohort --------------------------------------------------
p_aligned_reads <- alignment_by_aligners %>%
 filter(trimmed == "Trimmed") %>%
 group_by(cohort, sample_name, collection) %>%
 summarise(aligned_reads = mean(aligned_reads)) %>%
 ggplot(aes(cohort,
            aligned_reads,
            fill = cohort)) +
 geom_quasirandom(shape = 21,
                  colour = "black",
                  size = 2) +
 scale_fill_d3() +
 labs(x = "Cohort",
      y = "Aligned Reads") +
 theme(legend.position = "none") +
 ggtitle("Library Size")


# aligned reads - outliers ------------------------------------------------

t_aligned_reads_outliers_high <- alignment_by_aligners %>%
 group_by(cohort, sample_name, collection) %>%
 summarise(aligned_reads = mean(aligned_reads)) %>%
 group_by(cohort) %>%
 filter(outlier_iqr(aligned_reads, low_high = "high"))

outlier_tracking <-
 update_outlier_tracking(t_aligned_reads_outliers_high,
                         reason = "High aligned read count")

t_aligned_reads_outliers_low <- alignment_by_aligners %>%
 group_by(cohort, sample_name, collection) %>%
 summarise(aligned_reads = mean(aligned_reads)) %>%
 group_by(cohort) %>%
 filter(outlier_iqr(aligned_reads, low_high = "low"))

outlier_tracking <-
 update_outlier_tracking(t_aligned_reads_outliers_low,
                         reason = "Low aligned read count")

t_aligned_reads_outliers <- outlier_tracking %>%
  filter(str_detect(reason, "aligned")) %>%
  rename_with(remove_underscores) %>%
  rename_with(str_to_title)
  
# alignment feature -------------------------------------------------------

# multiqc_rnaseqc <-
#  read_delim(
#    "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/qc/multiqc/multiqc_data/multiqc_rna_seqc.txt",
#    delim = "\t"
#  ) %>%
#  mutate(fastq_id = str_extract(Sample, "[0-9|\\_]+(?=\\_csorted)")) %>%
#  inner_join(metadata) %>%
#  group_by(sample_name, cohort, fraction, region, age, genotype) %>% # average lanes
#  summarise(across(everything(), mean)) %>%
#  suppressMessages()

data <- multiqc_rnaseqc %>%
  pivot_longer(
    cols = c(
      "Exonic Rate",
      "Intronic Rate",
      "Intergenic Rate",
      "Ambiguous Alignment Rate"
    ),
    names_to = "alignment_type",
    values_to = "Proportion"
  ) %>%
  mutate(Percentage = Proportion * 100) 

data %>%
  group_by(alignment_type) %>%
  levene_test(formula = Percentage ~ fraction)

stat.test <- data %>%
  group_by(alignment_type) %>%
  wilcox_test(formula = Percentage ~ fraction) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  pval_asterisks(padj)

p_alignment_feature <- data %>%
  ggboxplot(x = "fraction", y = "Percentage", fill = "fraction", 
            facet.by = "alignment_type", scales = "free_y", 
            outlier.shape = NA) +
 geom_quasirandom(shape = 21,
                  colour = "black",
                  size = 0.5) +
 scale_fill_d3() +
 theme(legend.position = "none") +
 labs(x = "Fraction",
      y = "Percentage (%)") +
 scale_y_continuous(expand = c(0, 1)) +
  stat_pvalue_manual(stat.test, 
                     label = "ast", 
                     size = 5, 
                     y.position = c(7, 93, 7, 11))


multiqc_rnaseqc %>%
  pivot_longer(
    cols = c(
      "Exonic Rate",
      "Intronic Rate",
      "Intergenic Rate",
      "Ambiguous Alignment Rate"
    ),
    names_to = "Alignment Type",
    values_to = "Proportion"
  ) %>%
  # filter(cohort == "C1") %>%
  # ggplot(aes(x = fraction,
  #            y = Proportion)) +
  # geom_quasirandom() +
  # facet_wrap(vars(`Alignment Type`), scales = "free")
  group_by(cohort, `Alignment Type`) %>%
  do(tidy(lm(Proportion~fraction, data = .))) %>%
  filter(p.value < 0.05)







# salmon ------------------------------------------------------------------

# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/salmon"
files <-
 list.files(
   "input/data/salmon",
   pattern = "quant.sf",
   recursive =  T,
   full.names = T
 )

all(file.exists(files)) # everything exists

names(files) <-
 str_extract(files, "(?<=salmon\\/)[0-9|\\_]+(?=\\/quant.sf)") #name the files by the fastq_id from the enclosing folder

files <- files[names(files) %in% metadata$fastq_id] # remove blanks
files <-
 files[order(match(names(files), metadata$fastq_id))] #reorder based on metadata

metadata$files <- files
metadata$names <- metadata$fastq_id

if(!exists("se")){
  # se <- tximeta::tximeta(as.data.frame(metadata))
  se <- readRDS(file = "R/objects/se.rds")
}
# saveRDS(se, file = "R/objects/se.rds")

gse <- tximeta::summarizeToGene(se)

all(colnames(gse) == metadata$fastq_id)

dds_salmon <- DESeqDataSet(gse, design = ~ 1)

rownames(dds_salmon) <- gsub("\\..*$", "", rownames(dds_salmon))
rowData(dds_salmon)$gene_id <-
 gsub("\\..*$", "", rowData(dds_salmon)$gene_id)

# anno_salmon <- get_anno(rownames(dds_salmon))
anno_salmon <- readRDS("R/objects/anno_salmon.rds")

dds_salmon <- collapseReplicates(dds_salmon,
                                groupby = dds_salmon$sample_name,
                                run = dds_salmon$lane_code)

dds_salmon <- dds_salmon[rowSums(counts(dds_salmon)) > 0, ]

dds_salmon <- DESeq(
 dds_salmon,
 fitType = "local",
 sfType = "poscounts",
 parallel = TRUE,
 minReplicatesForReplace = Inf
)

# plotDispEsts(dds_salmon[rowData(dds_salmon)$dispGeneEst > 1e-2,])

vsd_salmon <- vst(dds_salmon)

# kallisto ----------------------------------------------------------------

# formerly /zfs/analysis/trap/active/snakemake/thesis_snakemake/input/index/references/annotation.gtf"
# txdb <- makeTxDbFromGFF(file = "input/index/references/annotation.gtf",
                        # dataSource="GENCODE", organism = "Mus musculus")
# txdb <- makeTxDbFromEnsembl(organism = "Mus musculus")
# tx2gene <- AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")

# formerly "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/kallisto"
files <- list.files("input/data/kallisto",
                    pattern = "abundance.h5",
                    recursive =  T,
                    full.names = T)

all(file.exists(files)) # everything exists

names(files) <- str_extract(files, "(?<=kallisto\\/)[0-9|\\_]+(?=\\/abundance.h5)") #name the files by the fastq_id from the enclosing folder

files <- files[names(files) %in% metadata$fastq_id] # remove blanks
files <- files[order(match(names(files), metadata$fastq_id))] #reorder based on metadata

metadata$files_kallisto <- files
metadata$names <- metadata$fastq_id

if(!exists("txi")){
  # txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
  txi <- readRDS(file = "R/objects/txi.rds")
}
# saveRDS(txi, file = "R/objects/txi.rds")


dds_kallisto <- DESeqDataSetFromTximport(txi,
                                         metadata,
                                         design = ~1)

rownames(dds_kallisto) <- gsub("\\..*$", "", rownames(dds_kallisto))

dds_kallisto <- collapseReplicates(dds_kallisto,
                                 groupby = dds_kallisto$sample_name,
                                 run = dds_kallisto$lane_code)

dds_kallisto <- dds_kallisto[rowSums(counts(dds_kallisto)) > 0,]

dds_kallisto <- DESeq(dds_kallisto,
                      fitType = "local",
                      sfType = "poscounts",
                      parallel = TRUE,
                      minReplicatesForReplace = Inf)

vsd_kallisto <- vst(dds_kallisto)

# STAR --------------------------------------------------------------------
# formerly /zfs/analysis/trap/active/snakemake/thesis_snakemake/output/star_pass_2/all_counts.txt"
star_counts0 <- vroom("input/data/star_pass_2/all_counts.txt",
                          delim = "\t")

star_counts <- star_counts0[,colnames(star_counts0) %in% metadata$fastq_id]
star_counts <- star_counts[,order(match(colnames(star_counts), metadata$fastq_id))]
star_counts <- as.matrix(star_counts)
rownames(star_counts) <- star_counts0$gene_id

all(colnames(star_counts) == metadata$fastq_id)

dds_star <- DESeqDataSetFromMatrix(star_counts,
                                   metadata,
                                   design = ~1)

rownames(dds_star) <- gsub("\\..*$", "", rownames(dds_star))

dds_star <- collapseReplicates(dds_star,
                                 groupby = dds_star$sample_name,
                                 run = dds_star$lane_code)
dim(dds_star)
dds_star <- dds_star[rowSums(counts(dds_star)) > 0,]

dds_star <- DESeq(dds_star,
                  fitType = "local",
                  sfType = "poscounts",
                  parallel = TRUE,
                  minReplicatesForReplace = Inf)

vsd_star <- vst(dds_star)

# aligner count comparison ------------------------------------------------

cts_aligners <- bind_rows(
  {
    counts(dds_salmon) %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      pivot_longer(-ensembl_gene_id, names_to = "sample_name", values_to = "count") %>%
      mutate(aligner = "salmon")
  },
  {
    counts(dds_kallisto) %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      pivot_longer(-ensembl_gene_id, names_to = "sample_name", values_to = "count") %>%
      mutate(aligner = "kallisto")
  },
  {
    counts(dds_star) %>%
      as_tibble(rownames = "ensembl_gene_id") %>%
      pivot_longer(-ensembl_gene_id, names_to = "sample_name", values_to = "count") %>%
      mutate(aligner = "star")
  }
) %>%
  inner_join(colData(dds_salmon), copy = TRUE) %>%
  group_by(ensembl_gene_id, aligner) %>%
  summarise(count = mean(count)) %>%
  pivot_wider(names_from = "aligner", values_from = "count")

p_cts_salmon_kallisto <- cts_aligners %>%
  ggplot(aes(x = log2(salmon+1),
             y = log2(kallisto+1))) +
  geom_point(colour = "#1f77b4", alpha = 0.25) +
  stat_cor() +
  labs(x = "Salmon Counts",
       y = "Kallisto Counts")

p_cts_salmon_star <- cts_aligners %>%
  # slice_sample(prop = 0.03) %>%
  ggplot(aes(x = log2(salmon+1),
             y = log2(star+1))) +
  geom_point(colour = "#ff7f0e", alpha = 0.25) +
  stat_cor() +
  labs(x = "Salmon Counts",
       y = "STAR Counts")

p_cts_kallisto_star <- cts_aligners %>%
  # slice_sample(prop = 0.03) %>%
  ggplot(aes(x = log2(kallisto+1),
             y = log2(star+1))) +
  geom_point(colour = "#2ca02c", alpha = 0.25) +
  stat_cor() +
  labs(x = "Kallisto Counts",
       y = "STAR Counts")


# dropout -----------------------------------------------------------------

# extract counts (as I could not set 0 values to NA in the DESeq object)
dropouts <- counts(dds_salmon) %>%
  as_tibble(rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample_name", values_to = "count") %>%
  mutate(count = na_if(count, 0)) %>%
  inner_join(colData(dds_salmon), copy = TRUE) %>%
  group_by(cohort, fraction, gene_id) %>%
  summarise(mean = mean(count, na.rm = T),
            dropouts = sum(is.na(count)),
            dropout_p = sum(is.na(count))/n()) %>%
  group_by(cohort, fraction) %>%
  mutate(gene_rank = rank(-mean),
         fraction = ifelse(fraction == "IP", "TRAP", "TOTAL"))

p_dropouts <- dropouts %>%
  ggplot(aes(x = mean,
             y = dropout_p,
             colour = fraction)) +
  geom_smooth() +
  scale_color_d3() +
  facet_wrap(vars(cohort), nrow = 3) +
  scale_x_log10(limits = c(1, 1e4)) +
  labs(x = "Mean Gene Count",
       y = "Dropout Proportion",
       colour = "Fraction") +
  theme(legend.position = "top") +
  panel_border()

# library sizes -----------------------------------------------------------
p_library_size <- colSums(counts(dds_salmon)) %>%
 as_tibble(rownames = "sample_name") %>%
 inner_join(colData(dds_salmon), copy = TRUE) %>%
 ggplot(aes(x = cohort,
            y = value,
            fill = cohort)) +
 geom_quasirandom(shape = 21,
                  size = 2,
                  colour = "black") +
 scale_y_continuous(breaks = seq(1e7, 5e7, 1e7)) +
 scale_fill_d3() +
 labs(x = "Cohort",
      y = "Number of Reads") +
 theme(legend.position = "none")

# number of genes detected with all reads - fraction labelled ---------------------------------
p_n_genes_detected_cohort_fraction <-
 colSums(counts(dds_salmon) > 0) %>%
 as_tibble(rownames = "sample_name") %>%
 rename(n_genes = value) %>%
 inner_join(colData(dds_salmon), copy = TRUE) %>%
 ggplot(aes(x = cohort,
            y = n_genes,
            fill = fraction)) +
 geom_quasirandom(shape = 21,
                  size = 2,
                  colour = "black") +
 scale_fill_d3() +
 labs(x = "Cohort", y = "Number of Genes Detected", fill = "Fraction") +
 theme(legend.position = "top")

# number of genes detected in IP by cohort and region ---------------------
# NB. in C3 there is no difference in n genes detected between regions
# even though RNA yield of striatal samples was low

p_n_genes_detected_cohort_region_ip <-
 colSums(counts(dds_salmon) > 0) %>%
 as_tibble(rownames = "sample_name") %>%
 rename(n_genes = value) %>%
 inner_join(colData(dds_salmon), copy = TRUE) %>%
 filter(fraction == "IP") %>%
 ggplot(aes(x = region,
            y = n_genes,
            fill = region)) +
 geom_quasirandom(shape = 21,
                  size = 2,
                  colour = "black") +
 scale_fill_d3() +
 labs(x = "region", y = "Number of Genes Detected", fill = "Region") +
 theme(legend.position = "top") +
 facet_wrap(vars(cohort)) +
 stat_compare_means(label = "p.signif",
                    label.x = 2)

p_n_genes_detected_cohort_region_ip
# or age/genotype

p_n_genes_detected_collection_ip <-
  colSums(counts(dds_salmon) > 0) %>%
  as_tibble(rownames = "sample_name") %>%
  rename(n_genes = value) %>%
  inner_join(colData(dds_salmon), copy = TRUE) %>%
  filter(cohort == "C3" &
           fraction == "IP") %>%
  mutate(
    # outlier = str_detect(collection, c("C3.3|C3.8")),
    collection = factor(collection, levels = paste0("C3.", seq(1, 10, 1)))
  ) %>%
  ggplot(aes(x = n_genes,
             y = collection,
             fill = collection)) +
  geom_density_ridges() +
  scale_fill_d3() +
 labs(x = "Number of Genes Detected", y = "Collection Batch") +
 theme(legend.position = "none")

p_n_genes_detected_collection_ip

# library size vs n genes detected ----------------------------------------

p_lib_size_vs_n_genes <-
 tibble(
   "sample_name" = colnames(dds_salmon),
   "library_size" = colSums(counts(dds_salmon)),
   n_genes_detected = colSums(counts(dds_salmon) > 0)
 ) %>%
 inner_join(colData(dds_salmon), copy = TRUE) %>%
 # filter(sample_name != "Vwithin) %>%
 ggplot(aes(x = library_size,
            y = n_genes_detected,
            fill = fraction,
            colour = fraction)) +
 geom_point(shape = 21,
            size = 2,
            colour = "black") +
  geom_smooth(method = "lm") +
 # stat_cor(method = "pearson", label.y = c(2500, 0)) +
 facet_wrap(vars(cohort), scales = "free_x", nrow = 3) +
 labs(
   x = "Library Size",
   y = "Number of Genes Detected",
   fill = "Fraction",
   colour = "Fraction"
 ) +
 scale_fill_d3() +
  scale_color_d3() +
  # scale_y_continuous(limits = c(0, NA)) +
  theme(legend.position = "top")

p_lib_size_vs_n_genes

# downsampling by fraction ------------------------------------------------
downsampled <- readRDS("R/objects/downsampled.rds")

if(!exists("downsampled")){
  depths <- seq(1e6, 11e6, 1e6)
  downsampled <- lapply(depths, function(depth){
    prop_vector <- depth/colSums(counts(dds_salmon))
    # print(depth)
    cts_downsampled <- downsampleMatrix(counts(dds_salmon), prop = prop_vector, bycol = TRUE)
    output <- colSums(cts_downsampled > 10) %>%
      as_tibble(rownames = "sample_name") %>%
      rename(value = "n_genes") %>%
      mutate(depth = depth)
    return(output)
  })
}

# p_downsampling_fraction <- downsampled %>%
#  bind_rows %>%
#  inner_join(colData(dds_salmon), copy = TRUE) %>%
#   filter(depth != 1.1e7) %>%
#  # filter(!collection %in% c("C3.3", "C3.8")) %>%
#  ggplot(aes(x = depth,
#             y = n_genes)) +
#  # geom_line(aes(group = sample_name,
#  #               colour = fraction),
#  #           alpha = 0.5) +
#  geom_smooth(aes(colour = fraction)) +
#  facet_wrap(vars(cohort), nrow = 3) +
#  scale_color_d3() +
#  scale_x_continuous(breaks = seq(1e6, 11e6, 2e6)) +
#  labs(x = "Library Size",
#       y = "Number of Genes Detected",
#       colour = "Fraction") +
#  theme(legend.position = "top",
#        axis.title.y = element_blank())

# show the gain in complexity per change in depth
p_downsampling_fraction <-
downsampled %>%
  bind_rows %>%
  inner_join(colData(dds_salmon), copy = TRUE) %>%
  # filter(depth %in% c(1e7, 9e6)) %>%
  select(sample_name, cohort, fraction, depth, n_genes) %>%
  group_by(sample_name) %>%
  arrange(depth, .by_group = TRUE) %>%
  mutate(diff = n_genes - lag(n_genes, default = dplyr::first(n_genes))) %>%
  filter(diff != 0) %>%
  ggplot(aes(x = depth,
             y = diff,
             colour = fraction)) +
  geom_smooth() +
  scale_color_d3() +
  facet_wrap(vars(cohort), nrow = 3) +
  labs(x = "Library Size",
       y = "Gain in Genes Detected") +
  scale_x_continuous(breaks = seq(2e6, 1e7, 2e6)) +
  theme(legend.position = "top") +
  labs(colour = "Fraction")


# pivot_wider(names_from = "depth",
#               values_from = "n_genes") %>%
#   mutate(diff = `1e+07` - `9e+06`) %>%
#   group_by(fraction) %>%
#   do(tidy(lm(diff~cohort, data = .)))

# # diversity - WIP --------------------------------------------------------------
#
# cts_salmon <- counts(dds_salmon) %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   pivot_longer(-ensembl_gene_id,
#                names_to = "sample_name",
#                values_to = "count") %>%
#   inner_join(colData(dds_salmon), copy = TRUE)
#
# prop_vector <- 1e7/colSums(counts(dds_salmon))
# cts_downsampled <- downsampleMatrix(counts(dds_salmon), prop = prop_vector, bycol = TRUE)
# cts_downsampled <- cts_downsampled %>%
#   as.matrix %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   pivot_longer(-ensembl_gene_id,
#                names_to = "sample_name",
#                values_to = "count") %>%
#   inner_join(colData(dds_salmon), copy = TRUE)
#
# diversity <- cts_downsampled %>%
#   group_by(sample_name) %>%
#   arrange(sample_name, desc(count)) %>%
#   mutate(n_genes = row_number(),
#          cum_sum = cumsum(count),
#          cum_prop = cum_sum/sum(count))
#
# diversity_cohort <- diversity %>%
#   group_by(cohort, n_genes) %>%
#   summarise(cum_prop = mean(cum_prop))
#
# diversity_cohort_fraction <- diversity %>%
#   group_by(cohort, fraction, n_genes) %>%
#   summarise(cum_prop = mean(cum_prop))
#
# p_diversity_cohort_fraction <- diversity_cohort_fraction %>%
#   ungroup() %>%
#   ggplot(aes(x = n_genes,
#              y = cum_prop,
#              colour = cohort)) +
#   geom_line() +
#   scale_x_log10() +
#   scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#   facet_wrap(vars(fraction), nrow = 3) +
#   scale_color_d3()
#
# p_diversity_cohort_fraction
#
# p_diversity_individual <- diversity %>%
#   filter(fraction == "IP") %>%
#   ggplot(aes(x = n_genes,
#              y = cum_prop,
#              colour = region)) +
#   geom_line(aes(group = sample_name)) +
#   scale_x_log10() +
#   scale_y_continuous(breaks = seq(0, 1, 0.2)) +
#   facet_wrap(vars(cohort, region), nrow = 3) +
#   scale_color_d3()
#
# p_diversity_individual
#
# # counts(dds_salmon) %>%
# #   as_tibble(rownames = "ensembl_gene_id") %>%
# #   pivot_longer(-ensembl_gene_id,
# #                names_to = "sample_name",
# #                values_to = "count") %>%
# #   inner_join(colData(dds_salmon), copy = TRUE) %>%
# #   inner_join(anno_salmon) %>%
# #   filter(external_gene_name == 'Th') %>% View
#
#
# # # count distribution - cohxort and fraction - WIP --------------------------------
# #
# # counts(dds_salmon) %>%
# #   as_tibble(rownames = "gene_id") %>%
# #   pivot_longer(-gene_id,
# #                names_to = "sample_name",
# #                values_to = "count") %>%
# #   inner_join(colData(dds_salmon), copy = TRUE) %>%
# #   # filter(fraction == "IP") %>%
# #   ggplot(aes(x = count,
# #              colour = fraction)) +
# #   geom_density() +
# #   scale_x_log10() +
# #   facet_wrap(vars(cohort), nrow = 3)
# #
# #
# # normalisation factors - cohort density - WIP ----------------------------------
#
#
# as_tibble(normalizationFactors(dds_salmon), rownames = "gene_id") %>%
#   pivot_longer(-gene_id, names_to = "sample_name", values_to = "normFactor") %>%
#   inner_join(colData(dds_salmon), copy = TRUE) %>%
#   ggplot(aes(x = normFactor,
#              fill = cohort)) +
#   geom_density() +
#   scale_x_continuous(limits = c(0, 2)) +
#   facet_wrap(vars(cohort),
#              nrow = 3) +
#   scale_fill_d3()
#
#
# as_tibble(normalizationFactors(dds_salmon), rownames = "gene_id") %>%
#   pivot_longer(-gene_id, names_to = "sample_name", values_to = "normFactor") %>%
#   inner_join(colData(dds_salmon), copy = TRUE) %>%
#   inner_join(rowData(dds_salmon), copy = TRUE) %>%
#   group_by(cohort, gene_id) %>%
#   summarise(baseMean = mean(baseMean),
#             normFactor = mean(normFactor)) %>%
#   slice_sample(prop = 0.1) %>%
#   ggplot(aes(x = baseMean,
#              y = normFactor)) +
#   geom_point() +
#   scale_y_log10() +
#   facet_wrap(vars(cohort)) +
#   scale_x_continuous(limits = c(0, 1e4))
#
#   group_by(sample_name) %>%
#   summarise(normFactor = mean(normFactor)) %>%
#   inner_join(colData(dds_salmon), copy = TRUE) %>%
#   ggplot(aes(x = cohort,
#              y = normFactor,
#              colour = fraction)) +
#   geom_quasirandom()
#



# gene body coverage ------------------------------------------------------

# bams <-
#   list.files(
#     "/zfs/analysis/trap/active/snakemake/thesis_snakemake/output/star_pass_2",
#     pattern = "csorted.bam$",
#     recursive =  T,
#     full.names = T
#   )

# txdb <- makeTxDbFromGFF("/zfs/analysis/trap/active/references/gencode_m25_data/gencode.vM25.annotation.gtf.gz")

# bfs <- Rsamtools::BamFileList(bams)
# bf <- Rsamtools::BamFile(bam)
# seqinfo(bfs)

# top_gene_ids <- rowSums(counts(dds_salmon)) %>%
#   sort(decreasing = TRUE) %>%
#   head(100) %>%
#   as_tibble(rownames = "gene_id") %>%
#   pull(gene_id)

# genes <- genes(txdb)[gsub("\\..*$", "", genes(txdb)$gene_id) %in% top_gene_ids,]

# covSigs <- sapply(bams, function(x){bamCoverage(x, genes, verbose = TRUE)})
# saveRDS(covSigs, file = "R/objects/covSigs.rds")

# covSigs <- lapply(covSigs, function(w) {
#   lapply(w, function(x) {
#     x[sum(x) > length(x) & length(x) >= 100]
#   })
# })

# covSigs <- readRDS("R/objects/covSigs.rds")
# covSigs <- lapply(covSigs, as.list)




# covSigs_summarised <- lapply(covSigs, function(x){
#   sapply(x, function(y){
#     y <- split(y, cut(seq_along(y), 100, labels = FALSE))
#     sapply(y, mean)
#   }) %>%
#     as.matrix %>%
#     apply(., 2, rescale) %>%
#     apply(., 1, mean) %>%
#     rescale
# })
#
# names(covSigs_summarised) <- str_extract(names(covSigs_summarised), "(?<=star_pass_2\\/\\/)[0-9|\\_]+(?=_csorted.bam)")

covSigs_summarised <- readRDS("R/objects/covSigs_summarised.rds")

bind_rows(covSigs_summarised) %>%
  mutate(.before = everything(), fastq_id = names(covSigs_summarised)) %>%
  pivot_longer(-fastq_id, names_to = "index", values_to = "proportion") %>%
  inner_join(metadata) %>%
  group_by(sample_name, cohort, index) %>%
  summarise(proportion = mean(proportion)) %>%
  mutate(index = as.numeric(index)) %>%
  ggplot(aes(x = index,
             y = proportion,
             colour = cohort)) +
  geom_smooth() +
  scale_color_d3()

# rm(covSigs)
# saveRDS(covSigs_summarised,
#         "R/objects/covSigs_summarised.rds")


# select aligner ----------------------------------------------------------

dds <- dds_salmon
colData(dds)$collection_outlier <- str_detect(colnames(dds), "C3.3|C3.8")


# C1 all ------------------------------------------------------------------

dds_C1 <- dds[,str_detect(colnames(dds), "C1")]
dds_C1 <- DESeq(dds_C1,
                # fitType = "local",
                # sfType = "poscounts",
                parallel = TRUE,
                minReplicatesForReplace = Inf)
vsd_C1 <- vst(dds_C1)
pca_C1 <- make_pca(vsd_C1, dds_object = dds_C1)

p_pca_c1_fraction <- pca_C1$x %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
  ggplot(aes(PC1, PC2, colour = fraction)) +
  geom_point() +
  scale_color_d3() +
  labs(colour = "Fraction") +
  theme(legend.position = "top") +
  # coord_fixed(ratio = pca_C1$var_explained$var[2]*6/pca_C1$var_explained$var[1]) +
  labs(x = paste("PC1:", pca_C1$var_explained$var[1], "%"),
       y = paste("PC2:", pca_C1$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))

p_pca_c1_age <- pca_C1$x %>%
  ggplot(aes(PC1, PC2, colour = age)) +
  geom_point() +
  scale_color_d3() +
  labs(colour = "Age") +
  theme(legend.position = "top") +
  # coord_fixed(ratio = pca_C1$var_explained$var[2]*6/pca_C1$var_explained$var[1]) +
  labs(x = paste("PC1:", pca_C1$var_explained$var[1], "%"),
       y = paste("PC2:", pca_C1$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))


# C1 ip -------------------------------------------------------------------
dds_C1_ip <- dds[,str_detect(colnames(dds), "C1") & !str_detect(colnames(dds), "TOTAL") &
                   !str_detect(colnames(dds), "POOL")]
dds_C1_ip <- DESeq(dds_C1_ip,
                   # fitType = "local",
                   # sfType = "poscounts",
                   parallel = TRUE,
                   minReplicatesForReplace = Inf)
vsd_C1_ip <- vst(dds_C1_ip)
pca_C1_ip <- make_pca(vsd_C1_ip, dds_object = dds_C1_ip)


# C2 all ------------------------------------------------------------------

dds_C2 <- dds[,str_detect(colnames(dds), "C2")]
dds_C2 <- DESeq(dds_C2,
                # fitType = "local",
                # sfType = "poscounts",
                parallel = TRUE,
                minReplicatesForReplace = Inf)
vsd_C2 <- vst(dds_C2)
pca_C2 <- make_pca(vsd_C2, dds_object = dds_C2)

p_pca_c2_region <- pca_C2$x %>%
  ggplot(aes(PC1, PC2, colour = region)) +
  geom_point() +
  geom_point(data = subset(pca_C2$x, str_detect(sample_name, "C2.3") & PC1 > -20 &
                             str_detect(sample_name, "DS|VS")),
             aes(PC1, PC2),
             shape = 21, colour = "black", stroke = 1) +
  scale_color_d3() +
  labs(colour = "Region") +
  theme(legend.position = "top") +
  # coord_fixed(ratio = pca_C2$var_explained$var[2]/pca_C2$var_explained$var[1]) +
  labs(x = paste("PC1:", pca_C2$var_explained$var[1], "%"),
       y = paste("PC2:", pca_C2$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))

p_pca_c2_fraction <- pca_C2$x %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
  ggplot(aes(PC1, PC2, colour = fraction)) +
  geom_point() +
  geom_point(data = subset(pca_C2$x, str_detect(sample_name, "C2.3") & PC1 > -20 &
                             str_detect(sample_name, "DS|VS")),
             aes(PC1, PC2),
             shape = 21, colour = "black", stroke = 1) +
  scale_color_d3() +
  labs(colour = "Fraction") +
  theme(legend.position = "top") +
  # coord_fixed(ratio = pca_C2$var_explained$var[2]/pca_C2$var_explained$var[1]) +
  labs(x = paste("PC1:", pca_C2$var_explained$var[1], "%"),
       y = paste("PC2:", pca_C2$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))

# C2 ip -------------------------------------------------------------------

dds_C2_ip <- dds[,str_detect(colnames(dds), "C2") & !str_detect(colnames(dds), "TOTAL") &
                   !str_detect(colnames(dds), "POOL")]
dds_C2_ip <- DESeq(dds_C2_ip,
                   # fitType = "local",
                   # sfType = "poscounts",
                   parallel = TRUE,
                   minReplicatesForReplace = Inf)
vsd_C2_ip <- vst(dds_C2_ip)
pca_C2_ip <- make_pca(vsd_C2_ip, dds_object = dds_C2_ip)

# C3 all ------------------------------------------------------------------

dds_C3 <- dds[,str_detect(colnames(dds), "C3")]
dds_C3 <- DESeq(dds_C3,
                fitType = "local",
                sfType = "poscounts",
                parallel = TRUE,
                minReplicatesForReplace = Inf)
vsd_C3 <- vst(dds_C3)
pca_C3 <- make_pca(vsd_C3, dds_object = dds_C3)

p_pca_c3_region <- pca_C3$x %>%
  ggplot(aes(PC1, PC2, colour = region)) +
  geom_point() +
  geom_point(data = subset(pca_C3$x, str_detect(sample_name, "C3.3|C3.8")),
             aes(PC1, PC2),
             shape = 21, colour = "black", stroke = 1) +
  scale_color_d3() +
  labs(colour = "Region") +
  theme(legend.position = "top") +
  # coord_fixed(ratio = pca_C3$var_explained$var[2]/pca_C3$var_explained$var[1]) +
  labs(x = paste("PC1:", pca_C3$var_explained$var[1], "%"),
       y = paste("PC2:", pca_C3$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))

p_pca_c3_fraction <- pca_C3$x %>%
  mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
  ggplot(aes(PC1, PC2, colour = fraction)) +
  geom_point() +
  geom_point(data = subset(pca_C3$x, str_detect(sample_name, "C3.3|C3.8")),
             aes(PC1, PC2),
             shape = 21, colour = "black", stroke = 1) +
  scale_color_d3() +
  labs(colour = "Fraction") +
  theme(legend.position = "top") +
  # coord_fixed(ratio = pca_C3$var_explained$var[2]/pca_C3$var_explained$var[1]) +
  labs(x = paste("PC1:", pca_C3$var_explained$var[1], "%"),
       y = paste("PC2:", pca_C3$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))

pca_C3$x %>%
  pivot_longer(paste0("PC", seq(1, 8, 1)), names_to = "PC", values_to = "coord") %>%
  group_by(PC, fraction, compartment) %>%
  do(tidy(lm(coord~age, data = .))) %>%
  filter(p.value < 0.05 & !str_detect(term, "Intercept"))

# C3 ip -------------------------------------------------------------------
dds_C3_ip <- dds[,str_detect(colnames(dds), "C3") & !str_detect(colnames(dds), "TOTAL") &
                   !str_detect(colnames(dds), "POOL")]
dds_C3_ip <- DESeq(dds_C3_ip,
                   fitType = "local",
                   sfType = "poscounts",
                   parallel = TRUE,
                   minReplicatesForReplace = Inf)
vsd_C3_ip <- vst(dds_C3_ip)
pca_C3_ip <- make_pca(vsd_C3_ip, dds_object = dds_C3_ip)

# pca_C3_ip$x %>%
#   pivot_longer(paste0("PC", seq(1, 8, 1)), names_to = "PC", values_to = "coord") %>%
#   ggplot(aes(PC, coord, colour = age)) +
#   geom_quasirandom() +
#   scale_color_d3() +
#   facet_wrap(vars(compartment, fraction))

# pca_C3_ip$x %>%
#   pivot_longer(paste0("PC", seq(1, 8, 1)), names_to = "PC", values_to = "coord") %>%
#   group_by(PC, compartment) %>%
#   do(tidy(lm(coord~collection+age, data = .))) %>%
#   filter(p.value < 0.05 & !str_detect(term, "Intercept")) %>% View

# pca_C3$x %>%
#   pivot_longer(paste0("PC", seq(1, 8, 1)), names_to = "PC", values_to = "coord") %>%
#   ggplot(aes(PC, coord, colour = genotype)) +
#   geom_quasirandom() +
#   scale_color_d3() +
#   facet_wrap(vars(region, fraction))






# Dropout further ----

# N of genes with dropout per cohort
plot_data <- bind_rows(
  {
    enframe(rowMeans(counts(dds_C1_ip))) %>%
      mutate(entire = rowSums(counts(dds_C1_ip) == 0) == 0,
             cohort = "C1")
  },
  {
    enframe(rowMeans(counts(dds_C2_ip))) %>%
      mutate(entire = rowSums(counts(dds_C2_ip) == 0) == 0,
             cohort = "C2")
  },
  {
    enframe(rowMeans(counts(dds_C3_ip))) %>%
      mutate(entire = rowSums(counts(dds_C3_ip) == 0) == 0,
             cohort = "C3")
  }
) %>%
  filter(value != 0) %>%
  mutate(dropouts = ifelse(entire, "No Dropouts",
                           "Dropouts")) %>%
  group_by(cohort, dropouts) %>%
  summarise(n = n()) %>%
  ungroup()

p_dropouts_n_genes <- ggplot(plot_data,
                             aes(x = cohort,
                                 y = n,
                                 fill = dropouts)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  stat_identity(geom = "text", colour = "white", size = 3.5,
                fontface = 2,
                aes(label = ifelse(n > 2500, n, "")),
                position=position_stack(vjust=0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "top") +
  labs(x = "Cohort",
       y = "Percentage of Genes",
       fill = "Dropout \nGenes")

# Cohort dropout library proportion
plot_data <- bind_rows(
  {
    enframe(rowMeans(counts(dds_C1_ip))) %>%
      mutate(entire = rowSums(counts(dds_C1_ip) == 0) == 0,
             cohort = "C1")
  },
  {
    enframe(rowMeans(counts(dds_C2_ip))) %>%
      mutate(entire = rowSums(counts(dds_C2_ip) == 0) == 0,
             cohort = "C2")
  },
  {
    enframe(rowMeans(counts(dds_C3_ip))) %>%
      mutate(entire = rowSums(counts(dds_C3_ip) == 0) == 0,
             cohort = "C3")
  }
) %>%
  group_by(cohort, entire) %>%
  summarise(count = sum(value)) %>%
  mutate(prop = count/sum(count),
         dropouts = ifelse(entire, "No Dropouts",
                           "Dropouts"))

p_dropout_cohort_percentage <- ggplot(plot_data,
       aes(x = cohort,
           y = prop,
           fill = dropouts)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  stat_identity(data = subset(plot_data, prop > 0.1),
                geom = "text", colour = "white", size = 3.5,
                fontface = 2,
                aes(label = scales::percent(prop)),
                position=position_stack(vjust=0.5)) +
  labs(x = "Cohort",
       y = "Percentage of Library",
       fill = "Genes with \nDropouts") +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "top")


# % of library representing genes with dropout in C3
entire_genes_C3 <- rownames(dds_C3_ip)[rowSums(counts(dds_C3_ip) == 0) == 0]

dropout_genes_C3 <- rownames(dds_C3_ip)[rowSums(counts(dds_C3_ip) == 0) != 0]

plot_data <- bind_rows(
  {
    enframe(rowMeans(counts(dds_C3_ip))) %>%
      mutate(entire = name %in% entire_genes_C3,
             cohort = "C3")
  },
  {
    enframe(rowMeans(counts(dds_C2_ip))) %>%
      mutate(entire = name %in% entire_genes_C3,
             cohort = "C2")
  },
  {
    enframe(rowMeans(counts(dds_C1_ip))) %>%
      mutate(entire = name %in% entire_genes_C3,
             cohort = "C1")
  }
) %>%
  group_by(cohort, entire) %>%
  summarise(count = sum(value)) %>%
  mutate(prop = count/sum(count),
         gene_type = ifelse(entire, "No Dropouts",
                            "Dropouts"))

p_dropout_c3_composition <- ggplot(plot_data,
         aes(x = cohort,
             y = prop,
             fill = gene_type)) +
  geom_col(position = "fill",
           colour = "black") +
  scale_fill_d3() +
  labs(x = "Cohort",
       y = "Percentage of Library",
       fill = "C3 Dropout \nGenes") +
  stat_identity(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = scales::percent(prop)),
             position=position_stack(vjust=0.5)) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "top")




# variation explained all ------------------------------------------------
var_explained_combined <- bind_rows(
  {pca_C3$var_explained %>%
      filter(PC %in% paste0("PC", seq(1, 8, 1))) %>%
      mutate(cohort = "C3",
             label = c("Region & Soma Age IP", "Fraction", "Soma Genotype TOTAL", "Axon Age IP", NA, "Soma Age IP & TOTAL", NA, NA))},
  {pca_C2$var_explained %>%
      filter(PC %in% paste0("PC", seq(1, 8, 1))) %>%
      mutate(cohort = "C2",
             label = c("Region", "Fraction & Axon Age IP", rep(NA, 5), "Axon Age IP"))},
  {pca_C1$var_explained %>%
      filter(PC %in% paste0("PC", seq(1, 8, 1))) %>%
      mutate(cohort = "C1",
             label = c("Fraction", "Age", NA, NA, "Soma Genotype IP", NA, NA, NA))}
)

p_var_explained <- var_explained_combined %>%
  ggplot(aes(PC, var, fill = PC, label = label)) +
  geom_col(colour = "black") +
  geom_text(angle = 90, hjust = -0.1, fontface = "bold") +
  scale_fill_d3() +
  facet_wrap(vars(cohort), nrow = 1) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Principal Component",
       y = "Variation Explained") +
  scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55, 5),
                     expand = expansion(mult = c(0, 0.1))) +
  panel_border()

# pca_C2$x %>%
#   pivot_longer(paste0("PC", seq(1, 8, 1)), names_to = "PC", values_to = "coord") %>%
#   group_by(PC, region, fraction) %>%
#   do(tidy(lm(coord~age, data = .))) %>%
#   filter(p.value < 0.05 & !str_detect(term, "Intercept")) %>% View


# cohorts combined --------------------------------------------------------

vsd <- vst(dds)
pca <- make_pca(vsd)

n_cohort_var_percent <- pca$var_explained$var[1]
n_region_var_percent <- pca$var_explained$var[2]

# everything
p_pca_cohort <- pca$x %>%
  ggplot(aes(PC1, PC2, colour = cohort)) +
  geom_point() +
  scale_color_d3() +
  labs(colour = "Cohort") +
  theme(legend.position = "top") +
  labs(x = paste("PC1:", pca$var_explained$var[1], "%"),
       y = paste("PC2:", pca$var_explained$var[2], "%")) +
  theme(legend.title = element_text(face = "bold"))

p_pca_cohort

p_pca_fraction <- plotPCA(vsd, intgroup = "fraction", ntop = Inf, returnData = TRUE) %>%
  ggplot(aes(PC1, PC2, colour = group)) +
  geom_point() +
  scale_color_d3() +
  labs(colour = "Fraction",
       x = paste("PC1:", pca$var_explained$var[1], "%"),
       y = paste("PC2:", pca$var_explained$var[2], "%")) +
  theme(legend.position = "top") +
  theme(legend.title = element_text(face = "bold"))

p_pca_fraction

p_pca_region <- plotPCA(vsd, intgroup = "region", ntop = Inf, returnData = TRUE) %>%
  ggplot(aes(PC1, PC2, colour = group)) +
  geom_point() +
  scale_color_d3() +
  labs(colour = "Region",
       x = paste("PC1:", pca$var_explained$var[1], "%"),
       y = paste("PC2:", pca$var_explained$var[2], "%")) +
  theme(legend.position = "top") +
  theme(legend.title = element_text(face = "bold"))

p_pca_region


p_pca_collection_outlier <- plotPCA(vsd, intgroup = "collection_outlier",
                                    ntop = Inf, returnData = TRUE) %>%
  ggplot(aes(PC1, PC2, colour = group)) +
  geom_point() +
  scale_color_d3() +
  labs(colour = "Collection C3.3 & C3.8",
       x = paste("PC1:", pca$var_explained$var[1], "%"),
       y = paste("PC2:", pca$var_explained$var[2], "%")) +
  theme(legend.position = "top")

p_pca_collection_outlier




# remove large objects ----
# rm(gse)
# rm(txi)
# rm(se)
# rm(list = ls(pattern = "dds"))
# rm(list = ls(pattern = "vsd"))
# rm(list = ls(pattern = "pca_C"))
# rm(pca)
# ----
# ----
# ----

# Save output ----

# minor edits
n_alignment_percentage_cohort_mean <- signif(n_alignment_percentage_cohort_mean, 3)

save(list = c(ls(pattern = "^p_"),
              ls(pattern = "^t_"),
              ls(pattern = "^n_"), 
              ls(pattern = "^outlier_tr")),
     file = "R/objects/03-results-objects.RData")

save.image("R/objects/RESULTS_1.RData")

