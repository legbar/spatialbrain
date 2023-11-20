# SETUP ----
# set seed for reproducibility
set.seed(821196)
# load packages and save package versions
source("R/renv.R")
# load sample metadata
source("R/metadata.R")
# load outlier functions
#Interquartile range outlier detection function (IROF)
outlier_iqr <- function(column, low_high = "both"){
  if(low_high == "both"){
    ifelse(column < quantile(column, 0.25) - 1.5*IQR(column) |
             column > quantile(column, 0.75) + 1.5*IQR(column),
           TRUE, FALSE)
  } else if(low_high == "low"){
    ifelse(column < quantile(column, 0.25) - 1.5*IQR(column),
           TRUE, FALSE)
  } else if(low_high == "high"){
    ifelse(column > quantile(column, 0.75) + 1.5*IQR(column),
           TRUE, FALSE)
  }
}

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
# RESULTS 2 ----
# ----
# ----
# ----
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

# # edgeR object ----
# y <- DGEList(counts = counts(dds_C1_MB_IP)[rownames(dds_C1_MB_IP) %in% MB_FRACTION_ENRICHED_GENES,],
#              group = colData(dds_C1_MB_IP)$age)
# y <- DGEList(counts = counts(dds_C1_MB_IP),
#              group = colData(dds_C1_MB_IP)$age)
# nrow(y)
# keep <- filterByExpr(y)
# sum(keep)
# y <- y[keep, , keep.lib.sizes = F]
# y <- calcNormFactors(y)
# y$samples$norm.factors
# sum(keep)
# design <- model.matrix(~y$samples$group)
#
# y <- estimateDisp(y, design)
# fit <- glmQLFit(y, design)
# qlf <- glmQLFTest(fit, coef = 2)
# topTags(qlf)
# qlf$table %>% as_tibble(rownames = "ensembl_gene_id") %>% left_join(anno) %>%
#   mutate(padj = p.adjust(PValue, method = "fdr")) %>%
#   filter(padj < 0.01) %>% View

# low count filter function -------------------------------------------------------

# Don't apply this filter to the DDS before running a FRACTION comparison!
# It will select for the enriched genes, particularly in C3
# and I think that skews the normalisation process, so enriched genes end up being classed as depleted!
# Filter after running the FRACTION comparison on the entire dataset.
#
#
# creating a list of genes that are detected in at least one condition
# allowing the definition of the proportion of samples that must express it (prop)
# per condition and the minimum number of counts (min)



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

# create dds objects for each cohort and fraction ---------------------------------------------------------------------

# core dds
dds <- dds_salmon
rm(dds_salmon)

# relevel fraction
colData(dds)$fraction <-
  relevel(colData(dds)$fraction, ref = "TOTAL")

# remove rRNA 106106
dds <- dds[rownames(dds) != "ENSMUSG00000106106", ]

# get anno
# anno <- get_anno(rownames(dds))
anno <- readRDS("R/objects/anno.rds")

# get biotype anno
# anno_biotype <- getBM(
#   attributes = c(
#     "ensembl_gene_id",
#     "external_gene_name",
#     "description",
#     "gene_biotype"
#   ),
#   filters = "ensembl_gene_id",
#   values = rownames(dds),
#   mart = ensembl,
#   uniqueRows = TRUE
# )

# remove pooled C1 IP samples
dds <- dds[, colData(dds)$collection != "C1.IPPOOL"]

# test: remove C3.3 and C3.8
dds <- dds[, !colData(dds)$collection %in% c("C3.3", "C3.8")]

# remove cohort 2 outlier samples (based on PCA and Cre count in AXON samples)
dds <- dds[, !colData(dds)$sample_name %in% c("DS_OLD_WT_MIXED_C2.3",
                                              "VS_YOUNG_WT_MIXED_C2.3",
                                              "MB_OLD_WT_MIXED_C2.6",
                                              "MB_YOUNG_WT_MIXED_C2.7",
                                              "DS_OLD_WT_MIXED_C2.TOTAL1-3")]

# filter zeros function
filter_zeros <- function(dds) {
  dds <- dds[rowSums(counts(dds)) > 0,]
}

# MB subset
dds_MB <- dds[, colData(dds)$compartment == "MB"] %>% filter_zeros()
# AXON subset
dds_AXON <-
  dds[, colData(dds)$compartment == "AXON"] %>% filter_zeros()

# IP subset
# dds_IP <- dds[, colData(dds)$fraction == "IP"] %>% filter_zeros()
# dds_MB_IP <- dds_IP[, colData(dds_IP)$compartment == "MB"] %>% filter_zeros()

# create cohort subsets
dds_C1 <- dds[, colData(dds)$cohort == "C1"] %>% filter_zeros()
dds_C2 <- dds[, colData(dds)$cohort == "C2"] %>% filter_zeros()
dds_C3 <- dds[, colData(dds)$cohort == "C3"] %>% filter_zeros()

# create cohort compartment subsets
dds_C1_MB <-
  dds_C1[, colData(dds_C1)$compartment == "MB"] %>% filter_zeros()
dds_C2_MB <-
  dds_C2[, colData(dds_C2)$compartment == "MB"] %>% filter_zeros()
dds_C3_MB <-
  dds_C3[, colData(dds_C3)$compartment == "MB"] %>% filter_zeros()
dds_C2_AXON <-
  dds_C2[, colData(dds_C2)$compartment == "AXON"] %>% filter_zeros()
dds_C3_AXON <-
  dds_C3[, colData(dds_C3)$compartment == "AXON"] %>% filter_zeros()

# # create cohort fraction subsets
# dds_C1_IP <-
#   dds[, colData(dds)$cohort == "C1" &
#         colData(dds)$fraction == "IP"] %>% filter_zeros()
dds_C2_IP <-
  dds[, colData(dds)$cohort == "C2" &
        colData(dds)$fraction == "IP"] %>% filter_zeros()
dds_C3_IP <-
  dds[, colData(dds)$cohort == "C3" &
        colData(dds)$fraction == "IP"] %>% filter_zeros()
# dds_C1_TOTAL <- dds[, colData(dds)$cohort == "C1" & colData(dds)$fraction == "TOTAL"] %>% filter_zeros()
dds_C2_TOTAL <- dds[, colData(dds)$cohort == "C2" & colData(dds)$fraction == "TOTAL"] %>% filter_zeros()
dds_C3_TOTAL <- dds[, colData(dds)$cohort == "C3" & colData(dds)$fraction == "TOTAL"] %>% filter_zeros()

# # create cohort compartment fraction subsets
dds_C1_MB_IP <-
  dds_C1_MB[, colData(dds_C1_MB)$fraction == "IP"] %>% filter_zeros()
dds_C2_MB_IP <-
  dds_C2_MB[, colData(dds_C2_MB)$fraction == "IP"] %>% filter_zeros()
dds_C3_MB_IP <-
  dds_C3_MB[, colData(dds_C3_MB)$fraction == "IP"] %>% filter_zeros()
dds_C1_MB_TOTAL <-
  dds_C1_MB[, colData(dds_C1_MB)$fraction == "TOTAL"] %>% filter_zeros()
dds_C2_MB_TOTAL <-
  dds_C2_MB[, colData(dds_C2_MB)$fraction == "TOTAL"] %>% filter_zeros()
dds_C3_MB_TOTAL <-
  dds_C3_MB[, colData(dds_C3_MB)$fraction == "TOTAL"] %>% filter_zeros()
dds_C2_AXON_IP <-
  dds_C2_AXON[, colData(dds_C2_AXON)$fraction == "IP"] %>% filter_zeros()
dds_C3_AXON_IP <-
  dds_C3_AXON[, colData(dds_C3_AXON)$fraction == "IP"] %>% filter_zeros()
dds_C2_AXON_TOTAL <-
  dds_C2_AXON[, colData(dds_C2_AXON)$fraction == "TOTAL"] %>% filter_zeros()
dds_C3_AXON_TOTAL <-
  dds_C3_AXON[, colData(dds_C3_AXON)$fraction == "TOTAL"] %>% filter_zeros()

# create young and old dds objects for genotype comparisons
dds_C1_MB_IP_YOUNG <- dds_C1_MB_IP[,colData(dds_C1_MB_IP)$age == "YOUNG"] %>% filter_zeros()
dds_C1_MB_IP_OLD <- dds_C1_MB_IP[,colData(dds_C1_MB_IP)$age == "OLD"] %>% filter_zeros()

dds_C3_MB_IP_YOUNG <- dds_C3_MB_IP[,colData(dds_C3_MB_IP)$age == "YOUNG"] %>% filter_zeros()
dds_C3_MB_IP_OLD <- dds_C3_MB_IP[,colData(dds_C3_MB_IP)$age == "OLD"] %>% filter_zeros()

# create MB TOTAL C1 C2 object to identify MB IP uniquely DEGs in AGE and GENOTYPE comparisons
dds_C1_C2_MB_TOTAL <- dds[, colData(dds)$compartment == "MB" &
                            colData(dds)$cohort %in% c("C1", "C2") &
                            colData(dds)$fraction == "TOTAL"] %>% filter_zeros()

# view counts function ----------------------------------------------------

view_counts <- function(dds, gene) {
  counts(dds) %>%
    as_tibble(rownames = "ensembl_gene_id") %>%
    pivot_longer(-ensembl_gene_id,
                 names_to = "sample_name",
                 values_to = "count") %>%
    left_join(anno) %>%
    filter(external_gene_name == gene)
}

view_vsd <- function(dds, gene) {
  assay(vst(dds)) %>%
    as_tibble(rownames = "ensembl_gene_id") %>%
    pivot_longer(-ensembl_gene_id,
                 names_to = "sample_name",
                 values_to = "count") %>%
    left_join(anno) %>%
    filter(external_gene_name == gene)
}

# long counts function ----------------------------------------------------

long_counts <- function(dds) {
  counts(dds) %>%
    as_tibble(rownames = "ensembl_gene_id") %>%
    pivot_longer(-ensembl_gene_id,
                 names_to = "sample_name",
                 values_to = "count")
}

# human anno --------------------------------------------------------------

# anno_human <- get_anno(rownames(dds), get_human = T)
anno_human <- readRDS("R/objects/anno_human.rds")

anno_human[anno_human$external_gene_name == "Fgf13",]$hsapiens_homolog_ensembl_gene <- "ENSG00000129682"
anno_human[anno_human$external_gene_name == "Fgf13",]$hsapiens_homolog_associated_gene_name <- "FGF13"

anno_human[anno_human$external_gene_name == "A2ml1",]$hsapiens_homolog_ensembl_gene <- "ENSG00000166535"
anno_human[anno_human$external_gene_name == "A2ml1",]$hsapiens_homolog_associated_gene_name <- "A2ML1"

# set up marts for human genes and getting SNP coordinates
ensembl_human <- useMart("ensembl",
                         dataset="hsapiens_gene_ensembl")
ensembl_variation <- useMart("ENSEMBL_MART_SNP",
                             dataset = "hsapiens_snp")

# convert gene names functions ---------------------------------------------

get_ensembl_gene_id <- function(external_gene_name) {
  anno[anno$external_gene_name == external_gene_name, ]$ensembl_gene_id
}
get_external_gene_name <- function(ensembl_gene_id) {
  anno[anno$ensembl_gene_id == ensembl_gene_id, ]$external_gene_name
}
get_human_gene_name <- function(ensembl_gene_id) {
  anno_human[anno_human$ensembl_gene_id == ensembl_gene_id, ]$hsapiens_homolog_associated_gene_name
}




# build STRING Interaction database ---------------------------------------

# setting to human database, for human interactions
string_db <- STRINGdb$new(version = "11",
                          # species = 10090, # mus musculus
                          species = 9606,
                          score_threshold = 400, # medium confidence
                          input_directory = "")

# setting to human database, for human interactions
string_db_mouse <- STRINGdb$new(version = "11",
                          species = 10090, # mus musculus
                          score_threshold = 400, # medium confidence
                          input_directory = "")

string_genes_df <- data.frame(gene = anno_human %>%
                                filter(hsapiens_homolog_ensembl_gene != "") %>%
                                pull(unique(hsapiens_homolog_ensembl_gene))) %>%
  left_join(anno_human, by = c("gene" = "hsapiens_homolog_ensembl_gene")) %>%
  distinct()

#String can't map non-protein-coding genes (and many others!)
string <- string_db$map(string_genes_df,
                        "hsapiens_homolog_associated_gene_name",
                        removeUnmappedRows = TRUE)

# STRING neighbour plots --------------------------------------

get_string_neighbours <- function(external_gene_name,
                                  top_n = 10) {
  # gene <- string_genes_df[string_genes_df$external_gene_name == external_gene_name, ]$hsapiens_homolog_associated_gene_name
  gene <- string_db$mp(external_gene_name)

  top_n <- tibble(to = string_db$get_neighbors(c(gene))) %>%
    mutate(.before = to,
           from = gene) %>%
    mutate(combined_score = sapply(to, function(x) {
      string_db$get_interactions(c(from, x))[1, 3]})) %>%
    # mutate(from = external_gene_name,
    #        to = string$external_gene_name[match(to, string$STRING_id)]) %>%
    drop_na() %>%
    slice_max(order_by = combined_score,
              n = top_n) %>%
    pull(to)

  neighbours <- string_db$get_interactions(c(gene, top_n)) %>%
    mutate(from = string$hsapiens_homolog_associated_gene_name[match(from, string$STRING_id)],
           to = string$hsapiens_homolog_associated_gene_name[match(to, string$STRING_id)]) %>%
    drop_na() %>%
    distinct()
}

# plot neighbours
plot_neighbours <- function(gene, neighbours){
  igr <- graph_from_data_frame(neighbours,
                               directed = F)

  V(igr)$size <- 5
  # V(igr)$frame.color <- "white"
  # V(igr)$color <- "orange"
  # V(igr)$label <- ""
  # E(igr)$arrow.mode <- 0
  E(igr)$width <- 2

  inc.edges <- incident(igr,  gene, mode="all")
  # Set colors to plot the selected edges.
  ecol <- rep("gray80", ecount(igr))
  ecol[inc.edges] <- "orange"
  E(igr)[inc.edges]$width <- 3
  vcol <- rep("lightblue", vcount(igr))
  vcol[names(V(igr))==gene] <- "gold"
  vfont <- rep(1, vcount(igr))
  vfont[names(V(igr)) == gene] <- 2
  vcol[names(V(igr)) %in% sapply(MB_FRACTION_ENRICHED_GENES, get_human_gene_name)] <- "green"

  vframecol <- rep("grey10", vcount(igr))
  # vframecol[names(V(igr)) %in% sapply(MB_FRACTION_ENRICHED_GENES, get_external_gene_name)] <- "#2CA02CFF"

  V(igr)$shape <- ifelse(names(V(igr)) == gene, "csquare", "circle")

  plot(igr,
       vertex.color=vcol,
       vertex.frame.color = vframecol,
       vertex.label.font = vfont,
       vertex.label.color = "black",
       vertex.label.family = "Helvetica",
       vertex.label.cex = 0.8,
       vertex.label.dist = 3,
       vertex.label.degree = -pi/2,
       # edge.curved = 0.1,
       edge.color=ecol,
       layout = layout.circle,
       # layout = layout.fruchterman.reingold
  )
  # title(gene)

}


# Allen Brain Atlas function ----------------------------------------------

get_aba <- function(external_gene_name,
                    region = "snpc",
                    experiment = 1,
                    downsample = 3,
                    plane = "coronal",
                    orientation = "antisense"){

  # set region
  if(is.numeric(region)){
    region_id = region
  }
  if(region == "snpc"){
    region_id = 374
  }
  if(region == "vta"){
    region_id = 749
  }
  if(region == "snpr"){
    region_id = 381
  }
  if(region == "snpl"){
    region_id = 615
  }
  if(region == "midbrain"){
    region_id = 313
  }
  if(region == "midbrain_motor"){
    region_id = 323
  }
  if(region == "midbrain_sensory"){
    region_id = 339
  }
  if(region == "midbrain_behaviour"){
    region_id = 348
  }

  # get the dataset id
  dataset_id = getGeneDatasets(
    gene = external_gene_name,
    planeOfSection = plane,
    probeOrientation = orientation
  )[experiment]

  # get the image id
  image_id = structureToImage(datasetID = dataset_id,
                              regionIDs = region_id)


  # atlas_id = imageToAtlas(image_id$section.image.id,
  #                         image_id$x,
  #                         image_id$y,
  #                         planeOfSection = 'coronal')

  # download the image
  downloadImage(
    imageID = image_id$section.image.id,
    view = "projection",
    outputFile = paste0('R/objects/ABA/', external_gene_name, '.jpg'),
    downsample = downsample)

  return(image_id)
}

# centering function
center_aba <- function(external_gene_name,
                       image_id,
                       x = c(0.1, 0.1),
                       y = c(0.1, 0.1)){
  centerImage(
    image = paste0('R/objects/ABA/', external_gene_name, '.jpg'),
    x = image_id$x,
    y = image_id$y,
    xProportions = x,
    yProportions = y,
    outputFile = paste0('R/objects/ABA/', external_gene_name, '_centered.jpg'),
    downsample = 3
  )
}



# mSigDB ------------

# msigdbr_collections() %>% View
cp_kegg_gene_sets <-
  msigdbr(species = "Mus musculus",
          category = "C2",
          subcategory = "CP:KEGG") %>%
  select("external_gene_name" = gene_symbol,
         gs_name)

# cp_reactome_gene_sets <-
#   msigdbr(species = "Mus musculus",
#           category = "C2",
#           subcategory = "CP:REACTOME") %>%
#   select("external_gene_name" = gene_symbol,
#          gs_name)
#
# cp_gocc_gene_sets <-
#   msigdbr(species = "Mus musculus",
#           category = "C5",
#           subcategory = "GO:CC") %>%
#   select("external_gene_name" = gene_symbol,
#          gs_name)
#
# cp_gobp_gene_sets <-
#   msigdbr(species = "Mus musculus",
#           category = "C5",
#           subcategory = "GO:BP") %>%
#   select("external_gene_name" = gene_symbol,
#          gs_name)

# keggs_of_interest <- c(
#   "KEGG_ALZHEIMERS_DISEASE",
#   "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
#   "KEGG_AXON_GUIDANCE",
#   "KEGG_CALCIUM_SIGNALING_PATHWAY",
#   "KEGG_HUNTINGTONS_DISEASE",
#   "KEGG_LYSOSOME",
#   "KEGG_PARKINSONS_DISEASE",
#   "KEGG_REGULATION_OF_AUTOPHAGY",
#   "KEGG_PROTEASOME",
#   "KEGG_SNARE_INTERACTIONS_IN_VESICULAR_TRANSPORT",
#   "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS"
# )

cp_kegg_pd <- cp_kegg_gene_sets %>%
  filter(gs_name == "KEGG_PARKINSONS_DISEASE")
#
# cp_kegg_interest <- cp_kegg_gene_sets %>%
#   filter(gs_name %in% keggs_of_interest) %>%
#   arrange(external_gene_name)
#
# jensen_knowledge <-
#   read_delim(
#     "R/objects/human_disease_knowledge_filtered.tsv",
#     delim = "\t",
#     col_names = FALSE
#   ) %>%
#   select(
#     "hsapiens_homolog_associated_gene_name" = X2,
#     "disease" = X4,
#     "disease_confidence" = X7
#   )

#
#
# MB FRACTION
# dGIDB ------
# dgidb <- read_delim("R/objects/dgidb_interactions.tsv",
#                     delim = "\t") %>%
#   inner_join(anno_human,
#              by = c("gene_name" = "hsapiens_homolog_associated_gene_name")) %>%
#   mutate(link = paste0("https://www.dgidb.org/genes/", gene_name))

#
#
# transcript translator ----

# anno_tx <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "ensembl_transcript_id",
#                  "external_transcript_name"),
#   filters = "ensembl_gene_id",
#   values = rownames(dds),
#   mart = ensembl,
#   uniqueRows = TRUE
# )

anno_tx <- readRDS("R/objects/anno_tx.rds")

# PUBMED ----

# set_entrez_key("73e832468a0883d3294f0fb858d6b70f3f08")

# publications_axon <- sapply(AXON_VS_MB_META_SIGNIF_PRIORITY$external_gene_name, function(x){try(get_pubmed_axon(x)$count)})
# publications_axon <- enframe(publications_axon) %>%
#   select(external_gene_name = name,
#          axon_pub = value)

publications_all_genes <- readRDS("R/objects/publications_all_genes.rds") %>%
  enframe %>%
  select(ensembl_gene_id = name,
         pub_pd = value)

publications_all_genes_axon <- readRDS("R/objects/publications_all_genes_axon.rds") %>%
  enframe %>%
  select(ensembl_gene_id = name,
         pub_axon = value)

# ----
# ----
# ----
# END OF SETUP ----
# ----
# ----
# ----
# Load ONT Data & QC ----

# Load counts
ont_tx_counts <-
  vroom(
    "/zfs/analysis/trap/active/testing_env/dtu_ont/merge_counts_output/merged_counts.tsv",
    delim = "\t"
  ) %>%
  separate(
    Reference,
    into = c(
      "ensembl_transcript_id",
      "ensembl_gene_id",
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    sep = "([|])"
  )

# The reference column has 8 components
# all(unlist(lapply(strsplit(ont_tx_counts$Reference, "|", fixed = T), length)) == 8)
# strsplit(ont_tx_counts$Reference[1], "|", fixed = T)
ont_metadata <- tibble(flowcell_barcode = colnames(ont_tx_counts)[-c(1, 2)],
                       barcode = str_extract(flowcell_barcode, "(?<=_)barcode[:digit:]{2}"),
                       flowcell = str_extract(flowcell_barcode, "[:alnum:]+(?=_)"))

ont_metadata0 <- vroom("input/metadata/trap_2020_ont_metadata_mb.csv",
                       delim = ",") %>%
  select("sample_name0" = sample_id,
         barcode)

ont_metadata <- ont_metadata %>%
  left_join(ont_metadata0,
            by  = "barcode") %>%
  left_join(
    colData(dds),
    copy = T
  ) %>%
  select(sample_name,
         barcode,
         flowcell,
         flowcell_barcode,
         cohort,
         region,
         age,
         genotype,
         sex,
         collection,
         fraction,
         compartment,
         sample_code)

# collapse transcripts
gene_id <- as.factor(ont_tx_counts$ensembl_gene_id)
length(gene_id) == nrow(ont_tx_counts)
sp <- split(seq(along = gene_id), gene_id)
countData <- as.matrix(ont_tx_counts[,-c(1,2)])
rownames(countData) <- ont_tx_counts$ensembl_transcript_id
countData <- sapply(sp, function(i) colSums(countData[i,, drop = FALSE]))
dim(countData)
countData <- t(countData)

dds_ONT <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = ont_metadata,
                                  design = ~1)
# remove flongle data
dds_ONT <- dds_ONT[,colData(dds_ONT)$flowcell != "AEL855"]

# Total counts per flowcell
p_ONT_counts_per_flowcell <- enframe(colSums(counts(dds_ONT))) %>%
  rename(flowcell_barcode = name,
         total_counts = value) %>%
  left_join(colData(dds_ONT), copy = T) %>%
  mutate(flowcell = factor(flowcell,
                           levels = c("FAL15042",
                                      "FAL02252",
                                      "FAL16168",
                                      "FAL02011"))) %>%
  group_by(flowcell) %>%
  mutate(sample = as.factor(row_number())) %>%
  ggplot(aes(x = sample,
             y = total_counts,
             fill = flowcell)) +
  geom_col(colour = "black") +
  facet_wrap(vars(flowcell),
             scales = "free_x",
             nrow = 1) +
  labs(x = "Sample",
       y = "Total Counts") +
  scale_fill_d3() +
  theme(legend.position = "none")

# calculate vsd
vsd_ONT <- vst(dds_ONT)

# plot flowcell pca
plotPCA(vsd_ONT,
        intgroup = "flowcell") +
  scale_color_d3()

# collapse flowcell replicates
dds_ONT <- collapseReplicates(dds_ONT,
                              groupby = dds_ONT$sample_name,
                              run = dds_ONT$flowcell)

# recalculate vsd
vsd_ONT <- vst(dds_ONT)

# plot age and genotype pca
plotPCA(vsd_ONT,
        intgroup = c("age", "genotype")) +
  scale_color_d3()

# remove version number from dds_ONT rownames
rownames(dds_ONT) <- substr(rownames(dds_ONT), 1, 18)

# Select (randomly for C1 and C2) the samples to compare to long read data
set.seed(821196)
ont_comparison_colnames <- colData(dds) %>%
  as_tibble() %>%
  filter(cohort != "C3" &
           compartment == "MB" &
           fraction == "IP") %>%
  group_by(cohort,
           age,
           genotype) %>%
  slice_sample(n = 3) %>%
  pull(sample_name) %>%
  c(colnames(dds_ONT))

# make the short dds subset
dds_ONT_short_equiv <- dds[,colnames(dds) %in% ont_comparison_colnames]
# make the "long equiv" dds subset
dds_ONT_long_equiv <- dds_ONT

# get the list of genes above threshold in both datasets
ont_comparison_rownames <- intersect(rownames(dds_ONT_short_equiv)[rowSums(counts(dds_ONT_short_equiv) > 3) == ncol(dds_ONT_short_equiv)],
                                     rownames(dds_ONT)[rowSums(counts(dds_ONT) > 3) == ncol(dds_ONT)])
# filter genes
dds_ONT_short_equiv <- dds_ONT_short_equiv[rownames(dds_ONT_short_equiv) %in% ont_comparison_rownames,]
dds_ONT_long_equiv <- dds_ONT_long_equiv[rownames(dds_ONT_long_equiv) %in% ont_comparison_rownames,]
# check the same genes are in each dds
all(rownames(dds_ONT_short_equiv) == rownames(dds_ONT_long_equiv))

# get number of genes used for comparison
n_dds_ONT_long_equiv <- nrow(dds_ONT_long_equiv)

# make cohort short, and long vsd objects
vsd_ONT_short_equiv_C1 <- vst(dds_ONT_short_equiv[,colData(dds_ONT_short_equiv)$cohort == "C1"])
vsd_ONT_short_equiv_C2 <- vst(dds_ONT_short_equiv[,colData(dds_ONT_short_equiv)$cohort == "C2"])
vsd_ONT_short_equiv_C3 <- vst(dds_ONT_short_equiv[,colData(dds_ONT_short_equiv)$cohort == "C3"])
vsd_ONT_long_equiv <- vst(dds_ONT_long_equiv)

# make matrix
short_long_matrix <- cbind(
  assay(vsd_ONT_short_equiv_C1),
  assay(vsd_ONT_short_equiv_C2),
  assay(vsd_ONT_short_equiv_C3),
  assay(vsd_ONT_long_equiv)
)

# uniquely name the columns
colnames(short_long_matrix) <- c(paste0(colnames(short_long_matrix)[1:30], "_short"),
                                 paste0(colnames(short_long_matrix)[31:42], "_long"))

# make an annotation df
annot <- data.frame(Cohort = c(rep("Cohort 1", 12),
                               rep("Cohort 2", 6),
                               rep("Cohort 3", 12),
                               rep("Cohort 3", 12)),
                    "Read Length" = c(rep("Short", 30),
                                      rep("Long", 12)))
rownames(annot) <- colnames(short_long_matrix)
mat_colors <- list(Cohort = pal_d3()(10)[c(1, 2, 4)],
                   Read.Length = pal_d3()(10)[c(8, 9)])
names(mat_colors$Cohort) <- unique(annot$Cohort)
names(mat_colors$Read.Length) <- unique(annot$Read.Length)

# calculate cor
cor_ONT_short_long <- cor(short_long_matrix)

# plot the correlation of each sample
p_ONT_cor_short_long <- as.grob(
  ~ pheatmap(
    mat = cor_ONT_short_long,
    treeheight_row = 0,
    # treeheight_col = 0,
    show_rownames = F,
    show_colnames = F,
    annotation_col = annot,
    annotation_colors = mat_colors,
    color = mako(20),
    border_color = NA,
    cellwidth = 6,
    cellheight = 6
  )
)

# Guppy QC

guppy_summary <- vroom("/zfs/analysis/trap/active/snakemake/guppy_basecalling/output/sequencing_summary_sample.txt",
                       delim = "\t") %>%
  mutate(flowcell = str_extract(filename, "[:alnum:]+(?=_)")) %>%
  select(flowcell,
         barcode = barcode_arrangement,
         sequence_length_template,
         passes_filtering,
         alignment_coverage)

# colnames(guppy_summary)
# guppy_summary %>%
#   slice_head(n = 10) %>%
#   View


plot_data <- guppy_summary %>%
  filter(!outlier_iqr(sequence_length_template) &
           flowcell != "AEL855" &
           passes_filtering) %>%
  inner_join(ont_metadata) %>%
  mutate(flowcell = factor(flowcell, levels = c("FAL15042",
                                                "FAL02252",
                                                "FAL16168",
                                                "FAL02011")))

p_ONT_read_length_flowcell <- ggplot(plot_data,
       aes(x = sequence_length_template,
             y = ..scaled..,
             colour = flowcell)) +
  geom_density() +
  scale_color_d3() +
  theme(legend.position = "top") +
  labs(x = "Read Length",
       y = "Scaled Proportion",
       color = "Flowcell") +
  guides(color = guide_legend(nrow = 2,
                             byrow=TRUE))

p_ONT_read_length_condition <- plot_data %>%
  unite(age, genotype, col = "age_genotype", sep = " ") %>%
  ggplot(aes(x = sequence_length_template,
             # y = ..scaled..,
             colour = age_genotype)) +
  stat_ecdf() +
  scale_color_d3() +
  theme(legend.position = "top") +
  labs(x = "Read Length",
       y = "Cumulative Proportion",
       color = "Condition") +
  guides(colour = guide_legend(nrow = 2,
                             byrow=TRUE))


plot_data %>%
  unite(age, genotype, col = "age_genotype", sep = " ") %>%
  ggplot(aes(x = alignment_coverage,
             # y = ..scaled..,
             colour = age_genotype)) +
  geom_density() +
  scale_color_d3() +
  theme(legend.position = "top") +
  labs(x = "Read Length",
       y = "Cumulative Proportion",
       color = "Condition") +
  guides(colour = guide_legend(nrow = 2,
                               byrow=TRUE)) +
  facet_wrap(vars(flowcell))




# ----
# ----
# ----
# AXON FRACTION: ABA data ----
# ABA
# # Get structure ids
# res <- GET("http://api.brain-map.org/api/v2/data/Structure/query.json?&num_rows=13709&start_row=0")
# stop_for_status(res)
# res <- fromJSON(toJSON(content(res)))
# aba_structures <- res$msg %>% flatten()

# Get mouse dataset ids
# res <- GET("http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse'],treatments[name$eq'ISH'],genes,plane_of_section,rma::options,[tabular$eq'plane_of_sections.name+as+plane','genes.acronym+as+gene','data_sets.id+as+section_data_set_id'],[order$eq'plane_of_sections.name,genes.acronym,data_sets.id']&num_rows=25000&start_row=0")
# stop_for_status(res)
# res <- fromJSON(toJSON(content(res)))
# res0 <- GET("http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse'],treatments[name$eq'ISH'],genes,plane_of_section,rma::options,[tabular$eq'plane_of_sections.name+as+plane','genes.acronym+as+gene','data_sets.id+as+section_data_set_id'],[order$eq'plane_of_sections.name,genes.acronym,data_sets.id']&num_rows=1079&start_row=25001")
# stop_for_status(res0)
# res0 <- fromJSON(toJSON(content(res0)))
#
# aba_datasets <- bind_rows(res$msg, res0$msg) %>%
#   unnest(cols = c(plane, gene, section_data_set_id))

aba_datasets <- readRDS("R/objects/aba_datasets.rds")

# # get expression summary
# res <- GET(
#   "http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,rma::include,structure_unionizes(structure[id$eq477]),genes(gene_aliases(gene))&num_rows=20000&start_row=0"
# )
# stop_for_status(res)
# res <- fromJSON(toJSON(content(res)))
# res$total_rows
# res0 <- GET(
#   "http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,rma::include,structure_unionizes(structure[id$eq477]),genes(gene_aliases(gene))&num_rows=10435&start_row=20000"
# )
# stop_for_status(res0)
# res0 <- fromJSON(toJSON(content(res0)))
# res0

# aba_expression <- bind_rows(
#   {
#     res$msg$structure_unionizes %>% bind_rows %>% flatten
#   }, {
#     res0$msg$structure_unionizes %>% bind_rows %>% flatten
#   }
# )

# saveRDS(aba_expression, "R/objects/aba_expression.rds")

aba_expression <- readRDS("R/objects/aba_expression.rds")

aba_expression_axon <- aba_expression %>%
  select(section_data_set_id,
         expression_energy) %>%
  unnest(c(section_data_set_id,
           expression_energy)) %>%
  mutate(section_data_set_id = as.character(section_data_set_id)) %>%
  inner_join(aba_datasets) %>% # filter for AXON ENRICHED AGNOSTIC GENES
  select(external_gene_name = gene,
         everything()) %>%
  group_by(external_gene_name) %>%
  summarise(expression_energy = mean(expression_energy))
# ----
# ----
# ----
# MB fraction subset common genes -----------------------------------------


# filter genes for common detected
# Not applying this for fraction comparison, to avoid conflicts,
# where gene might be enriched in all cohorts but detected in few cohort 3 samples, so is excluded from analysis.
# dds_C1_MB_FILTER <- filter_genes(dds_C1_MB,
#                                  grouping = c("fraction", "gene_id"))
# dds_C2_MB_FILTER <- filter_genes(dds_C2_MB,
#                                  grouping = c("fraction", "gene_id"))
# dds_C3_MB_FILTER <- filter_genes(dds_C3_MB,
#                                  grouping = c("fraction", "gene_id"))

# genes_MB <-
#   intersect(dds_C1_MB_FILTER, intersect(dds_C2_MB_FILTER,
#                                         dds_C3_MB_FILTER))

# find common genes in all MB IP and TOTAL samples
genes_MB <-
  intersect(rownames(dds_C1_MB), intersect(rownames(dds_C2_MB),
                                           rownames(dds_C3_MB)))

dds_C1_MB <-
  dds_C1_MB[genes_MB, ]
dds_C2_MB <-
  dds_C2_MB[genes_MB, ]
dds_C3_MB <-
  dds_C3_MB[genes_MB, ]

rm(genes_MB)

n_MB_FRACTION <- nrow(dds_C1_MB)

# MB fraction DESeq2 -----------------------------------------------------
# this segment generates a list of genes enriched in midbrain by TRAP
# It uses DESeq2 with a filtered number of genes that have to be expressed in at least
# 3/5ths of all samples in either IP or TOTAL
# Each cohort is analysed separately
# One sided p values are used to calculate a meta p value with stouffer's method
# fdr adjustment is applied and sumz pvalues less than 0.01 are considered significant

dds_C1_MB@design <- ~ fraction
dds_C2_MB@design <- ~ fraction
dds_C3_MB@design <- ~ fraction

dds_C1_MB <- DESeq(
  dds_C1_MB,
  # sfType = "poscounts",
  # fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

dds_C2_MB <- DESeq(
  dds_C2_MB,
  # sfType = "poscounts",
  # fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

dds_C3_MB <- DESeq(
  dds_C3_MB,
  sfType = "poscounts",
  fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

# plotDispEsts(dds_C1_MB)
# plotDispEsts(dds_C2_MB)
# plotDispEsts(dds_C3_MB)

res_C1_MB <- DESeq2::results(
  dds_C1_MB,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C1_MB <- lfcShrink(
  dds_C1_MB,
  res = res_C1_MB,
  contrast = c("fraction", "IP", "Input"),
  type = "ashr",
  parallel = T
)

res_C2_MB <- DESeq2::results(
  dds_C2_MB,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_MB <- lfcShrink(
  dds_C2_MB,
  res = res_C2_MB,
  contrast = c("fraction", "IP", "Input"),
  type = "ashr",
  parallel = T
)

res_C3_MB <- DESeq2::results(
  dds_C3_MB,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_MB <- lfcShrink(
  dds_C3_MB,
  res = res_C3_MB,
  contrast = c("fraction", "IP", "Input"),
  type = "ashr",
  parallel = T
)

# summary(res_C1_MB)
# summary(res_C2_MB)
# summary(res_C3_MB)

MB_FRACTION_META <-
  as_tibble(res_C1_MB, rownames = "ensembl_gene_id") %>%
  left_join(
    as_tibble(res_C2_MB, rownames = "ensembl_gene_id"),
    by = "ensembl_gene_id",
    suffix = c("_C1", "_C2")
  ) %>%
  left_join({
    as_tibble(res_C3_MB, rownames = "ensembl_gene_id") %>%
      rename_with(~ paste0(.x, "_C3"), starts_with(c("base", "log", "lfc", "pva", "padj")))
  },
  by = "ensembl_gene_id") %>%
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
  mutate(
    pvalue_C1 = ifelse(
      # correct pvalues where there is a conflict in log2foldchange between groups
      sign(log2FoldChange_C1) == sign(# if the sign of a cohort l2fc is the same as the median, keep the pvalue, otherwise invert it
        median(
          c(log2FoldChange_C1,
            log2FoldChange_C2,
            log2FoldChange_C3)
        )),
      pvalue_C1,
      1 - pvalue_C1
    ),
    pvalue_C2 = ifelse(
      sign(log2FoldChange_C2) == sign(median(
        c(log2FoldChange_C1,
          log2FoldChange_C2,
          log2FoldChange_C3)
      )),
      pvalue_C2,
      1 - pvalue_C2
    ),
    pvalue_C3 = ifelse(
      sign(log2FoldChange_C3) == sign(median(
        c(log2FoldChange_C1,
          log2FoldChange_C2,
          log2FoldChange_C3)
      )),
      pvalue_C3,
      1 - pvalue_C3
    )
  ) %>%
  mutate(across(starts_with("pvalue"),
                ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
  mutate(
    # sumlog = sumlog(c_across(pvalue_C1:pvalue_C3))$p,
    # calculate sumlog
    sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
    # calculate stouffers
    log2FoldChange = median(c(
      log2FoldChange_C1, # calculate median l2fc
      log2FoldChange_C2,
      log2FoldChange_C3
    ))
  ) %>%
  ungroup() %>%
  mutate(
    # sumlog_adj = p.adjust(sumlog, method = "fdr"),
    # correction for multiple comparisons
    sumz_adj = p.adjust(sumz, method = "fdr"),
    # conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C2) |
    #   sign(log2FoldChange_C2) != sign(log2FoldChange_C3), # previous approach
    conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C2) # Only report conflict if C1 and C2 disagree: C3 disagreements can be put down to high dispersion and dropouts (see Casr as an example)
    ) %>% # state whether there is a conflict in l2fc direction
  mutate(score = -log10(sumz_adj) * log2FoldChange) %>%
  arrange(desc(score), desc(log2FoldChange))

# make a simple MB FRACTION META table
MB_FRACTION_META_SIMPLE <- MB_FRACTION_META %>%
  select(external_gene_name,
         description,
         "baseMean" = baseMean_C1,
         log2FoldChange,
         sumz_adj,
         conflict,
         score) %>%
  mutate(rank_enrichment = row_number())

# sumlog finds a greater proportion of conflicting genes significant
# sumz is more conservative
# MB_FRACTION_META %>%
#   group_by(conflict) %>%
#   summarise(n = n(),
#             n_sig = sum(sumz_adj < 0.01),
#             n_sig_sumlog = sum(sumlog_adj < 0.01))

# get list of enriched genes
MB_FRACTION_ENRICHED_GENES <- MB_FRACTION_META %>%
  filter(sumz_adj < 0.01) %>%
  rowwise() %>%
  filter(log2FoldChange > 0) %>%
  pull(ensembl_gene_id)

# get list of depleted genes
MB_FRACTION_DEPLETED_GENES <- MB_FRACTION_META %>%
  filter(sumz_adj < 0.01) %>%
  rowwise() %>%
  filter(log2FoldChange < 0) %>%
  pull(ensembl_gene_id)

n_MB_FRACTION_ENRICHED_GENES <- length(MB_FRACTION_ENRICHED_GENES)
n_MB_FRACTION_DEPLETED_GENES <- length(MB_FRACTION_DEPLETED_GENES)


# TABLE 1: MB fraction table summary -----------------------------------------------
t_mb_ip_enrichment_summary <- tibble(Metric = c("Number of Genes", "Percentage of Library"),
                                     Enriched = c(as.integer(length(MB_FRACTION_ENRICHED_GENES)),
                                                  paste(signif(sum((rowSums(counts(dds_C1_MB_IP))/sum(colSums(counts(dds_C1_MB_IP))))[MB_FRACTION_ENRICHED_GENES]*100, na.rm = T), 3), "%")),
                                     Depleted = c(as.integer(length(MB_FRACTION_DEPLETED_GENES)),
                                                  paste(signif(sum((rowSums(counts(dds_C1_MB_IP))/sum(colSums(counts(dds_C1_MB_IP))))[MB_FRACTION_DEPLETED_GENES]*100, na.rm = T), 3), "%")),
                                     Unchanged = c(as.integer(nrow(MB_FRACTION_META)-length(MB_FRACTION_ENRICHED_GENES)-length(MB_FRACTION_DEPLETED_GENES)),
                                                   paste(signif(sum((rowSums(counts(dds_C1_MB_IP))/sum(colSums(counts(dds_C1_MB_IP))))[!rownames(dds_C1_MB_IP) %in% c(MB_FRACTION_ENRICHED_GENES, MB_FRACTION_DEPLETED_GENES)] * 100, na.rm = T), 3), "%")))


# FIGURE 1: MB enrichment and cell type validation
# FIGURE A: MB fraction plot counts -------------------------------------------------

# MA-plot style
plot_data <- MB_FRACTION_META %>%
  # filter(abs(log2FoldChange) < 10) %>%
  mutate(enrichment = ifelse(sumz_adj > 0.01, "Unchanged",
                             ifelse(log2FoldChange > 0, "Enriched", "Depleted")))

# highlight markers
highlight_markers <- plot_data %>%
  left_join(anno) %>%
  filter(external_gene_name %in% c(
    c("Slc6a3", "Th", "Ddc", "Slc18a2"),
    c("Gfap", "Gad2", "S100b", "Gad1", "Aldh1l1")
  ))

p_mb_ip_enrichment_count_plot <- plot_data %>%
  ggplot(aes(x = baseMean_C1,
             y = log2FoldChange,
             colour = enrichment)) +
  geom_point(size = 0.5,
             alpha = 0.5) +
  geom_point(data = highlight_markers,
             aes(fill = enrichment),
             size = 2,
             shape = 21,
             colour = "black") +
  geom_label_repel(data = highlight_markers,
                   aes(label = external_gene_name),
                   colour = "black",
                   arrow = arrow(length = unit(0.02, "npc")),
                   box.padding = 1,
                   # nudge_x = -1
  ) +
  scale_x_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)) +
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
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ,
         fill = FALSE)

# FIGURE B: MB fraction panglao --------------------------------------------------

# load panglao DB
panglao <- read_delim("R/objects/PanglaoDB_markers_27_Mar_2020.tsv",
                      delim = "\t") %>%
  # filter(ifelse(str_detect(`official gene symbol`, "TH|SLC6A3|SLC18A2|DDC") & !str_detect(`cell type`, "drenergic|otonergic|nterneuron"),
  #               TRUE, ifelse(!str_detect(`official gene symbol`, "TH|SLC6A3|SLC18A2|DDC"), TRUE, FALSE ))) %>%
  group_by(`cell type`) %>%
  mutate(group_size = length(unique(`official gene symbol`)))

# internal
# score for each cell type
# adrenergic/nor-/sero- also enriched,
# but this is because of Slc6a3/Slc18a2/Th
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
  distinct

# table of dopaminergic marker enrichment
t_mb_ip_enrichment_dopamine_markers <- MB_FRACTION_META %>%
  left_join(anno_human) %>%
  left_join(panglao,
            by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
  mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
         score = -log10(sumz_adj) * log2FoldChange_C1) %>%
  group_by(`cell type`) %>%
  filter(`cell type` == "Dopaminergic neurons" &
           log2FoldChange_C1 > 0 &
           sumz_adj < 0.01) %>%
  rowwise() %>%
  mutate(log2FoldChange = mean(log2FoldChange_C1,
                               log2FoldChange_C2,
                               log2FoldChange_C3)) %>%
  ungroup() %>%
  arrange(sumz_adj, desc(log2FoldChange)) %>%
  select("Gene" = external_gene_name,
         "Log2 Fold Change" = log2FoldChange,
         "Adjusted P Value" = sumz_adj)

n_mb_ip_enrichment_da_markers <- nrow(t_mb_ip_enrichment_dopamine_markers)

# MB_FRACTION_META %>%
#   left_join(anno_human) %>%
#   left_join(panglao,
#             by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#   mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
#          score = -log10(sumz_adj) * log2FoldChange_C1) %>%
#   group_by(`cell type`) %>%
#   filter(sumz_adj < 0.01) %>%
#   rowwise() %>%
#   mutate(log2FoldChange = mean(log2FoldChange_C1,
#                                log2FoldChange_C2,
#                                log2FoldChange_C3)) %>%
#   ungroup() %>%
#   # filter(`cell type` %in% c("Astrocytes",
#   #                           "Cholinergic neurons",
#   #                           "Dopaminergic neurons",
#   #                           "Endothelial cells",
#   #                           "GABAergic neurons",
#   #                           "Glutaminergic neurons",
#   #                           "Interneurons",
#   #                           "Microglia",
#   #                           "Oligodendrocytes",
#   #                           "Oligodendrocyte progenitor cells"
#   # )) %>%
# select("Cell Type" = `cell type`,
#        "Gene" = external_gene_name,
#        "Log2 Fold Change" = log2FoldChange,
#        "Adjusted P Value" = sumz_adj) %>%
#   group_by(`Cell Type`) %>%
#   summarise(lfc = mean(`Log2 Fold Change`)) %>% View

# plot cell type enrichment
p_mb_celltypes <- MB_FRACTION_META %>%
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
                            "Oligodendrocytes",
                            "Oligodendrocyte progenitor cells"
  )) %>%
  mutate(`cell type` = factor(`cell type`)) %>%
  mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
  ggplot(aes(x = prop_score,
             y = `cell type`,
             fill = `cell type`)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  labs(x = "Enrichment Score",
       y = "Cell Type") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.y = element_blank())



#
#
#
#
#
#
# TABLE 2: MB fraction GO search ------------------------------------------------

# do a GO search
# using all MB fraction tested genes as background
# and the enriched genes ordered by p value
# as input
MB_FRACTION_GO <- gost(MB_FRACTION_ENRICHED_GENES,
                       organism = "mmusculus",
                       # exclude_iea = T,
                       custom_bg = MB_FRACTION_META$ensembl_gene_id,
                       ordered_query = T)

# view the GO result
# ordered by pvalue
# MB_FRACTION_GO$result %>%
#   arrange(p_value) %>%
#   filter(source != "TF") %>%
#   View

# term_filter <- sapply(MB_FRACTION_GO$result$term_id, function(x){
#   str_detect(unlist(MB_FRACTION_GO$result$parents), x)
# })
#
#
# MB_FRACTION_GO$result[!rowSums(term_filter),] %>% View

# Extract the top 10 entries for KEGG,
# GO:BP and GO:CC
t_mb_ip_enrichment_go <- MB_FRACTION_GO$result %>%
  mutate(score = -log10(p_value)) %>%
  filter(source %in% c("KEGG", "GO:BP", "GO:CC")) %>%
  group_by(source) %>%
  slice_max(order_by = score, n = 10) %>%
  mutate(term_name = str_to_sentence(term_name)) %>%
  select(
    "Source" = source,
    "Term name" = term_name,
    "Term size" = term_size,
    "Intersection size" = intersection_size,
    "P value" = p_value)

# ggplot(aes(x = score,
#            y = term_name,
#            fill = source)) +
# geom_col(colour = "black") +
# facet_wrap(vars(source),
#            scales = "free_y",
#            nrow = 3) +
# scale_fill_d3() +
# labs(x = "-log10 P value") +
# theme(axis.title.y = element_blank(),
#       legend.position = "none",
#       strip.text = element_text(size = 14)
#       # axis.text = element_text(size = 12)
# )


#
#
# FIGURE 2
# MB FRACTION fgsea ----

# get entrez IDs
# anno_fgsea <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "entrezgene_id"),
#   filters = "ensembl_gene_id",
#   values = MB_FRACTION_META$ensembl_gene_id,
#   mart = ensembl,
#   uniqueRows = TRUE
# )

anno_fgsea <- readRDS("R/objects/anno_fgsea.rds")

# load gmt

# pathways <- fgsea::gmtPathways("R/objects/msigdb.v7.4.symbols.gmt.txt")

# pathways <- fgsea::gmtPathways("R/objects/c2.cp.kegg.v7.4.symbols.gmt.txt")
pathways <- fgsea::gmtPathways("R/objects/c5.go.v7.4.symbols.gmt.txt")
# pathways <- fgsea::gmtPathways("R/objects/c5.go.bp.v7.4.symbols.gmt.txt")
# pathways <- fgsea::gmtPathways("R/objects/c5.go.mf.v7.4.symbols.gmt.txt")

ranks <- MB_FRACTION_META %>%
  # filter(!conflict
  #          # ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES
  #        ) %>%
  mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
         score = ifelse(score == Inf, max(score[score != Inf])+log2FoldChange*100, score)) %>%
  mutate(rank = score) %>%
  inner_join(anno_human) %>%
  select(hsapiens_homolog_associated_gene_name,
         rank) %>%
  filter(hsapiens_homolog_associated_gene_name != "") %>%
  distinct(across(-rank), .keep_all = T) %>%
  deframe

MB_FRACTION_FGSEA <- fgsea(pathways,
                           ranks,
                           eps = 0,
                           minSize = 10,
                           maxSize = 500,
                           nPermSimple = 10000,
                           BPPARAM = MulticoreParam(),
                           nproc = 22)

p_mb_ip_enrichment_fgsea <- ggarrange(
  {plotEnrichment(pathways[["GOCC_RIBOSOME"]],
                 ranks) + labs(title="GO:CC Ribosome",
                               y = "Enrichment Score") +
    theme(axis.title.x = element_blank())},
  {plotEnrichment(pathways[["GOBP_MITOCHONDRIAL_GENE_EXPRESSION"]],
                 ranks) + labs(title="GO:BP Mitochondrial Gene\nExpression") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())},
  plotEnrichment(pathways[["GOBP_ENSHEATHMENT_OF_NEURONS"]],
                 ranks) + labs(title="GO:BP Ensheathment of Neurons",
                               y = "Enrichment Score",
                               x = "Rank"),
  plotEnrichment(pathways[["GOBP_NEGATIVE_REGULATION_OF_GLIOGENESIS"]],
                 ranks) + labs(title="GO:BP Regulation of Gliogenesis",
                               x = "Rank") +
    theme(axis.title.y = element_blank())
)





# MB FRACTION WGCNA ----

# This code could be used for an appendix figure about how the
# translated genes in somal TRAP samples are VERY correlated
# with one another and that correlation doesn't differ greatly
# so modules are similar to one another.


# vsd_mb_ip_old_transl <- vst(dds_C1_MB_IP_OLD[rownames(dds_C1_MB_IP_OLD) %in% MB_FRACTION_TRANSLATED_GENES,])
#
# library(WGCNA)
#
# options(stringsAsFactors = FALSE)
# enableWGCNAThreads()
#
# # prepare data
# datExpr0 <- as.data.frame(t(assay(vsd_mb_ip_old_transl)))
# rownames(datExpr0) <- colnames(vsd_mb_ip_old_transl)
# gsg <- goodSamplesGenes(datExpr0, verbose = 3)
# gsg$allOK
#
# sampleTree <- hclust(dist(datExpr0), method = "average")
# # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# # The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(12,9)
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
#
# datExpr <- datExpr0
#
# # Choose a set of soft-thresholding powers
# powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# # Call the network topology analysis function
# sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#
# net = blockwiseModules(datExpr, power = 5,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "MB_OLD_IP_TRANSL_C1",
#                        verbose = 3)
#
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)






# MB FRACTION PubMed Query ----
# entrez_dbs()
# entrez_db_searchable("pubmed")
# entrez_db_searchable("gene")

# # get pubmed record
# get_pubmed <- function(input_gene){
#
#   gene <- entrez_search(db = "gene", term = paste0(input_gene, "[PREF] AND 9606[TID]"))
#   gene$ids
#   if(length(gene$ids) == 0){
#     gene <- entrez_search(db = "gene", term = paste0(input_gene, "[PREF] AND 10090[TID]"))
#   }
#   gene$ids
#   if(length(gene$ids) == 0){
#     gene$name <- input_gene
#     query <- paste0(c(gene$name, "(dopamin* OR parkinson*)"), collapse = " AND ")
#   } else {
#     gene <- entrez_summary(db = "gene", id = gene$ids[1])
#     gene$name <- ifelse(nchar(gene$name) <= 2, gene$description, gene$name)
#
#     if(gene$otheraliases != ""){
#       gene$otheraliases <- sapply(unlist(strsplit(gene$otheraliases, ", ")), function(x){ifelse(nchar(x) <= 2, "", x)})
#       # gene$otheraliases <- sapply(gene$otheraliases, function(x){ifelse(str_count(x, "[0-9]") >= 2 & nchar(x) == 3, "", x)})
#       gene$otheraliases <- sapply(gene$otheraliases, function(x){ifelse(str_detect(x, "[\\[]"), "", x)})
#       gene$otheraliases <- ifelse(str_detect(gene$name, "Rik"), "", gene$otheraliases)
#
#       query <- paste0(c(paste0(c(gene$name, unlist(strsplit(gene$otheraliases, ", "))
#                                  # unlist(strsplit(gene$description, ", "))
#       ), collapse = " OR "), "(dopamin* OR parkinson*)"), collapse = " AND ")
#     } else {
#       query <- paste0(c(gene$name, "(dopamin* OR parkinson*)"), collapse = " AND ")
#     }
#   }
#
#   print(query)
#
#   r_search <- entrez_search(db = "pubmed", term = query)
#
#   print(r_search)
#   return(r_search)
# }
#
# # get number of publications for multiple records
# get_pubmed_n <- function(pubmed){
#   sapply(pubmed, function(x){x[["count"]]}) %>% unname
# }
#
# # get titles of publications
# get_pubmed_titles <- function(pubmed){
#   summary <- entrez_summary(db = "pubmed",
#                             id = pubmed$ids)
#   if(length(pubmed$ids) > 1){
#     titles <- sapply(summary, function(x){x[["title"]]}) %>% unname
#   } else {
#     titles <- summary[["title"]]
#   }
#   return(titles)
# }

# get number of publications for individual records
# get_pubmed("S100a10")$count
# get_pubmed("9030622O22Rik")$count
# get_pubmed("F2r")$count
# get_pubmed("Pbcd1")$count
# get_pubmed("Cpne7")$count
#
# get_pubmed("Fam210b") %>%
#   get_pubmed_titles()
#
# get_pubmed("Syt17") %>%
#   get_pubmed_titles()
#
#
#
#
# get_pubmed("Cpne7")$count


# retrieve publication number for list of genes (takes a long time)
# MB_FRACTION_META_SIMPLE_PUBMED <- MB_FRACTION_META_SIMPLE %>%
#   filter(external_gene_name %in% sapply(MB_FRACTION_ENRICHED_GENES, get_external_gene_name) &
#            log2FoldChange >= 2) %>%
#   # left_join(as_tibble(pubmed_n,
#   #                     rownames = "external_gene_name")) %>%
#   arrange(desc(score), desc(log2FoldChange)) %>%
#   mutate(score = ifelse(score == Inf, max(score[score != Inf])+log2FoldChange*100, score)) %>% # rank by p value and l2fc
#   arrange(desc(score)) %>%
#   mutate(rank_enrichment = row_number()) %>%
#   mutate(publications = sapply(external_gene_name, get_pubmed_n)) %>%
#   arrange(desc(publications)) %>%
#   mutate(rank_publications = row_number())

# saveRDS(MB_FRACTION_META_SIMPLE_PUBMED,
#         "R/objects/MB_FRACTION_META_SIMPLE_PUBMED.rds")
MB_FRACTION_META_SIMPLE_PUBMED <- readRDS("R/objects/MB_FRACTION_META_SIMPLE_PUBMED.rds") %>%
  mutate(publications = ifelse(external_gene_name == "Fam210b", publications-4, publications),
         publications = ifelse(external_gene_name == "Syt17", publications-1, publications)) # Remove spurious papers


plot_data <- MB_FRACTION_META_SIMPLE_PUBMED %>%
  filter(external_gene_name %in% sapply(MB_FRACTION_ENRICHED_GENES, get_external_gene_name)) %>%
  arrange(desc(rank_enrichment)) %>%
  # filter(rank_enrichment < 50) %>%
  mutate(external_gene_name = factor(external_gene_name),
         low_info = publications <= 5)

low_info_genes <- plot_data %>%
  filter(low_info) %>%
  pull(external_gene_name)

# FIGURE A: MB fraction ranking by enrichment magnitude -----------------------------

# plot enrichment cutoff
# create the data (add a row number as a rank of lfc)
plot_data <- MB_FRACTION_META %>%
  filter(baseMean_C1 > 100) %>%
  filter(sumz_adj < 0.01 &
           # !conflict &
           log2FoldChange > 0) %>%
  arrange(desc(score), desc(log2FoldChange)) %>%
  mutate(score = ifelse(score == Inf, max(score[score != Inf])+log2FoldChange*100, score)) %>% # rank by p value and l2fc
  arrange(desc(score)) %>%
  mutate(rn = row_number(),
         low_info = external_gene_name %in% low_info_genes)

# count the number of genes with l2fc >= 2
# plot_nr <- plot_data %>%
#   filter(log2FoldChange >= 2) %>%
#   summarise(n = n()) %>%
#   mutate(cutoff = 2)

# plot the l2fc/rank
# pos <- position_jitter(width = 0.1, seed = 1, height = 0)

p_mb_ip_enrichment_rank_cutoff <- plot_data %>%
  arrange(log2FoldChange) %>%
  filter(rn < 75) %>%
  mutate(external_gene_name = ifelse(external_gene_name %in% c("Th",
                                                               "Slc6a3",
                                                               "Cpne7",
                                                               "Syt17",
                                                               "Tagln3",
                                                               "Fgf20",
                                                               "Snca",
                                                               "Sncg",
                                                               "Fam210b"),
                                     external_gene_name,
                                     "")) %>%
  ggplot(aes(x = rn,
             y = score,
             color = low_info)) +
  geom_point(size = 1.5) +
              # shape = 21,
             # position = pos) +
  # geom_text(data = plot_nr,
  #           aes(x = 15,
  #               y = 2,
  #               vjust = -0.5,
  #               label = paste0(n, " genes > 4-fold enriched")),
  #           size = 6) +
  geom_label_repel(aes(x = rn, y = score,
                       label = external_gene_name),
                   # position = pos,
                   colour = "black",
                   # arrow = arrow(length = unit(0.02, "npc")),
                   # nudge_x = 0.15,
                   box.padding = 1,
                   max.overlaps = 100,
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20) +
  scale_color_d3() +
  # scale_x_log10() +
  # scale_y_continuous(breaks = seq(0, 10, 2)) +
  # geom_hline(yintercept = plot_nr$cutoff,
  #            linetype = "dotted",
  #            colour = "black") +
  # geom_vline(xintercept = plot_nr$n,
  #            linetype = "dotted",
  #            colour = "black") +
  labs(x = "Rank",
       y = "Enrichment Score",
       color = "Limited Information") +
  theme(legend.position = "top")

n_MB_FRACTION_ENRICHED_HIGH <- MB_FRACTION_META %>% filter(log2FoldChange >= 2 & sumz_adj < 0.01) %>% nrow

# the greater the magnitude of median l2fc, the lower the proportion of conflicts

# FIGURE B: MB fraction Pubmed number of publications -------------------------------------------

plot_data <- MB_FRACTION_META_SIMPLE_PUBMED %>%
  filter(external_gene_name %in% sapply(MB_FRACTION_ENRICHED_GENES, get_external_gene_name)) %>%
  arrange(desc(rank_enrichment)) %>%
  filter(rank_enrichment < 50) %>%
  mutate(external_gene_name = factor(external_gene_name),
         low_info = publications <= 5)

# plot the publications for top 20 enriched genes
p_mb_ip_enrichment_pubmed_n <- plot_data %>%
  ggplot(aes(y = fct_reorder(external_gene_name, -rank_enrichment),
             x = publications,
             fill = rank_enrichment)) +
  geom_col(colour = "black") +
  scale_x_log10(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Number of Publications",
       y = "Gene (Ordered by enrichment)") +
  scale_fill_distiller(palette = "Blues") +
  theme(legend.position = "none",
        axis.text.y = element_text(face = ifelse(plot_data$low_info == TRUE, "bold", "plain"),
                                   colour = ifelse(plot_data$low_info == TRUE, pal_d3()(2)[2], pal_d3()(1)),
                                   hjust = 1,
                                   size = 8)) +
  geom_text(aes(label = publications),
            hjust = -0.3,
            size = 2.5)

# ----
# ----
# ----
# MB fraction ranking by enrichment magnitude - GO ------------------------

# # get GO anno + biotype
# anno_go <- getBM(
#   attributes = c(
#     "ensembl_gene_id",
#     "external_gene_name",
#     "description",
#     "gene_biotype",
#     "name_1006",
#     "definition_1006"
#
#   ),
#   filters = "ensembl_gene_id",
#   values = rownames(dds),
#   mart = ensembl,
#   uniqueRows = TRUE
# )

anno_go <- readRDS("R/objects/anno_go.rds")

# Rank the genes by fold enrichment and assign GO terms
# tally up the GO terms and calculate the mean enrichment per term
t_mb_ip_enrichment_max_go_rank <- MB_FRACTION_META %>%
  filter(sumz_adj < 0.01) %>%
  filter(log2FoldChange >= 2) %>%
  left_join(anno_go) %>%
  group_by(name_1006) %>%
  mutate(n = n()) %>%
  filter(name_1006 != "") %>%
  arrange(name_1006) %>%
  select(ensembl_gene_id,
         external_gene_name,
         description,
         name_1006,
         definition_1006,
         log2FoldChange) %>%
  group_by(name_1006) %>%
  summarise(total = n(),
            enrichment = mean(log2FoldChange)) %>%
  arrange(desc(enrichment)) %>%
  filter(total >= 5) %>%
  mutate(name_1006 = str_to_sentence(name_1006),
         name_1006 = str_replace(name_1006, "Dna|dna", "DNA"),
         name_1006 = str_replace(name_1006, "rna", "RNA"),
         name_1006 = str_replace(name_1006, "ii", "II")) %>%
  rename("GO term" = name_1006,
         "Number of Genes" = total,
         "Mean Fold-Enrichment (Log2)" = enrichment)

# see the full term-joined list
# N.B. this list excludes genes that have no GO term
MB_FRACTION_META_GO_INDIV <- MB_FRACTION_META_SIMPLE %>%
  filter(sumz_adj < 0.01) %>%
  filter(log2FoldChange >= 2) %>%
  left_join(anno_go) %>%
  group_by(name_1006) %>%
  mutate(n = n()) %>%
  filter(name_1006 != "") %>%
  arrange(name_1006) %>%
  select(ensembl_gene_id,
         external_gene_name,
         description,
         name_1006,
         definition_1006,
         log2FoldChange,
         score) %>%
  arrange(desc(score), desc(log2FoldChange))

# view the DNA binding elements
# MB_FRACTION_META_GO_INDIV %>%
#   filter(str_detect(name_1006, "DNA")) %>%
#   ungroup() %>%
#   select(-c(name_1006,
#             definition_1006)) %>%
#   distinct() %>%
#   View
#
# FIGURE 2: MB enrichment rank and publication number
#
#
# FIGURE 3
# MB fraction ISH and STRING -----------------------
# Slc10a4 - No longer - too much is known about it!
# Tagln3
# Cpne7
# Syt17
# Atp6v1g1
# ######## Slc10a4
# # ABA
# aba_slc10a4 <- get_aba("Slc10a4",
#                        region = "snpc",
#                        experiment = 1,
#                        downsample = 3,
#                        plane = "coronal",
#                        orientation = "antisense")
# center_aba("Slc10a4",
#            image_id = aba_slc10a4,
#            x = c(.35, .1),
#            y = c(0.2, 0.15))
#
# img <- image_read("R/objects/ABA/Slc10a4_centered.jpg")
# p_slc10a4 <- ggdraw() +
#   draw_image(img)
#
# # STRING
# slc10a4_neighbours <- get_string_neighbours("Slc10a4")
# p_slc10a4_string <- base2grob( ~ plot_neighbours("SLC10A4", slc10a4_neighbours))

######## Tagln3
# ABA
aba_tagln3 <- get_aba("Tagln3",
                       region = "snpc",
                       experiment = 1,
                       downsample = 3,
                       plane = "coronal",
                       orientation = "antisense")
center_aba("Tagln3",
           image_id = aba_tagln3,
           x = c(.32, .1),
           y = c(0.1, 0.25))

img <- image_read("R/objects/ABA/Tagln3_centered.jpg")
p_tagln3 <- ggdraw() +
  draw_image(img)

# STRING
tagln3_neighbours <- get_string_neighbours("Tagln3")
p_tagln3_string <- base2grob( ~ plot_neighbours("TAGLN3", tagln3_neighbours))

######## Cpne7
# ABA
aba_cpne7 <- get_aba("Cpne7",
                     region = "snpc",
                     experiment = 1,
                     downsample = 3,
                     plane = "coronal",
                     orientation = "antisense")
center_aba("Cpne7",
           image_id = aba_cpne7,
           x = c(.35, .1),
           y = c(0.15, 0.2))

img <- image_read("R/objects/ABA/Cpne7_centered.jpg")
p_cpne7 <- ggdraw() +
  draw_image(img)

# STRING
cpne7_neighbours <- get_string_neighbours("Cpne7")
p_cpne7_string <- base2grob( ~ plot_neighbours("CPNE7", cpne7_neighbours))

######## Syt17
# ABA
aba_syt17 <- get_aba("Syt17",
                     region = "snpc",
                     experiment = 1,
                     downsample = 3,
                     plane = "coronal",
                     orientation = "antisense")
center_aba("Syt17",
           image_id = aba_syt17,
           x = c(.35, .1),
           y = c(0.1, 0.25))

img <- image_read("R/objects/ABA/Syt17_centered.jpg")
p_syt17 <- ggdraw() +
  draw_image(img)

# STRING - All interactors are just putative association with Syt1 - these are false positives - discuss.
syt17_neighbours <- get_string_neighbours("Syt17")
# syt17_neighbours <- syt17_neighbours %>%
#   filter(combined_score > 750)
p_syt17_string <- base2grob( ~ plot_neighbours("SYT17", syt17_neighbours))

# p_syt17_string <- ggplot() +
#   annotate("text", x = 1, y = 1, size = 8, label = paste("No Interactors")) +
#   theme_void()

######## Fam210b
# ABA
aba_fam210b <- get_aba("Fam210b",
                       region = "midbrain",
                       experiment = 1,
                       downsample = 3,
                       plane = "sagittal",
                       orientation = "antisense")
center_aba("Fam210b",
           image_id = aba_fam210b,
           x = c(0.1, .07),
           y = c(-0.1, 0.5))

img <- image_read("R/objects/ABA/Fam210b_centered.jpg")
p_fam210b <- ggdraw() +
  draw_image(img)

# STRING
fam210b_neighbours <- get_string_neighbours("Fam210b")
p_fam210b_string <- base2grob( ~ plot_neighbours("FAM210B", fam210b_neighbours))


# Other genes
{

  ######## Atp6v1g1
  # ABA
  # aba_atp6v1g1 <- get_aba("Atp6v1g1",
  #                         region = 749,
  #                         experiment = 1,
  #                         downsample = 3,
  #                         plane = "sagittal",
  #                         orientation = "antisense")
  # center_aba("Atp6v1g1",
  #            image_id = aba_atp6v1g1,
  #            x = c(0, 0.15),
  #            y = c(0.2, .15))
  #
  # img <- image_read("R/objects/ABA/Atp6v1g1_centered.jpg")
  # p_atp6v1g1 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # atp6v1g1_neighbours <- get_string_neighbours("Atp6v1g1")
  # p_atp6v1g1_string <- base2grob( ~ plot_neighbours("ATP6V1G1", atp6v1g1_neighbours))

  # ######## Lix1
  # # ABA
  # aba_lix1 <- get_aba("Lix1",
  #                     region = 321,
  #                     experiment = 1,
  #                     downsample = 3,
  #                     plane = "coronal",
  #                     orientation = "antisense")
  # center_aba("Lix1",
  #            image_id = aba_lix1,
  #            x = c(0.5, 0.05),
  #            y = c(0, 0.3))
  #
  # img <- image_read("R/objects/ABA/Lix1_centered.jpg")
  # p_lix1 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # lix1_neighbours <- get_string_neighbours("Lix1")
  # p_lix1_string <- base2grob( ~ plot_neighbours("LIX1", lix1_neighbours))
  #
  # ######## Fgf13
  # # ABA
  # aba_fgf13 <- get_aba("Fgf13",
  #                      region = 330,
  #                      experiment = 1,
  #                      downsample = 3,
  #                      plane = "coronal",
  #                      orientation = "antisense")
  # center_aba("Fgf13",
  #            image_id = aba_fgf13,
  #            x = c(0.32, 0.2),
  #            y = c(0.9, -0.2))
  #
  # img <- image_read("R/objects/ABA/Fgf13_centered.jpg")
  # p_fgf13 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # fgf13_neighbours <- get_string_neighbours("Fgf13")
  # p_fgf13_string <- base2grob( ~ plot_neighbours("FGF13", fgf13_neighbours))
  #
  # ######## Tagln3
  # # ABA
  # aba_tagln3 <- get_aba("Tagln3",
  #                       region = 374,
  #                       experiment = 1,
  #                       downsample = 3,
  #                       plane = "coronal",
  #                       orientation = "antisense")
  # center_aba("Tagln3",
  #            image_id = aba_tagln3,
  #            x = c(0.35, 0.15),
  #            y = c(0.07, .25))
  #
  # img <- image_read("R/objects/ABA/Tagln3_centered.jpg")
  # p_tagln3 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # tagln3_neighbours <- get_string_neighbours("Tagln3")
  # p_tagln3_string <- base2grob( ~ plot_neighbours("TAGLN3", tagln3_neighbours))
  #
  # ######## Bex2 - May not present as high background, but
  # ######## definitely discuss as plays and important
  # ######## role in autophagy
  # # ABA
  # aba_bex2 <- get_aba("Bex2",
  #                     region = 330,
  #                     experiment = 1,
  #                     downsample = 3,
  #                     plane = "sagittal",
  #                     orientation = "antisense")
  # center_aba("Bex2",
  #            image_id = aba_bex2,
  #            x = c(0, 0.2),
  #            y = c(0.7, -0.1))
  #
  # img <- image_read("R/objects/ABA/Bex2_centered.jpg")
  # p_bex2 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # bex2_neighbours <- get_string_neighbours("Bex2")
  # p_bex2_string <- base2grob( ~ plot_neighbours("BEX2", bex2_neighbours))
  #
  # ######## Fam167a
  # # ABA
  # aba_fam167a <- get_aba("Fam167a",
  #                        region = 749,
  #                        experiment = 1,
  #                        downsample = 3,
  #                        plane = "coronal",
  #                        orientation = "antisense")
  # center_aba("Fam167a",
  #            image_id = aba_fam167a,
  #            x = c(0.32, 0.2),
  #            y = c(0.13, .15))
  #
  # img <- image_read("R/objects/ABA/Fam167a_centered.jpg")
  # p_fam167a <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # fam167a_neighbours <- get_string_neighbours("Fam167a")
  # p_fam167a_string <- base2grob( ~ plot_neighbours("FAM167A", fam167a_neighbours))
  #
  # ######## Anxa1
  # # ABA
  # aba_anxa1 <- get_aba("Anxa1",
  #                      region = 331,
  #                      experiment = 1,
  #                      downsample = 3,
  #                      plane = "coronal",
  #                      orientation = "antisense")
  # center_aba("Anxa1",
  #            image_id = aba_anxa1,
  #            x = c(0.25, 0.2),
  #            y = c(0.25, .15))
  #
  # img <- image_read("R/objects/ABA/Anxa1_centered.jpg")
  # p_anxa1 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # anxa1_neighbours <- get_string_neighbours("Anxa1")
  # p_anxa1_string <- base2grob( ~ plot_neighbours("ANXA1", anxa1_neighbours))
  #
  # ######## C130021I20Rik
  # # ABA
  # aba_c130021I20Rik <- get_aba("C130021I20Rik",
  #                              region = 749,
  #                              experiment = 1,
  #                              downsample = 3,
  #                              plane = "coronal",
  #                              orientation = "antisense")
  # center_aba("C130021I20Rik",
  #            image_id = aba_c130021I20Rik,
  #            x = c(0.3, 0.15),
  #            y = c(0.1, .2))
  #
  # img <- image_read("R/objects/ABA/C130021I20Rik_centered.jpg")
  # p_c130021I20Rik <- ggdraw() +
  #   draw_image(img)
  #
  # # NO STRING
  # p_c130021I20Rik_string <- ggdraw() +
  #   draw_text("ncRNA: \nNo STRING Interactors")
  # # https://pubmed.ncbi.nlm.nih.gov/24066094/
  # # Interaction with Lmx1a
  #
  # ######## A2ml1
  # ######## No ABA data
  # p_a2ml1 <- ggdraw() +
  #   draw_text("No ABA Data")
  # # STRING
  # a2ml1_neighbours <- get_string_neighbours("A2ml1")
  # p_a2ml1_string <- base2grob( ~ plot_neighbours("A2ML1", a2ml1_neighbours))
  #
  # ######## Kcns3
  # # ABA
  # aba_kcns3 <- get_aba("Kcns3",
  #                      region = 749,
  #                      experiment = 1,
  #                      downsample = 3,
  #                      plane = "coronal",
  #                      orientation = "antisense")
  # center_aba("Kcns3",
  #            image_id = aba_kcns3,
  #            x = c(0.3, 0.2),
  #            y = c(0.1, .2))
  #
  # img <- image_read("R/objects/ABA/Kcns3_centered.jpg")
  # p_kcns3 <- ggdraw() +
  #   draw_image(img)
  # # STRING
  # kcns3_neighbours <- get_string_neighbours("Kcns3")
  # p_kcns3_string <- base2grob( ~ plot_neighbours("KCNS3", kcns3_neighbours))
  #
  #
  # ######## Gng3
  # # ABA
  # aba_gng3 <- get_aba("Gng3",
  #                     region = 321,
  #                     experiment = 1,
  #                     downsample = 3,
  #                     plane = "coronal",
  #                     orientation = "antisense")
  # center_aba("Gng3",
  #            image_id = aba_gng3,
  #            x = c(0.45, 0),
  #            y = c(0, .3))
  #
  # img <- image_read("R/objects/ABA/Gng3_centered.jpg")
  # p_gng3 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # gng3_neighbours <- get_string_neighbours("Gng3")
  # p_gng3_string <- base2grob( ~ plot_neighbours("GNG3", gng3_neighbours))
  #
  # ######## Prmt2
  # # ABA
  # aba_prmt2 <- get_aba("Prmt2",
  #                      region = 331,
  #                      experiment = 1,
  #                      downsample = 3,
  #                      plane = "coronal",
  #                      orientation = "antisense")
  # center_aba("Prmt2",
  #            image_id = aba_gng3,
  #            x = c(0.45, 0),
  #            y = c(0.05, 0.5))
  #
  # img <- image_read("R/objects/ABA/Prmt2_centered.jpg")
  # p_prmt2 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # prmt2_neighbours <- get_string_neighbours("Prmt2")
  # p_prmt2_string <- base2grob( ~ plot_neighbours("PRMT2", prmt2_neighbours))
  #
  # ######## Lbh
  # # ABA
  # aba_lbh <- get_aba("Lbh",
  #                    region = 313,
  #                    experiment = 1,
  #                    downsample = 3,
  #                    plane = "sagittal",
  #                    orientation = "antisense")
  # center_aba("Lbh",
  #            image_id = aba_gng3,
  #            x = c(-0.2, 0.4),
  #            y = c(-0.15, 0.5))
  #
  # img <- image_read("R/objects/ABA/Lbh_centered.jpg")
  # p_lbh <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # lbh_neighbours <- get_string_neighbours("Lbh")
  # p_lbh_string <- base2grob( ~ plot_neighbours("LBH", lbh_neighbours))
  #
  # ######## S100a10
  # # ABA
  # aba_s100a10 <- get_aba("S100a10",
  #                        region = 374,
  #                        experiment = 1,
  #                        downsample = 3,
  #                        plane = "coronal",
  #                        orientation = "antisense")
  # center_aba("S100a10",
  #            image_id = aba_s100a10,
  #            x = c(0.35, 0.1),
  #            y = c(0.1, 0.5))
  #
  # img <- image_read("R/objects/ABA/S100a10_centered.jpg")
  # p_s100a10 <- ggdraw() +
  #   draw_image(img)
  #
  # # STRING
  # s100a10_neighbours <- get_string_neighbours("S100a10")
  # p_s100a10_string <- base2grob( ~ plot_neighbours("S100A10", s100a10_neighbours))
  #
  # # MB fraction high enrichment low publication high background -------------
  #
  # # Bex2 - autophagy
  # # Tuba1a - cytoskeletal
  # # 6330403K07Rik - ncRNA?
  #
  #
  # # MB fraction high enrichment low expression ------------------------------
  #
  # # Sult5a1
  # # Ctxn2 - No ABA signal
  #
  # # MB fraction high enrichment medium publications -------------------------
  #
  # ######## Cadps2 - a lot is known within a few papers
  # get_pubmed("Cadps2") %>%
  #   get_pubmed_titles()
  #
  # pubmed <- get_pubmed("Gng3")
  # entrez_summary(db = "pubmed",
  #                id = pubmed$ids[1])$title
  # ######## Pbx3 - 1 paper, but Pbx1 very published
  # ######## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5553812/
  #
  #
  # # MB fraction high enrichment low publications no information -------------
  #
  # # Pbdc1
  # # ABA data available
  #
  #
  # # MB fraction example reserve ---------------------------
  #
  #
  # ```{r mb-ip-enrichment-second, fig.asp=1.3}
  # ggarrange(
  #   ggarrange(p_fam210b,
  #             p_fam210b_string,
  #             nrow = 1),
  #   ggarrange(p_lix1,
  #             p_lix1_string,
  #             nrow = 1),
  #   ggarrange(p_fgf13,
  #             p_fgf13_string,
  #             nrow = 1),
  #   ggarrange(p_tagln3,
  #             p_tagln3_string,
  #             nrow = 1),
  #   nrow = 4,
  #   labels = c("Fam210b",
  #              "Lix1",
  #              "Fgf13",
  #              "Tagln3"))
  # ```
  #
  # ```{r mb-ip-enrichment-third, fig.asp=1.3}
  # ggarrange(
  #   ggarrange(p_fam167a,
  #             p_fam167a_string,
  #             nrow = 1),
  #   ggarrange(p_anxa1,
  #             p_anxa1_string,
  #             nrow = 1),
  #   ggarrange(p_c130021I20Rik,
  #             p_c130021I20Rik_string,
  #             nrow = 1),
  #   ggarrange(p_kcns3,
  #             p_kcns3_string,
  #             nrow = 1),
  #   nrow = 4,
  #   labels = c("Fam167a",
  #              "Anxa1",
  #              "C130021I20Rik",
  #              "Kcns3"))
  #
  # # A2ml1 removed because no ABA data
  # ```
  #
  # ```{r mb-ip-enrichment-fourth, fig.asp=1.3}
  # ggarrange(
  #   ggarrange(p_gng3,
  #             p_gng3_string,
  #             nrow = 1),
  #   ggarrange(p_prmt2,
  #             p_prmt2_string,
  #             nrow = 1),
  #   ggarrange(p_lbh,
  #             p_lbh_string,
  #             nrow = 1),
  #   ggarrange(p_s100a10,
  #             p_s100a10_string,
  #             nrow = 1),
  #   nrow = 4,
  #   labels = c("Gng3",
  #              "Prmt2",
  #              "Lbh",
  #              "S100a10"))
  #
  # ```
}
#
#
#
#


# -------------
# -------------
# ----
# FIND THE DOMINANT ISOFORMS: MB ----

# Use the ONT data for this: Will have more accurate transcript-level estimates
# But have to use short read data for FRACTION DTU analysis, as don't have TOTAL long-read yet

plot_data <- tibble(
  ensembl_gene_id = substr(ont_tx_counts$ensembl_gene_id, 1, 18),
  ensembl_transcript_id = substr(ont_tx_counts$ensembl_transcript_id, 1, 18),
  count = rowSums(ont_tx_counts[,3:ncol(ont_tx_counts)])) %>%
  group_by(ensembl_gene_id) %>%
  mutate(prop = count/sum(count),
         n = n()) %>%
  filter(n > 1 &
           sum(count > 10)) %>%
  arrange(n,
          ensembl_gene_id,
          desc(prop)) %>%
  filter(n <= 5) %>%
  mutate(rank = factor(row_number()),
         n = paste(n, "transcripts"))

p_mb_fraction_isoform_diversity <- plot_data %>%
  ggplot(aes(x = rank,
             y = prop,
             fill = factor(n))) +
  geom_violin(draw_quantiles = c(0.5),
              scale = "width") +
  facet_wrap(vars(n),
             scales = "free_x",
             nrow = 1) +
  scale_fill_manual(values = pal_d3()(5)[c(1, 2, 3, 5)]) +
  theme(legend.position = "none") +
  labs(x = "Transcript",
       y = "Proportion of \nTotal Expression") +
  panel_border()

n_median_primary_transcript_proportion <- plot_data %>%
  filter(rank == 1) %>%
  pull(prop) %>%
  mean*100
n_median_primary_transcript_proportion <- signif(n_median_primary_transcript_proportion, 3)

# MB Fraction: DTU ----


# Keep this code, but start from the readRDS function (the commented code has already been run)
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
#   str_extract(files, "(?<=salmon\\/)[0-9|\\_]+(?=\\/quant.sf)")
#
# files <- files[names(files) %in% metadata$fastq_id] # remove blanks
# files <-
#   files[order(match(names(files), metadata$fastq_id))] #reorder based on metadata
#
# metadata$files_salmon <- files
# metadata$names <- metadata$fastq_id
#
# txdb <-
#   makeTxDbFromGFF(file = "/zfs/analysis/trap/active/snakemake/thesis_snakemake/input/index/references/annotation.gtf",
#                   dataSource = "GENCODE",
#                   organism = "Mus musculus")
#
# tx2gene <-
#   AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")

# # import counts
# txi <- tximport(
#   files,
#   txOut = TRUE,
#   type = "salmon",
#   tx2gene = tx2gene,
#   ignoreTxVersion = F,
#   countsFromAbundance = "dtuScaledTPM"
# )

# # extract counts
# cts <- txi$counts
#
# # remove txi
# rm(txi)
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
# # remove duplicated features
# length(rownames(cts)) == length(unique(rownames(cts)))

# IF YOU WANT TO START FROM HERE, YOU CAN
# cts <- readRDS("R/objects/cts.rds")
#
# # create a data.frame for DRIMSeq
# counts <- data.frame(gene_id = tx2gene$GENEID,
#                      feature_id = tx2gene$TXNAME,
#                      cts,
#                      check.names = F) # it changes hyphens to full stops! C2 TOTAL names get changed
#
# # filter counts to remove depleted (matching the substring name - version number removed temporarily)
# counts <- counts[!substr(counts$gene_id, 1, 18) %in% MB_FRACTION_DEPLETED_GENES,]
#
# # filter samples
# # MB_YOUNG_WT_MIXED_C2.TOTAL1.5
# # MB_YOUNG_WT_MIXED_C2.TOTAL1-5
#
# # Use cohort 1 for IP vs TOTAL comparison
# counts <-
#   counts[, colnames(counts) %in% c("gene_id",
#                                    "feature_id",
#                                    (
#                                      metadata %>% filter(
#                                        !collection %in% c("C1.IPPOOL") &
#                                          cohort == "C1")
#                                      ) %>% pull(sample_name) %>% unique
#                                    )]
#
# # check dimensions
# dim(counts)
#
# # create sample metadata for DRIMSeq
# samps <- metadata %>%
#   filter(sample_name %in% colnames(counts)) %>%
#   select("sample_id" = sample_name,
#          fraction) %>%
#   distinct
# samps <- as.data.frame(samps)
# samps$fraction <- relevel(samps$fraction, ref = "TOTAL")
#
# # create drimseq object
# d <- dmDSdata(counts = counts,
#               samples = samps)
#
# methods(class = class(d))
# #
# substr(counts(d)$gene_id, 1, 18) %>% unique
#
# # DRIMSeq filter
# n <- nrow(samps)
#
# table(DRIMSeq::samples(d)$fraction)
# n.small <- sum(samps$fraction == "TOTAL")
# d <- dmFilter(
#   d,
  # min_samps_gene_expr = n * 0.75,
  # min_gene_expr = 10,
  # min_samps_feature_expr = n.small,
  # min_feature_expr = 10,
  # min_samps_feature_prop = n.small,
  # min_feature_prop = 0.1,
# )

# saveRDS(d,
#         "R/objects/d_MB_FRACTION.rds")

# code to run on isolated environment
# d <- readRDS("R/objects/d_MB_FRACTION.rds")
#
# design_full <-
#   model.matrix( ~ fraction, data = DRIMSeq::samples(d))
# colnames(design_full)
#
# set.seed(821196)
# system.time({
#   d <-
#     dmPrecision(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
#   d <-
#     dmFit(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
#   d <-
#     dmTest(d, coef = "fractionIP", BPPARAM = BiocParallel::MulticoreParam())
# })
#
# saveRDS(d,
#         "R/objects/d_MB_FRACTION.rds")

# load the drimseq processed object: MB FRACTION (Cohort 1)
d_MB_FRACTION <- readRDS("R/objects/d_MB_FRACTION.rds")
# plot the number of transcripts per gene
counts_MB_FRACTION_DTU <- readRDS("R/objects/counts_MB_FRACTION_DTU.rds")

single_isoform_genes <- substr(names(table(counts_MB_FRACTION_DTU$gene_id)[table(counts_MB_FRACTION_DTU$gene_id) == 1]), 1, 18)
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

p_mb_fraction_dtu_enrichment <- plot_data %>%
  ggplot(aes(x = enriched,
             fill = signif)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  labs(x = "Gene-level Enrichment",
       y = "Number of Genes",
       fill = "DTU") +
  theme(legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  stat_count(geom = "text", colour = "white", size = 3.5,
             fontface = 2,
             aes(label = ..count..),
             position=position_stack(vjust=0.5))

n_dtu_not_enriched <- plot_data %>%
  filter(enriched == "Not enriched" & signif) %>%
  nrow

n_dtu_mb_fraction <- plot_data %>%
  filter(signif) %>%
  nrow

MB_FRACTION_DTU_NOT_ENRICHED_GENES <- plot_data %>%
  filter(enriched == "Not enriched" & signif) %>%
  pull(ensembl_gene_id)

MB_FRACTION_DTU_GENES <- plot_data %>%
  filter(signif) %>%
  pull(ensembl_gene_id)

# Plot a DTU example: Cdc42

data <- counts(d_MB_FRACTION) %>%
  as_tibble %>%
  filter(str_detect(gene_id, "ENSMUSG00000006699")) %>%
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

# STAGER
# get annotation
anno_drimseq <- get_anno(res_MB_FRACTION$gene_id_simple, get_human = TRUE)

# stageR
# create a vector of gene pvalues, stripped of version number
pScreen_MB_FRACTION <- res_MB_FRACTION$pvalue
strp <- function(x) substr(x,1,18)
names(pScreen_MB_FRACTION) <- strp(res_MB_FRACTION$gene_id)

# create a vector of transcript pvalues
pConfirmation_MB_FRACTION <- matrix(res.txp_MB_FRACTION$pvalue, ncol=1)
rownames(pConfirmation_MB_FRACTION) <- strp(res.txp_MB_FRACTION$feature_id)

# create a df with transcript and gene ids (tx2gene)
tx2gene <- res.txp_MB_FRACTION[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# run stageR at alpha of 0.05
stageRObj_MB_FRACTION <- stageRTx(pScreen=pScreen_MB_FRACTION, pConfirmation=pConfirmation_MB_FRACTION,
                                  pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj_MB_FRACTION <- stageWiseAdjustment(stageRObj_MB_FRACTION, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj_MB_FRACTION <- getAdjustedPValues(stageRObj_MB_FRACTION, order=FALSE,
                                              onlySignificantGenes=TRUE)
})
drim.padj_MB_FRACTION <- drim.padj_MB_FRACTION %>%
  as_tibble %>%
  left_join(anno,
            by = c("geneID" = "ensembl_gene_id"))

stat.test <- drim.padj_MB_FRACTION %>%
  filter(external_gene_name == "Cdc42") %>%
  left_join(anno_tx, by = c("txID" = "ensembl_transcript_id")) %>%
  select(external_transcript_name, gene, transcript)

stat.test <- stat.test %>%
  mutate(group1 = "TOTAL", group2 = "TRAP") %>%
  arrange(external_transcript_name) %>%
  mutate(gene = as.character(signif(gene, 3)), 
         transcript = as.character(signif(transcript, 3))) %>%
  pval_asterisks(transcript)

p_mb_fraction_dtu_cdc42 <- ggboxplot(data, 
                                     x = "fraction", 
                                     y = "count", 
                                     facet.by = "external_transcript_name", 
                                     fill = "fraction") +
  geom_quasirandom(size = 0.5, shape = 21) +
  stat_pvalue_manual(stat.test, label = "ast", y.position = c(1, 1)) +
  labs(x = "Fraction", y = "Proportion of \nTotal Expression", title = "Cdc42") +
  theme_PK(font_size = 11) +
  theme(legend.position = "none") +
  panel_border() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_d3()

# p_mb_fraction_dtu_cdc42 <- counts(d_MB_FRACTION) %>%
#   as_tibble %>%
#   filter(str_detect(gene_id, "ENSMUSG00000006699")) %>%
#   select(ensembl_gene_id = gene_id,
#          ensembl_transcript_id = feature_id,
#          everything()) %>%
#   mutate(ensembl_gene_id = str_extract(ensembl_gene_id,
#                                        "[:alnum:]+(?=\\.[:digit:]+)"),
#          ensembl_transcript_id = str_extract(ensembl_transcript_id,
#                                              "[:alnum:]+(?=\\.[:digit:]+)")) %>%
#   left_join(anno_tx) %>%
#   # mutate(feature_id = paste0("Transcript ", seq(1, nrow(.), 1))) %>%
#   pivot_longer(-c(ensembl_gene_id, ensembl_transcript_id, external_transcript_name),
#                names_to = "sample_name",
#                values_to = "count") %>%
#   group_by(sample_name) %>%
#   filter(median(count) > 0) %>%
#   mutate(count = count / sum(count)) %>%
#   inner_join(colData(dds), copy = TRUE) %>%
#   ggplot(aes(x = fraction,
#              y = count,
#              fill = fraction)) +
#   geom_violin() +
#   geom_quasirandom(shape = 21,
#                    size = 3,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Fraction",
#        y = "Proportion of \nTotal Expression",
#        title = "Cdc42") +
#   theme(legend.position = "none") +
#   facet_wrap(vars(ensembl_transcript_id),
#              nrow = 1,
#              strip.position = "top") +
#   panel_border() +
#   scale_y_continuous(limits = c(0, 1))


# MB Fraction: DEU ----
# Can't find any results that are not just a case of enriched genes being detected

# Keep this code, but don't run it: Load from saved RDS

# files <-
#   list.files(
#     "/zfs/analysis/trap/active/testing_env/dexseq/output",
#     pattern = ".txt",
#     recursive =  T,
#     full.names = T
#   )
#
# all(file.exists(files)) # everything exists
# names(files) <-
#   str_extract(files, "(?<=output\\/)[0-9]+(?=.txt)")
#
# files <- files[names(files) %in% colData(dds)$sample_code]
# names(files[order(match(names(files), colData(dds)$sample_code))]) == colData(dds)$sample_code
# files <- files[order(match(names(files), colData(dds)$sample_code))]
#
# dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(
#   files,
#   sampleData = as.data.frame(colData(dds)),
#   design = ~ sample + exon + fraction:exon,
#   flattenedfile = list.files("/zfs/analysis/trap/active/testing_env/dexseq",
#                              pattern = "dexseq.gtf",
#                              full.names = T)
# )
#
# # filter dxd for cohort 1 IP only
# dxd <- dxd[,colData(dxd)$cohort == "C1"]
# dim(dxd)
# saveRDS(dxd, "R/objects/dxd_MB_FRACTION.rds")

# dxd <- readRDS("R/objects/dxd_MB_FRACTION.rds")
# dxd <- estimateSizeFactors(dxd)
# dxd <- estimateDispersions(dxd,
#                            BPPARAM = MulticoreParam())
#
# dxd <- DEXSeq::testForDEU(dxd,
#                   BPPARAM = MulticoreParam())
#
# dxd = DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar="compartment",
#                                        BPPARAM = MulticoreParam())
# saveRDS(dxd, "R/objects/dxd_MB_FRACTION.rds")

# dxd_MB_FRACTION <- readRDS("R/objects/dxd_MB_FRACTION.rds")
#
# dxr1_MB_FRACTION = DEXSeqResults( dxd_MB_FRACTION )
# # dxr1 %>% as.data.frame %>% View
#
# # number of genes considered for DEU
# dxr1_MB_FRACTION$groupID %>% unique %>% length
#
# # number of genes with signif DEU at alpha of 0.01
# dxr1_MB_FRACTION %>% as_tibble %>% filter(padj < 0.01) %>% pull(groupID) %>% unique %>% length
#
# dxr1_MB_FRACTION %>%
#   as_tibble %>%
#   filter(padj < 0.01) %>%
#   mutate(.before = groupID,
#          ensembl_gene_id = substr(groupID, 1, 18)) %>%
#   select(ensembl_gene_id,
#          groupID,
#          featureID,
#          exonBaseMean,
#          padj,
#          TOTAL,
#          IP,
#          lfc = log2fold_IP_TOTAL,
#          start = genomicData.start,
#          end = genomicData.end,
#          width = genomicData.width,
#          strand = genomicData.strand,
#          chr = genomicData.seqnames) %>%
#   mutate(translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES) %>%
#   # filter(translated) %>%
#   left_join(anno) %>%
#   arrange(padj) %>%
#   filter(exonBaseMean > 500) %>%
#   group_by(ensembl_gene_id) %>% View
#   filter(n() <= 2) %>% View
#
# # RCAN2
#
# plotDEXSeq(dxr1_MB_FRACTION,
#             "ENSMUSG00000019370.10",
#             fitExpToVar = "fraction",
#             legend = T,
#             displayTranscripts = F,
#             names = T,
#             color = pal_d3()(2),
#             # splicing = T,
#             transcriptDb = txdb,
#             sub = NULL)

# ----
# ----
# ----
# CREATE THE MB TRANSLATED LIST ----

MB_FRACTION_TRANSLATED_GENES <- c(MB_FRACTION_ENRICHED_GENES,
  MB_FRACTION_DTU_NOT_ENRICHED_GENES)

# ----
# ----
# ----
# AXON fraction subset common genes -----------------------------------------

# find common genes in all AXON IP and TOTAL samples
# Decided not to subset for only MB translated genes at this stage, because it is
# useful to know how many non-MB translated genes are enriched in Axons
genes_AXON <- intersect(
  # intersect(
  rownames(
    dds_C2_AXON
  ),
  rownames(
    dds_C3_AXON
  )
  # )
  # MB_FRACTION_TRANSLATED_GENES
)

dds_C2_AXON <-
  dds_C2_AXON[genes_AXON, ]
dds_C3_AXON <-
  dds_C3_AXON[genes_AXON, ]

# AXON fraction DESeq2 -----------------------------------------------------
# this segment generates a list of genes enriched in axon by TRAP

# Each cohort is analysed separately
# One sided p values are used to calculate a meta p value with stouffer's method
# fdr adjustment is applied and sumz pvalues less than 0.01 are considered significant

dds_C2_AXON@design <- ~ fraction
dds_C3_AXON@design <- ~ fraction

dds_C2_AXON <- DESeq(
  dds_C2_AXON,
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

dds_C3_AXON <- DESeq(
  dds_C3_AXON,
  sfType = "poscounts",
  fitType = "local",
  minReplicatesForReplace = Inf,
  parallel = TRUE
)

res_C2_AXON <- DESeq2::results(
  dds_C2_AXON,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C2_AXON <- lfcShrink(
  dds_C2_AXON,
  res = res_C2_AXON,
  contrast = c("fraction", "IP", "Input"),
  type = "ashr",
  parallel = T
)

res_C3_AXON <- DESeq2::results(
  dds_C3_AXON,
  cooksCutoff = Inf,
  filterFun = ihw,
  parallel = T
)
res_C3_AXON <- lfcShrink(
  dds_C3_AXON,
  res = res_C3_AXON,
  contrast = c("fraction", "IP", "Input"),
  type = "ashr",
  parallel = T
)

AXON_FRACTION_META <-
  as_tibble(res_C2_AXON, rownames = "ensembl_gene_id") %>%
  left_join(
    as_tibble(res_C3_AXON, rownames = "ensembl_gene_id"),
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
  mutate(
    pvalue_C3 = ifelse(
      sign(log2FoldChange_C2) == sign(log2FoldChange_C3),
      pvalue_C3,
      1 - pvalue_C3
    )
  ) %>%
  mutate(across(starts_with("pvalue"),
                ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
  mutate(
    sumz = sumz(c_across(pvalue_C2:pvalue_C3))$p,
    log2FoldChange = log2FoldChange_C2
  ) %>%
  ungroup() %>%
  mutate(
    sumz_adj = p.adjust(sumz, method = "fdr"),
    conflict = sign(log2FoldChange_C2) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
  mutate(score = -log10(sumz_adj) * log2FoldChange) %>%
  arrange(desc(score), desc(log2FoldChange)) %>%
  mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES)

# AXON FRACTION enriched lists ----
# get list of enriched genes that are not necessarily MB translated
AXON_FRACTION_ENRICHED_AGNOSTIC_GENES <- AXON_FRACTION_META %>%
  filter(sumz_adj < 0.01
         # !conflict # this only removes about 28 genes. Keep them in and just label by conflict status
  ) %>%
  rowwise() %>%
  filter(log2FoldChange > 0) %>%
  pull(ensembl_gene_id)

n_AXON_FRACTION_ENRICHED_AGNOSTIC_GENES <- AXON_FRACTION_ENRICHED_AGNOSTIC_GENES %>% length

AXON_FRACTION_ENRICHED_GENES <- AXON_FRACTION_META %>%
  filter(sumz_adj < 0.01 & mb_translated & log2FoldChange > 0) %>%
  pull(ensembl_gene_id)

n_AXON_FRACTION_ENRICHED_GENES <- length(AXON_FRACTION_ENRICHED_GENES)

# AXON FRACTION depleted ----

# get list of depleted genes
AXON_FRACTION_DEPLETED_GENES <- AXON_FRACTION_META %>%
  filter(sumz_adj < 0.01) %>%
  rowwise() %>%
  filter(log2FoldChange < 0) %>%
  pull(ensembl_gene_id)

n_AXON_FRACTION_DEPLETED_GENES <- AXON_FRACTION_DEPLETED_GENES %>% length

# FIGURE A: AXON fraction plot counts -------------------------------------------------

# MA-plot style
plot_data <- AXON_FRACTION_META %>%
  # filter(abs(log2FoldChange) < 10) %>%
  mutate(enrichment = ifelse(sumz_adj > 0.01, "Unchanged",
                             ifelse(log2FoldChange > 0 & mb_translated, "Enriched",
                                    ifelse(log2FoldChange > 0 & !mb_translated, "Filtered",
                                           "Depleted"))))

# highlight markers
highlight_markers <- plot_data %>%
  left_join(anno) %>%
  filter(external_gene_name %in% c(
    c("Slc6a3", "Th", "Ddc", "Slc18a2"),
    c("Gfap", "Gad2", "S100b", "Gad1", "Aldh1l1")
  ))

p_axon_ip_enrichment_count_plot <- plot_data %>%
  ggplot(aes(x = baseMean_C2,
             y = log2FoldChange,
             colour = enrichment)) +
  geom_point(size = 0.5,
             alpha = 0.5) +
  geom_point(data = highlight_markers,
             aes(fill = enrichment),
             size = 2,
             shape = 21,
             colour = "black") +
  geom_label_repel(data = highlight_markers,
                   aes(label = external_gene_name),
                   colour = "black",
                   arrow = arrow(length = unit(0.02, "npc")),
                   box.padding = 1,
                   # nudge_x = -1
  ) +
  scale_x_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)) +
  scale_y_continuous(limits = c(-6, 6),
                     breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
                     oob = scales::squish) +
  scale_color_manual(values = c("#1F77B4FF",
                                "#2CA02CFF",
                                "#8C564BFF",
                                "#C1C1C1")) +
  scale_fill_manual(values = c("#1F77B4FF",
                               "#2CA02CFF",
                               "#8C564BFF",
                               "#C1C1C1")) +
  labs(x = "Mean Counts",
       y = "Log2 Fold Change",
       colour = "Enrichment") +
  theme(legend.position = "top",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ) ,
         fill = FALSE)

# Keep this commented out code: This is the plot when a hard filter is placed, based
# on lfc
# # MA-plot style
# plot_data <- AXON_FRACTION_META_SIMPLE %>%
#   filter(abs(log2FoldChange) < 10 &
#            !conflict) %>%
#   mutate(enrichment = ifelse(sumz_adj > 0.01, "Unchanged",
#                              ifelse(log2FoldChange > 0 & mb_translated, "Enriched", "Depleted")))
# # obtain slope for filtering based on dopamine markers
# slope_data <- plot_data %>%
#   filter(external_gene_name %in% c("Slc18a2",
#                                    "Ddc"))
# model <- lm(log2FoldChange ~ log10(baseMean), data = slope_data)
#
# # predict lfc requirement based on dopamine markers
# plot_data <- plot_data %>%
#   rowwise() %>%
#   mutate(lfc_required = model$coefficients[1] + log10(baseMean)*model$coefficients[2],
#          enrichment_stringent = ifelse(enrichment == "Enriched",
#                                        ifelse(log2FoldChange >= lfc_required*0.999,
#                                               "Enriched", "Filtered"), enrichment),
#          enrichment_stringent = ifelse(sumz_adj < 0.01,
#                                        ifelse(log2FoldChange > 0,
#                                               ifelse(sapply(external_gene_name, get_ensembl_gene_id) %in% MB_FRACTION_TRANSLATED_GENES,
#                                                      "Enriched", "Filtered"),
#                                               "Depleted"),
#                                        "Unchanged"),
#          enrichment_stringent = factor(enrichment_stringent, levels = c("Depleted", "Enriched", "Filtered", "Unchanged")))
#
# # apply prediction to entire AXON table
# AXON_FRACTION_META <- AXON_FRACTION_META %>%
#   rowwise() %>%
#   mutate(lfc_required = model$coefficients[1] + log10(baseMean_C2)*model$coefficients[2],
#          enrichment_stringent = ifelse(sumz_adj < 0.01,
#                                        ifelse(log2FoldChange >= lfc_required*0.999,
#                                               "Enriched", "Filtered"), "Unchanged"),
#          enrichment_stringent = ifelse(sumz_adj < 0.01,
#                                        ifelse(log2FoldChange > 0,
#                                               ifelse(sapply(external_gene_name, get_ensembl_gene_id) %in% MB_FRACTION_TRANSLATED_GENES,
#                                                      "Enriched", "Filtered"),
#                                               "Depleted"),
#                                        "Unchanged"))
#
# AXON_FRACTION_META_SIMPLE <- AXON_FRACTION_META_SIMPLE %>%
#   rowwise() %>%
#   mutate(lfc_required = model$coefficients[1] + log10(baseMean)*model$coefficients[2],
#          enrichment_stringent = ifelse(sumz_adj < 0.01,
#                                        ifelse(log2FoldChange >= lfc_required*0.999,
#                                               "Enriched", "Filtered"), "Unchanged"),
#          enrichment_stringent = ifelse(sumz_adj < 0.01,
#                                        ifelse(log2FoldChange > 0,
#                                               ifelse(sapply(external_gene_name, get_ensembl_gene_id) %in% MB_FRACTION_TRANSLATED_GENES,
#                                                      "Enriched", "Filtered"),
#                                               "Depleted"),
#                                        "Unchanged"))
#
# # subset plot_data for line, so that the line doesn't extend below statistically enriched genes
# plot_line_data <- plot_data %>%
#   filter(baseMean > 2.5e1)
#
# # highlight markers
# highlight_markers <- plot_data %>%
#   left_join(anno) %>%
#   filter(external_gene_name %in% c(
#     c("Slc6a3", "Th", "Ddc", "Slc18a2")
#     # c("Gfap", "Gad2", "S100b", "Gad1", "Aldh1l1", "Chat")
#   )
#   )
#
# p_axon_ip_enrichment_count_plot <- plot_data %>%
#   ggplot(aes(x = baseMean,
#              y = log2FoldChange,
#              colour = enrichment_stringent)) +
#   geom_point(size = 0.5,
#              alpha = 0.5) +
#   geom_point(data = highlight_markers,
#              aes(fill = enrichment_stringent),
#              size = 2,
#              shape = 21,
#              colour = "black") +
#   geom_label_repel(data = highlight_markers,
#                    aes(label = external_gene_name),
#                    colour = "black",
#                    arrow = arrow(length = unit(0.02, "npc")),
#                    box.padding = 1,
#                    # nudge_x = -1
#   ) +
#   # geom_line(data = plot_line_data,
#   #           aes(x = baseMean,
#   #               y = lfc_required),
#   #           colour = "black",
#   #           linetype = "dotted") +
#   scale_x_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)) +
#   scale_y_continuous(limits = c(-6, 6),
#                      breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
#                      oob = scales::squish) +
#   scale_color_manual(values = c("#1F77B4FF",
#                                 "#2CA02CFF",
#                                 "#808080",
#                                 "#C1C1C1")) +
#   scale_fill_manual(values = c("#2CA02CFF",
#                                "#808080",
#                                "#C1C1C1")) +
#   labs(x = "Mean Counts",
#        y = "Log2 Fold Change",
#        colour = "Enrichment") +
#   theme(legend.position = "top",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14)) +
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1) ),
#          fill = FALSE)

# FIGURE B: AXON fraction panglao --------------------------------------------------

# internal
# score for each cell type
# adrenergic/nor-/sero- also enriched,
# but this is because of Slc6a3/Slc18a2/Th
# AXON_FRACTION_META %>%
#   left_join(anno_human) %>%
#   left_join(panglao,
#             by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#   mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
#          score = -log10(sumz_adj) * log2FoldChange) %>%
#   group_by(`cell type`) %>%
#   summarise(score = mean(score),
#             group_proportion = n()/group_size,
#             group_size = group_size,
#             prop_score = score * group_proportion) %>%
#   arrange(desc(prop_score)) %>%
#   distinct
#
# AXON_FRACTION_META_SIMPLE %>%
#   filter(enrichment_stringent != "Filtered") %>%
#   left_join(anno_human) %>%
#   left_join(panglao,
#             by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#   mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
#          score = -log10(sumz_adj) * log2FoldChange) %>%
#   group_by(`cell type`) %>%
#   summarise(score = mean(score),
#             group_proportion = n()/group_size,
#             group_size = group_size,
#             prop_score = score * group_proportion) %>%
#   arrange(desc(prop_score)) %>%
#   distinct

# table of dopaminergic marker enrichment
t_axon_ip_enrichment_dopamine_markers <- AXON_FRACTION_META %>%
  left_join(anno_human) %>%
  left_join(panglao,
            by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
  mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
         score = -log10(sumz_adj) * log2FoldChange) %>%
  group_by(`cell type`) %>%
  filter(`cell type` == "Dopaminergic neurons" &
           log2FoldChange > 0 &
           sumz_adj < 0.01) %>%
  select("Gene" = external_gene_name,
         "Log2 Fold Change" = log2FoldChange,
         "Adjusted P Value" = sumz_adj)

# AXON_FRACTION_META %>%
#   left_join(anno_human) %>%
#   left_join(panglao,
#             by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#   mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
#          score = -log10(sumz_adj) * log2FoldChange)

# # plot panglao scores before and after filtering on lfc
# p_axon_celltypes <- bind_rows({
#   AXON_FRACTION_META %>%
#     left_join(anno_human) %>%
#     left_join(panglao,
#               by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#     mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
#            score = -log10(sumz_adj) * log2FoldChange) %>%
#     group_by(`cell type`) %>%
#     summarise(score = mean(score),
#               group_proportion = n()/group_size,
#               group_size = group_size,
#               prop_score = score * group_proportion) %>%
#     arrange(desc(prop_score)) %>%
#     distinct() %>%
#     filter(`cell type` %in% c("Astrocytes",
#                               "Cholinergic neurons",
#                               "Dopaminergic neurons",
#                               "Endothelial cells",
#                               "GABAergic neurons",
#                               "Glutaminergic neurons",
#                               "Interneurons",
#                               "Microglia",
#                               "Oligodendrocytes",
#                               "Oligodendrocyte progenitor cells"
#     )) %>%
#     mutate(`cell type` = factor(`cell type`)) %>%
#     mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
#     mutate(filtering = "Unfiltered")
# },
# {
#   AXON_FRACTION_META %>%
#     filter(enrichment_stringent != "Filtered") %>%
#     left_join(anno_human) %>%
#     left_join(panglao,
#               by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#     mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
#            score = -log10(sumz_adj) * log2FoldChange) %>%
#     group_by(`cell type`) %>%
#     summarise(score = mean(score),
#               group_proportion = n()/group_size,
#               group_size = group_size,
#               prop_score = score * group_proportion) %>%
#     arrange(desc(prop_score)) %>%
#     distinct() %>%
#     filter(`cell type` %in% c("Astrocytes",
#                               "Cholinergic neurons",
#                               "Dopaminergic neurons",
#                               "Endothelial cells",
#                               "GABAergic neurons",
#                               "Glutaminergic neurons",
#                               "Interneurons",
#                               "Microglia",
#                               "Oligodendrocytes",
#                               "Oligodendrocyte progenitor cells"
#     )) %>%
#     mutate(`cell type` = factor(`cell type`)) %>%
#     mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
#     mutate(filtering = "Filtered")
# }) %>%
#   mutate(filtering = factor(filtering, levels = c("Unfiltered", "Filtered"))) %>%
#   ggplot(aes(x = prop_score,
#              y = `cell type`,
#              fill = `cell type`)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Enrichment Score",
#        y = "Cell Type") +
#   theme(legend.position = "none",
#         axis.text = element_text(size = 12),
#         axis.title.y = element_blank()) +
#   facet_wrap(vars(filtering))


# plot cell type enrichment
p_axon_celltypes <- AXON_FRACTION_META %>%
  left_join(anno_human) %>%
  left_join(panglao,
            by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
  mutate(sumz_adj = ifelse(sumz_adj == 0, min(sumz_adj[sumz_adj > 0]), sumz_adj),
         score = -log10(sumz_adj) * log2FoldChange_C2) %>%
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
                            "Oligodendrocytes",
                            "Oligodendrocyte progenitor cells"
  )) %>%
  mutate(`cell type` = factor(`cell type`)) %>%
  mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
  ggplot(aes(x = prop_score,
             y = `cell type`,
             fill = `cell type`)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  labs(x = "Enrichment Score",
       y = "Cell Type") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.y = element_blank())




#
#
# FIGURE C: AXON-MB fraction venn ----------
d <- tibble(value = unique(c(MB_FRACTION_TRANSLATED_GENES,
                             AXON_FRACTION_ENRICHED_AGNOSTIC_GENES))) %>%
  mutate(`Cell Body` = value %in% MB_FRACTION_TRANSLATED_GENES,
         `Axon` = value %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)

p_axon_ip_enrichment_venn <- ggvenn(d,
                                    fill_color = pal_d3()(4),
                                    fill_alpha = 0.55,
                                    text_size = 3)

# FIGURE D: AXON Cre Counts ---------------------------------------------------------
files <-
  list.files(
    "/zfs/analysis/trap/active/testing_env/cre/kallisto",
    pattern = "abundance.tsv",
    recursive =  T,
    full.names = T
  )

names(files) <-
  str_extract(files, "(?<=kallisto\\/)[0-9|\\_]+(?=\\/abundance.tsv)")

p_axon_cre_counts <- lapply(files, function(x){
  read_tsv(x) %>%
    mutate(tpm = as.double(tpm))
           } 
  ) %>%
  bind_rows %>%
  mutate(sample_code = names(files)) %>%
  left_join(colData(dds), copy = T) %>%
  filter(!is.na(compartment) &
           cohort == "C2") %>%
  ggplot(aes(x = fraction,
             y = log2(est_counts + 1),
             fill = fraction)) +
  geom_quasirandom(colour = "black",
                   shape = 21,
                   size = 3) +
  scale_fill_d3() +
  facet_wrap(vars(compartment)) +
  labs(x = "Fraction",
       y = "Log2 Expression Count") +
  theme(legend.position = "none")


#
#
#
#

# TABLE 5: AXON fraction Cibersortx -----
t_cibersort_sc <- read_delim("R/objects/CIBERSORTx_SC.txt",
                             delim = "\t") %>%
  inner_join(colData(dds),
             by = c("Mixture" = "sample_name0"),
             copy = T) %>%
  select(Neuron:NSC, region, fraction) %>%
  group_by(region, fraction) %>%
  summarise(across(Neuron:NSC, mean)) %>%
  rename("Region" = region,
         "Fraction" = fraction,
         "Astrocyte" = Astro,
         "Oligodendrocyte" = Oligo) %>%
  select(-c(Vascular,
            `Ependy-Sec`,
            # Microglia,
            OPC)) %>%
  mutate(across(Neuron:NSC, function(x){x*100}),
         Fraction = ifelse(Fraction == "IP", "TRAP", "TOTAL")) %>%
  select(Region,
         Fraction,
         sort(colnames(.)))

t_cibersort_chen <- read_delim("R/objects/CIBERSORTx_ChenTRAP.txt",
                               delim = "\t") %>%
  inner_join(colData(dds),
             by = c("Mixture" = "sample_name0"),
             copy = T) %>%
  select(Endothelial:Astrocyte, region, fraction) %>%
  group_by(region, fraction) %>%
  summarise(across(Endothelial:Astrocyte, mean)) %>%
  rename("Region" = region,
         "Fraction" = fraction) %>%
  select(-c(Endothelial,
            MO)) %>%
  mutate(across(Microglia:Astrocyte, function(x){x*100}),
         Fraction = ifelse(Fraction == "IP", "TRAP", "TOTAL")) %>%
  select(Region,
         Fraction,
         sort(colnames(.)))

t_cibersort_mixedTRAP <- read_delim("R/objects/CIBERSORTx_mixedTRAP.txt",
                                    delim = "\t") %>%
  inner_join(colData(dds),
             by = c("Mixture" = "sample_name0"),
             copy = T) %>%
  select(astrocytestriatal2year:astrocytestriatal10week, region, fraction) %>%
  group_by(region, fraction) %>%
  summarise(across(astrocytestriatal2year:astrocytestriatal10week, mean)) %>%
  rowwise() %>%
  mutate(.before = astrocytestriatal2year,
         Astrocyte = sum(astrocytestriatal10week,
                         astrocytestriatal2year)) %>%
  mutate(Oligodendrocyte = sum(oligodendrocytesControl,
                               oligodendrocytesfed)) %>%
  rename("Region" = region,
         "Fraction" = fraction,
         "D1 MSN" = D1,
         "D2 MSN" = D2,
         "Midbrain Dopamine Neuron" = danbrichta) %>%
  select(-c(astrocytestriatal10week,
            astrocytestriatal2year,
            oligodendrocytesControl,
            oligodendrocytesfed)) %>%
  mutate(across(Astrocyte:Oligodendrocyte, function(x){x*100}),
         Fraction = ifelse(Fraction == "IP", "TRAP", "TOTAL")) %>%
  select(Region,
         Fraction,
         sort(colnames(.)))

#
#
#
#
#
#
#
#
#
#
# ----
# ----
# ----
# # AXON AGE - LIBRARY COMPOSITION CHECK ----
#
#
# C2_AXON_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
#       select(age, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
#       select("Th" = count)
#   }, {
#     #   plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
#     #     select("Ddc" = count)
#     # }, {
#     plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-age,
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C2")
#
#
# C3_AXON_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
#       select(age, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
#       select("Th" = count)
#   }, {
#     #   plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
#     #     select("Ddc" = count)
#     # }, {
#     plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C3_AXON_IP, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-age,
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C3")
#
#
#
# p_axon_ip_age_library_composition <- bind_rows(
#   C2_AXON_IP_MARKERS,
#   C3_AXON_IP_MARKERS) %>%
#   mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Slc18a2"),
#                          "Dopaminergic", "Other")) %>%
#   ggplot(aes(x = age,
#              y = log2(count + 1),
#              colour = marker,
#              group = gene)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(vars(cohort), nrow = 1) +
#   panel_border() +
#   scale_color_d3() +
#   labs(x = "Cohort and Age",
#        y = "Log2 Expression Count",
#        colour = "Cell Type Marker") +
#   theme(legend.position = "top")
#
#
#
#
#
#
# C2_AXON_TOTAL_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
#       select(age, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
#       select("Th" = count)
#   }, {
#   #   plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
#   #     select("Ddc" = count)
#   # }, {
#     plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C2_AXON_TOTAL, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-age,
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C2")
#
#
# C3_AXON_TOTAL_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Slc6a3"), "age", returnData = T) %>%
#       select(age, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Th"), "age", returnData = T) %>%
#       select("Th" = count)
#   }, {
#   #   plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Ddc"), "age", returnData = T) %>%
#   #     select("Ddc" = count)
#   # }, {
#     plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Slc18a2"), "age", returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Gfap"), "age", returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("Gad2"), "age", returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C3_AXON_TOTAL, get_ensembl_gene_id("S100b"), "age", returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-age,
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C3")
#
#
# p_axon_total_age_library_composition <- bind_rows(C2_AXON_TOTAL_MARKERS,
#                                                   C3_AXON_TOTAL_MARKERS) %>%
#   mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
#                          "Dopaminergic", "Other")) %>%
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
#
#
#
#
#
#
#
# # AXON AGE - ALTERNATIVE COMPOSITION CHECK ----
# # This code is commented because I don't know which genes are reliably axonal in striatum:
# # If I look at the proportion of library that is enriched, that will include
# # other cell types, so the proportion doesnt change so much
# # The first composition check (plotting DA markers) is informative.
# # LIBRARY_PROPORTIONS <- bind_rows(
# #   {
# #     as_tibble(counts(dds_C2_AXON_IP), rownames = "ensembl_gene_id") %>%
# #       mutate(cohort = "C2")
# #   }, {
# #     as_tibble(counts(dds_C3_AXON_IP), rownames = "ensembl_gene_id") %>%
# #       mutate(cohort = "C3")
# #   }
# # ) %>%
# #   pivot_longer(-c(ensembl_gene_id, cohort),
# #                names_to = "sample_name",
# #                values_to = "count") %>%
# #   mutate(status = ifelse(ensembl_gene_id %in% AXON_FRACTION_ENRICHED_GENES,
# #                          "Enriched",
# #                          ifelse(ensembl_gene_id %in% AXON_FRACTION_ENRICHED_GENES,
# #                                 "Depleted", "Unchanged")),
# #          status = factor(status, levels = c("Depleted", "Unchanged", "Enriched"))) %>%
# #   mutate(age = ifelse(str_detect(sample_name, "YOUNG"), "YOUNG", "OLD"),
# #          age = factor(age, levels = c("YOUNG", "OLD"))) %>%
# #   group_by(cohort, age, status) %>%
# #   summarise(count = sum(count, na.rm = T)) %>%
# #   mutate(prop = count/sum(count))
# #
# # # DON'T KNOW WHICH TO USE: Don't know which genes are definitely axon specific in striatum.
# # # AXON_FRACTION_ENRICHED_GENES
# # # AXON_FRACTION_ENRICHED_ABACONF_GENES
# # # AXON_FRACTION_ENRICHED_AGNOSTIC_GENES
# # # MB_FRACTION_ENRICHED_GENES
# #
# # p_age_mb_library_proportions <- LIBRARY_PROPORTIONS %>%
# #   ggplot(aes(x = age,
# #              y = prop,
# #              fill = status)) +
# #   geom_col() +
# #   facet_wrap(vars(cohort), nrow = 1)  +
# #   theme(legend.position = "top") +
# #   geom_text(aes(y = prop, label = ifelse(prop > 0.1, scales::percent(prop, accuracy = 1), "")),
# #             position = position_stack(vjust = 0.5),
# #             show.legend = FALSE,
# #             colour = "white",
# #             fontface = 2) +
# #   scale_fill_manual(values = pal_d3()(4)[c(4, 1, 3)]) +
# #   scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
# #                      labels = scales::percent) +
# #   panel_border() +
# #   labs(x = "Age and Cohort",
# #        y = "Percentage of library",
# #        fill = "Enrichment")
#
# # AXON AGE TOTAL DESeq2 -----------------------------------------------------
#
# # Using C2 and C3 for total to maximise number of axonal samples.
# # Will check whether combining worsens detection power like TRAP samples
#
# # Filter for genes that are reliably detected in both cohorts
# dds_C2_AXON_TOTAL_AGE_FILTER <- filter_genes(dds_C2_AXON_TOTAL,
#                                     grouping = c("age", "gene_id"))
# dds_C3_AXON_TOTAL_AGE_FILTER <- filter_genes(dds_C3_AXON_TOTAL,
#                                     grouping = c("age", "gene_id"))
#
# genes_AXON_TOTAL_AGE <- intersect(dds_C2_AXON_TOTAL_AGE_FILTER,
#                          dds_C3_AXON_TOTAL_AGE_FILTER)
#
# dds_C2_AXON_TOTAL_AGE <- dds_C2_AXON_TOTAL[rownames(dds_C2_AXON_TOTAL) %in% genes_AXON_TOTAL_AGE,]
# dds_C3_AXON_TOTAL_AGE <- dds_C3_AXON_TOTAL[rownames(dds_C3_AXON_TOTAL) %in% genes_AXON_TOTAL_AGE,]
#
# dds_C2_AXON_TOTAL_AGE@design <- ~ age
# dds_C3_AXON_TOTAL_AGE@design <- ~ age
#
# colData(dds_C2_AXON_TOTAL_AGE) <- droplevels(colData(dds_C2_AXON_TOTAL_AGE))
# colData(dds_C3_AXON_TOTAL_AGE) <- droplevels(colData(dds_C3_AXON_TOTAL_AGE))
#
# dds_C2_AXON_TOTAL_AGE <- DESeq(
#   dds_C2_AXON_TOTAL_AGE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# dds_C3_AXON_TOTAL_AGE <- DESeq(
#   dds_C3_AXON_TOTAL_AGE,
#   fitType = "local",
#   sfType = "poscounts",
#   minReplicatesForReplace = Inf,
#   parallel = TRUE,
# )
#
# # dds_C3_AXON_TOTAL_AGE <- dds_C3_AXON_TOTAL_AGE[rowData(dds_C3_AXON_TOTAL_AGE)$baseMean > 500 ,]
#
# # results
#
# res_C2_AXON_TOTAL_AGE <- DESeq2::results(
#   dds_C2_AXON_TOTAL_AGE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C2_AXON_TOTAL_AGE <- lfcShrink(
#   dds_C2_AXON_TOTAL_AGE,
#   res = res_C2_AXON_TOTAL_AGE,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# res_C3_AXON_TOTAL_AGE <- DESeq2::results(
#   dds_C3_AXON_TOTAL_AGE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C3_AXON_TOTAL_AGE <- lfcShrink(
#   dds_C3_AXON_TOTAL_AGE,
#   res = res_C3_AXON_TOTAL_AGE,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# AXON_AGE_TOTAL_META <-
#   as_tibble(res_C2_AXON_TOTAL_AGE, rownames = "ensembl_gene_id") %>%
#   left_join(
#     as_tibble(res_C3_AXON_TOTAL_AGE, rownames = "ensembl_gene_id"),
#     by = "ensembl_gene_id",
#     suffix = c("_C2", "_C3")
#   ) %>%
#   left_join(anno) %>%
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
#   # filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
#   drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
#   rowwise() %>%
#   mutate(
#     pvalue_C3 = ifelse(
#       sign(log2FoldChange_C3) == sign(log2FoldChange_C2),
#       pvalue_C3,
#       1 - pvalue_C3)
#   ) %>%
#   mutate(across(starts_with("pvalue"),
#                 ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
#   mutate(
#     # sumlog = sumlog(c_across(pvalue_C2:pvalue_C3))$p,
#     # calculate sumlog
#     sumz = sumz(c_across(pvalue_C2:pvalue_C3))$p,
#     # calculate stouffers
#     # log2FoldChange = mean(c(
#     #   log2FoldChange_C2,
#     #   log2FoldChange_C3
#     # ))
#     log2FoldChange = log2FoldChange_C2
#   ) %>%
#   ungroup() %>%
#   mutate(
#     # sumlog_adj = p.adjust(sumlog, method = "fdr"),
#     # correction for multiple comparisons
#     sumz_adj = p.adjust(sumz, method = "fdr"),
#     conflict = sign(log2FoldChange_C2) != sign(log2FoldChange_C3)) %>%
#   mutate(outcome = ifelse(sumz_adj < 0.1,
#                           ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
#                           "Unchanged"))
#
# # AXON AGE TOTAL OUTCOMES ----
# AXON_AGE_TOTAL_OUTCOMES <- AXON_AGE_TOTAL_META %>%
#   select(ensembl_gene_id,
#          outcome)
#
#
#
# # AXON AGE TRAP DESeq2 -----------------------------------------------------
#
# # Just using Cohort 2: The dispersion in cohort 3 is just too large.
# # Also, only 7000 genes are measurable in both
# # versus 17000 in C2 alone
# # Combining cohort 2 and 3 by meta analysis results in fewer detections,
# # and the results of cohort 2 actually looking realistic (DA markers downregulated in old)
#
# dds_C2_AXON_IP_AGE_FILTER <- filter_genes(dds_C2_AXON_IP,
#                                     grouping = c("age", "gene_id"))
#
# dds_C2_AXON_IP_AGE <- dds_C2_AXON_IP[dds_C2_AXON_IP_AGE_FILTER,]
#
# dds_C2_AXON_IP_AGE@design <- ~ collection + age
#
# colData(dds_C2_AXON_IP_AGE) <- droplevels(colData(dds_C2_AXON_IP_AGE))
#
# dds_C2_AXON_IP_AGE <- DESeq(
#   dds_C2_AXON_IP_AGE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_C2_AXON_IP_AGE <- DESeq2::results(
#   dds_C2_AXON_IP_AGE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C2_AXON_IP_AGE <- lfcShrink(
#   dds_C2_AXON_IP_AGE,
#   res = res_C2_AXON_IP_AGE,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# # AXON AGE: Use depleted genes to add to AXON_FRACTION_ENRICHED_GENES list ----
# AXON_FRACTION_ENRICHED_GENES <- res_C2_AXON_IP_AGE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   select(-c(lfcSE, pvalue)) %>%
#   filter(padj < 0.01) %>%
#   arrange(log2FoldChange) %>%
#   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
#          axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
#   filter(axon_translated &
#            !mb_translated &
#            log2FoldChange < 0) %>%
#   pull(ensembl_gene_id)
#

# ----
# ----
# ----
# AXON FRACTION: Update META with PubMed, ABA and MB enrichment data ----

AXON_FRACTION_META <- AXON_FRACTION_META %>%
  select(ensembl_gene_id,
         external_gene_name,
         description,
         baseMean = baseMean_C2,
         log2FoldChange,
         sumz_adj,
         conflict,
         score,
         mb_translated) %>%
  left_join(aba_expression_axon) %>%
  left_join(publications_all_genes_axon) %>%
  left_join(publications_all_genes)

# AXON FRACTION enriched and MB enriched ----
# get list of enriched genes



AXON_FRACTION_ENRICHED_GENES <- c(AXON_FRACTION_ENRICHED_GENES,
                                  {AXON_FRACTION_META %>%
                                      filter(sumz_adj < 0.01 & mb_translated) %>%
                                      rowwise() %>%
                                      filter(log2FoldChange > 0) %>%
                                      pull(ensembl_gene_id)})

# number of genes assayed
n_AXON_FRACTION <- AXON_FRACTION_META %>% nrow


# AXON FRACTION enriched and MB enriched and ABA confirmed ----
AXON_FRACTION_ENRICHED_ABACONF_GENES <- AXON_FRACTION_META %>%
  filter(mb_translated & sumz_adj < 0.01 & log2FoldChange > 0 & expression_energy < 0.5) %>%
  pull(ensembl_gene_id)

n_AXON_ENRICHED_ABACONF <- length(AXON_FRACTION_ENRICHED_ABACONF_GENES)

# AXON FRACTION enriched agnostic and ABA confirmed ----
# AXON_FRACTION_ENRICHED_AGNOSTIC_ABACONF_GENES <- AXON_FRACTION_META %>%
#   filter(expression_energy < 0.5) %>%
#   pull(ensembl_gene_id)



# ----
# ----
# ----
# TABLE 6: AXON fraction top 30 with ABA check ----
t_axon_fraction_top_enriched <- AXON_FRACTION_META %>%
  filter(!conflict) %>%
  filter(expression_energy < 0.5 &
           mb_translated) %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(description = str_to_sentence(str_extract(description, ".+(?=\\[)"))) %>%
  filter(!str_detect(description, "Predicted|Riken|pseudogene")) %>%
  slice_max(order_by = log2FoldChange,
            n = 30) %>%
  select(external_gene_name,
         description,
         log2FoldChange,
         expression_energy)

# TABLE 6: AXON fraction table summary -----------------------------------------------
t_axon_ip_enrichment_summary <- tibble(Metric = c("Number of Genes", "Percentage of Library"),
                                       Enriched = c(length(AXON_FRACTION_ENRICHED_AGNOSTIC_GENES),
                                                    paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[AXON_FRACTION_ENRICHED_AGNOSTIC_GENES]*100, na.rm = T), 3), "%")),
                                       Depleted = c(length(AXON_FRACTION_DEPLETED_GENES),
                                                    paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[AXON_FRACTION_DEPLETED_GENES]*100, na.rm = T), 3), "%")),
                                       Unchanged = c(nrow(AXON_FRACTION_META)-length(AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)-length(AXON_FRACTION_DEPLETED_GENES),
                                                     paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[!rownames(dds_C2_AXON_IP) %in% c(AXON_FRACTION_ENRICHED_AGNOSTIC_GENES, AXON_FRACTION_DEPLETED_GENES)] * 100, na.rm = T), 3), "%")))

AXON_FRACTION_ENRICHED_FILTERED_GENES <- AXON_FRACTION_ENRICHED_AGNOSTIC_GENES[!AXON_FRACTION_ENRICHED_AGNOSTIC_GENES %in% c(AXON_FRACTION_ENRICHED_GENES)]

t_axon_ip_enrichment_summary_filtered <- tibble(Metric = c("Number of Genes", "Percentage of Library"),
                                                "Enriched: Retained" = c(length(AXON_FRACTION_ENRICHED_GENES),
                                                                         paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[AXON_FRACTION_ENRICHED_GENES]*100, na.rm = T), 3), "%")),
                                                "Enriched: Filtered" = c(length(AXON_FRACTION_ENRICHED_FILTERED_GENES),
                                                                         paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[AXON_FRACTION_ENRICHED_FILTERED_GENES]*100, na.rm = T), 3), "%")),
                                                Depleted = c(length(AXON_FRACTION_DEPLETED_GENES),
                                                             paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[AXON_FRACTION_DEPLETED_GENES]*100, na.rm = T), 3), "%")),
                                                Unchanged = c(nrow(AXON_FRACTION_META)-length(AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)-length(AXON_FRACTION_DEPLETED_GENES),
                                                              paste(signif(sum((rowSums(counts(dds_C2_AXON_IP))/sum(colSums(counts(dds_C2_AXON_IP))))[!rownames(dds_C2_AXON_IP) %in% c(AXON_FRACTION_ENRICHED_AGNOSTIC_GENES, AXON_FRACTION_DEPLETED_GENES)] * 100, na.rm = T), 3), "%")))
#
#


# ----


# ----
# ----


# AXON Age: This is here so that age-depleted genes can be added to the
# AXON FRACTION ENRICHED GENES list
# GWAS prioritisation - Hook 2018 --------------------------------------------------------------------

# # start with Hook, 2018 validation
# # load Hook 2018 GWAS scores
# hook18_scores <- read_delim("R/objects/park17_gwas_chang_scores.txt",
#                             delim = "\t") %>%
#   group_by(snp) %>%
#   arrange(desc(score.pLI)) %>%
#   filter(!is.na(MouseSymbol)) %>%
#   mutate(rank_park = row_number())
#
# # plot agreement between Hook and TRAP
# # Discarding pLI information, because that is not part of the single cell/enrichment method
# # and unfairly distinguishes equally ranked genes from one another in the Hook dataset
#
# gwas_hook18_trap <- tibble(hsapiens_homolog_associated_gene_name = hook18_scores$HumanSymbol,
#                            snp = hook18_scores$snp,
#                            score_park = hook18_scores$score,
#                            score_park.pli = hook18_scores$score.pLI) %>%
#   left_join(anno_human) %>%
#   left_join(
#     {
#       MB_FRACTION_META %>%
#         select(ensembl_gene_id,
#                log2FoldChange,
#                sumz_adj,
#                score_trap = score,
#                conflict)
#     }
#   ) %>%
#   filter(!is.na(ensembl_gene_id) &
#            !is.na(log2FoldChange)) %>%
#   group_by(snp) %>%
#   arrange(snp, desc(score_park)) %>%
#   select(snp,
#          external_gene_name,
#          ensembl_gene_id,
#          score_park,
#          score_trap) %>%
#   mutate(rank_park = rank(-score_park,
#                           ties.method = "min"),
#          rank_trap = rank(-score_trap))
#
# gwas_hook18_trap_top_park <- gwas_hook18_trap %>%
#   filter(rank_park == 1) %>%
#   slice_min(order_by = rank_trap,
#             n = 1)
#
# # plot agreement between park and TRAP
# p_hook18_trap_comparison <- gwas_hook18_trap_top_park %>%
#   ggplot(aes(x = rank_trap,
#              fill = factor(rank_trap,
#                            levels = seq(max(rank_trap), 1, -1)))) +
#   geom_bar(colour = "black") +
#   scale_fill_brewer() +
#   geom_text(aes(label = ..count..,
#                 x = rank_trap,
#                 y = ..count..),
#             stat = "count",
#             vjust = -0.5) +
#   scale_x_continuous(breaks = seq(1, 15, 1)) +
#   theme(panel.grid.major.x = element_blank(),
#         legend.position = "none") +
#   labs(x = "Ranking by Somal Enrichment",
#        y = "Number of Genes") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
#
# # Get top TRAP gene for disagreeing SNPs
# gwas_hook18_trap_top_trap_filtered <- gwas_hook18_trap %>%
#   slice_min(order_by = rank_trap,
#             n = 1) %>%
#   select(snp, external_gene_name, rank_park) %>%
#   filter(rank_park != 1)
#
# # Final list for presentation
# t_gwas_hook18_trap <- gwas_hook18_trap_top_park %>%
#   left_join(gwas_hook18_trap_top_trap_filtered,
#             by  = c("snp"))  %>%
#   replace_na(list(external_gene_name.y = "",
#                   rank_park.y = "")) %>%
#   select(SNP = snp,
#          "Top Gene in Park" = external_gene_name.x,
#          # "Rank (Park)" = rank_park.x,
#          "Rank in TRAP" = rank_trap,
#          "Top Gene in TRAP" = external_gene_name.y,
#          "Rank in Park" = rank_park.y)
#
# t_gwas_hook18_trap %>% View
#
#
#
# # get coordinates of Park 2018 snps - Not doing this code anymore: The above code is final.
# hook18_snps <- getBM(attributes = c("refsnp_id",
#                                     "chr_name",
#                                     "chrom_start",
#                                     "associated_gene",
#                                     "ensembl_gene_name"),
#                      filters = c("snp_filter"),
#                      values = unique(hook18_scores$snp),
#                      mart = ensembl_variation,
#                      uniqueRows = TRUE)
#
#
#
# filter(!str_detect(chr_name, "CHR")) %>%
#   select("variant_id" = refsnp_id,
#          "chromosome_name" = chr_name,
#          "chromosome_position" = chrom_start)

# # get variants for Park, 2018 study
# pd_gwas_variants <- park18_snps
#
# # get +-1Mb position for each variant
# pd_gwas_variants_ranges <- tibble(
#   variant = pd_gwas_variants$variant_id,
#   chr = pd_gwas_variants$chromosome_name,
#   start = pd_gwas_variants$chromosome_position-1e6,
#   end = pd_gwas_variants$chromosome_position+1e6
# )
#
# # get all human genes that fall within those boundaries
# # retrieve genes within ranges including mouse homologs
# pd_gwas_variants_genes_range <- apply(pd_gwas_variants_ranges, 1, function(x){
#   getBM(attributes = c("ensembl_gene_id",
#                        "external_gene_name",
#                        "mmusculus_homolog_associated_gene_name"),
#         filters = c("chromosome_name", "start", "end"),
#         values = list(x[2],
#                       x[3],
#                       x[4]),
#         mart = ensembl_human,
#         uniqueRows = TRUE) %>%
#     mutate(variant = x[1])
# }) %>%
#   bind_rows()

# join META FRACTION to Park 2018 scores

# # Trying to see if age and genotype improve the correlation - doesn't look like it
# # but then the Park dataset is not the gold standard, of course!
# could propose a revised table at the end incorporating age and genotype and ds vs vs info
# park18_scores %>%
#   left_join(MB_FRACTION_META_SIMPLE,
#             by = c("MouseSymbol" = "external_gene_name")) %>%
#   arrange(snp, desc(log2FoldChange)) %>%
#   mutate(rank_trap = row_number()) %>%
#   select(MouseSymbol,
#          snp,
#          locus,
#          score = score.x,
#          score.pLI,
#          rank_park,
#          score_enrichment = score.y, rank_trap) %>%
#   left_join(anno,
#             by = c("MouseSymbol" = "external_gene_name")) %>%
#   filter(!is.na(MouseSymbol)) %>%
#   mutate(ageing = ensembl_gene_id %in% MB_AGE_DE_GENES,
#          ovx = ensembl_gene_id %in% ) %>% View

# park18_scores_vs_trap_cor <- park18_scores %>%
#   left_join(MB_FRACTION_META_SIMPLE,
#             by = c("MouseSymbol" = "external_gene_name")) %>%
#   arrange(snp, desc(log2FoldChange)) %>%
#   mutate(rank_trap = row_number()) %>%
#   group_by(snp) %>%
#   summarise(cor = cor(rank_park, rank_trap, method = "spearman"))
#
# park18_scores_vs_trap$park_top <- park18_scores %>%
#   group_by(snp) %>%
#   slice_max(order_by = score.pLI,
#             n = 1) %>%
#   pull(MouseSymbol)
#
# park18_scores %>%
#   group_by(snp) %>%
#   slice_max(order_by = score.pLI,
#             n = 3)
#
# park18_scores_vs_trap$trap_top <- park18_scores %>%
#   left_join(MB_FRACTION_META_SIMPLE,
#             by = c("MouseSymbol" = "external_gene_name")) %>%
#   group_by(snp) %>%
#   slice_max(order_by = log2FoldChange,
#             n = 3) %>%
#   select(snp, MouseSymbol, rank_trap)
#
#   deframe()
#   pull(MouseSymbol)
#
# t_park18_scores_vs_trap <- park18_scores_vs_trap %>%
#   rename("Variant" = snp,
#          "Rank Correlation" = cor,
#          "Park Top Candidate" = park_top,
#          "TRAP Top Candidate" = trap_top)


#
#

# ----
# ----
# ----
# GWAS - Nalls 2019 ----

# find all PD studies
# pd_gwas <- get_studies(efo_id = 'EFO_0002508')

# get variants
# pd_gwas_variants <- get_variants(study_id = "GCST009325") # Nalls, 2019
# pd_gwas_variants <- get_variants(study_id = "GCST004902") # Chang, 2017
# pd_gwas_variants <- pd_gwas_variants@variants
# pd_gwas_variants %>% View

# load Nalls summary statistics
nalls <- read_excel("R/objects/Table S2. Detailed summary statistics on all nominated risk variants, known and novel_.xlsx")
nalls <- nalls %>%
  filter(str_detect(SNP, "^rs") &
           !str_detect(CHR, "Signif|Hardy|Average"))

# # GET EQTL TISSUES
# ext <- "/eqtl/tissue/homo_sapiens?"
#
# r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
#
# stop_for_status(r)
#
# sort(names(unlist(fromJSON(toJSON(content(r))))))


# # GET VARIANTS AND GENES WITHIN 500kb WINDOW WITH LD +0.6
# get_variant_genes_and_consequences <- function(lead_variant){
#
#   # SET SERVER
#   server <- "https://rest.ensembl.org"
#
#   # GET LD VARIANTS +0.6 LD
#    ext <- paste0("/ld/human/",
#                   lead_variant,
#                   "/1000GENOMES:phase_3:GBR?r2=0.6")
#     r <-
#       GET(paste(server, ext, sep = ""),
#           content_type("application/json"))
#     stop_for_status(r)
#
#     linked_variants <- c(lead_variant,
#       unlist(fromJSON(toJSON(content(
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
#       variants_consequences <- fromJSON(toJSON(content(r)))
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
#
#
# # GET EQTL GENES PER VARIANT
#
# # get +-250kb position for each variant
# lead_variant_range <- tibble(
#   variant = pd_gwas_variants$variant_id,
#   chr = pd_gwas_variants$chromosome_name,
#   start = pd_gwas_variants$chromosome_position-2.5e5,
#   end = pd_gwas_variants$chromosome_position+2.5e5
# )
#
# ext <- paste0("/overlap/region/human/",
#               paste0(lead_variant_range$chr,
#                      ":",
#                      lead_variant_range$start,
#                      "-",
#                      lead_variant_range$end),
#               "?feature=gene")
#
# eqtl_genes <- lapply(ext, function(x){
#   r <- GET(paste(server, x, sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   fromJSON(toJSON(content(r))) %>%
#     pull(id) %>%
#     unlist()
# })
#
# names(eqtl_genes) <- lead_variant_range$variant

# get Nalls genes and homologs ----
# load nalls genes within 1 MB r^2 of loci
nalls_genes <- read_excel("R/objects/nalls_genes.xlsx") %>%
  select("external_gene_name" = "Nearest gene")

# add mouse homologs

homologs_ncbi <- nalls_genes %>%
  left_join({
    homologene(nalls_genes$external_gene_name, inTax = 9606, outTax = 10090) %>%
      select(external_gene_name = "9606",
             mmusculus_homolog_associated_gene_name = "10090")
  })

# homologs_ensembl <- getBM(
#   attributes = c(
#     "ensembl_gene_id",
#     "external_gene_name",
#     "mmusculus_homolog_associated_gene_name"
#   ),
#   filters = "external_gene_name",
#   values = nalls_genes$external_gene_name,
#   mart = ensembl_human,
#   uniqueRows = TRUE
# )

homologs_ensembl <- readRDS("R/objects/homologs_ensembl.rds")

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

# these humans genes are not in the database/data
nalls_genes$external_gene_name[!nalls_genes$external_gene_name %in% homologs$external_gene_name]

# add TRAP to nalls -----
# add MB and AXON FRACTION results
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
#
# nalls_consequences <- fromJSON(toJSON(content(r)))

nalls_consequences <- readRDS("R/objects/nalls_consequences.rds")

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

# GET GENES WITHIN +-1MB OF NALLS SNPS
ext <- paste0("/overlap/region/human/",
              paste0(nalls_loci$chr,
                     ":",
                     nalls_loci$start,
                     "-",
                     nalls_loci$end,
                     "?feature=gene"))

# nalls_all_genes <- lapply(ext, function(x){
#
#   r <- GET(paste(server, x, sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   fromJSON(toJSON(content(r))) %>%
#     pull(id) %>%
#     unlist()
#
# })

# names(nalls_all_genes) <- nalls_loci$SNP

# nalls_all_genes <- nalls_all_genes %>%
#   enframe %>%
#   unnest(value) %>%
#   select(SNP = name,
#          hsapiens_homolog_ensembl_gene = value)

nalls_all_genes <- readRDS("R/objects/nalls_all_genes.rds")

# make some summary numbers ----
# total number of snps
n_nalls_snps <- length(unique(nalls$SNP))
# total number of qtls selected by nalls
n_nalls_qtl_selected <- length(unique(nalls$`QTL Nominated Gene (nearest QTL)`[!is.na(nalls$`QTL Nominated Gene (nearest QTL)`)]))
# total number of genes selected by nalls
n_nalls_genes_selected <- length(unique(nalls$`Nearest Gene`))
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

n_nalls_genes_biotypes <- readRDS("R/objects/n_nalls_genes_biotypes.rds")

n_nalls_genes <- length(unique(n_nalls_genes_biotypes$external_gene_name))

# make a gene summary table ----

# summary tibble
plot_data <- tibble(
  external_gene_name = nalls_genes$external_gene_name,
  "Nalls, 2019" = TRUE,
  "Mouse\nHomolog" = external_gene_name %in% homologs$external_gene_name,
  "Detected\nin TRAP" = external_gene_name %in% nalls_trap$external_gene_name,
  "Closest Gene" = external_gene_name %in% nalls$`Nearest Gene`,
  "QTL" = external_gene_name %in% nalls$`QTL Nominated Gene (nearest QTL)`
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



# simplify the nalls results ----
nalls <- nalls %>%
  select(SNP,
         nearest_select = "Nearest Gene",
         qtl_select = "QTL Nominated Gene (nearest QTL)")

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

nalls_trap_axon_specific <- AXON_FRACTION_META %>%
  left_join(aba_expression_axon) %>%
  filter(external_gene_name %in% nalls_trap_combined$mmusculus_homolog_associated_gene_name) %>%
  select(external_gene_name,
         expression_energy) %>%
  filter(expression_energy < 0.5) %>%
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

# Categorise Nalls selections by dopa status----

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



# ----
# ----
# ----
# GWAS Nalls - old --------------------------------------------------

# number of human genes implicated
# pd_gwas_variants_genes_range$external_gene_name %>% unique %>% length
# proportion with matching homolog
# (pd_gwas_variants_genes_range$mmusculus_homolog_associated_gene_name %>%
#    unique %>% length())/(pd_gwas_variants_genes_range$external_gene_name %>% unique %>% length())

# remove human genes with no mouse homolog
# pd_gwas_variants_genes_range <- pd_gwas_variants_genes_range %>%
#   filter(mmusculus_homolog_associated_gene_name != "") %>%
#   arrange(variant)

# find all PD studies
# pd_gwas <- get_studies(efo_id = 'EFO_0002508')

# get variants from Nalls (need to think about whether to get all ~500? or Chang's)
# pd_gwas_variants <- get_variants(study_id = "GCST009325")@variants # Nalls, 2019
# pd_gwas_variants <- get_variants(study_id = "GCST004902")@variants # Chang, 2017

# # get +-1Mb position for each variant
# pd_gwas_variants_ranges <- tibble(
#   variant = pd_gwas_variants$variant_id,
#   chr = pd_gwas_variants$chromosome_name,
#   start = pd_gwas_variants$chromosome_position-1e6,
#   end = pd_gwas_variants$chromosome_position+1e6
# )
#
# # get all human genes that fall within those boundaries
# # retrieve genes within ranges including mouse homologs
# pd_gwas_variants_genes_range <- apply(pd_gwas_variants_ranges, 1, function(x){
#   getBM(attributes = c("ensembl_gene_id",
#                        "external_gene_name",
#                        "mmusculus_homolog_associated_gene_name"),
#         filters = c("chromosome_name", "start", "end"),
#         values = list(x[2],
#                       x[3],
#                       x[4]),
#         mart = ensembl_human,
#         uniqueRows = TRUE) %>%
#     mutate(variant = x[1])
# }) %>%
#   bind_rows()
#
# t_gwas_nalls <- pd_gwas_variants_genes_range %>%
#   left_join(anno_human,
#             by = c("ensembl_gene_id" = "hsapiens_homolog_ensembl_gene")) %>%
#   left_join(MB_FRACTION_META_SIMPLE, by = c("external_gene_name.y" = "external_gene_name")) %>%
#   filter(!is.na(baseMean)) %>%
#   arrange(variant, desc(score), desc(log2FoldChange)) %>%
#   group_by(variant) %>%
#   mutate(rank_trap = row_number()) %>%
#   slice_min(order_by = rank_trap,
#             n = 1) %>%
#   mutate(description.x = str_extract(description.x, ".+(?=\\[)")) %>%
#   select("Variant" = variant,
#          "TRAP Top Candidate" = mmusculus_homolog_associated_gene_name,
#          "Description" = description.x)
# ----
# ----
# ----
# TALON Novel isoform detection -------

# Load read QC file
talon_qc <- vroom("/zfs/analysis/trap/active/testing_env/talon/TALON_QC.log",
                  delim = "\t",
                  skip = 5) %>%
  mutate(passed_QC = as.factor(ifelse(passed_QC == 1, "Pass", "Fail")),
         primary_mapped = as.factor(ifelse(primary_mapped == 1, "Mapped", "Unmapped")))



# load the transcript/gene annotation file
talon_annot <- vroom("/zfs/analysis/trap/active/testing_env/talon/TALON_talon_read_annot.tsv",
                     delim = "\t")
# load the transcript/gene whitelist
talon_whitelist <- vroom("/zfs/analysis/trap/active/testing_env/talon/TALON_whitelist",
                         delim = ",",
                         col_names = F) %>%
  rename(gene_ID = X1,
         transcript_ID = X2)
# filter passed QC transcripts and genes
talon_annot <- talon_annot %>%
  filter(gene_ID %in% talon_whitelist$gene_ID &
           transcript_ID %in% talon_whitelist$transcript_ID)
# load the CAGE support data
talon_CAGE <- read_delim("/zfs/analysis/trap/active/testing_env/talon/tx_with_CAGE_support.txt",
                         delim = "\t",
                         col_names = F) %>%
  rename(annot_transcript_id = X1)
# add the CAGE support column to talon annotation
talon_annot <- talon_annot %>%
  mutate(cage = annot_transcript_id %in% talon_CAGE$annot_transcript_id)
# load the poly-A support data
talon_polyA <- vroom("/zfs/analysis/trap/active/testing_env/talon/tx_polyA_distance.txt",
                     delim = "\t",
                     col_names = F) %>%
  select("tx" = X4,
         "distance" = X18)
# add the poly-A distance data
talon_annot <- talon_annot %>%
  left_join(talon_polyA,
            by = c("annot_transcript_id" = "tx"))
# categorise whether polya within 20 bp up or downstream
talon_annot <- talon_annot %>%
  mutate(polya = distance < 20)
# categorise cage + polya support
talon_annot <- talon_annot %>%
  mutate(cage_polya = ifelse(cage & polya, "CAGE & Poly-A",
                             ifelse(cage, "CAGE",
                                    ifelse(polya, "Poly-A",
                                           "No Support"))),
         cage_polya = factor(cage_polya, levels = c("CAGE",
                                                    "Poly-A",
                                                    "CAGE & Poly-A",
                                                    "No Support")))

# Transdecoder: ORF support
talon_orfs <- vroom("/zfs/analysis/trap/active/testing_env/transdecoder/single_best_only/tx_with_orfs.txt",
                    delim = "\t",
                    col_names = "transcript_id")

talon_annot <- talon_annot %>%
  mutate(orf = annot_transcript_id %in% talon_orfs$transcript_id)

# Transdecoder:AGAT Start stop codon list
talon_start_stop_codons <- vroom("/zfs/analysis/trap/active/testing_env/transdecoder/single_best_only/agat_start_stop_transcript_entries.txt",
                                 delim = "\t",
                                 col_names = "transcript_id") %>%
  group_by(transcript_id) %>%
  summarise(n = n()) %>%
  mutate(coding_potential = ifelse(n > 1, "Coding", "Non-coding")) %>%
  select(-n)

talon_annot <- talon_annot %>%
  left_join(talon_start_stop_codons,
            by = c("annot_transcript_id" = "transcript_id")) %>%
  mutate(coding_potential = ifelse(is.na(coding_potential), "No ORF", coding_potential))

n_talon_reads <- nrow(talon_qc)
n_talon_reads_passed <- talon_qc %>%
  filter(passed_QC == "Pass") %>%
  nrow
n_talon_reads_failed <- talon_qc %>%
  filter(passed_QC == "Fail") %>%
  nrow

# talon_annot %>%
#   filter(annot_gene_name == "Syt17") %>%
#   select(annot_transcript_name,
#          gene_novelty,
#          transcript_novelty,
#          ISM_subtype) %>%
#   distinct() %>% View



# FIGURE A: Reads pass/fail ----
# the number of pass and failed reads
p_talon_pass_fail <- talon_qc %>%
  slice_sample(prop = 0.01) %>%
  ggplot(aes(x = passed_QC,
             fill = passed_QC)) +
  geom_bar(aes(y = (..count..)/sum(..count..)),
           colour = "black") +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.15))) +
  geom_text(aes( label = scales::percent((..count..)/sum(..count..)),
                 y= (..count..)/sum(..count..) ), stat= "count", vjust = -.5) +
  scale_fill_d3() +
  labs(x = "Passed QC",
       y = "Percentage of Reads",
       fill = "Passed QC") +
  theme(legend.position = "none")

# FIGURE B: Read length ----
# read length of passed and failed reads
p_talon_read_length_qc <- talon_qc %>%
  slice_sample(prop = 0.01) %>%
  ggplot(aes(x = read_length,
             fill = factor(passed_QC))) +
  geom_density(colour = "black",
               alpha = 0.75,
               aes(y = ..count..)) +
  scale_fill_d3() +
  scale_x_continuous(limits = c(0, 5000)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  # theme(legend.position = "top") +
  labs(x = "Read Length (bp)",
       y = "Number of Reads",
       fill = "Passed QC")

# the proportion of passed/failed reads with a primary mapping
# p_talon_primary_mapping <- talon_qc %>%
#   slice_sample(prop = 0.01) %>%
#   ggplot(aes(x = passed_QC,
#              fill = primary_mapped)) +
#   geom_bar(aes(y = (..count..)/sum(..count..)),
#            colour = "black") +
#   scale_y_continuous(labels = scales::percent) +
#   geom_text(aes( label = scales::percent((..count..)/sum(..count..)),
#                  y= (..count..)/sum(..count..) ), stat= "count", vjust = -.5,
#             colour = "white",
#             position=position_stack(vjust=0.5)) +
#   scale_fill_d3() +
#   labs(x = "Passed QC",
#        y = "Percentage of Reads",
#        fill = "Primary Mapping") +
#   theme(legend.position = "top")
# p_talon_primary_mapping

# the alignment fraction of passed and failed reads
# I think this means the fraction of the read that aligned?
# p_talon_alignment_fraction <- talon_qc %>%
#   slice_sample(prop = 0.01) %>%
#   ggplot(aes(x = fraction_aligned,
#              fill = factor(passed_QC))) +
#   geom_density(colour = "black",
#                alpha = 0.75,
#                aes(y = ..count..)) +
#   scale_fill_d3() +
#   scale_x_continuous() +
#   theme(legend.position = "top") +
#   labs(x = "Fraction of Alignment",
#        y = "Count",
#        fill = "Passed QC")

# FIGURE C: Transcript models by number ----
# make plot data
plot_data <- talon_annot %>%
  mutate(group = 1,
         transcript_novelty = factor(transcript_novelty,
                                     levels = c("Known",
                                                "ISM",
                                                "NIC",
                                                "NNC",
                                                "Antisense"))) %>%
  select(annot_transcript_id,
         gene_novelty,
         transcript_novelty,
         ISM_subtype,
         cage,
         polya,
         cage_polya,
         orf,
         coding_potential,
         group) %>%
  distinct() %>%
  group_by(transcript_novelty) %>%
  mutate(cage_prop = sum(cage)/n(),
         n_models = n()) %>%
  ungroup %>%
  mutate(num_prop = n_models/length(n_models))

# talon categories by number of models
p_talon_categories_numbers <- ggplot(plot_data,
                                     aes(x = transcript_novelty,
                                         fill = transcript_novelty)) +
  geom_bar(
    stat="count",
    colour = "black") +
  scale_fill_d3() +
  theme(legend.position = "none") +
  labs(x = "Transcript Novelty",
       y = "Number of Models") +
  geom_text(aes( label = scales::percent(num_prop),
                 y= n_models ), vjust = -.5,
            colour = "black",
            check_overlap = T) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# FIGURE D: Transcript models by library percentage ----
# talon categories by percentage of library
p_talon_categories_representation <- talon_annot %>%
  mutate(group = 1,
         transcript_novelty = factor(transcript_novelty,
                                     levels = c("Known",
                                                "ISM",
                                                "NIC",
                                                "NNC",
                                                "Antisense"))) %>%
  group_by(transcript_novelty) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = transcript_novelty,
             y = prop,
             fill = transcript_novelty)) +
  geom_col(colour = "black") +
  scale_fill_d3() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank()) +
  labs(x = "Transcript Novelty",
       y = "Percentage of Library")  +
  geom_text(aes( label = scales::percent(prop) ), vjust = -.5) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.15)))

# FIGURE E: CAGE/Poly-A Support by model ----
# percentage of transcript models with CAGE and/or poly-A support
p_talon_cage_polya <- plot_data %>%
  ggplot(aes(x = transcript_novelty,
             fill = cage_polya)) +
  geom_bar(position = "fill",
           colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Transcript Novelty",
       y = "Percentage",
       fill = "Model Support") +
  theme(legend.position = "none")

# CAGE support
# p_talon_cage <- ggplot(plot_data,
#        aes(x = transcript_novelty,
#              # fill = transcript_novelty,
#              group = group)) +
#   geom_bar(
#     # aes(y = ..prop..),
#            stat="count",
#            colour = "black",
#            fill = NA) +
#   geom_bar(data = (plot_data %>% filter(cage)),
#            aes(fill = factor(..x..)),
#            stat="count",
#            colour = "black") +
#   scale_fill_d3() +
#   theme(legend.position = "none") +
#   labs(x = "Transcript Novelty",
#        y = "Number of Models") +
#   geom_text(aes( label = scales::percent(cage_prop),
#                  y= n_models ), vjust = -.5, check_overlap = T)


# distance to nearest poly-a site (probably won't include)
# talon_annot %>%
#   slice_sample(prop = 0.01) %>%
#   ggplot(aes(x = distance+1,
#              y = transcript_novelty,
#              fill = transcript_novelty)) +
#   geom_density_ridges(alpha = 0.75,
#                       scale = 0.9) +
#   # scale_x_continuous(limits = c(0, 1000)) +
#   scale_x_log10() +
#   scale_fill_d3() +
#   labs(x = "Distance to Nearest Poly-A Signal",
#        y = "Transcript Novelty") +
#   theme(legend.position = "none")

# transcript novelty category filled by CAGE and poly-A support categories
# talon_annot %>%
#   ggplot(aes(x = transcript_novelty,
#              group = cage_polya)) +
#   geom_bar(aes(y = ..prop..,
#                fill = cage_polya),
#            stat="count",
#            colour = "black") +
#   scale_fill_d3()

# FIGURE F: ISM subtypes and CAGE/polyA support ----
# ISM subtype distribution: More Suffix models, because of 3' bias reasons:
# RNA degradation, incomplete RT, pore blocking
p_talon_cage_polya_ism <- plot_data %>%
  filter(transcript_novelty == "ISM") %>%
  ggplot(aes(x = ISM_subtype,
             fill = cage_polya)) +
  geom_bar(colour = "black") +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  # theme(legend.position = "none") +
  labs(x = "ISM Subtype",
       y = "Number of Models",
       fill = "Model Support")

# FIGURE G: Coding potential by transcript support category ----
# coding potential (by ORF and codons) per model support category (CAGE/polyA)
p_talon_coding_potential_by_cage_polya <- plot_data %>%
  mutate(coding_potential = factor(coding_potential,
                                   levels = c("Coding",
                                              "Non-coding",
                                              "No ORF")),
         cage_polya = factor(cage_polya,
                             levels = c("No Support",
                                        "CAGE",
                                        "Poly-A",
                                        "CAGE & Poly-A"))) %>%
  filter(transcript_novelty != "Known") %>%
  ggplot(aes(x = cage_polya,
             fill = coding_potential)) +
  geom_bar(position = "fill",
           colour = "black") +
  scale_fill_manual(values = pal_d3()(4)[c(3, 1, 4)]) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Model Support",
       y = "Percentage of Models",
       fill = "Coding Potential") +
  theme(legend.position = "top")

# FIGURE H: Transcript models library percentage after filtering for CAGE & polyA support ----
# after filtering transcripts with CAGE and poly-A support
# plot the proportion of library represented by novelty category
p_talon_categories_representation_filtered <- talon_annot %>%
  filter(cage_polya == "CAGE & Poly-A" &
           coding_potential != "No ORF") %>%
  group_by(transcript_novelty,
           coding_potential) %>%
  summarise(n = n()) %>%
  ungroup %>%
  mutate(prop = n / sum(n)) %>%
  group_by(transcript_novelty) %>%
  mutate(prop_transcript_novelty = sum(prop)) %>%
  ggplot(aes(x = transcript_novelty,
             fill = coding_potential)) +
  geom_col(colour = "black",
           aes(y = prop)) +
  geom_text(aes(y = prop_transcript_novelty,
                label = scales::percent(prop_transcript_novelty)),
            check_overlap = T,
            vjust = -0.5) +
  scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2)]) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Transcript Novelty",
       y = "Percentage of Library",
       fill = "Coding Potential") +
  theme(legend.position = "top")

# plot_data %>%
#   filter(cage_polya == "CAGE & Poly-A" &
#            coding_potential != "No ORF") %>%
#   group_by(transcript_novelty,
#            coding_potential) %>%
#   summarise(n = n()) %>%
#   ungroup %>%
#   mutate(prop = n / sum(n)) %>%
#   group_by(transcript_novelty) %>%
#   mutate(prop_transcript_novelty = sum(prop)) %>%
#   ggplot(aes(x = transcript_novelty,
#              fill = coding_potential)) +
#   geom_col(colour = "black",
#            aes(y = prop)) +
#   geom_text(aes(y = prop_transcript_novelty,
#                 label = scales::percent(prop_transcript_novelty)),
#             check_overlap = T,
#             vjust = -0.5) +
#   scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2)]) +
#   scale_y_continuous(labels = scales::percent,
#                      limits = c(0, 1)) +
#   labs(x = "Transcript Novelty",
#        y = "Percentage of Library",
#        fill = "Coding Potential") +
#   theme(legend.position = "top")




# Finding an interesting novel transcript - see comments ----

# talon_annot %>%
#   filter(orf &
#            cage_polya == "CAGE & Poly-A" &
#            transcript_novelty != "Known") %>%
#   pull(annot_transcript_id) %>%
#   unique %>%
#   length
# colnames(talon_annot)
# talon_annot %>%
#   slice_head(n = 10) %>%
#   View
# talon_annot %>%
#   filter(transcript_novelty == "Known" |
#            orf &
#            cage_polya == "CAGE & Poly-A"
#          ) %>%
#   # select(annot_gene_id,
#   #        annot_transcript_id,
#   #        annot_gene_name,
#   #        annot_transcript_name,
#   #        transcript_novelty,
#   #        ISM_subtype,
#   #        coding_potential) %>%
#   group_by(annot_gene_id,
#            annot_transcript_id,
#            annot_gene_name,
#            annot_transcript_name,
#            transcript_novelty,
#            ISM_subtype,
#            coding_potential) %>%
#   summarise(n = n()) %>%
#   group_by(annot_gene_id,
#            annot_gene_name) %>%
#   mutate(prop = paste0(signif((100 * n / sum(n)), 3), "%")) %>%
#   group_by(annot_gene_id) %>%
#   filter(sum(transcript_novelty != "Known") >= 1) %>% View
#
# talon_annot %>%
#   filter(transcript_novelty == "Known" |
#            cage_polya != "No Support"
#   ) %>%
#   group_by(annot_gene_id,
#            annot_transcript_id,
#            annot_gene_name,
#            annot_transcript_name,
#            transcript_novelty,
#            ISM_subtype,
#            coding_potential) %>%
#   summarise(n = n()) %>%
#   group_by(annot_gene_id,
#            annot_gene_name) %>%
#   mutate(prop = paste0(signif((100 * n / sum(n)), 3), "%")) %>%
#   group_by(annot_gene_id) %>%
#   filter(sum(transcript_novelty != "Known") >= 1) %>% View

# There are many cases of exon boundaries being 1-9 bp longer or shorter than existing reference transcripts
# Genes where this happens: Fdft1, Macroh2a1, Ubxn6
# Short read data does not support these novel boundaries
# This is likely an error of ONT sequencing/alignment, which should be discussed.
# This needs to be written about: ONT sequencing or alignment error around exon boundaries
# Also potential skipping of small exons, such at in Syt17: Almost no short-read support for this skipping
# Variable 5' and 3' UTRs are found with TALON, but these are not primary interest,
# and also there is scope for just premature UTR termination
# As a solution, I'm now trying the reverse: rMATS to look for skipped exons primarily.
# Then manually see if these skipped exons are in the long read data

# ----
# ----
# ----
# # LIBRARY METRICS ----
# # Make and enrichment tibble ----
# # Make an enrichment status tibble
# enrichment_status <- tibble(
#   ensembl_gene_id = rownames(dds_C2_IP),
#   enrichment = ifelse(
#     ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#     ifelse(
#       ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
#       "Soma and Axon",
#       "Soma only"
#     ),
#     ifelse(
#       ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
#       "Axon only",
#       "No enrichment"
#     )
#   ),
# ) %>%
#   mutate(
#     enrichment = factor(
#       enrichment,
#       levels = c("Soma only",
#                  "Axon only",
#                  "Soma and Axon",
#                  "No enrichment")
#     )
#   )
#
# # Make a library proportion tibble - by compartment ----
# # Make a library proportion tibble - split by compartment
# library_proportion_info <- enrichment_status %>%
#   left_join(enframe(rowMeans(counts(dds_C2_IP)[, colData(dds_C2_IP)$compartment == "MB"]),
#                     name = "ensembl_gene_id",
#                     value = "mb")) %>%
#   left_join(enframe(rowMeans(counts(dds_C2_IP)[, colData(dds_C2_IP)$compartment == "AXON"]),
#                     name = "ensembl_gene_id",
#                     value = "axon"))
#
#
# library_proportion_info_total <- enrichment_status %>%
#   left_join(enframe(rowMeans(counts(dds_C2_TOTAL)[, colData(dds_C2_TOTAL)$compartment == "MB"]),
#                     name = "ensembl_gene_id",
#                     value = "mb")) %>%
#   left_join(enframe(rowMeans(counts(dds_C2_TOTAL)[, colData(dds_C2_TOTAL)$compartment == "AXON"]),
#                     name = "ensembl_gene_id",
#                     value = "axon"))
#
# # Plot just the proportion of library in each enrichment category by sample type (axon or mb) ----
#
#
# p_library_proportion_compartment <- bind_rows({
#   library_proportion_info %>%
#     pivot_longer(c(mb, axon),
#                  names_to = "compartment",
#                  values_to = "count") %>%
#     mutate(compartment = ifelse(compartment == "mb", "Soma", "Axon")) %>%
#     group_by(compartment,
#              enrichment) %>%
#     summarise(count = sum(count)) %>%
#     mutate(prop = count / sum(count),
#            fraction = "TRAP")
# }, {
#   library_proportion_info_total %>%
#     drop_na %>%
#     pivot_longer(c(mb, axon),
#                  names_to = "compartment",
#                  values_to = "count") %>%
#     mutate(compartment = ifelse(compartment == "mb", "Soma", "Axon")) %>%
#     group_by(compartment,
#              enrichment) %>%
#     summarise(count = sum(count)) %>%
#     mutate(prop = count / sum(count),
#            fraction = "TOTAL")
# }) %>%
#   ggplot(aes(x = compartment,
#              y = prop,
#              fill = enrichment)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   scale_y_continuous(
#     labels = scales::percent,
#     limits = c(0, 1),
#     expand = expansion(mult = c(0, 0.05))
#   ) +
#   labs(x = "TRAP Fraction and Compartment",
#        y = "Percentage of Library",
#        fill = "Enrichment") +
#   geom_text(
#     aes(
#       x = compartment,
#       y = prop,
#       label = ifelse(prop > 0.05,
#                      scales::percent(prop,
#                                      accuracy = 1),
#                      "")
#     ),
#     position = position_stack(vjust = 0.5),
#     colour = "white",
#     fontface = 2
#   ) +
#   theme(legend.position = "top") +
#   guides(fill = guide_legend(
#     nrow = 1,
#     byrow = TRUE,
#     title.position = "top",
#     title.hjust = 0.5
#   ))  +
#   facet_wrap(vars(fraction)) +
#   panel_border()
#
#
# # Plot the percentage of soma enrichment per cohort and age in soma TRAP samples ----------
#
#
# p_mb_ip_age_library_composition <- bind_rows({
#   enrichment_status %>%
#     left_join(enframe(rowMeans(counts(dds_C1_MB_IP)[, colData(dds_C1_MB_IP)$age == "YOUNG"]),
#                       name = "ensembl_gene_id",
#                       value = "Young")) %>%
#     left_join(enframe(rowMeans(counts(dds_C1_MB_IP)[, colData(dds_C1_MB_IP)$age == "OLD"]),
#                       name = "ensembl_gene_id",
#                       value = "Old")) %>%
#
#     pivot_longer(c(Young, Old),
#                  names_to = "age",
#                  values_to = "count") %>%
#     mutate(cohort = "C1")
# }, {
#   enrichment_status %>%
#     left_join(enframe(rowMeans(counts(dds_C2_MB_IP)[, colData(dds_C2_MB_IP)$age == "YOUNG"]),
#                       name = "ensembl_gene_id",
#                       value = "Young")) %>%
#     left_join(enframe(rowMeans(counts(dds_C2_MB_IP)[, colData(dds_C2_MB_IP)$age == "OLD"]),
#                       name = "ensembl_gene_id",
#                       value = "Old")) %>%
#
#     pivot_longer(c(Young, Old),
#                  names_to = "age",
#                  values_to = "count") %>%
#     mutate(cohort = "C2")
# }, {
#   enrichment_status %>%
#     left_join(enframe(rowMeans(counts(dds_C3_MB_IP)[, colData(dds_C3_MB_IP)$age == "YOUNG"]),
#                       name = "ensembl_gene_id",
#                       value = "Young")) %>%
#     left_join(enframe(rowMeans(counts(dds_C3_MB_IP)[, colData(dds_C3_MB_IP)$age == "OLD"]),
#                       name = "ensembl_gene_id",
#                       value = "Old")) %>%
#
#     pivot_longer(c(Young, Old),
#                  names_to = "age",
#                  values_to = "count") %>%
#     mutate(cohort = "C3")
# }) %>%
#   mutate(age = factor(age, levels = c("Young", "Old")),
#          enrichment = ifelse(str_detect(enrichment, "oma"), "Enriched", "Not enriched"),
#          enrichment = factor(enrichment, levels = c("Not enriched", "Enriched"))) %>%
#   drop_na %>%
#   group_by(age,
#            cohort,
#            enrichment) %>%
#   summarise(count = sum(count)) %>%
#   mutate(prop = count / sum(count)) %>%
#   ggplot(aes(x = age,
#              y = prop,
#              fill = enrichment)) +
#   geom_col(colour = "black") +
#   scale_fill_manual(values = pal_d3()(2)[c(2, 1)]) +
#   scale_y_continuous(
#     labels = scales::percent,
#     limits = c(0, 1),
#     expand = expansion(mult = c(0, 0.05))
#   ) +
#   labs(x = "TRAP Cohort and Age",
#        y = "Percentage of Library",
#        fill = "Enrichment") +
#   geom_text(
#     aes(
#       x = age,
#       y = prop,
#       label = ifelse(prop > 0.05,
#                      scales::percent(prop,
#                                      accuracy = 1),
#                      "")
#     ),
#     position = position_stack(vjust = 0.5),
#     colour = "white",
#     fontface = 2
#   ) +
#   theme(legend.position = "top") +
#   guides(fill = guide_legend(
#     nrow = 1,
#     byrow = TRUE,
#     title.position = "top",
#     title.hjust = 0.5
#   )) +
#   facet_wrap(vars(cohort))
#
#
#
# # The number of genes in each category ----------
# p_library_composition_enrichment_n_genes <- library_proportion_info %>%
#   filter(enrichment != "No enrichment") %>%
#   group_by(enrichment) %>%
#   tally %>%
#   mutate(prop = n/sum(n)) %>%
#   ggplot(aes(x = enrichment,
#              y = n,
#              fill = enrichment)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   geom_text(aes(label = n),
#             vjust = -0.5,
#             fontface = 2) +
#   geom_text(aes(label = scales::percent(prop),
#                 y = n/2),
#             colour = "white",
#             fontface = 2) +
#   labs(x = "Enrichment",
#        y = "Number of Genes") +
#   theme(legend.position = "none")
#
#
#
#
# # ----
# # ----
# # ----
# # Basic comparison of AXON VS MB ----
#
# # Plot the counts of Soma only genes by fraction and compartment: Are axonal TRAP samples depleted of
# # some DA transcripts?
# # bind_rows(
# #   {library_proportion_info_total %>% mutate(fraction = "TOTAL")},
# #   {library_proportion_info %>% mutate(fraction = "TRAP")}
# # ) %>%
# #   # filter(enrichment == "Soma only") %>%
# #   pivot_longer(c(mb, axon),
# #                names_to = "compartment",
# #                values_to = "count") %>%
# #   ggplot(aes(x = fraction,
# #              y = count)) +
# #   geom_boxplot() +
# #   scale_y_log10() +
# #   facet_wrap(vars(compartment,
# #                   enrichment))
# #
# #
# library_proportion_info %>%
#   group_by(enrichment) %>%
#   summarise(mean(mb))
#   # filter(enrichment == "Axon only") %>%
#   pivot_longer(c(mb, axon),
#                names_to = "compartment",
#                values_to = "count") %>%
#   ggplot(aes(x = compartment,
#              y = count)) +
#   geom_boxplot() +
#   scale_y_log10()
# #
# # library_proportion_info %>%
# #   filter(enrichment == "Soma and Axon") %>%
# #   pivot_longer(c(mb, axon),
# #                names_to = "compartment",
# #                values_to = "count") %>%
# #   ggplot(aes(x = count,
# #              colour = compartment)) +
# #   stat_ecdf() +
# #   scale_x_log10()
# #
# #
# # axon_vs_mb_trap <- library_proportion_info %>%
# #   filter(enrichment != "No enrichment") %>%
# #   mutate(diff = log2(axon/mb)) %>%
# #   arrange(desc(diff)) %>%
# #   mutate(compartment = ifelse(diff > 1, "axon", "mb"),
# #          dopaminergic = ensembl_gene_id %in% putative_dopaminergic) %>%
# #   left_join(anno) %>%
# #   # filter(mb > 20 | axon > 20) %>%
# #   group_by(enrichment) %>%
# #   arrange(desc(diff)) %>%
# #   mutate(rank = row_number()/dplyr::n()) %>%
# #   ggplot(aes(x = rank,
# #              y = diff,
# #              colour = enrichment)) +
# #   geom_line() +
# #   geom_hline(yintercept = 0) +
# #   geom_vline(xintercept = 0.5)
# #
# #
# #
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
# ggplot(aes(x = rank,
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
#
# #
# #
# #   axon_vs_mb_trap <- library_proportion_info %>%
# #     # filter(enrichment != "No enrichment") %>%
# #     mutate(diff = log2(axon/mb)) %>%
# #     arrange(desc(diff)) %>%
# #     mutate(compartment = ifelse(diff > 1, "axon", "mb"),
# #            dopaminergic = ensembl_gene_id %in% putative_dopaminergic,
# #            non_dopaminergic = ensembl_gene_id %in% putative_nondopaminergic) %>%
# #     left_join(anno) %>%
# #     left_join(publications_all_genes) %>%
# #     left_join(publications_all_genes_axon) %>%
# #     # filter(mb > 20 | axon > 20) %>%
# #     group_by(enrichment) %>%
# #     arrange(desc(diff)) %>%
# #     mutate(fraction = "TRAP") %>%
# #     left_join(axon_vs_mb_trap_total)
# #
# #
# #
# #   axon_vs_mb_trap %>%
# #     filter(!enrichment %in% c("No enrichment", "Axon only")) %>%
# #     filter(!is.infinite(diff)) %>%
# #     filter(diff_total < 2 &
# #              diff_total > -2)
# #
# #
# #   ggplot(aes(x = diff,
# #                y = diff_total)) +
# #     geom_point(alpha = 0.1) +
# #     stat_cor() +
# #     geom_hline(yintercept = -2) +
# #     geom_hline(yintercept = 2)
# #
# #
#   axon_vs_mb_trap_total <- library_proportion_info_total %>%
#     filter(enrichment != "No enrichment") %>%
#     mutate(diff = log2(axon / mb)) %>%
#     arrange(desc(diff)) %>%
#     mutate(
#       compartment = ifelse(diff > 1, "axon", "mb"),
#       dopaminergic = ensembl_gene_id %in% putative_dopaminergic,
#       non_dopaminergic = ensembl_gene_id %in% putative_nondopaminergic
#     ) %>%
#     left_join(anno) %>%
#     left_join(publications_all_genes) %>%
#     left_join(publications_all_genes_axon) %>%
#     # filter(mb > 20 | axon > 20) %>%
#     group_by(enrichment) %>%
#     arrange(desc(diff)) %>%
#     mutate(fraction = "TRAP") %>%
#     ungroup() %>%
#     select(ensembl_gene_id,
#            diff_total = diff)
#
#
#
#
# # # get entrez IDs
# # anno_fgsea <- getBM(
# #   attributes = c("ensembl_gene_id",
# #                  "hsapiens_homolog_associated_gene_name"),
# #   filters = "ensembl_gene_id",
# #   values = axon_vs_mb_trap$ensembl_gene_id,
# #   mart = ensembl,
# #   uniqueRows = TRUE
# # )
# #
# # temp <- gost(axon_vs_mb_trap$ensembl_gene_id,
# #      organism = "mmusculus",
# #      ordered_query = T,
# #      evcodes = T)
# # temp$result %>% filter(source != "TF") %>% View
# #
# #
# # pd_genes <- temp$result %>%
# #   filter(term_name == "Parkinson disease") %>%
# #   pull(intersection) %>%
# #   strsplit(",") %>%
# #   unlist
# #
# # axon_vs_mb_trap %>%
# #   filter(ensembl_gene_id %in% pd_genes) %>% View
# #   sapply(., get_external_gene_name) %>%
# #   unique
# #
# #
# #
# # # load gmt
# # # pathways <- fgsea::gmtPathways("R/objects/msigdb.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c2.cp.kegg.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.bp.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.mf.v7.4.symbols.gmt.txt")
# #
# # ranks <- axon_vs_mb_trap %>%
# #   left_join(anno_fgsea) %>%
# #   select(hsapiens_homolog_associated_gene_name,
# #          rank = diff) %>%
# #   filter(hsapiens_homolog_associated_gene_name != "") %>%
# #   distinct(across(-rank), .keep_all = T) %>%
# #   deframe
# #
# # names(ranks)
# #
# # AXON_VS_MB_FGSEA <- fgsea(pathways,
# #                            ranks,
# #                            eps = 0,
# #                            minSize = 10,
# #                            maxSize = 500,
# #                            nPermSimple = 10000,
# #                            BPPARAM = MulticoreParam(),
# #                            nproc = 22)
# #
# # plot(ranks)
# #
# #
# # # Soma only genes are abundant in soma samples and very low in axon samples
# # # Soma and axon genes are either shared DA neuron genes, with greater expression in axonal compartment
# # # or are genes that are less specific to DA neurons compared to soma only genes
# # # Can check that by comparing the log2foldchange
# #
# # # Axon only genes take up a tiny proportion of soma samples: These genes are likely
# # # not DAergic at all
# #
# # # Not enriched genes are of no interest and take up less than 10% of the library,
# # # so dont plot them in the next plots
# #
# # # Soma and axon genes are less specific: the log2foldchange is lower
# #
# # p_enrichment_specificity_lfc <- MB_FRACTION_META %>%
# #   filter(sumz_adj < 0.01) %>%
# #   select(ensembl_gene_id,
# #          log2FoldChange) %>%
# #   inner_join(enrichment_status) %>%
# #   filter(enrichment != "No enrichment") %>%
# #   ggplot(aes(x = log2FoldChange,
# #              colour = enrichment)) +
# #   geom_density() +
# #   scale_color_d3() +
# #   scale_x_continuous(limits = c(-3, 3))
# #
# #
# # AXON_FRACTION_META %>%
# #   filter(sumz_adj < 0.01) %>%
# #   select(ensembl_gene_id,
# #          log2FoldChange) %>%
# #   inner_join(enrichment_status) %>%
# #   filter(enrichment == "Soma and Axon") %>%
# #   select(ensembl_gene_id,
# #          log2FoldChange) %>%
# #   left_join(
# #     {MB_FRACTION_META %>%
# #         filter(sumz_adj < 0.01) %>%
# #         select(ensembl_gene_id,
# #                log2FoldChange) %>%
# #         inner_join(enrichment_status) %>%
# #         filter(enrichment == "Soma and Axon") %>%
# #         select(ensembl_gene_id,
# #                log2FoldChange)},
# #     by = "ensembl_gene_id",
# #     suffix = c("_axon", "_soma")
# #   ) %>%
# #   left_join(anno) %>% View
# #
# #
# # # Plot the number of genes in each enrichment group by compartment preference
# # library_proportion_info %>%
# #   filter(enrichment != "No enrichment") %>%
# #   inner_join(AXON_VS_MB_META %>% filter(padj < 0.01) %>% select(ensembl_gene_id, outcome = compartment)) %>%
# #   group_by(outcome,
# #            enrichment) %>%
# #   tally %>%
# #   mutate(prop = n/sum(n),
# #          sum = sum(n)) %>%
# #   ggplot(aes(x = outcome,
# #              y = prop,
# #              fill = enrichment)) +
# #   geom_col(colour = "black") +
# #   scale_fill_d3() +
# #   geom_text(aes(label = paste("Total:", sum),
# #                 y = 1),
# #             check_overlap = T,
# #             vjust = -0.5,
# #             fontface = 2,
# #             size = 4) +
# #   scale_y_continuous(labels = scales::percent,
# #                      limits = c(0, 1),
# #                      expand = expansion(mult = c(0, 0.1))) +
# #   labs(x = "Change in Age",
# #        y = "Percentage of Genes",
# #        fill = "Enrichment") +
# #   geom_text(
# #     aes(
# #       x = outcome,
# #       y = prop,
# #       label = ifelse(prop > 0.05,
# #                      scales::percent(prop,
# #                                      accuracy = 1),
# #                      "")
# #     ),
# #     position = position_stack(vjust = 0.5),
# #     colour = "white",
# #     fontface = 2
# #   ) +
# #   theme(
# #     legend.position = "top"
# #   ) +
# #   guides(fill = guide_legend(nrow = 2,
# #                              byrow=TRUE,
# #                              title.position="top",
# #                              title.hjust = 0.5))
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # # Plot the number of genes in each enrichment group by age direction
# # p_axon_age_library_composition <- library_proportion_info %>%
# #   inner_join(AXON_AGE_META %>% filter(padj < 0.01) %>% select(ensembl_gene_id, outcome = outcome_TRAP)) %>%
# #   group_by(outcome,
# #            enrichment) %>%
# #   tally %>%
# #   mutate(prop = n/sum(n),
# #          sum = sum(n)) %>%
# #   ggplot(aes(x = outcome,
# #              y = prop,
# #              fill = enrichment)) +
# #   geom_col(colour = "black") +
# #   scale_fill_d3() +
# #   geom_text(aes(label = paste("Total:", sum),
# #                 y = 1),
# #             check_overlap = T,
# #             vjust = -0.5,
# #             fontface = 2,
# #             size = 4) +
# #   scale_y_continuous(labels = scales::percent,
# #                      limits = c(0, 1),
# #                      expand = expansion(mult = c(0, 0.1))) +
# #   labs(x = "Change in Age",
# #        y = "Percentage of Genes",
# #        fill = "Enrichment") +
# #   geom_text(
# #     aes(
# #       x = outcome,
# #       y = prop,
# #       label = ifelse(prop > 0.05,
# #                      scales::percent(prop,
# #                                      accuracy = 1),
# #                      "")
# #     ),
# #     position = position_stack(vjust = 0.5),
# #     colour = "white",
# #     fontface = 2
# #   ) +
# #   theme(
# #     legend.position = "top"
# #   ) +
# #   guides(fill = guide_legend(nrow = 2,
# #                              byrow=TRUE,
# #                              title.position="top",
# #                              title.hjust = 0.5))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # ----
# # ----
# # ----
# # AXON VS MB - LIBRARY COMPOSITION CHECK ----
#
# C2_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C2_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype", "region", "compartment"), returnData = T) %>%
#       select(age, genotype, region, compartment, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C2_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#     # }, {
#     #   plotCounts(dds_C2_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#     #     select("Ddc" = count)
#   }, {
#     plotCounts(dds_C2_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C2_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C2_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C2_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype, region, compartment),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C2")
#
#
# C3_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C3_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype", "region", "compartment"), returnData = T) %>%
#       select(age, genotype, region, compartment, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C3_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#     # }, {
#     #   plotCounts(dds_C3_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#     #     select("Ddc" = count)
#   }, {
#     plotCounts(dds_C3_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C3_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C3_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C3_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype, region, compartment),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C3")
#
#
# p_axon_vs_mb_ip_library_composition <- bind_rows(C2_IP_MARKERS,
#                                                  C3_IP_MARKERS) %>%
#   mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
#                          "Dopaminergic", "Other"),
#          compartment = ifelse(compartment == "MB", "Soma", "Axon")) %>%
#   filter(count > 1) %>%
#   # group_by(age, gene) %>%
#   # summarise(count = mean(count)) %>%
#   ggplot(aes(x = compartment,
#              y = log2(count + 1),
#              colour = marker,
#              group = gene)) +
#   geom_quasirandom() +
#   # geom_smooth(method = "lm") +
#   # facet_wrap(vars(cohort), nrow = 1) +
#   panel_border() +
#   scale_color_d3() +
#   labs(x = "Compartment",
#        y = expression(Log[2] ~ Expression ~ Count),
#        colour = "Cell Type Marker") +
#   theme(legend.position = "top") +
#   guides(colour = guide_legend(title.position="top", title.hjust = 0.5))
#
# # AXON VS MB - ALTERNATIVE COMPOSITION CHECK ----
#
# LIBRARY_PROPORTIONS <- bind_rows(
#   {
#     as_tibble(counts(dds_C2_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C2")
#   }, {
#     as_tibble(counts(dds_C3_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C3")
#   }
# ) %>%
#   pivot_longer(-c(ensembl_gene_id, cohort),
#                names_to = "sample_name",
#                values_to = "count") %>%
#   left_join(colData(dds), copy = T) %>%
#   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
#          enrichment = ifelse(mb_translated & axon_translated,
#                            "Soma and Axon",
#                            ifelse(mb_translated,
#                                   "Soma only",
#                                   ifelse(axon_translated,
#                                          "Axon only",
#                                          "No enrichment"))),
#          enrichment = factor(enrichment,
#                            levels = c("Soma only",
#                                       "Axon only",
#                                       "Soma and Axon",
#                                       "No enrichment"
#                            )),
#          compartment = ifelse(compartment == "MB", "Soma", "Axon")) %>%
#   group_by(compartment, enrichment) %>%
#   summarise(count = sum(count, na.rm = T)) %>%
#   mutate(prop = count/sum(count)) %>%
#   drop_na()
#
#
#
# # LIBRARY_PROPORTIONS %>%
# #   ggplot(aes(x = compartment,
# #              y = prop,
# #              fill = enrichment)) +
# #   geom_col(colour = "black") +
# #   # facet_wrap(vars(cohort), nrow = 1)  +
# #   theme(legend.position = "top") +
# #   geom_text(aes(y = prop, label = ifelse(prop > 0.1, scales::percent(prop, accuracy = 1), "")),
# #             position = position_stack(vjust = 0.5),
# #             show.legend = FALSE,
# #             colour = "white",
# #             fontface = 2) +
# #   scale_fill_d3() +
# #   # scale_fill_manual(values = pal_d3()(4)[c(4, 1, 3)]) +
# #   scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
# #                      labels = scales::percent) +
# #   # panel_border() +
# #   labs(x = "Compartment",
# #        y = "Percentage of library",
# #        fill = "Enrichment")  +
# #   guides(fill = guide_legend(title.position="top",
# #                              title.hjust = 0.5,
# #                              nrow = 2,
# #                              byrow = T))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # AXON VS MB: THIS IS COMPLETE ----------------------
# # AXON VS MB: TOTAL DESeq2 ----
#
# dds_C2_TOTAL_COMPARTMENT_FILTER <- filter_genes(dds_C2_TOTAL,
#                                              grouping = c("compartment", "gene_id"))
# dds_C3_TOTAL_COMPARTMENT_FILTER <- filter_genes(dds_C3_TOTAL,
#                                              grouping = c("compartment", "gene_id"))
#
# genes_TOTAL_COMPARTMENT <- intersect(dds_C2_TOTAL_COMPARTMENT_FILTER,
#                                   dds_C3_TOTAL_COMPARTMENT_FILTER)
#
# dds_C2_TOTAL_COMPARTMENT <- dds_C2_TOTAL[rownames(dds_C2_TOTAL) %in% genes_TOTAL_COMPARTMENT,]
# dds_C3_TOTAL_COMPARTMENT <- dds_C3_TOTAL[rownames(dds_C3_TOTAL) %in% genes_TOTAL_COMPARTMENT,]
#
# dds_C2_TOTAL_COMPARTMENT@design <- ~ compartment
# dds_C3_TOTAL_COMPARTMENT@design <- ~ compartment
#
# colData(dds_C2_TOTAL_COMPARTMENT) <- droplevels(colData(dds_C2_TOTAL_COMPARTMENT))
# colData(dds_C3_TOTAL_COMPARTMENT) <- droplevels(colData(dds_C3_TOTAL_COMPARTMENT))
#
# dds_C2_TOTAL_COMPARTMENT <- DESeq(
#   dds_C2_TOTAL_COMPARTMENT,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# dds_C3_TOTAL_COMPARTMENT <- DESeq(
#   dds_C3_TOTAL_COMPARTMENT,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_C2_TOTAL_COMPARTMENT <- DESeq2::results(
#   dds_C2_TOTAL_COMPARTMENT,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C2_TOTAL_COMPARTMENT <- lfcShrink(
#   dds_C2_TOTAL_COMPARTMENT,
#   res = res_C2_TOTAL_COMPARTMENT,
#   contrast = c("compartment", "AXON", "MB"),
#   type = "ashr",
#   parallel = T
# )
#
# res_C3_TOTAL_COMPARTMENT <- DESeq2::results(
#   dds_C3_TOTAL_COMPARTMENT,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C3_TOTAL_COMPARTMENT <- lfcShrink(
#   dds_C3_TOTAL_COMPARTMENT,
#   res = res_C3_TOTAL_COMPARTMENT,
#   contrast = c("compartment", "AXON", "MB"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_C2_TOTAL_COMPARTMENT)
# summary(res_C3_TOTAL_COMPARTMENT)
#
# AXON_VS_MB_TOTAL_META <-
#   as_tibble(res_C2_TOTAL_COMPARTMENT, rownames = "ensembl_gene_id") %>%
#   left_join(
#     as_tibble(res_C3_TOTAL_COMPARTMENT, rownames = "ensembl_gene_id"),
#     by = "ensembl_gene_id",
#     suffix = c("_C2", "_C3")
#   ) %>%
#   left_join(anno) %>%
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
#   mutate(across(starts_with("pvalue"),
#                 ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
#   mutate(
#     sumz = sumz(c_across(pvalue_C2:pvalue_C3))$p,
#     log2FoldChange = log2FoldChange_C2
#   ) %>%
#   ungroup() %>%
#   mutate(
#     sumz_adj = p.adjust(sumz, method = "fdr"),
#     conflict = sign(log2FoldChange_C2) != sign(log2FoldChange_C3)) # state whether there is a conflict in l2fc direction
#
# # AXON VS MB: TOTAL outcomes plot ----
#
# AXON_VS_MB_TOTAL_OUTCOMES <- AXON_VS_MB_TOTAL_META %>%
#   mutate(outcome = ifelse(sumz_adj > 0.01,
#                           "Unchanged",
#                           ifelse(log2FoldChange > 0,
#                                  "Axon", "Soma"))) %>%
#   select(external_gene_name,
#          outcome_total = outcome)
#
# p_axon_vs_mb_total_outcome <- AXON_VS_MB_TOTAL_OUTCOMES %>%
#   filter(outcome_total != "Unchanged") %>%
#   ggplot(aes(x = outcome_total,
#              fill = outcome_total)) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2, 4)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   labs(x = "TOTAL Region",
#        y = "Number of Genes") +
#   geom_text(aes(label = paste("Total:", ..count..),
#                 y = ..count..),
#             stat = "count",
#             check_overlap = T,
#             vjust = -0.5,
#             fontface = 2,
#             size = 4) +
#   theme(legend.position = "none")
#
# # AXON vs MB DESeq2 -------------------------------------------------------
#
# dds_C2_IP@design <- ~ collection + age + compartment
#
# colData(dds_C2_IP) <- droplevels(colData(dds_C2_IP))
#
# dds_C2_IP_COMPARTMENT <- dds_C2_IP %>% filter_zeros()
#
# dds_C2_IP_COMPARTMENT <- DESeq(
#   dds_C2_IP_COMPARTMENT,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_C2_IP_COMPARTMENT <- DESeq2::results(
#   dds_C2_IP_COMPARTMENT,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C2_IP_COMPARTMENT <- lfcShrink(
#   dds_C2_IP_COMPARTMENT,
#   res = res_C2_IP_COMPARTMENT,
#   contrast = c("compartment", "AXON", "MB"),
#   type = "ashr",
#   parallel = T
# )
#
# # summary(res_C2_IP_COMPARTMENT)
#
# AXON_VS_MB_META <-
# res_C2_IP_COMPARTMENT %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   select(-c(lfcSE, pvalue)) %>%
#   arrange(log2FoldChange) %>%
#   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
#   left_join((AXON_FRACTION_META %>%
#                select("external_gene_name",
#                       expression_energy)),
#             by = c("external_gene_name")) %>%
#   mutate(compartment = ifelse(log2FoldChange > 0, "Axon", "Soma")
#
#          ) %>%
#   left_join(enrichment_status)
#
#          # evidence = ifelse(axon_translated,
#          #                   ifelse(aba_specific & mb_translated,
#          #                          "Axon, Soma and ABA",
#          #                          ifelse(aba_specific,
#          #                                 "Axon and ABA",
#          #                                 ifelse(mb_translated,
#          #                                        "Axon and Soma",
#          #                                        "Axon only"))),
#          #                   ifelse(aba_specific & mb_translated,
#          #                          "Soma and ABA",
#          #                          "Soma only")),
#          # evidence = ifelse(!mb_translated & evidence == "Soma only",
#          #                   "No support", evidence),
#          # evidence = factor(evidence, levels = c("Axon, Soma and ABA",
#          #                                        "Axon and Soma",
#          #                                        "Axon and ABA",
#          #                                        "Axon only",
#          #                                        "Soma and ABA",
#          #                                        "Soma only",
#          #                                        "No support")))
#
#
# # AXON vs MB: Confidence stratified (PLOT) ----
#
# plot_data <- AXON_VS_MB_META %>%
#   filter(padj < 0.01)
#   # filter(mb_translated | axon_translated) # Commented because I am plotting all the evidence; even no support genes
#
# p_axon_vs_mb_enrichment <- plot_data %>%
#   left_join(enrichment_status) %>%
# group_by(compartment, enrichment) %>%
#   tally %>%
#   mutate(prop = n/sum(n),
#          sum = sum(n)) %>%
#   ggplot(aes(x = compartment,
#              y = prop,
#              fill = enrichment)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   scale_y_continuous(labels = scales::percent,
#                      limits = c(0, 1),
#                      expand = expansion(mult = c(0, 0.1))) +
#   labs(x = "TRAP Region",
#        y = "Percentage of Genes",
#        fill = "Axonal Enrichment") +
#   geom_text(aes(x = compartment,
#                 y = 1,
#                 label = paste("Total:", sum)),
#             check_overlap = T,
#             vjust = -0.5) +
#   geom_text(aes(x = compartment,
#                 y = prop,
#                 label = ifelse(prop > 0.05,
#                                scales::percent(prop,
#                                                accuracy = 0.1),
#                                "")),
#             position = position_stack(vjust = 0.5),
#             colour = "white",
#             fontface = 2) +
#   theme(legend.position = "top") +
#   guides(fill = guide_legend(nrow = 2,
#                       byrow = T,
#                       title.position = "top",
#                       title.hjust = 0.5))
#
#
#
# # AXON VS MB: Number of DEGS ----
#
# AXON_VS_MB_META_SIGNIF_SUMMARY <- AXON_VS_MB_META %>%
#   filter(padj < 0.01) %>%
#   filter(mb_translated)
#
# n_AXON_VS_MB_SIGNIF_GENES <- AXON_VS_MB_META %>%
#   filter(padj < 0.01) %>%
#   nrow
#
# # require that genes be soma enriched.
# n_AXON_VS_MB_SIGNIF_ENRICHED_GENES <- AXON_VS_MB_META %>%
#   filter(padj < 0.01) %>%
#   filter(mb_translated) %>%
#   nrow
#
#
# # AXON VS MB: AXON NUMBERS ----
#
# n_AXON_VS_MB_AXON_ENRICHED <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(log2FoldChange > 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# n_AXON_VS_MB_AXON_ENRICHED_LOW_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(enrichment == "Axon only" &
#            log2FoldChange > 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# n_AXON_VS_MB_AXON_ENRICHED_MEDIUM_EVIDENCE_SOMA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(enrichment == "Soma and Axon" &
#            log2FoldChange > 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# # n_AXON_VS_MB_AXON_ENRICHED_MEDIUM_EVIDENCE_ABA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
# #   filter(enrichment == "Axon and ABA" &
# #            log2FoldChange > 0) %>%
# #   pull(external_gene_name) %>%
# #   length
#
# n_AXON_VS_MB_AXON_ENRICHED_SOMA_ONLY <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(enrichment == "Soma only" &
#            log2FoldChange > 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# # n_AXON_VS_MB_AXON_ENRICHED_HIGH_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
# #   filter(enrichment == "Axon, Soma and ABA" &
# #            log2FoldChange > 0) %>%
# #   pull(external_gene_name) %>%
# #   length
#
# # AXON VS MB: SOMA NUMBERS ----
# n_AXON_VS_MB_SOMA_ENRICHED <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# n_AXON_VS_MB_SOMA_ENRICHED_LOW_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon only" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# n_AXON_VS_MB_SOMA_ENRICHED_MEDIUM_EVIDENCE_SOMA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon and Soma" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# n_AXON_VS_MB_SOMA_ENRICHED_MEDIUM_EVIDENCE_ABA <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon and ABA" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length
#
# n_AXON_VS_MB_SOMA_ENRICHED_HIGH_EVIDENCE <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   filter(evidence == "Axon, Soma and ABA" &
#            log2FoldChange < 0) %>%
#   pull(external_gene_name) %>%
#   length
#
#
# # AXON VS MB: Prioritising ----
#
# AXON_VS_MB_META_SIGNIF_PRIORITY <- AXON_VS_MB_META_SIGNIF_SUMMARY %>%
#   left_join(enrichment_status) %>%
#   filter(enrichment %in% c(
#                            "Soma and Axon",
#                          "Soma only") &
#            log2FoldChange > 0) %>%
#   left_join(publications_all_genes) %>%
#   left_join(publications_all_genes_axon) %>%
#   # filter(expression_energy < 1) %>%
#   select(external_gene_name,
#          description,
#          log2FoldChange,
#          padj,
#          expression_energy,
#          enrichment,
#          pub_pd,
#          pub_axon,
#          ensembl_gene_id) %>%
#   distinct() %>%
#   mutate(outcome_trap = ifelse(padj > 0.01,
#                                "Unchanged",
#                                ifelse(log2FoldChange > 0,
#                                       "Axon", "Soma"))) %>%
#   left_join(axon_vs_mb_trap_total) %>%
#   filter(diff_total < 0.5) %>%
#   left_join(AXON_VS_MB_TOTAL_OUTCOMES) %>%
#   # mutate(relationship = ifelse(outcome_trap == outcome_total,
#   #                              "Shared",
#   #                              ifelse(outcome_trap == "Axon",
#   #                                     ifelse(outcome_total == "Unchanged",
#   #                                            "Axon specific",
#   #                                            "Axon opposing"),
#   #                                     ifelse(outcome_total == "Unchanged",
#   #                                            "Soma specific",
#   #                                            "Soma opposing")))) %>%
#   distinct
#
#
# n_AXON_VS_MB_AXON_ENRICHED_PRIORITY <- nrow(AXON_VS_MB_META_SIGNIF_PRIORITY)
#
# n_AXON_VS_MB_AXON_ENRICHED_HIGH_EVIDENCE_TOTAL_OPPOSING <- AXON_VS_MB_META_SIGNIF_PRIORITY %>%
#   filter(relationship == "Axon opposing") %>%
#   nrow
#
# # AXON VS MB: TOTAL outcome comparison with HIGHEST CONFIDENCE TRAP GENES ----
#
# p_axon_vs_mb_total_comparison <- AXON_VS_MB_META_SIGNIF_PRIORITY %>%
#   filter(!is.na(relationship)) %>%
#   mutate(relationship = ifelse(relationship == "Axon opposing",
#                                "Axon\nopposing", relationship)) %>%
# ggplot(aes(x = relationship,
#              fill = relationship)) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(4)[c(3, 1, 4)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   geom_text(aes(label = ..count..,
#                 y = ..count..),
#             vjust = -0.5,
#             stat = "count",
#             fontface = 2) +
#   theme(legend.position = "none") +
#   labs(x = "Relationship to TOTAL Samples",
#        y = "Number of Genes")
#
# n_AXON_VS_MB_TOTAL_BREAKDOWN <- AXON_VS_MB_META_SIGNIF_PRIORITY %>%
#   filter(!is.na(relationship)) %>%
#   group_by(relationship) %>%
#   tally() %>%
#   mutate(prop = n / sum(n)) %>%
#   pull(n)
#
# # AXON VS MB: Example Atg7 ----
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
#
# # AXON VS MB: Example Snx1 ----
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
#
# # AXON VS MB: Example Cul3 ----
# p_axon_vs_mb_cul3 <- plotCounts(dds,
#                                  get_ensembl_gene_id("Cul3"),
#                                  c("cohort", "fraction", "compartment", "region", "age", "genotype"),
#                                  returnData = T,
#                                  normalized = F) %>%
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
#
# # FIGURE: AXON vs MB DEG Volcano - OLD ----------------
# # # plot AXON vs MB DEGs
# # plot_data <- AXON_VS_MB_META_SIMPLE %>%
# #   mutate(location = ifelse(log2FoldChange > 0 &
# #                              sumz_adj < 0.01, "Axonal",
# #                            ifelse(log2FoldChange < 0 &
# #                                     sumz_adj < 0.01, "Somal",
# #                                   "Balanced")))
# # highlight_markers <- plot_data %>%
# #   filter(external_gene_name %in% c("Tubb4a", "Chn1", "Kalrn", "Gch1"))
# #
# # highlight_markers <- plot_data %>%
# #   filter(external_gene_name %in% cp_kegg_pd$external_gene_name)
# #
# # highlight_labels <- plot_data %>%
# #   filter(external_gene_name %in% cp_kegg_pd$external_gene_name &
# #            -log10(sumz_adj) > 30)
# #
# # axonal_label <- plot_data %>%
# #   filter(location == "Axonal") %>%
# #   summarise(n = n())
# # somal_label <- plot_data %>%
# #   filter(location == "Somal") %>%
# #   summarise(n = n())
# # balanced_label <- plot_data %>%
# #   filter(location == "Balanced") %>%
# #   summarise(n = n())
# #
# # p_axon_vs_mb_deg <- ggplot(plot_data,
# #                            aes(x = log2FoldChange,
# #                                y = -log10(sumz_adj))) +
# #   geom_point(aes(fill = location,
# #                  size = baseMean *1e-3),
# #              shape = 21,
# #              colour = "black",
# #              alpha = 0.75) +
# #   scale_fill_manual(values = c(pal_d3()(1),
# #                                "#C1C1C1",
# #                                pal_d3()(3)[3])) +
# #   scale_x_continuous(limits = c(-8, 8),
# #                      oob = squish) +
# #   labs(fill = "Location") +
# #   guides(size = "none") +
# #   theme(legend.position = "top") +
# #   guides(fill = guide_legend(override.aes = list(size=5))) +
# #   geom_vline(xintercept = 0,
# #              linetype = "dotted",
# #              colour = "black") +
# #   geom_point(data = highlight_markers,
# #              aes(size = baseMean *1e-3),
# #              shape = 21,
# #              colour = "black",
# #              fill = "red",
# #              alpha = 0.75) +
# #   geom_label_repel(data = highlight_labels,
# #                    aes(label = external_gene_name),
# #                    # colour = "black",
# #                    arrow = arrow(length = unit(0.02, "npc")),
# #                    box.padding = 1,
# #                    # nudge_x = 1,
# #                    colour = "black",
# #                    # fill = "red",
# #                    segment.colour = "black",
# #                    max.overlaps = 60
# #   ) +
# #   geom_text(data = axonal_label,
# #             aes(x = 5,
# #                 y = 100,
# #                 vjust = -0.5,
# #                 fontface = 2,
# #                 label = paste0(n, " Axonal")),
# #             size = 5) +
# #   geom_text(data = balanced_label,
# #             aes(x = 0,
# #                 y = 125,
# #                 vjust = -0.5,
# #                 fontface = 2,
# #                 label = paste0(n, " Balanced")),
# #             size = 4) +
# #   geom_text(data = somal_label,
# #             aes(x = -5,
# #                 y = 100,
# #                 vjust = -0.5,
# #                 fontface = 2,
# #                 label = paste0(n, " Somal")),
# #             size = 5)
#
#
#
# # AXON vs MB fgsea OLD ---------------
#
# # # get entrez IDs
# # anno_fgsea <- getBM(
# #   attributes = c("ensembl_gene_id",
# #                  "entrezgene_id"),
# #   filters = "ensembl_gene_id",
# #   values = AXON_VS_MB_META$ensembl_gene_id,
# #   mart = ensembl,
# #   uniqueRows = TRUE
# # )
# #
# # # load gmt
# #
# # pathways <- fgsea::gmtPathways("R/objects/msigdb.v7.4.symbols.gmt.txt")
# #
# # pathways <- fgsea::gmtPathways("R/objects/c2.cp.kegg.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.bp.v7.4.symbols.gmt.txt")
# # pathways <- fgsea::gmtPathways("R/objects/c5.go.mf.v7.4.symbols.gmt.txt")
# #
# # ranks <- AXON_VS_MB_META %>%
# #   filter(mb_translated) %>%
# #   mutate(rank = -log10(padj) * log2FoldChange) %>%
# #   inner_join(anno_human) %>%
# #   select(hsapiens_homolog_associated_gene_name,
# #          rank) %>%
# #   filter(hsapiens_homolog_associated_gene_name != "") %>%
# #   distinct(across(-rank), .keep_all = T) %>%
# #   deframe
# # names(ranks)
# # AXON_VS_MB_FGSEA <- fgsea(pathways,
# #                            ranks,
# #                            eps = 0,
# #                            minSize = 10,
# #                            maxSize = 500,
# #                            nPermSimple = 10000,
# #                            BPPARAM = MulticoreParam(),
# #                            nproc = 22)
# #
# # AXON_VS_MB_FGSEA %>% View
# #
# # AXON_VS_MB_FGSEA$leadingEdge[[1156]]
# # ----
# # ----
# # ----
# # AXON vs MB DTU: Counts for IP and TOTAL SETUP ------------
#
# # Keep this code, but start from the readRDS function (the commented code has already been run)
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
#   str_extract(files, "(?<=salmon\\/)[0-9|\\_]+(?=\\/quant.sf)")
#
# files <- files[names(files) %in% metadata$fastq_id] # remove blanks
# files <-
#   files[order(match(names(files), metadata$fastq_id))] #reorder based on metadata
#
# metadata$files_salmon <- files
# metadata$names <- metadata$fastq_id
#
# txdb <-
#   makeTxDbFromGFF(file = "/zfs/analysis/trap/active/snakemake/thesis_snakemake/input/index/references/annotation.gtf",
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
#                                          fraction == "IP" &
#                                          cohort %in% c("C1", "C2") &
#                                            collection != "C1.IPPOOL" &
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
#                                          fraction == "TOTAL" &
#                                            cohort %in% c("C1", "C2") &
#                                            !sample_name %in% c("DS_OLD_WT_MIXED_C2.3",
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
#
# counts_IP <- readRDS("R/objects/counts_IP.rds")
# counts_TOTAL <- readRDS("R/objects/counts_TOTAL.rds")
#
# # AXON VS MB DTU: TOTAL DRIMSeq -----
#
# counts <- counts_TOTAL
#
# # create sample metadata for DRIMSeq
# samps <- metadata %>%
#   filter(sample_name %in% colnames(counts)) %>%
#   select("sample_id" = sample_name,
#          cohort,
#          compartment) %>%
#   distinct
# samps <- as.data.frame(samps)
# samps$cohort <- relevel(samps$cohort, ref = "C2")
# samps$compartment <- relevel(samps$compartment, ref = "MB")
#
# # create drimseq object
# d <- dmDSdata(counts = counts,
#               samples = samps)
#
# # DRIMSeq filter
# n <- nrow(samps)
#
# table(samples(d)$cohort, samples(d)$compartment)
# # n.small <- sum(samps$compartment == "MB" & samps$cohort == "C3") # all the cohort 2 MB samples multiplied by 2
# # n.small <- 44
# n.small <- sum(samps$compartment == "AXON")
# d <- dmFilter(
#   d,
#   min_samps_gene_expr = n * 0.75,
#   min_gene_expr = 10,
#   min_samps_feature_expr = n.small,
#   min_feature_expr = 10,
#   min_samps_feature_prop = n.small,
#   min_feature_prop = 0.1,
# )
#
# saveRDS(d, "R/objects/d_TOTAL_COMPARTMENT.rds")
#
# # code to run in isolate env
# # d <- readRDS("R/objects/d_TOTAL_COMPARTMENT.rds")
# #
# # design_full <-
# #   model.matrix( ~ cohort + compartment, data = DRIMSeq::samples(d))
# # colnames(design_full)
# #
# # set.seed(821196)
# # system.time({
# #   d <-
# #     dmPrecision(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
# #   d <-
# #     dmFit(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
# #   d <-
# #     dmTest(d, coef = "compartmentAXON", BPPARAM = BiocParallel::MulticoreParam())
# # })
# #
# # saveRDS(d,
# #         "R/objects/d_TOTAL_COMPARTMENT_PROCESSED.rds")
#
# d_TOTAL_COMPARTMENT <- readRDS("R/objects/d_TOTAL_COMPARTMENT_PROCESSED.rds")
#
# # get results
# res_TOTAL_COMPARTMENT <- DRIMSeq::results(d_TOTAL_COMPARTMENT)
# res.txp_TOTAL_COMPARTMENT <- DRIMSeq::results(d_TOTAL_COMPARTMENT, level = "feature")
#
# # replace NA values with 1
# no.na <- function(x){
#   ifelse(is.na(x), 1, x)
# }
#
# res_TOTAL_COMPARTMENT$pvalue <- no.na(res_TOTAL_COMPARTMENT$pvalue)
# res.txp_TOTAL_COMPARTMENT$pvalue <- no.na(res.txp_TOTAL_COMPARTMENT$pvalue)
#
# # remove gene id version number
# res_TOTAL_COMPARTMENT$gene_id_simple <- substr(res_TOTAL_COMPARTMENT$gene_id, 1, 18)
#
#
# # AXON VS MB DTU: IP DRIMSeq -----
#
# counts <- counts_IP
#
# # create sample metadata for DRIMSeq
# samps <- metadata %>%
#   filter(sample_name %in% colnames(counts)) %>%
#   select("sample_id" = sample_name,
#          cohort,
#          compartment) %>%
#   distinct
# samps <- as.data.frame(samps)
# samps$cohort <- relevel(samps$cohort, ref = "C2")
# samps$compartment <- relevel(samps$compartment, ref = "MB")
#
# # create drimseq object
# d <- dmDSdata(counts = counts,
#               samples = samps)
#
# methods(class = class(d))
# #
# substr(counts(d)$gene_id, 1, 18) %>% unique
#
# # DRIMSeq filter
# n <- nrow(samps)
#
# table(samples(d)$cohort, samples(d)$compartment)
# # n.small <- sum(samps$compartment == "MB" & samps$cohort == "C3") # all the cohort 2 MB samples multiplied by 2
# # n.small <- 44
# n.small <- sum(samps$compartment == "AXON")
# d <- dmFilter(
#   d,
#   min_samps_gene_expr = n * 0.75,
#   min_gene_expr = 10,
#   min_samps_feature_expr = n.small,
#   min_feature_expr = 10,
#   min_samps_feature_prop = n.small,
#   min_feature_prop = 0.1,
# )
#
# saveRDS(d, "R/objects/d_IP_COMPARTMENT.rds")
#
# # code to run in isolate env
# # d <- readRDS("R/objects/d_IP_COMPARTMENT.rds")
# #
# # design_full <-
# #   model.matrix( ~ cohort + compartment, data = DRIMSeq::samples(d))
# # colnames(design_full)
# #
# # set.seed(821196)
# # system.time({
# #   d <-
# #     dmPrecision(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
# #   d <-
# #     dmFit(d, design = design_full, BPPARAM = BiocParallel::MulticoreParam())
# #   d <-
# #     dmTest(d, coef = "compartmentAXON", BPPARAM = BiocParallel::MulticoreParam())
# # })
# #
# # saveRDS(d,
# #         "R/objects/d_IP_COMPARTMENT_PROCESSED.rds")
#
# d_IP_COMPARTMENT <- readRDS("R/objects/d_IP_COMPARTMENT_PROCESSED.rds")
#
# # get results
# res_IP_COMPARTMENT <- DRIMSeq::results(d_IP_COMPARTMENT)
# res.txp_IP_COMPARTMENT <- DRIMSeq::results(d_IP_COMPARTMENT, level = "feature")
#
# # replace NA values with 1
# no.na <- function(x){
#   ifelse(is.na(x), 1, x)
# }
#
# res_IP_COMPARTMENT$pvalue <- no.na(res_IP_COMPARTMENT$pvalue)
# res.txp_IP_COMPARTMENT$pvalue <- no.na(res.txp_IP_COMPARTMENT$pvalue)
#
# # remove gene id version number
# res_IP_COMPARTMENT$gene_id_simple <- substr(res_IP_COMPARTMENT$gene_id, 1, 18)
#
# # NOT DOING STAGE-R FOR THESIS
# # get annotation
# # anno_drimseq <- get_anno(res_compartment$gene_id_simple, get_human = TRUE)
#
# # # stageR
# # # create a vector of gene pvalues, stripped of version number
# # pScreen_compartment <- res_compartment$pvalue
# # strp <- function(x) substr(x,1,18)
# # names(pScreen_compartment) <- strp(res_compartment$gene_id)
# #
# # # create a vector of transcript pvalues
# # pConfirmation_compartment <- matrix(res.txp_compartment$pvalue, ncol=1)
# # rownames(pConfirmation_compartment) <- strp(res.txp_compartment$feature_id)
# #
# # # create a df with transcript and gene ids (tx2gene)
# # tx2gene <- res.txp_compartment[,c("feature_id", "gene_id")]
# # for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
# #
# # # run stageR at alpha of 0.05
# # stageRObj_compartment <- stageRTx(pScreen=pScreen_compartment, pConfirmation=pConfirmation_compartment,
# #                                   pScreenAdjusted=FALSE, tx2gene=tx2gene)
# # stageRObj_compartment <- stageWiseAdjustment(stageRObj_compartment, method="dtu", alpha=0.05)
# # suppressWarnings({
# #   drim.padj_compartment <- getAdjustedPValues(stageRObj_compartment, order=FALSE,
# #                                               onlySignificantGenes=TRUE)
# # })
# # drim.padj_compartment <- drim.padj_compartment %>%
# #   as_tibble %>%
# #   left_join(anno,
# #             by = c("geneID" = "ensembl_gene_id"))
#
# # post-hoc filtering on the standard deviation in proportions
# # res.txp.filt <- DRIMSeq::results(d, level="feature")
# # smallProportionSD <- function(d, filter=0.1) {
# #   cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
# #   gene.cts <- rowsum(cts, counts(d)$gene_id)
# #   total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
# #   props <- cts/total.cts
# #   propSD <- sqrt(rowVars(props))
# #   propSD < filter
# # }
# # filt <- smallProportionSD(d)
# # res.txp.filt$pvalue[filt] <- 1
# # res.txp.filt$adj_pvalue[filt] <- 1
#
# # res_COMPARTMENT %>%
# #   as_tibble %>%
# #   arrange(adj_pvalue) %>%
# #   filter(adj_pvalue < 0.01) %>%
# #   inner_join(anno, by = c("gene_id_simple" = "ensembl_gene_id"))
#
# # mutate(gwas = hsapiens_homolog_ensembl_gene %in% pd_gwas_genes_human) %>%
# # filter(external_gene_name %in% t_gwas_nalls$`TRAP Top Candidate`) %>%
#
# # AXON VS MB DTU: PLOTTING FUNCTION ----
#
#
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
#                  y = count)) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = compartment),
#                    shape = 21,
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
#
#
#
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
#   geom_violin() +
#   geom_quasirandom(aes(fill = compartment),
#                    shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Compartment",
#        y = "Proportion of Total Expression",
#        title = "Oxr1") +
#   theme(legend.position = "none") +
#   facet_grid(rows = vars(fraction),
#              cols = vars(external_transcript_name)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   panel_border()
#
#
#
#
#
# # AXON VS MB DTU: TOTAL SUMMARY ----
#
# AXON_VS_MB_TOTAL_DTU <- res_TOTAL_COMPARTMENT %>%
#   as_tibble %>%
#   arrange(adj_pvalue) %>%
#   select(ensembl_gene_id = gene_id_simple, adj_pvalue) %>%
#   left_join(anno) %>%
#   left_join(aba_expression_axon) %>%
#   left_join(publications_all_genes) %>%
#   mutate(signif = adj_pvalue < 0.01,
#          mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
#          # axon_abundant = ensembl_gene_id %in% AXON_VS_MB_META_SIGNIF_PRIORITY$ensembl_gene_id,
#          aba_specific = expression_energy < 0.5,
#          aba_specific = replace_na(aba_specific, FALSE),
#          evidence = ifelse(axon_translated,
#                            ifelse(aba_specific & mb_translated,
#                                   "Axon, Soma and ABA",
#                                   ifelse(aba_specific,
#                                          "Axon and ABA",
#                                          ifelse(mb_translated,
#                                                 "Axon and Soma",
#                                                 "Axon only"))),
#                            ifelse(aba_specific & mb_translated,
#                                   "Soma and ABA",
#                                   "Soma only")),
#          evidence = ifelse(!mb_translated & evidence == "Soma only",
#                            "No support", evidence),
#          evidence = factor(evidence, levels = c("Axon, Soma and ABA",
#                                                 "Axon and Soma",
#                                                 "Axon and ABA",
#                                                 "Axon only",
#                                                 "Soma and ABA",
#                                                 "Soma only",
#                                                 "No support"))) %>%
#   filter(signif)
#
# AXON_VS_MB_TOTAL_DTU %>%
#   ggplot(aes(x = evidence)) +
#   geom_bar()
#
# # AXON VS MB DTU: IP SUMMARY ----
#
# AXON_VS_MB_TRAP_DTU <- res_IP_COMPARTMENT %>%
#   as_tibble %>%
#   arrange(adj_pvalue) %>%
#   select(ensembl_gene_id = gene_id_simple, adj_pvalue) %>%
#   left_join(anno) %>%
#   left_join(aba_expression_axon) %>%
#   left_join(publications_all_genes) %>%
#   left_join(publications_all_genes_axon) %>%
#   mutate(signif = adj_pvalue < 0.01,
#          mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
#          # axon_abundant = ensembl_gene_id %in% AXON_VS_MB_META_SIGNIF_PRIORITY$ensembl_gene_id,
#          aba_specific = expression_energy < 0.5,
#          aba_specific = replace_na(aba_specific, FALSE),
#          evidence = ifelse(axon_translated,
#                            ifelse(aba_specific & mb_translated,
#                                   "Axon, Soma and ABA",
#                                   ifelse(aba_specific,
#                                          "Axon and ABA",
#                                          ifelse(mb_translated,
#                                                 "Axon and Soma",
#                                                 "Axon only"))),
#                            ifelse(aba_specific & mb_translated,
#                                   "Soma and ABA",
#                                   "Soma only")),
#          evidence = ifelse(!mb_translated & evidence == "Soma only",
#                            "No support", evidence),
#          evidence = factor(evidence, levels = c("Axon, Soma and ABA",
#                                                 "Axon and Soma",
#                                                 "Axon and ABA",
#                                                 "Axon only",
#                                                 "Soma and ABA",
#                                                 "Soma only",
#                                                 "No support")),
#          specific = !ensembl_gene_id %in% AXON_VS_MB_TOTAL_DTU$ensembl_gene_id)
#
# # AXON VS MB DTU: NUMBERS ----
#
# # sort this out: how many genes did we really consider to begin with?
# n_AXON_VS_MB_DTU_CONSIDERED <- length(unique(counts_IP$gene_id))
# n_AXON_VS_MB_DTU_CONSIDERED <- nrow(AXON_VS_MB_META)
#
# n_AXON_VS_MB_DTU_RETAINED <- nrow(AXON_VS_MB_TRAP_DTU)
#
# n_AXON_VS_MB_DTU_FILTER_PERCENT <- 100*signif(n_AXON_VS_MB_DTU_RETAINED/n_AXON_VS_MB_DTU_CONSIDERED, 2)
#
# n_AXON_VS_MB_DTU_ANYCONFIDENCE <- AXON_VS_MB_TRAP_DTU %>%
#   filter(adj_pvalue < 0.01) %>%
#   nrow
#
# # AXON VS MB DTU: IP vs TOTAL PLOT ----
#
# # Filtered for genes enriched in either axon or soma
# p_axon_vs_mb_dtu_specificity <- AXON_VS_MB_TRAP_DTU %>%
#   filter(signif) %>%
#   filter(mb_translated) %>%
#   left_join(enrichment_status) %>%
#   select(-evidence) %>%
#   rename(evidence = enrichment) %>%
#   # mutate(evidence = ifelse(evidence == "Axon, Soma and ABA",
#   #                          "Axon, Soma\nand ABA",
#   #                          ifelse(evidence == "Axon and Soma",
#   #                                 "Axon and\nSoma",
#   #                                 ifelse(evidence == "Soma and ABA",
#   #                                        "Soma and\nABA",
#   #                                        "Soma only")))) %>%
#   group_by(evidence) %>%
#   mutate(n = dplyr::n()) %>%
#   mutate(specific = ifelse(ensembl_gene_id %in% {AXON_VS_MB_TOTAL_DTU %>%
#       filter(signif) %>%
#       pull(ensembl_gene_id)},
#                            "Shared", "Specific")) %>%
#   ggplot(aes(x = evidence,
#              fill = specific)) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   geom_text(aes(label = n,
#                 y = n),
#             check_overlap = T,
#             vjust = -0.5,
#             fontface = 2,
#             size = 4) +
#   stat_count(geom = "text", colour = "white", size = 3.5,
#              fontface = 2,
#              aes(label = ifelse(..count.. > 50, ..count.., "")),
#              position=position_stack(vjust=0.5)) +
#   theme(legend.position = "top",
#         axis.text.x = element_text(size = 9)) +
#   # guides(fill = guide_legend(nrow = 4,
#   #                            byrow=TRUE)) +
#   labs(x = "Enrichment Category",
#        y = "Number of Genes with DTU",
#        fill = "Specificity")
#
#
# # AXON VS MB DTU: Isoforms BAR CHART -----------------------
#
# # counts <- counts_IP
# #
# # # plot number of genes DTU by number of isoforms
# # n_tx <- tibble(gene_id = unique(counts$gene_id),
# #                n_tx = sapply(unique(counts$gene_id), function(x) sum(counts(d_IP_COMPARTMENT)$gene_id == x))) %>%
# #   mutate(signif = ifelse(gene_id %in% res_IP_COMPARTMENT$gene_id[res_IP_COMPARTMENT$adj_pvalue < 0.01], "DTU", "NS"))
# # n_tx$n_tx[table(counts$gene_id)[match(unique(counts$gene_id), names(table(counts$gene_id)))] == 1] <- 1
# # n_tx$n_tx[n_tx$n_tx == 0] <- "Filtered"
# # n_tx$n_tx <- factor(n_tx$n_tx, levels = c("Filtered", seq(1, 6, 1)))
#
# n_tx <- readRDS("R/objects/n_tx.rds")
#
# # Filtered for genes enriched in either axon or soma
# p_axon_vx_mb_dtu_bar <- n_tx %>%
#   filter(substr(gene_id, 1, 18) %in% unique(MB_FRACTION_TRANSLATED_GENES,
#                                             AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)) %>%
#   filter(!is.na(n_tx) &
#            n_tx != "6") %>%
#   ggplot(aes(x = n_tx,
#                                    fill = signif)) +
#   geom_bar(colour = "black") +
#   stat_count(geom = "text", colour = "white", size = 3.5, fontface = 2,
#              aes(label = ifelse(..count.. < 100, "", ..count..)),
#              position=position_stack(vjust=0.5)) +
#   scale_fill_manual(values = pal_d3()(2)[c(2, 1)]) +
#   labs(x = "Number of Isoforms",
#        y = "Number of Genes",
#        fill = "Outcome") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#   theme(legend.position = c(0.75, 0.75))
#
#
# # ----
# # ----
# # ----
# # FIGURE B: DEG vs DTU - OLD DO NOT USE -----------------------------
# # DEG vs DTU plot
#
# # plot_data <- AXON_VS_MB_META_SIMPLE %>%
# #   filter(ensembl_gene_id %in% substr(res_COMPARTMENT$gene_id, 1, 18)) %>%
# #   select(ensembl_gene_id, external_gene_name, sumz_adj) %>%
# #   left_join((res_COMPARTMENT %>%
# #                mutate(gene_id = substr(gene_id, 1, 18)) %>%
# #                select(gene_id, adj_pvalue) %>%
# #                distinct()),
# #             by = c("ensembl_gene_id" = "gene_id")) %>%
# #   mutate(signif = ifelse(sumz_adj < 0.01 &
# #                            adj_pvalue < 0.01,
# #                          "DEG & DTU",
# #                          ifelse(sumz_adj < 0.01,
# #                                 "DEG",
# #                                 ifelse(adj_pvalue < 0.01,
# #                                        "DTU",
# #                                        "NS")))) %>%
# #   filter(!is.na(signif))
# #
# # p_axon_vs_mb_deg_dtu <- ggplot(plot_data,
# #                                aes(x = -log10(sumz_adj),
# #                                    y = -log10(adj_pvalue))) +
# #   geom_point(colour = "black",
# #              shape = 21,
# #              size = 2,
# #              alpha = 1,
# #              aes(fill = signif)) +
# #   # scale_x_continuous(limits = c(0, 50),
# #   #                    oob = squish) +
# #   # scale_y_continuous(limits = c(0, 50),
# #   #                    oob = squish) +
# #   scale_fill_d3() +
# #   scale_x_log10(limits = c(0.5, 100)) +
# #   scale_y_log10(limits = c(0.5, 100)) +
# #   coord_fixed() +
# #   geom_vline(xintercept = 2,
# #              size = 1,
# #              linetype = "dashed",
# #              colour = "black") +
# #   geom_hline(yintercept = 2,
# #              size = 1,
# #              linetype = "dashed",
# #              colour = "black") +
# #   theme(legend.position = "top",
# #         legend.title = element_blank()) +
# #   labs(x = "DEG -log10 P",
# #        y = "DTU -log10 P")  +
# #   guides(fill = guide_legend(nrow = 2,
# #                              byrow=TRUE))
#
# # FIGURE C: DTU example - OLD ------------------------------
# # plot Snap25 isoform switching
# # p_axon_vs_mb_dtu_example <- counts(d_IP_COMPARTMENT) %>%
# #   as_tibble %>%
# #   rename_with(~ str_replace(.x, "TOTAL1.5", "TOTAL1-5")) %>%
# #   filter(str_detect(gene_id, "ENSMUSG00000027273")) %>%
# #   mutate(feature_id = paste0("Isoform ", seq(1, nrow(.), 1))) %>%
# #   pivot_longer(-c(gene_id, feature_id),
# #                names_to = "sample_name",
# #                values_to = "count") %>%
# #   group_by(sample_name) %>%
# #   filter(median(count) > 0) %>%
# #   mutate(count = count / sum(count)) %>%
# #   inner_join(colData(dds), copy = TRUE) %>%
# #   mutate(COMPARTMENT = ifelse(compartment == "MB", "CELL BODY", "AXON")) %>%
# #   ggplot(aes(x = COMPARTMENT,
# #              y = count,
# #              fill = COMPARTMENT)) +
# #   geom_violin() +
# #   geom_quasirandom(shape = 21,
# #                    size = 3,
# #                    colour = "black") +
# #   scale_fill_d3() +
# #   labs(y = "Proportion of Counts",
# #        title = "Snap25") +
# #   # scale_x_discrete(labels = c("Isoform 1", "Isoform 2", "Isoform 3")) +
# #   theme(legend.position = "none",
# #         axis.title.x = element_blank()) +
# #   facet_wrap(vars(feature_id),
# #              nrow = 1,
# #              strip.position = "top") +
# #   panel_border()
# # ggsave(filename = "output/images/dcun1d1_COMPARTMENT_switching.png")
#
# # plotProportions(d, res$gene_id[which(res$gene_id_simple == "ENSMUSG00000027708")], "compartment")
#
#
# # AXON vs MB DEU: DEXSeq - OLD ---------
# #
# # # Keep this code, but don't run it: Load from saved RDS
# #
# # files <-
# #   list.files(
# #     "/zfs/analysis/trap/active/testing_env/dexseq/output",
# #     pattern = ".txt",
# #     recursive =  T,
# #     full.names = T
# #   )
# #
# # all(file.exists(files)) # everything exists
# # names(files) <-
# #   str_extract(files, "(?<=output\\/)[0-9]+(?=.txt)")
# #
# # files <- files[names(files) %in% colData(dds)$sample_code]
# # names(files[order(match(names(files), colData(dds)$sample_code))]) == colData(dds)$sample_code
# # files <- files[order(match(names(files), colData(dds)$sample_code))]
# #
# # dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(
# #   files,
# #   sampleData = as.data.frame(colData(dds)),
# #   design = ~ sample + exon + compartment:exon,
# #   flattenedfile = list.files("/zfs/analysis/trap/active/testing_env/dexseq",
# #                              pattern = "dexseq.gtf",
# #                              full.names = T)
# # )
# #
# # # filter dxd for axon enriched genes (that must be mb translated by definition)
# # dxd <- dxd[substr(rownames(dxd), 1, 18) %in% AXON_FRACTION_ENRICHED_GENES,]
# # # filter dxd for cohort 2 IP only
# # dxd <- dxd[,colData(dxd)$fraction == "IP" &
# #              colData(dxd)$cohort == "C2"]
# # dim(dxd)
# # saveRDS(dxd, "R/objects/dxd.rds")
# #
# # # dxd <- readRDS("R/objects/dxd.rds")
# # # dxd <- estimateSizeFactors(dxd)
# # # dxd <- estimateDispersions(dxd,
# # #                            BPPARAM = MulticoreParam())
# # #
# # # dxd <- DEXSeq::testForDEU(dxd,
# # #                   BPPARAM = MulticoreParam())
# # #
# # # dxd = DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar="compartment",
# # #                                        BPPARAM = MulticoreParam())
# # # saveRDS(dxd, "R/objects/dxd_COMPARTMENT.rds")
# #
# # dxd_COMPARTMENT <- readRDS("R/objects/dxd_COMPARTMENT.rds")
# #
# # dxr1_COMPARTMENT = DEXSeqResults( dxd_COMPARTMENT )
# # # dxr1 %>% as.data.frame %>% View
# #
# # # number of genes considered for DEU
# # dxr1_COMPARTMENT$groupID %>% unique %>% length
# #
# # # number of genes with signif DEU at alpha of 0.01
# # dxr1_COMPARTMENT %>% as_tibble %>% filter(padj < 0.01) %>% pull(groupID) %>% unique %>% length
# #
# # # FIGURE A: DEU Venn - OLD ----------------
# #
# # # DEU DTU intersection at P = 0.05
# # AXON_VS_MB_DEU_GENES <- dxr1_COMPARTMENT %>%
# #   as_tibble %>%
# #   filter(padj < 0.05) %>%
# #   pull(groupID) %>%
# #   unique %>%
# #   substr(., 1, 18)
# #
# # AXON_VS_MB_DTU_GENES <- res_COMPARTMENT %>%
# #   as_tibble %>%
# #   filter(adj_pvalue < 0.05) %>%
# #   pull(gene_id) %>%
# #   unique %>%
# #   substr(., 1, 18)
# #
# # venn <- tibble(value = AXON_FRACTION_ENRICHED_GENES) %>%
# #   mutate(`DEG` = value %in% c(AXON_VS_MB_AXON_ENRICHED_GENES,
# #                               AXON_VS_MB_MB_ENRICHED_GENES),
# #          `DTU` = value %in% AXON_VS_MB_DTU_GENES,
# #          DEU = value %in% AXON_VS_MB_DEU_GENES)
# #
# # p_axon_vs_mb_deg_dtu_deu_venn <- ggvenn(venn,
# #                                         fill_color = pal_d3()(4),
# #                                         fill_alpha = 0.55,
# #                                         text_size = 2.5)
# #
# # # FIGURE B: DEU example - OLD ---------------
# #
# # # Keep the commented code! I have to edit the plotDEXSeq function to remove the
# # # ensembl gene title from the plot
# # # so have just saved the plot as an RDS to load for rendering.
# #
# # # dxr1 %>%
# # #   as_tibble(rownames = "exonID") %>%
# # #   # group_by(groupID) %>%
# # #   # slice_tail(n = 1) %>%
# # #   # ungroup() %>%
# # #   filter(padj < 0.01) %>%
# # #   arrange(padj) %>%
# # #   select(groupID,
# # #          exonID,
# # #          exonBaseMean,
# # #          padj,
# # #          log2fold_AXON_MB,
# # #          genomicData.seqnames,
# # #          genomicData.start,
# # #          genomicData.end,
# # #          genomicData.strand,
# # #          genomicData.width) %>%
# # #   mutate(ensembl_gene_id = substr(groupID, 1, 18)) %>%
# # #   left_join(anno) %>% View
# #
# # # dxr1 %>% as_tibble() %>%
# # #   mutate(ensembl_gene_id = strp(groupID)) %>%
# # #   group_by(groupID, ensembl_gene_id) %>%
# # #   summarise(baseMean = mean(exonBaseMean),
# # #             padj = min(padj, na.rm = T)) %>%
# # #   filter(padj < 0.1) %>%
# # #   arrange(padj) %>%
# # #   left_join(anno) %>% View
# #
# # ## edit the plotDEXSeq function
# # ## a more robust way is described here: https://stackoverflow.com/questions/24331690/modify-package-function
# # # trace(plotDEXSeq, edit = TRUE)
# #
# # # p_dexseq <- base2grob( ~plotDEXSeq(dxr1_COMPARTMENT,
# # #                                    "ENSMUSG00000027273.13",
# # #                                    fitExpToVar = "compartment",
# # #                                    legend = T,
# # #                                    displayTranscripts = F,
# # #                                    names = T,
# # #                                    color = pal_d3()(2),
# # #                                    # splicing = T,
# # #                                    transcriptDb = txdb,
# # #                                    sub = NULL))
# #
# # # saveRDS(p_dexseq, "R/objects/p_dexseq.rds")
# #
# # p_dexseq <- readRDS("R/objects/p_dexseq.rds")
# #
# # # plot exon coverage - not very useful for comparing because there are different y axes
# # # for MB and AXON
# # # set strict chromosome naming for loading BAMs
# # # options(ucscChromosomeNames = T)
# # # # load BAMs
# # # mb_track <- AlignmentsTrack(file.path("/zfs/analysis/trap/active/testing_env/stringtie/bam/281193.bam"),
# # #                             isPaired = T)
# # # names(mb_track) <- "Cell Body"
# # # axon_track <- AlignmentsTrack(file.path("/zfs/analysis/trap/active/testing_env/stringtie/bam/282181.bam"),
# # #                               isPaired = T)
# # # names(axon_track) <- "Axon"
# # #
# # # # set non-strict chromosome naming for loading gene region track
# # # options(ucscChromosomeNames = F)
# # #
# # # # plot exon coverage function
# # # plot_exon_cov <- function(exon,
# # #                           buffer_l = 100,
# # #                           buffer_r = 100){
# # #
# # #   # set non-strict chromosome naming for loading gene region track
# # #   options(ucscChromosomeNames = F)
# # #
# # #   coords <- as_tibble(dxr1,
# # #                       rownames = "exon_id") %>%
# # #     filter(exon_id == exon) %>%
# # #     select(c(chrom = genomicData.seqnames,
# # #              start = genomicData.start,
# # #              end = genomicData.end))
# # #
# # #   bmt <- GeneRegionTrack(txdb,
# # #                          chromosome = coords$chrom,
# # #                          start = coords$start - 100,
# # #                          end = coords$end + 100)
# # #   displayPars(bmt)$size <- 15 # enlarge the annotation track
# # #
# # #   plotTracks(c(bmt, mb_track, axon_track),
# # #              from = coords$start - buffer_l,
# # #              to = coords$end + buffer_r,
# # #              chromosome = coords$chrom,
# # #              type = "coverage")
# # # }
# # #
# # # plot_exon_cov("ENSMUSG00000027273.13:E009", 200, 200)
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # ----
# # ----
# # ----
# # AXON DS VS VS Subset common genes ----
# dds_C2_AXON_IP_REGION_FILTER <- filter_genes(dds_C2_AXON_IP,
#                                              grouping = c("region", "gene_id"))
#
# # AXON DS VS VS TOTAL DESEQ2 ----
#
# dds_C2_AXON_TOTAL_REGION <- dds_C2_AXON_TOTAL[rownames(dds_C2_AXON_TOTAL) %in% dds_C2_AXON_IP_REGION_FILTER,]
#
# dds_C2_AXON_TOTAL_REGION@design <- ~ region
#
# colData(dds_C2_AXON_TOTAL_REGION) <- droplevels(colData(dds_C2_AXON_TOTAL_REGION))
#
# colData(dds_C2_AXON_TOTAL_REGION)$region <- relevel(colData(dds_C2_AXON_TOTAL_REGION)$region, ref = "VS")
#
# dds_C2_AXON_TOTAL_REGION <- DESeq(
#   dds_C2_AXON_TOTAL_REGION,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_C2_AXON_TOTAL_REGION <- DESeq2::results(
#   dds_C2_AXON_TOTAL_REGION,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C2_AXON_TOTAL_REGION <- lfcShrink(
#   dds_C2_AXON_TOTAL_REGION,
#   res = res_C2_AXON_TOTAL_REGION,
#   contrast = c("region", "DS", "VS"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_C2_AXON_TOTAL_REGION)
#
# DS_VS_VS_TOTAL_META <- res_C2_AXON_TOTAL_REGION %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   mutate(outcome = ifelse(padj < 0.1,
#                           ifelse(log2FoldChange > 0, "Dorsal", "Ventral"),
#                           "Unchanged"))
#
# DS_VS_VS_TOTAL_OUTCOMES <- DS_VS_VS_TOTAL_META %>%
#   select(ensembl_gene_id,
#          outcome)
#
# # AXON DS VS VS TRAP DESEQ2 ----
#
# dds_C2_AXON_IP_REGION <- dds_C2_AXON_IP[dds_C2_AXON_IP_REGION_FILTER,]
#
# dds_C2_AXON_IP_REGION@design <- ~ collection + region
#
# colData(dds_C2_AXON_IP_REGION) <- droplevels(colData(dds_C2_AXON_IP_REGION))
#
# colData(dds_C2_AXON_IP_REGION)$region <- relevel(colData(dds_C2_AXON_IP_REGION)$region, ref = "VS")
#
# dds_C2_AXON_IP_REGION <- DESeq(
#   dds_C2_AXON_IP_REGION,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_C2_AXON_IP_REGION <- DESeq2::results(
#   dds_C2_AXON_IP_REGION,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
#
# res_C2_AXON_IP_REGION <- lfcShrink(
#   dds_C2_AXON_IP_REGION,
#   res = res_C2_AXON_IP_REGION,
#   contrast = c("region", "DS", "VS"),
#   type = "ashr",
#   parallel = T
# )
#
# DS_VS_VS_META <- res_C2_AXON_IP_REGION %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES) %>%
#   left_join(enrichment_status) %>%
#   mutate(outcome = ifelse(padj < 0.01,
#                           ifelse(log2FoldChange > 0, "Dorsal", "Ventral"),
#                           "Unchanged")) %>%
#   left_join(DS_VS_VS_TOTAL_OUTCOMES,
#             by = "ensembl_gene_id",
#             suffix = c("_TRAP", "_TOTAL")) %>%
#   mutate(specific = ifelse(outcome_TRAP == outcome_TOTAL,
#                            "Shared",
#                            ifelse(outcome_TRAP == "Dorsal",
#                                   ifelse(outcome_TOTAL == "Unchanged",
#                                          "Specific",
#                                          "Specific Opposing"),
#                                   ifelse(outcome_TOTAL == "Unchanged",
#                                          "Specific",
#                                          "Specific Opposing"))),
#          specific = replace_na(specific, "Specific")) %>%
#   left_join(publications_all_genes_axon) %>%
#   left_join(publications_all_genes) %>%
#   dplyr::select(ensembl_gene_id,
#                 external_gene_name,
#                 description,
#                 baseMean,
#                 log2FoldChange,
#                 padj,
#                 mb_translated,
#                 axon_translated,
#                 # mb_age_relationship,
#                 pub_pd,
#                 pub_axon,
#                 # expression_energy,
#                 # aba_specific,
#                 # evidence,
#                 specific,
#                 outcome_TRAP,
#                 outcome_TOTAL)
#
#   # left_join((AXON_FRACTION_META %>%
#   #              select("external_gene_name",
#   #                     expression_energy)),
#   #           by = c("external_gene_name")) %>%
#   # mutate(aba_specific = expression_energy < 0.5,
#   #        aba_specific = replace_na(aba_specific, FALSE),
#   #        # compartment = ifelse(log2FoldChange > 0, "DS", "VS"),
#   #        evidence = ifelse(axon_translated,
#   #                          ifelse(aba_specific & mb_translated,
#   #                                 "Axon, Soma and ABA",
#   #                                 ifelse(aba_specific,
#   #                                        "Axon and ABA",
#   #                                        ifelse(mb_translated,
#   #                                               "Axon and Soma",
#   #                                               "Axon only"))),
#   #                          ifelse(aba_specific & mb_translated,
#   #                                 "Soma and ABA",
#   #                                 "Soma only")),
#   #        evidence = ifelse(!mb_translated & evidence == "Soma only",
#   #                          "No support", evidence),
#   #        evidence = factor(evidence, levels = c("Axon, Soma and ABA",
#   #                                               "Axon and Soma",
#   #                                               "Axon and ABA",
#   #                                               "Axon only",
#   #                                               "Soma and ABA",
#   #                                               "Soma only",
#   #                                               "No support")),
#
#
#
# # AXON DS VS VS: EVIDENCE STRATIFIED (2 PLOTS) ----
#
# # FILTER FOR MB TRANSLATED GENES
#
# DS_VS_VS_META_SIGNIF_SUMMARY <- DS_VS_VS_META %>%
#   filter(padj < 0.01) %>%
#   mutate(region = outcome_TRAP)
#
# n_DS_VS_VS_DE_NOFILTER <- nrow(DS_VS_VS_META_SIGNIF_SUMMARY)
#
# # JUST THE NUMBER OF GENES FIRST
# p_ds_vs_vs_deg <- DS_VS_VS_META_SIGNIF_SUMMARY %>%
#   group_by(region) %>%
#   mutate(n = dplyr::n()) %>%
#   ggplot(aes(x = region,
#              fill = region
#   )) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2, 4)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   labs(x = "TRAP Region",
#        y = "Number of Genes") +
#   geom_text(aes(label = paste("Total:", n),
#                 y = n),
#             check_overlap = T,
#             vjust = -0.5,
#             fontface = 2,
#             size = 4) +
#   theme(legend.position = "none") +
#   guides(fill = guide_legend(nrow = 2,
#                              byrow=TRUE))
#
# # THEN EVIDENCE BREAKDOWN
# p_ds_vs_vs_deg_evidence <- DS_VS_VS_META_SIGNIF_SUMMARY %>%
#   left_join(enrichment_status) %>%
#   select(everything(),
#          evidence = enrichment) %>%
#   group_by(region, evidence) %>%
#   tally %>%
#   mutate(prop = n/sum(n),
#          sum = sum(n)) %>%
#   ggplot(aes(x = region,
#              y = prop,
#              fill = evidence)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   scale_y_continuous(labels = scales::percent,
#                      expand = expansion(mult = c(0, 0.1))) +
#   labs(x = "TRAP Region",
#        y = "Percentage of Genes",
#        fill = "Axonal Enrichment") +
#   geom_text(aes(x = region,
#                 y = 1.05,
#                 label = paste("Total:", sum)),
#             check_overlap = T) +
#   geom_text(aes(x = region,
#                 y = prop,
#                 label = ifelse(prop > 0.05,
#                                scales::percent(prop,
#                                                accuracy = 0.1),
#                                "")),
#             position = position_stack(vjust = 0.5),
#             colour = "white",
#             fontface = 2)
#
# # AXON DS VS VS: Prioritising ----
#
# # Selecting genes that are enriched in soma, and up in dorsal striatum.
# # There are still interesting genes that are soma enriched, not axon enriched
#
# DS_VS_VS_META_SIGNIF_PRIORITY <- DS_VS_VS_META_SIGNIF_SUMMARY %>%
#   filter(mb_translated) %>%
#   left_join(publications_all_genes_axon) %>%
#   left_join(publications_all_genes) %>%
#   # filter(expression_energy < 1) %>%
#   select(ensembl_gene_id,
#          external_gene_name,
#          description,
#          log2FoldChange,
#          padj,
#          axon_translated,
#          # expression_energy,
#          # evidence,
#          pub_pd,
#          pub_axon) %>%
#   distinct() %>%
#   filter(log2FoldChange > 0)
#
# n_DS_VS_VS_DE_RETAINED <- nrow(DS_VS_VS_META_SIGNIF_PRIORITY)
#
# # AXON DS VS VS: TOTAL outcomes plot ----
#
# DS_VS_VS_TOTAL_OUTCOMES <- DS_VS_VS_TOTAL_META %>%
#   mutate(outcome = ifelse(pvalue > 0.1,
#                           "Unchanged",
#                           ifelse(log2FoldChange > 0,
#                                  "Dorsal", "Ventral"))) %>%
#   select(external_gene_name,
#          outcome_total = outcome,
#          pvalue_TOTAL = pvalue)
#
# p_ds_vs_vs_total_outcome <- DS_VS_VS_TOTAL_OUTCOMES %>%
#   filter(outcome_total != "Unchanged") %>%
#   ggplot(aes(x = outcome_total,
#              fill = outcome_total)) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(4)[c(3, 1, 2, 4)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   labs(x = "TOTAL Region",
#        y = "Number of Genes") +
#   geom_text(aes(label = paste("Total:", ..count..),
#                 y = ..count..),
#             stat = "count",
#             check_overlap = T,
#             vjust = -0.5,
#             fontface = 2,
#             size = 4) +
#   theme(legend.position = "none")
#
# # AXON DS VS VS: TOTAL outcome comparison with HIGHEST CONFIDENCE TRAP GENES ----
#
# plot_data <- DS_VS_VS_META_SIGNIF_PRIORITY %>%
#   mutate(outcome_trap = ifelse(padj > 0.01,
#                                "Unchanged",
#                                ifelse(log2FoldChange > 0,
#                                       "Dorsal", "Ventral"))) %>%
#   left_join(DS_VS_VS_TOTAL_OUTCOMES) %>%
#   mutate(relationship = ifelse(outcome_trap == outcome_total,
#                                "Shared",
#                                ifelse(outcome_trap == "Dorsal",
#                                       ifelse(outcome_total == "Unchanged",
#                                              "Dorsal specific",
#                                              "Dorsal opposing"),
#                                       ifelse(outcome_total == "Unchanged",
#                                              "Ventral specific",
#                                              "Ventral opposing"))))
#
# DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON <- plot_data
#
# p_ds_vs_vs_total_comparison <- plot_data %>%
#   filter(!is.na(relationship)) %>%
#   ggplot(aes(x = relationship,
#              fill = relationship)) +
#   geom_bar(colour = "black") +
#   # scale_fill_manual(values = pal_d3()(4)[c(3, 1, 4)]) +
#   scale_fill_d3() +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   geom_text(aes(label = ..count..,
#                 y = ..count..),
#             vjust = -0.5,
#             stat = "count",
#             fontface = 2) +
#   theme(legend.position = "none") +
#   labs(x = "Relationship to TOTAL Samples",
#        y = "Number of Genes")
#
# n_DS_VS_VS_DE_PRIORITY_SPECIFIC <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific") %>%
#   nrow
#
# # n_DS_VS_VS_DE_PRIORITY_SPECIFIC_ABACONF <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
# #   filter(relationship == "Dorsal specific" &
# #            expression_energy < 0.5) %>%
# #   nrow
# #
# # n_DS_VS_VS_DE_PRIORITY_SPECIFIC_AXON_ENRICHED <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
# #   filter(relationship == "Dorsal specific" &
# #            axon_translated) %>%
# #   nrow
# #
# # n_DS_VS_VS_DE_PRIORITY_SPECIFIC_ABACONF_AXON_ENRICHED <- DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
# #   filter(relationship == "Dorsal specific" &
# #            expression_energy < 0.5 &
# #            axon_translated) %>%
# #   nrow
#
#
# # DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
# #   filter(relationship == "Dorsal specific" &
# #            expression_energy < 0.5) %>% View
#
# DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
#   filter(relationship == "Dorsal specific")
#
# # No GWAS targets
# # DS_VS_VS_META_SIGNIF_PRIORITY_HIGHEST_CONFIDENCE_TOTAL_COMPARISON %>%
# #   filter(relationship == "Dorsal specific") %>%
# #   mutate(gwas = external_gene_name %in% GWAS_GENES_BROAD) %>%
# #   filter(gwas)
#
# # plotCounts(dds_C2_AXON_IP, get_ensembl_gene_id("Chrna4"), "region")
#
# # AXON DS VS VS: Examples ----
# p_ds_vs_vs_ptgs2 <- view_counts(dds_C2,
#                                 "Ptgs2") %>%
#   left_join(colData(dds_C2), copy = T) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   filter(count != 0) %>%
#   filter(fraction == "IP") %>%
#   ggplot(aes(x = region,
#              y = log2(count+1))) +
#   geom_quasirandom(aes(fill = region),
#                    shape = 21,
#                    size = 2,
#                    colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Fraction and Region",
#        y = "Log2 Expression Count",
#        title = "Ptgs2 (COX-2)") +
#   panel_border() +
#   scale_y_continuous(limits = c(0, NA)) +
#   theme(legend.position = "none")
#
#
#
#
# p_ds_vs_vs_oxr1 <- view_counts(dds,
#             "Oxr1") %>%
#   left_join(colData(dds), copy = T) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   filter(count != 0) %>%
#   # filter(fraction == "IP") %>%
#   ggplot(aes(x = region,
#              y = log2(count+1))) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = region),
#                    shape = 21,
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
#
# p_ds_vs_vs_micu1 <- view_counts(dds,
#             "Prkaa2") %>%
#   inner_join(colData(dds), copy = T) %>%
#   mutate(fraction = ifelse(fraction == "IP", "TRAP", "TOTAL")) %>%
#   group_by(cohort, fraction, region) %>%
#   filter(!outlier_iqr(count)) %>%
#   ggplot(aes(x = region,
#              y = log2(count+1))) +
#   geom_violin() +
#   geom_quasirandom(aes(fill = region),
#                    shape = 21,
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
#
#
# # ----
# # ----
# # ----
# # FIGURE: DS vs VS DEGs - OLD ----
#
# # # plot AXON vs MB DEGs
# # plot_data <- DS_VS_VS_META_SIMPLE %>%
# #   mutate(location = ifelse(log2FoldChange > 0 &
# #                              sumz_adj < 0.01, "Dorsal",
# #                            ifelse(log2FoldChange < 0 &
# #                                     sumz_adj < 0.01, "Ventral",
# #                                   "Balanced")),
# #          location = factor(location, levels = c("Balanced", "Dorsal", "Ventral")))
# # plot_data$location <- droplevels(plot_data$location)
# #
# # highlight_markers <- plot_data %>%
# #   filter(external_gene_name %in% cp_kegg_pd$external_gene_name)
# #
# # highlight_labels <- plot_data %>%
# #   filter(external_gene_name %in% cp_kegg_pd$external_gene_name &
# #            -log10(sumz_adj) > 10)
# #
# # dorsal_label <- plot_data %>%
# #   filter(location == "Dorsal") %>%
# #   summarise(n = n())
# # ventral_label <- plot_data %>%
# #   filter(location == "Ventral") %>%
# #   summarise(n = n())
# # balanced_label <- plot_data %>%
# #   filter(location == "Balanced") %>%
# #   summarise(n = n())
# #
# # p_ds_vs_vs_deg <- plot_data %>%
# #   filter(!is.na(location)) %>%
# #   ggplot(aes(x = baseMean,
# #              y = log2FoldChange)) +
# #   geom_point(
# #     aes(fill = location,
# #         size = -log10(sumz_adj)),
# #     shape = 21,
# #     colour = "black",
# #     alpha = 1
# #   ) +
# #   scale_fill_manual(values = c("#C1C1C1",
# #                                pal_d3()(10)[c(1, 2)])) +
# #   scale_x_log10() +
# #   scale_y_continuous(limits = c(-2.5, 2.5),
# #                      oob = squish) +
# #   labs(x = "Mean Counts",
# #        y = "Log2 Fold Change",
# #        fill = "Location") +
# #   guides(size = "none") +
# #   theme(legend.position = "top") +
# #   guides(fill = guide_legend(override.aes = list(size = 5))) +
# #   geom_hline(yintercept = 0,
# #              linetype = "dotted",
# #              colour = "black")  +
# #   geom_text(data = dorsal_label,
# #             aes(x = 1e5,
# #                 y = 2,
# #                 vjust = -0.5,
# #                 fontface = 2,
# #                 label = paste0(n, " Dorsal")),
# #             size = 5) +
# #   geom_text(data = balanced_label,
# #             aes(x = 1,
# #                 y = 0.3,
# #                 vjust = 0,
# #                 fontface = 2,
# #                 label = paste0(n, "\nBalanced")),
# #             size = 4) +
# #   geom_text(data = ventral_label,
# #             aes(x = 1e5,
# #                 y = -2,
# #                 vjust = -0.5,
# #                 fontface = 2,
# #                 label = paste0(n, " Ventral")),
# #             size = 5)
#
#
# # ----
# # ----
# # ----
# # MB AGE - Plan ----
# # C3 ONT is now being used in place of C3 short read, based on correlation findings.
# # MB AGE - LIBRARY COMPOSITION CHECK ----
#
# C1_MB_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#     # }, {
#     #   plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#     #     select("Ddc" = count)
#   }, {
#     plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C1_MB_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C1")
#
#
# C2_MB_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#     # }, {
#     #   plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#     #     select("Ddc" = count)
#   }, {
#     plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C2_MB_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C2")
#
#
# C3_MB_IP_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#     # }, {
#     #   plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#     #     select("Ddc" = count)
#   }, {
#     plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_C3_MB_IP, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C3")
#
#
#
#
#
# C3_MB_IP_ONT_MARKERS <- bind_cols(
#   {
#     plotCounts(dds_ONT, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
#       select(age, genotype, "Slc6a3" = count)
#   }, {
#     plotCounts(dds_ONT, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
#       select("Th" = count)
#     # }, {
#     #   plotCounts(dds_ONT, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
#     #     select("Ddc" = count)
#   }, {
#     plotCounts(dds_ONT, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
#       select("Slc18a2" = count)
#   }, {
#     plotCounts(dds_ONT, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
#       select("Gfap" = count)
#   }, {
#     plotCounts(dds_ONT, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
#       select("Gad2" = count)
#   }, {
#     plotCounts(dds_ONT, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
#       select("S100b" = count)
#   }
# ) %>%
#   pivot_longer(-c(age, genotype),
#                names_to = "gene",
#                values_to = "count") %>%
#   mutate(cohort = "C3 ONT")
#
#
#
# # Take home points from this plot:
# #
# # Enrichment can be seen (and depletion of non-DA markers)
# #
# # In some cohorts (C1 and a little bit C2/3), DA markers go up
# # while in C3 ONT, markers go down
# # this creates a problem in each cohort individually, as other DA genes
# # may be classed as DE
# #
# # The different collection strategy of C1 shows: DA markers are definitely higher there
# # and non DA are lower (this is the only cohort where this is the case)
# #
# # By using meta-analysis, differences in DA markers and associated genes that follow
# # are screened, preserving only non-conflicting changes across cohorts
# #
# # Other point: non-DA markers are really going up in all cohorts apart from C1
# # This emphasises the importance of holding imited value in results from
# # DEPLETED genes.
# p_mb_ip_age_markers <- bind_rows(C1_MB_IP_MARKERS,
#                                              C2_MB_IP_MARKERS,
#                                              C3_MB_IP_MARKERS,
#                                              C3_MB_IP_ONT_MARKERS) %>%
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
#   labs(x = "Age and Cohort",
#        y = "Log2 Expression Count",
#        colour = "Cell Type Marker") +
#   theme(legend.position = "top")
#
#
#
#
#
#
# # C1_MB_TOTAL_MARKERS <- bind_cols(
# #   {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
# #       select(age, genotype, "Slc6a3" = count)
# #   }, {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
# #       select("Th" = count)
# #   }, {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
# #       select("Ddc" = count)
# #   }, {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
# #       select("Slc18a2" = count)
# #   }, {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
# #       select("Gfap" = count)
# #   }, {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
# #       select("Gad2" = count)
# #   }, {
# #     plotCounts(dds_C1_MB_TOTAL, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
# #       select("S100b" = count)
# #   }
# #
# # ) %>%
# #   pivot_longer(-c(age, genotype),
# #                names_to = "gene",
# #                values_to = "count") %>%
# #   mutate(cohort = "C1")
# #
# #
# # C2_MB_TOTAL_MARKERS <- bind_cols(
# #   {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
# #       select(age, genotype, "Slc6a3" = count)
# #   }, {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
# #       select("Th" = count)
# #   }, {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
# #       select("Ddc" = count)
# #   }, {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
# #       select("Slc18a2" = count)
# #   }, {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
# #       select("Gfap" = count)
# #   }, {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
# #       select("Gad2" = count)
# #   }, {
# #     plotCounts(dds_C2_MB_TOTAL, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
# #       select("S100b" = count)
# #   }
# # ) %>%
# #   pivot_longer(-c(age, genotype),
# #                names_to = "gene",
# #                values_to = "count") %>%
# #   mutate(cohort = "C2")
# #
# #
# # C3_MB_TOTAL_MARKERS <- bind_cols(
# #   {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Slc6a3"), c("age", "genotype"), returnData = T) %>%
# #       select(age, genotype, "Slc6a3" = count)
# #   }, {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Th"), c("age", "genotype"), returnData = T) %>%
# #       select("Th" = count)
# #   }, {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Ddc"), c("age", "genotype"), returnData = T) %>%
# #       select("Ddc" = count)
# #   }, {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Slc18a2"), c("age", "genotype"), returnData = T) %>%
# #       select("Slc18a2" = count)
# #   }, {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Gfap"), c("age", "genotype"), returnData = T) %>%
# #       select("Gfap" = count)
# #   }, {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("Gad2"), c("age", "genotype"), returnData = T) %>%
# #       select("Gad2" = count)
# #   }, {
# #     plotCounts(dds_C3_MB_TOTAL, get_ensembl_gene_id("S100b"), c("age", "genotype"), returnData = T) %>%
# #       select("S100b" = count)
# #   }
# # ) %>%
# #   pivot_longer(-c(age, genotype),
# #                names_to = "gene",
# #                values_to = "count") %>%
# #   mutate(cohort = "C3")
# #
# #
# # # Not very revealing: Or perhaps it is: DA and non DA markers with overlapping abundance,
# # # some up in age, some down, so not a clear picture.
# # # Probably not evidence of a difference in library composition here.
# # # The number of detected genes is no different, either (First results chapter).
# # # This is good!
# # p_mb_total_age_library_composition <- bind_rows(C1_MB_TOTAL_MARKERS,
# #                                                 C2_MB_TOTAL_MARKERS,
# #                                                 C3_MB_TOTAL_MARKERS) %>%
# #   mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
# #                          "Dopaminergic", "Other")) %>%
# #   # group_by(age, gene) %>%
# #   # summarise(count = mean(count)) %>%
# #   ggplot(aes(x = age,
# #              y = log2(count + 1),
# #              colour = marker,
# #              group = gene)) +
# #   geom_smooth(method = "lm") +
# #   facet_wrap(vars(cohort), nrow = 1) +
# #   panel_border() +
# #   scale_color_d3() +
# #   labs(x = "Age",
# #        y = "Log2 Expression Count",
# #        colour = "Cell Type Marker") +
# #   theme(legend.position = "top")
#
#
#
#
#
#
#
# # MB AGE - ALTERNATIVE COMPOSITION CHECK ----
#
# LIBRARY_PROPORTIONS <- bind_rows(
#   {
#     as_tibble(counts(dds_C1_MB_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C1")
#   }, {
#     as_tibble(counts(dds_C2_MB_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C2")
#   }, {
#     as_tibble(counts(dds_C3_MB_IP), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C3")
#   }, {
#     as_tibble(counts(dds_ONT), rownames = "ensembl_gene_id") %>%
#       mutate(cohort = "C3 ONT")
#   }
# ) %>%
#   pivot_longer(-c(ensembl_gene_id, cohort),
#                names_to = "sample_name",
#                values_to = "count") %>%
#   mutate(status = ifelse(ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
#                          "Enriched",
#                          ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                 "Depleted", "Unchanged")),
#          status = factor(status, levels = c("Depleted", "Unchanged", "Enriched"))) %>%
#   mutate(age = ifelse(str_detect(sample_name, "YOUNG"), "YOUNG", "OLD"),
#          age = factor(age, levels = c("YOUNG", "OLD"))) %>%
#   group_by(cohort, age, status) %>%
#   summarise(count = sum(count, na.rm = T)) %>%
#   mutate(prop = count/sum(count))
#
# p_age_mb_library_proportions <- LIBRARY_PROPORTIONS %>%
#   ggplot(aes(x = age,
#              y = prop,
#              fill = status)) +
#   geom_col(colour = "black") +
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
#
# # MB AGE TOTAL DESeq2 -----------------------------------------------------
# # TRAP will not have been performed on these samples, so there isnt a batch difference to worry about
#
# # C1 TOTAL MB
# dds_MB_TOTAL_AGE <- dds_C1_MB_TOTAL %>% filter_zeros()
#
# dds_MB_TOTAL_AGE@design <- ~ age
#
# colData(dds_MB_TOTAL_AGE) <- droplevels(colData(dds_MB_TOTAL_AGE))
#
# dds_MB_TOTAL_AGE <- DESeq(
#   dds_MB_TOTAL_AGE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_MB_TOTAL_AGE <- DESeq2::results(
#   dds_MB_TOTAL_AGE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_MB_TOTAL_AGE <- lfcShrink(
#   dds_MB_TOTAL_AGE,
#   res = res_MB_TOTAL_AGE,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_MB_TOTAL_AGE)
#
# # dds_C3_MB_TOTAL_AGE <- dds_C3_MB_TOTAL[rownames(dds_C3_MB_TOTAL) %in% genes_MB_IP_AGE,]
# #
# # dds_C3_MB_TOTAL_AGE@design <- ~ age
# #
# # colData(dds_C3_MB_TOTAL_AGE) <- droplevels(colData(dds_C3_MB_TOTAL_AGE))
# #
# # dds_C3_MB_TOTAL_AGE <- DESeq(
# #   dds_C3_MB_TOTAL_AGE,
# #   sfType = "poscounts",
# #   fitType = "local",
# #   minReplicatesForReplace = Inf,
# #   parallel = TRUE
# # )
# #
# # res_C3_MB_TOTAL_AGE <- DESeq2::results(
# #   dds_C3_MB_TOTAL_AGE,
# #   cooksCutoff = Inf,
# #   filterFun = ihw,
# #   parallel = T
# # )
# # res_C3_MB_TOTAL_AGE <- lfcShrink(
# #   dds_C3_MB_TOTAL_AGE,
# #   res = res_C3_MB_TOTAL_AGE,
# #   contrast = c("age", "OLD", "YOUNG"),
# #   type = "ashr",
# #   parallel = T
# # )
# #
# # summary(res_C3_MB_TOTAL_AGE)
#
#
# MB_AGE_TOTAL_META <-
#   as_tibble(res_MB_TOTAL_AGE, rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   select(
#     ensembl_gene_id,
#     external_gene_name,
#     description,
#     everything()
#   ) %>%
#   mutate(outcome = ifelse(log2FoldChange > 0,
#                           "Up", "Down"),
#          outcome = ifelse(padj < 0.1,
#                           outcome, "Unchanged"),
#          mb_translation = ifelse(ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#                                  "Translated",
#                                  ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                         "Depleted", "Unchanged")))
#
#
#
# MB_AGE_TOTAL_OUTCOMES <- MB_AGE_TOTAL_META %>%
#   select(ensembl_gene_id,
#          outcome)
#
#
# # IF ALSO USING C3 DATA
# # MB_AGE_TOTAL_META <-
# #   as_tibble(res_C1_MB_TOTAL_AGE, rownames = "ensembl_gene_id") %>%
# #   left_join(
# #     as_tibble(res_C3_MB_TOTAL_AGE, rownames = "ensembl_gene_id"),
# #     by = "ensembl_gene_id",
# #     suffix = c("_C1", "_C3")
# #   ) %>%
# #   left_join(anno) %>%
# #   select(
# #     ensembl_gene_id,
# #     external_gene_name,
# #     description,
# #     everything()
# #   ) %>%
# #   mutate(
# #     across(starts_with("pvalue"), ~ replace_na(.x, 1 - 1e-10)),
# #     # replace NA with *almost* 1 (it can't be perfectly 1)
# #     across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)),
# #     # replace 0 values with minimum non-zero of the column
# #     across(starts_with("pvalue"), ~ .x / 2),
# #     # divide 2-sided pvalue by 2 into 1-sided
# #     across(starts_with("pvalue"), ~ ifelse(.x == 0, min(.x[.x > 0]), .x)) # again replace 0 values with minimum non-zero of the column
# #   ) %>%
# #   select(sort(names(.))) %>%
# #   select(
# #     ensembl_gene_id,
# #     external_gene_name,
# #     description,
# #     everything(),
# #     -starts_with("padj")
# #   ) %>%
# #   filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
# #   drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
# #   rowwise() %>%
# #   mutate(
# #     pvalue_C3 = ifelse(
# #       # correct pvalues where there is a conflict in log2foldchange between groups
# #       sign(log2FoldChange_C1) == sign(log2FoldChange_C3),
# #       pvalue_C3,
# #       1 - pvalue_C3
# #     )) %>%
# #   mutate(across(starts_with("pvalue"),
# #                 ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
# #   mutate(
# #     sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
# #     log2FoldChange = log2FoldChange_C1
# #   ) %>%
# #   ungroup() %>%
# #   mutate(
# #     sumz_adj = p.adjust(sumz, method = "fdr"),
# #     conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
# #   arrange(desc(log2FoldChange)) %>%
# #   mutate(outcome = ifelse(log2FoldChange > 0,
# #                           "Up", "Down"),
# #          outcome = ifelse(sumz_adj < 0.1,
# #                           outcome, "Unchanged"),
# #          mb_translation = ifelse(ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
# #                                  "Translated",
# #                                  ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
# #                                         "Depleted", "Unchanged")))
#
#
# # MB AGE - C1 C2 C3_ONT MERGED DATASET ----
#
# dds_C1_C2_C3ONT_MB_IP <-
#   dds[rownames(dds) %in% rownames(dds_ONT),
#       colData(dds)$cohort != "C3" &
#         colData(dds)$fraction == "IP" &
#         colData(dds)$region == "MB"]
# dds_C1_C2_C3ONT_MB_IP_AGE_FILTER <-
#   filter_genes(dds_C1_C2_C3ONT_MB_IP,
#                grouping = c("cohort", "age", "gene_id"))
# dds_ONT_AGE_FILTER <- filter_genes(dds_ONT,
#                                    grouping = c("age", "gene_id"))
#
# genes_MB_IP_AGE <- intersect(
#   dds_C1_C2_C3ONT_MB_IP_AGE_FILTER,
#   dds_ONT_AGE_FILTER
# )
#
# counts_C1_C2_C3ONT_MB_IP <- counts(dds_C1_C2_C3ONT_MB_IP)[rownames(dds_C1_C2_C3ONT_MB_IP) %in% genes_MB_IP_AGE,]
# counts_ONT <- counts(dds_ONT)[rownames(dds_ONT) %in% genes_MB_IP_AGE,]
#
# all(rownames(counts_ONT) == rownames(counts_C1_C2_C3ONT_MB_IP))
#
# dim(counts_C1_C2_C3ONT_MB_IP)
# dim(counts_ONT)
#
# counts_C1_C2_C3ONT_MB_IP <- cbind(counts_C1_C2_C3ONT_MB_IP,
#                                   counts_ONT)
#
# coldata_C1_C2_C3ONT_MB_IP <- colData(dds_C1_C2_C3ONT_MB_IP)[,c("sample_name", "cohort", "collection", "age", "genotype")]
# coldata_ONT <- colData(dds_ONT)[,c("sample_name", "cohort", "collection", "age", "genotype")]
# coldata_C1_C2_C3ONT_MB_IP <- rbind(coldata_C1_C2_C3ONT_MB_IP, coldata_ONT)
#
# dds_C1_C2_C3ONT_MB_IP <- DESeqDataSetFromMatrix(counts_C1_C2_C3ONT_MB_IP,
#                                                 coldata_C1_C2_C3ONT_MB_IP,
#                                                 design = ~ collection + age)
#
# # dds_C1_C2_C3ONT_MB_IP <- dds_C1_C2_C3ONT_MB_IP[rownames(dds_C1_C2_C3ONT_MB_IP) %in% MB_FRACTION_TRANSLATED_GENES, ]
#
# dds_C1_C2_C3ONT_MB_IP <- DESeq(dds_C1_C2_C3ONT_MB_IP,
#                                minReplicatesForReplace = Inf,
#                                parallel = TRUE)
#
# res_C1_C2_C3ONT_MB_IP <- DESeq2::results(
#   dds_C1_C2_C3ONT_MB_IP,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   alpha = 0.01,
#   lfcThreshold = log2(1.05),
#   parallel = T
# )
# res_C1_C2_C3ONT_MB_IP <- lfcShrink(
#   dds_C1_C2_C3ONT_MB_IP,
#   res = res_C1_C2_C3ONT_MB_IP,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "apeglm",
#   lfcThreshold = log2(1.05),
#   parallel = T
# )
#
# summary(res_C1_C2_C3ONT_MB_IP)
#
# MB_AGE_META <- res_C1_C2_C3ONT_MB_IP %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   # filter(padj < 0.01) %>%
#   mutate(outcome = ifelse(log2FoldChange > 0,
#                           "Up", "Down"),
#          outcome = ifelse(padj < 0.01, outcome, "Unchanged"),
#          status = ifelse(ensembl_gene_id %in% MB_FRACTION_ENRICHED_GENES,
#                          "Enriched",
#                          ifelse(ensembl_gene_id %in% MB_FRACTION_DTU_NOT_ENRICHED_GENES &
#                                   !ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                 "Spliced",
#                                 ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                        "Depleted", "Unchanged")))) %>%
#   left_join(publications_all_genes) %>%
#   left_join(MB_AGE_TOTAL_OUTCOMES, by = "ensembl_gene_id", suffix = c("_TRAP", "_TOTAL")) %>%
#   mutate(specific = ifelse(outcome_TRAP == outcome_TOTAL,
#                            "Shared",
#                            ifelse(outcome_TRAP == "Up",
#                                   ifelse(outcome_TOTAL == "Unchanged",
#                                          "Specific",
#                                          "Specific Opposing"),
#                                   ifelse(outcome_TOTAL == "Unchanged",
#                                          "Specific",
#                                          "Specific Opposing"))),
#          specific = replace_na(specific, "Specific"))
#
# # MB AGE ENRICHMENT CATEGORY ----
#
# p_age_translation_status <- MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   left_join(enrichment_status) %>%
#   ggplot(aes(x = enrichment,
#              fill = enrichment)) +
#   geom_bar(colour = "black") +
#   scale_fill_d3() +
#   # scale_fill_manual(values = pal_d3()(5)[c(4, 3, 1, 5)]) +
#   labs(x = "Enrichment Category",
#        y = "Number of Genes") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 7)) +
#   geom_text(aes(label = ..count..,
#                 y = ..count..),
#             stat = "count",
#             fontface = 2,
#             vjust = -0.5) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#
#  # MB_AGE_META %>%
#  #  filter(padj < 0.01) %>%
#  #  ggplot(aes(x = status,
#  #             fill = status)) +
#  #  geom_bar(colour = "black") +
#  #  scale_fill_manual(values = pal_d3()(5)[c(4, 3, 1, 5)]) +
#  #  labs(x = "Translation Category",
#  #       y = "Number of Genes") +
#  #  theme(legend.position = "none",
#  #        axis.text.x = element_text(size = 7)) +
#  #  geom_text(aes(label = ..count..,
#  #                y = ..count..),
#  #            stat = "count",
#  #            fontface = 2,
#  #            vjust = -0.5) +
#  #  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#
# # MB AGE TRAP ENRICHMENT OUTCOME ----
#
# plot_data <- MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   left_join(enrichment_status)
#
# p_mb_age_enrichment_outcome <- plot_data %>%
#   ggplot(aes(x = outcome_TRAP)) +
#   geom_bar(position = "fill",
#            colour = "black",
#            aes(fill = enrichment)) +
#   scale_fill_d3() +
#   # scale_fill_manual(values = pal_d3()(5)[c(4, 3, 1, 5)]) +
#   labs(x = "Change in Age",
#        y = "Percentage of Genes",
#        fill = "TRAP\nSpecificity") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   theme(legend.position = "top",
#         legend.title = element_blank()) +
#   geom_text(data = . %>%
#               group_by(outcome_TRAP, enrichment) %>%
#               tally() %>%
#               mutate(p = n / sum(n)) %>%
#               ungroup(),
#             aes(y = p,
#                 group = enrichment,
#                 label = ifelse(p > 0.05, scales::percent(signif(p, 2)), "")),
#             position = position_stack(vjust = 0.5),
#             show.legend = FALSE,
#             colour = "white",
#             fontface = 2) +
#   geom_text(data = plot_data %>%
#               group_by(outcome_TRAP) %>%
#               tally(),
#             aes(y = 1,
#                 label = paste("Total:", n)),
#             vjust = -0.5) +
#   guides(fill = guide_legend(nrow = 2,
#                              byrow=TRUE))
#
# # MB AGE TRAP SPECIFICITY ----
#
# p_mb_age_specificity <- MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   filter(status %in% c("Enriched", "Spliced")) %>%
#   group_by(outcome_TRAP) %>%
#   mutate(n = dplyr::n()) %>%
#   ggplot(aes(x = outcome_TRAP,
#              fill = specific)) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(4)[c(4, 3, 1)]) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
#   geom_text(aes(label = n,
#                 y = n),
#             stat = "identity",
#             vjust = -0.5,
#             check_overlap = T,
#             fontface = 2) +
#   labs(x = "Change in Age",
#        y = "Number of Genes",
#        fill = "TRAP\nSpecificity") +
#   theme(legend.position = "top",
#         legend.title = element_blank()) +
#   guides(fill = guide_legend(nrow = 3,
#                              byrow=TRUE))
#
#
#
# # MB AGE - MAKE A LIST OF ALL CONFIDENT DEGS -------
#
# MB_AGE_DE_SPECIFIC_META <- MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   filter(specific != "Shared" &
#            status %in% c("Enriched", "Spliced")) %>%
#   select(ensembl_gene_id,
#          external_gene_name,
#          description,
#          baseMean,
#          log2FoldChange,
#          padj,
#          status,
#          pub_pd,
#          specific,
#          outcome_TRAP)
#
# MB_AGE_DE_SPECIFIC_UP <- MB_AGE_DE_SPECIFIC_META %>%
#   filter(log2FoldChange > 0)
#
# MB_AGE_DE_SPECIFIC_DOWN <- MB_AGE_DE_SPECIFIC_META %>%
#   filter(log2FoldChange < 0)
#
# # MB AGE Specificity Examples (PLOT) ----
# gene_ids <- bind_rows({
#   #   MB_AGE_DE_SPECIFIC_META %>%
#   #     filter(baseMean > 50 &
#   #              specific == "Shared") %>%
#   #     filter(external_gene_name == "Ass1") %>% # Or Tyrobp
#   #     select(specific,
#   #            ensembl_gene_id,
#   #            external_gene_name)
#   # }, {
#   MB_AGE_DE_SPECIFIC_META %>%
#     filter(baseMean > 50 &
#              specific == "Specific") %>%
#     filter(external_gene_name == "Mtfr1") %>%
#     select(specific,
#            ensembl_gene_id,
#            external_gene_name)
# }, {
#   MB_AGE_DE_SPECIFIC_META %>%
#     filter(baseMean > 50 &
#              specific == "Specific Opposing") %>%
#     filter(external_gene_name == "Scg2") %>%
#     select(specific,
#            ensembl_gene_id,
#            external_gene_name)
# }
# )
#
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
#
#
#
#
# # MB AGE NUMBERS ----
#
# n_MB_AGE_DEG <- MB_AGE_META %>%
#   filter(padj < 0.01) %>%
#   nrow
#
# n_MB_AGE_DEG_SPECIFIC_OPPOSING_TRANSLATED <- MB_AGE_DE_SPECIFIC_META %>%
#   filter(specific == "Specific Opposing") %>%
#   nrow
#
# # MB age STRING interactions ------------------------------------
#
# # I want to visualise the genes with DE to see how they cluster
# # and whether clusters display similar changes in gene expression
# # (e.g. a mitchondrial cluster goes down)
# # retrieve the string IDs for MB age DEs
#
# # USING STRING HOMOLOGS: More genes are retained in the end (around 75 %)
#
# string_MB_AGE_HUMAN <- string_db$map(data.frame(ensembl_gene_id = MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
#                                                 external_gene_name = MB_AGE_DE_SPECIFIC_META$external_gene_name),
#                                      "external_gene_name",
#                                      removeUnmappedRows = FALSE)
#
# # string_MB_AGE_MOUSE <- string_db_mouse$map(data.frame(ensembl_gene_id = MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
# #                                                 external_gene_name = MB_AGE_DE_SPECIFIC_META$external_gene_name),
# #                                "ensembl_gene_id",
# #                                removeUnmappedRows = FALSE)
#
# input_strings <- string_MB_AGE_HUMAN %>%
#   drop_na %>%
#   pull(STRING_id)
#
# # % of AGE DE that are in STRING
# length(input_strings) / nrow(MB_AGE_DE_SPECIFIC_META)
#
# # retrieve the interactions
# interactions <- string_db$get_interactions(input_strings) %>%
#   mutate(from = string_MB_AGE_HUMAN$external_gene_name[match(from, string_MB_AGE_HUMAN$STRING_id)],
#          to = string_MB_AGE_HUMAN$external_gene_name[match(to, string_MB_AGE_HUMAN$STRING_id)]) %>%
#   drop_na() %>%
#   distinct()
#
# # % of STRING-mapped AGE DEs with known interactions
# length(unique(c(interactions$from, interactions$to))) / length(input_strings)
#
# # create the igraph
# igr <- graph_from_data_frame(interactions,
#                              directed = F)
#
# # add ensembl_gene_id
# V(igr)$ensembl_gene_id <- string_MB_AGE_HUMAN[match(V(igr)$name, string_MB_AGE_HUMAN$external_gene_name),]$ensembl_gene_id
#
# # rename the main name to ensembl_gene_id
# V(igr)$external_gene_name <- V(igr)$name
# V(igr)$name <- V(igr)$ensembl_gene_id
#
# # convert combined score to decimal
# E(igr)$combined_score <- E(igr)$combined_score/1000
# # plot the distribution of scores
# plot(hist(E(igr)$combined_score))
# # STRING says that a score of 0.5 means 1 in 2 are likely false positives
# # I choose to filter so that 1 in 5 are false positives (score >= 0.8)
# # filter edges for confidence and remove orphan nodes
# igr <- delete.edges(igr, which(E(igr)$combined_score < 0.8 ))
#
# # MB AGE DEG Degree (PLOT) ----
#
# plot_data <- enframe(degree(igr)) %>%
#   # left_join(anno_human,
#   #           by = c("name" = "hsapiens_homolog_associated_gene_name")) %>%
#   select(ensembl_gene_id = name,
#          value) %>%
#   filter(ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_META$ensembl_gene_id) %>%
#   right_join(MB_AGE_DE_SPECIFIC_META) %>%
#   # mutate(value = replace_na(value, 0)) %>%
#   select(ensembl_gene_id:description,
#          degree = value) %>%
#   mutate(degree_bin =
#            ifelse(is.na(degree),
#                   "No\nSTRING\nIdentifier",
#                   ifelse(degree == 0,
#                          "0",
#                          ifelse(degree < 5,
#                                 "1-4",
#                                 ifelse(degree < 10,
#                                        "5-9",
#                                        "10+")))),
#          degree_bin = factor(degree_bin,
#                              levels = c("No\nSTRING\nIdentifier",
#                                         "0",
#                                         "1-4",
#                                         "5-9",
#                                         "10+")))
#
# MB_AGE_DEGREE <- plot_data %>%
#   select(external_gene_name,
#          degree)
#
# p_mb_age_degree <- plot_data %>%
#   ggplot(aes(x = degree_bin, fill = degree_bin)) +
#   geom_bar(colour = "black") +
#   geom_text(aes( label = ..count..,
#                  y= ..count.. ), stat= "count", vjust = -.5) +
#   scale_fill_d3() +
#   labs(x = "Number of STRING\nInteractions",
#        y = "Number of Genes") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   theme(legend.position = "none",
#         panel.grid.major.x = element_blank(),
#         axis.text.x = element_text(size = 8))
#
# # MB AGE TABLE OF HIGHEST DEGREE GENES ----
#
# plot_data %>%
#   filter(degree_bin %in% c("5-9", "10+"))
#
# # MB AGE further STRING analysis ----
#
# # add lfc values
# V(igr)$log2FoldChange <- MB_AGE_DE_SPECIFIC_META[match(V(igr)$ensembl_gene_id, MB_AGE_DE_SPECIFIC_META$ensembl_gene_id),]$log2FoldChange
# # check for NA values
# sum(is.na(V(igr)$log2FoldChange))
#
# # check degree
# # degree(igr)
#
# # remove 0 degree vertices
# igr <- induced_subgraph(igr, degree(igr) > 0)
#
# # remove quadruplets and smaller
# igr <- induced_subgraph(igr, components(igr)$membership %in% which(components(igr)$csize > 2))
#
# # define a layout
# layout <- layout_nicely(igr)
# # layout <- layout.fruchterman.reingold(igr)
# # layout <- layout.kamada.kawai(igr)
#
# plot(igr,
#      layout = layout,
#      edge.color = "grey70",
#      vertex.label = NA,
#      vertex.size = 5,
#      vertex.color = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
# )
#
# # simplify igr
# # igr <- simplify(igr)
#
# # cluster with fast greedy
# clusters <- cluster_fast_greedy(igr,
#                                 weights = E(igr)$combined_score)
#
# plot(
#   clusters,
#   igr,
#   layout = layout,
#   edge.color = "grey70",
#   vertex.label = NA,
#   vertex.size = 4,
#   # col = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
# )
#
# # percentage of edges that cross communities (I am deleting these)
# crossing(clusters, igr) %>% sum / length(E(igr))
#
# # delete crossing edges
# igr <- delete.edges(igr, which(crossing(clusters, igr)))
#
# # redo layout
# layout <- layout_nicely(igr)
#
# # MB AGE STRING Clusters ----
# # plot
# p_mb_age_deg_igraph <- as.grob(
#   ~ plot(
#     clusters,
#     igr,
#     layout = layout,
#     edge.color = "grey70",
#     vertex.label = NA,
#     vertex.size = 3,
#     vertex.color = ifelse(V(igr)$log2FoldChange < 0, "blue", "red"),
#   )
# )
#
# cluster_order <- membership(clusters)[order(membership(clusters))]
#
# sapply(split(names(cluster_order), cluster_order), function(x){sapply(x, get_external_gene_name)})
#
# # Trying to automate GO search and labelling is a nightmare
# # Opting to label group with the gene of highest degree per group instead
# # temp <- gost(query = split(names(cluster_order), cluster_order),
# #      organism = "mmusculus",
# #      multi_query = T,
# #      custom_bg = NULL)
# # temp$result %>%
# #   mutate(significant = str_remove_all(significant, "c\\(|\\)")) %>%
# #   separate(col = significant,
# #            into = paste0("signif_", unique(cluster_order)),
# #            sep = ", ") %>%
# #   pivot_longer(starts_with("signif"),
# #                names_to = "group",
# #                values_to = "signif",
# #                names_prefix = "signif_") %>%
# #   mutate(signif = as.logical(signif)) %>%
# #   filter(signif) %>%
# #   arrange(group) %>%
# #   filter(source == "REAC") %>% View
#
# # MB AGE Clusters Heatmap ----
#
# # get genes with highest degree
# gene_labels <- sapply(split(degree(igr), membership(clusters)), function(x) names(sort(x, decreasing = T)[1]) )
# # Convert to readable names
# gene_labels <- sapply(gene_labels, get_external_gene_name)
#
# # expand to the number of genes
# gene_labels <- gene_labels[membership(clusters)]
#
# # plot an ageing heatmap
# plot_data <- t(scale(t(assay(vst(dds_C1_MB_IP))[match(names(membership(clusters)), rownames(dds_C1_MB_IP)),order(colData(dds_C1_MB_IP)$age)])))
# plot_data_quantiles <- quantile(plot_data, probs = seq(0, 1, length.out = 10))
#
# # make an annotation column df
# annot_cols <- data.frame(Age = str_extract(colnames(plot_data), pattern = "(?<=MB_)[:alpha:]+"))
# rownames(annot_cols) <- colnames(plot_data)
# # make an annotation row df
# annot_rows <- data.frame(Group = gene_labels)
# # correct non-unique rownames
# rownames(plot_data) <- make.unique(rownames(plot_data))
# rownames(annot_rows) <- rownames(plot_data)
#
# mat_colors <- list(Age = pal_d3()(10)[c(3, 2)],
#                    Group = pal_d3(palette = "category20")(20)[c(1:length(unique(gene_labels)))])
# names(mat_colors$Age) <- unique(annot_cols$Age)
# names(mat_colors$Group) <- unique(gene_labels)
#
# p_mb_age_deg_heatmap <- as.grob(
#   ~ pheatmap(
#     mat = plot_data[order(names(gene_labels)),],
#     scale = "none",
#     treeheight_row = 0,
#     show_rownames = F,
#     show_colnames = F,
#     annotation_col = annot_cols,
#     annotation_row = annot_rows,
#     annotation_colors = mat_colors,
#     # color = brewer.pal(n = length(plot_data_quantiles) - 1,
#     #                    name = "PiYG"),
#     color = colorspace::diverging_hcl(11,
#                                       h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7)),
#     # color = colorspace::diverging_hcl(11,
#     #                                   h = c(180, 50), c = 80, l = c(20, 95), power = c(0.7, 1.3)),
#     # color = brewer.pal(n = 8,
#     #                    name = "RdBu"),
#     breaks = plot_data_quantiles,
#     border_color = NA,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     gaps_row = which(unlist(
#       lapply(seq(2:length(names(gene_labels)[order(names(gene_labels))])), function(n)
#         names(gene_labels)[order(names(gene_labels))][n] > names(gene_labels)[order(names(gene_labels))][n - 1])
#     )),
#     gaps_col = which(unlist(
#       lapply(seq(2:ncol(plot_data)), function(n)
#         colData(dds_C1_MB_IP)$age[n] != colData(dds_C1_MB_IP)$age[n-1])
#     ))
#     # cellwidth = 6,
#     # cellheight = 2
#   )
# )
# # ----
# # ----
# # ----
# # MB AGE DE SPECIFIC SIMPLE VIEWING TABLE ----
# MB_AGE_DE_SPECIFIC_META
#
# # MB AGE CLUSTER TABLE ----
#
# MB_AGE_TRAP_OUTCOMES <- MB_AGE_DE_SPECIFIC_META %>%
#   select(ensembl_gene_id,
#          status,
#          specific,
#          outcome_TRAP,
#          pub_pd)
#
# t_mb_age_cluster_membership <- enframe(membership(clusters), name = "ensembl_gene_id", value = "cluster") %>%
#   mutate(external_gene_name = sapply(ensembl_gene_id, get_external_gene_name),
#          cluster = as.character(cluster)) %>%
#   left_join(
#     enframe(gene_labels, name = "cluster", value = "lead") %>% distinct()) %>%
#   left_join(MB_AGE_DEGREE) %>%
#   left_join(MB_AGE_TRAP_OUTCOMES) %>%
#   left_join(anno) %>%
#   mutate(description = str_extract(description, ".+(?=\\[)"),
#          description = ifelse(str_detect(description,
#                                          "^[:upper:]"),
#                               description,
#                               str_to_sentence(description)),
#          cluster = factor(paste0("Cluster ", cluster), levels = paste0("Cluster ", seq(1, length(unique(cluster)))))) %>%
#   arrange(cluster, desc(degree))
#
# # t_mb_age_cluster_membership %>%
# #   mutate(GWAS = external_gene_name %in% GWAS_GENES_BROAD) %>%
# #   filter(GWAS) %>% View
#
# # ----
# # ----
# # ----
# # OLD
# # AXON AGE CONTINUED - MAKE THE METADATA TABLE ----
#
#
# AXON_AGE_META <- res_C2_AXON_IP_AGE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   left_join((
#     AXON_FRACTION_META %>%
#       select("external_gene_name",
#              expression_energy)
#   ),
#   by = c("external_gene_name")) %>%
#   mutate(
#     mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#     axon_translated = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES,
#     aba_specific = expression_energy < 0.5,
#     aba_specific = replace_na(aba_specific, FALSE),
#     # evidence = ifelse(
#     #   axon_translated,
#     #   ifelse(
#     #     aba_specific & mb_translated,
#     #     "Axon, Soma and ABA",
#     #     ifelse(
#     #       aba_specific,
#     #       "Axon and ABA",
#     #       ifelse(mb_translated,
#     #              "Axon and Soma",
#     #              "Axon only")
#     #     )
#     #   ),
#     #   ifelse(aba_specific & mb_translated,
#     #          "Soma and ABA",
#     #          "Soma only")
#     # ),
#     # evidence = ifelse(!mb_translated & evidence == "Soma only",
#     #                   "No support", evidence),
#     # evidence = factor(
#     #   evidence,
#     #   levels = c(
#     #     "Axon, Soma and ABA",
#     #     "Axon and Soma",
#     #     "Axon and ABA",
#     #     "Axon only",
#     #     "Soma and ABA",
#     #     "Soma only",
#     #     "No support"
#     #   )
#     # ),
#     # mb_age = ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
#     # mb_age_relationship = ifelse(
#     #   ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id,
#     #   ifelse(log2FoldChange > 0, "Common Upregulation",
#     #          "Opposing"),
#     #   ifelse(
#     #     ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_DOWN$ensembl_gene_id,
#     #     ifelse(log2FoldChange < 0, "Common Downregulation",
#     #            "Opposing"),
#     #     "Axon Specific"
#     #   )
#     # ),
#     # mb_age_relationship = ifelse(padj < 0.01, mb_age_relationship, "NS in Axon"),
#     outcome = ifelse(
#       padj < 0.01,
#       ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
#       "Unchanged"
#     )
#   ) %>%
#   left_join(AXON_AGE_TOTAL_OUTCOMES,
#             by = "ensembl_gene_id",
#             suffix = c("_TRAP", "_TOTAL")) %>%
#   mutate(
#     specific = ifelse(
#       outcome_TRAP == outcome_TOTAL,
#       "Shared",
#       ifelse(
#         outcome_TRAP == "Up",
#         ifelse(outcome_TOTAL == "Unchanged",
#                "Specific",
#                "Specific Opposing"),
#         ifelse(outcome_TOTAL == "Unchanged",
#                "Specific",
#                "Specific Opposing")
#       )
#     ),
#     specific = replace_na(specific, "Specific")
#   ) %>%
#   left_join(publications_all_genes_axon) %>%
#   left_join(publications_all_genes) %>%
#   dplyr::select(
#     ensembl_gene_id,
#     external_gene_name,
#     description,
#     # baseMean,
#     log2FoldChange,
#     padj,
#     mb_translated,
#     axon_translated,
#     # mb_age_relationship,
#     pub_pd,
#     pub_axon,
#     # expression_energy,
#     aba_specific,
#     # evidence,
#     specific,
#     outcome_TRAP,
#     outcome_TOTAL
#   ) %>%
#   left_join(enrichment_status)
#
# # AXON AGE - Number of genes in enrichment groups ----
#
# p_axon_age_enrichment_groups <- AXON_AGE_META %>%
#   filter(padj < 0.01 & log2FoldChange > 0) %>%
#   ggplot(aes(x = enrichment,
#              fill = enrichment %in% c("Soma only", "Soma and Axon"))) +
#   geom_bar(colour = "black") +
#   theme(legend.position = "none") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   scale_fill_manual(values = c("grey", pal_d3()(1)[1])) +
#   geom_text(aes(y = ..count..,
#                 label = ..count..),
#             stat = "count",
#             vjust = -0.5) +
#   labs(x = "Enrichment Category",
#        y = "Number of Genes")
#
# # AXON AGE TRAP SPECIFICITY (PLOT) ----
#
# p_axon_age_trap_specificity <- AXON_AGE_META %>%
#   filter(padj < 0.01 & log2FoldChange > 0 & enrichment %in% c("Soma only", "Soma and Axon")) %>%
#   ggplot(aes(x = specific,
#              fill = specific)) +
#   geom_bar(colour = "black") +
#   scale_fill_manual(values = pal_d3()(4)[c(4, 1)]) +
#   labs(x = "Enrichment Category",
#        y = "Number of Genes",
#        fill = "TRAP\nSpecificity") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   geom_text(aes(y = ..count..,
#                 label = ..count..),
#             stat = "count",
#             vjust = -0.5)
#
# AXON_AGE_SHARED_DE_META <- AXON_AGE_META %>%
#   filter(specific == "Shared") %>%
#   nrow
#
#
# AXON_AGE_META %>% filter(padj < 0.01 &
#                            log2FoldChange > 0 &
#                            enrichment %in% c("Soma only", "Soma and Axon")) %>%
#   filter(ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id) %>% View
#
# # AXON AGE EVIDENCE (PLOT) ----
#
# p_axon_age_evidence <- AXON_AGE_META %>%
#   left_join(enrichment_status) %>%
#   filter(padj < 0.1) %>%
#   filter(enrichment != "No enrichment") %>%
#   mutate(outcome_TRAP = ifelse(log2FoldChange > 0, "Increased\nabundance", "Decreased\nabundance")) %>%
#   # filter(outcome_TRAP != "Unchanged") %>%
#   # mutate(evidence = ifelse(mb_translated & axon_translated,
#   #                          "Soma and Axon",
#   #                          ifelse(mb_translated,
#   #                                 "Soma only",
#   #                                 ifelse(axon_translated,
#   #                                        "Axon only",
#   #                                        "No enrichment"))),
#   #        evidence = factor(evidence,
#   #                          levels = c("Soma only",
#   #                                     "Axon only",
# #                                     "Soma and Axon",
# #                                     "No enrichment"
# #                                     ))) %>%
# group_by(outcome_TRAP, enrichment) %>%
#   tally %>%
#   mutate(prop = n / sum(n),
#          sum = sum(n)) %>%
#   ggplot(aes(x = outcome_TRAP,
#              y = prop,
#              fill = enrichment)) +
#   geom_col(colour = "black") +
#   scale_fill_d3() +
#   geom_text(aes(label = paste("Total:", sum),
#                 y = 1),
#             check_overlap = T,
#             vjust = -0.5,
#             fontface = 2,
#             size = 4) +
#   scale_y_continuous(labels = scales::percent,
#                      limits = c(0, 1),
#                      expand = expansion(mult = c(0, 0.1))) +
#   labs(x = "Change in Age",
#        y = "Percentage of Genes",
#        fill = "Enrichment") +
#   geom_text(
#     aes(
#       x = outcome_TRAP,
#       y = prop,
#       label = ifelse(prop > 0.1,
#                      scales::percent(prop,
#                                      accuracy = 1),
#                      "")
#     ),
#     position = position_stack(vjust = 0.5),
#     colour = "white",
#     fontface = 2
#   ) +
#   theme(
#     legend.position = "top"
#   ) +
#   guides(fill = guide_legend(nrow = 2,
#                              byrow=TRUE,
#                              title.position="top",
#                              title.hjust = 0.5))
#
# # Axon AGE: Plot log fold changes and putative dopaminergic selection ----
#
#
# p_axon_age_putative_dopaminergic <- AXON_AGE_META %>%
#   left_join(enrichment_status) %>%
#   filter(padj < 0.01) %>%
#   mutate(dopaminergic = ifelse(enrichment == "Soma only" |
#                                  enrichment == "Soma and Axon" &
#                                  log2FoldChange < 0 |
#                                  enrichment == "Axon only" &
#                                  log2FoldChange < 0 |
#                                  enrichment == "No enrichment" &
#                                  log2FoldChange < 0, "Dopaminergic", "Other Cell Types"),
#          dopaminergic = ifelse(enrichment == "Soma and Axon" &
#                                  log2FoldChange > 0,
#                                "Potentially\nDopaminergic", dopaminergic)) %>%
#   ggplot(aes(x = enrichment,
#              y = log2FoldChange,
#              colour = dopaminergic)) +
#   geom_quasirandom() +
#   scale_y_continuous(limits = c(-2, 2)) +
#   scale_color_d3() +
#   theme(legend.position = "top",
#         legend.title = element_blank()) +
#   labs(x = "Enrichment Category",
#        y = expression(Log[2] ~ Fold ~ Change)) +
#   geom_hline(yintercept = 0,
#              linetype = "dotted")
#
#
#
# # AXON AGE: Assign putative dopaminergic based on aged axon TRAP downregulation ----
#
#
# putative_dopaminergic <- AXON_AGE_META %>%
#   left_join(enrichment_status) %>%
#   filter(padj < 0.1) %>%
#   mutate(dopaminergic = ifelse(enrichment == "Soma only" |
#                                  enrichment == "Soma and Axon" &
#                                  log2FoldChange < 0 |
#                                  enrichment == "Axon only" &
#                                  log2FoldChange < 0 |
#                                  enrichment == "No enrichment" &
#                                  log2FoldChange < 0, TRUE, FALSE)) %>%
#   filter(dopaminergic) %>%
#   pull(ensembl_gene_id)
#
#
#
# putative_nondopaminergic <- AXON_AGE_META %>%
#   left_join(enrichment_status) %>%
#   filter(padj < 0.1) %>%
#   mutate(dopaminergic = ifelse(enrichment == "Soma and Axon" &
#                                  log2FoldChange > 0 |
#                                  enrichment == "Axon only" &
#                                  log2FoldChange > 0 |
#                                  enrichment == "No enrichment" &
#                                  log2FoldChange > 0, TRUE, FALSE)) %>%
#   filter(dopaminergic) %>%
#   pull(ensembl_gene_id)
#
#
# # Using ABA specificity is no better than flipping a coin
# # AXON_AGE_META %>%
# #   left_join(enrichment_status) %>%
# #   filter(padj < 0.1) %>%
# #   # filter(enrichment != "No enrichment") %>%
# #   mutate(outcome = ifelse(log2FoldChange > 0, "Increased\nabundance", "Decreased\nabundance")) %>%
# #   group_by(outcome, enrichment, aba_specific) %>%
# #   tally %>%
# #   mutate(aba_specific = aba_specific,
# #          prop = n / sum(n),
# #          sum = sum(n)) %>%
# #   ggplot(aes(x = enrichment,
# #              y = prop,
# #              fill = aba_specific)) +
# #   geom_col(colour = "black") +
# #   scale_fill_d3() +
# #   facet_wrap(vars(outcome),
# #              nrow = 1)
#
#
#
# # ----
# # ----
# # ----
# # AXON AGE DA MARKER DOWNREG EXAMPLE (PLOT) ----
#
# p_axon_age_ds_vs <- view_counts(dds_C2_IP,
#                                gene = "Slc6a3") %>%
#  left_join(colData(dds_C2_IP),
#            copy = T) %>%
#  # filter(region != "MB") %>%
#  filter(fraction == "IP") %>%
#  group_by(region) %>%
#  filter(!outlier_iqr(count, "high")) %>%
#  ggplot(aes(x = age,
#             y = count)) +
#  geom_violin(colour = "black") +
#  geom_quasirandom(shape = 21,
#                   size = 3,
#                   colour = "black",
#                   aes(fill = region)) +
#  facet_wrap(vars(region),
#             scales = "free_y") +
#  scale_fill_d3() +
#  labs(x = "Region and Age",
#       y = "Expression Count",
#       title = "Slc6a3") +
#  theme(legend.position = "none")
#
#
# # AXON AGE ATG7 EXAMPLE (PLOT) ----
#
# p_axon_age_ds_vs_atg7 <- view_counts(dds_C2_IP,
#                                     gene = "Atg7") %>%
#  left_join(colData(dds_C2_IP),
#            copy = T) %>%
#  # filter(region != "MB") %>%
#  filter(fraction == "IP") %>%
#  group_by(region) %>%
#  filter(!outlier_iqr(count, "high")) %>%
#  ggplot(aes(x = age,
#             y = count)) +
#  geom_violin(colour = "black") +
#  geom_quasirandom(shape = 21,
#                   size = 3,
#                   colour = "black",
#                   aes(fill = region)) +
#  facet_wrap(vars(region),
#             scales = "free_y") +
#  scale_fill_d3() +
#  labs(x = "Region and Age",
#       y = "Expression Count",
#       title = "Atg7") +
#  theme(legend.position = "none")
#
#
#
#
# # AXON AGE PANGLAO CELL TYPES ----
#
# p_axon_age_panglao <- AXON_AGE_META %>%
#  left_join(anno_human) %>%
#  left_join(panglao,
#            by = c("hsapiens_homolog_associated_gene_name" = "official gene symbol")) %>%
#  mutate(padj = ifelse(padj == 0, min(padj[padj > 0]), padj),
#         score = -log10(padj) * log2FoldChange) %>%
#  group_by(`cell type`) %>%
#  summarise(score = mean(score),
#            group_proportion = dplyr::n()/group_size,
#            group_size = group_size,
#            prop_score = score * group_proportion) %>%
#  arrange(desc(prop_score)) %>%
#  distinct() %>%
#  filter(`cell type` %in% c("Astrocytes",
#                            "Dopaminergic neurons",
#                            "Microglia",
#                            "Oligodendrocytes",
#                            "Interneurons"
#  )) %>%
#  mutate(`cell type` = factor(`cell type`)) %>%
#  mutate(`cell type` = fct_relevel(`cell type`, "Dopaminergic neurons")) %>%
#  ggplot(aes(x = prop_score,
#             y = `cell type`,
#             fill = `cell type`)) +
#  geom_col(colour = "black") +
#  scale_fill_d3() +
#  labs(x = "Differential Expression\nScore",
#       y = "Cell Type") +
#  theme(legend.position = "none",
#        axis.text = element_text(size = 12),
#        axis.title.y = element_blank()) +
#   geom_vline(xintercept = 0,
#              linetype = "dotted")
#
# # AXON AGE GENES THAT GO UP ----
#
#
#
# n_AXON_AGE_DE <- AXON_AGE_META %>% filter(padj < 0.01) %>% nrow
#
# n_AXON_AGE_DE_UP_SOMA_ONLY <- AXON_AGE_META %>% filter(padj < 0.01 &
#                                                          mb_translated &
#                                                          !axon_translated &
#                                                          log2FoldChange > 0) %>%
#   nrow()
#
# n_AXON_AGE_DE_UP_SOMA_AND_AXON <- AXON_AGE_META %>% filter(padj < 0.01 &
#                                                          mb_translated &
#                                                          axon_translated &
#                                                          log2FoldChange > 0) %>%
#   nrow()
#
# n_AXON_AGE_DE_UP_SOMA_AND_AXON_MB_AGREE <- AXON_AGE_META %>% filter(padj < 0.01 &
#                            log2FoldChange > 0 &
#                            enrichment %in% c("Soma and Axon")) %>%
#   filter(ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id) %>%
#   nrow()
#
#
# # n_AXON_AGE_DE_UP_DOPAMINERGIC <- AXON_AGE_META %>%
# #   left_join(enrichment_status) %>%
# #  filter(padj < 0.01 &
# #           log2FoldChange > 0 &
# #           enrichment != "Axon only") %>%
# #  filter(mb_translated |
# #           axon_translated) %>%
# #  nrow
# #
# # n_AXON_AGE_DE_UP_DOPAMINERGIC_COMMON <- AXON_AGE_META %>%
# #  filter(padj < 0.01 &
# #           log2FoldChange > 0 &
# #           confidence != "Axon TRAP and ABA") %>%
# #  filter(mb_translated |
# #           axon_enriched & aba_specific) %>%
# #  filter(mb_age_relationship == "Common Upregulation") %>%
# #  nrow
#
#
# # AXON AGE SARM1 DS/VS ----
#
# dds_C2_DS_IP_AGE <- dds_C2_AXON_IP[,colData(dds_C2_AXON_IP)$region == "DS"] %>%
#  filter_zeros()
# dds_C2_DS_IP_AGE <- dds_C2_DS_IP_AGE[dds_C2_AXON_IP_AGE_FILTER,]
# dds_C2_DS_IP_AGE@design <- ~ collection + age
# colData(dds_C2_DS_IP_AGE) <- droplevels(colData(dds_C2_DS_IP_AGE))
# dds_C2_DS_IP_AGE <- DESeq(
#  dds_C2_DS_IP_AGE,
#  minReplicatesForReplace = Inf,
#  parallel = TRUE
# )
# res_C2_DS_IP_AGE <- DESeq2::results(
#  dds_C2_DS_IP_AGE,
#  cooksCutoff = Inf,
#  filterFun = ihw,
#  parallel = T
# )
# res_C2_DS_IP_AGE <- lfcShrink(
#  dds_C2_DS_IP_AGE,
#  res = res_C2_DS_IP_AGE,
#  contrast = c("age", "OLD", "YOUNG"),
#  type = "ashr",
#  parallel = T
# )
# DS_AGE_META <- as_tibble(res_C2_DS_IP_AGE, rownames = "ensembl_gene_id") %>%
#  mutate(outcome_DS_AGE = ifelse(padj < 0.01,
#                                 ifelse(log2FoldChange < 0,
#                                        "Downregulated", "Upregulated"),
#                                 "NS in DS")) %>%
#  left_join(anno) %>%
#  arrange(padj)
# # select(ensembl_gene_id, outcome_DS_AGE)
#
# # PLOT GENE TESTERS -----
#
#
#
# # AGE IP TESTER
# plotCounts(
#  dds,
#  get_ensembl_gene_id("Sarm1"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(
#    fraction == "IP"
#    # cohort != "C3"
#  ) %>%
#  ggplot(aes(x = age, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = age)) +
#  facet_grid(rows = vars(cohort), cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Age and Region",
#       y = "Expression Count",
#       title = "Cpne6") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
# # ENRICHMENT TESTER
# plotCounts(
#  dds,
#  get_ensembl_gene_id("Sarm1"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  # filter(fraction == "IP" &
#  #          cohort != "C3") %>%
#  ggplot(aes(x = fraction, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = fraction)) +
#  facet_grid(rows = vars(cohort), cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Age and Region",
#       y = "Expression Count",
#       title = "Sarm1") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
# # GENOTYPE IP TESTER
# plotCounts(
#  dds_C1_C2_C3ONT_MB_IP,
#  get_ensembl_gene_id("Caln1"),
#  c("cohort", "age", "genotype"),
#  returnData = T
# ) %>%
#  filter(
#    # fraction == "IP" &
#    cohort != "C2"
#  ) %>%
#  ggplot(aes(x = genotype, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = genotype)) +
#  facet_grid(rows = vars(cohort), cols = vars(age), scales = "free_y") +
#  panel_border() +
#  labs(x = "Genotype and Region",
#       y = "Expression Count",
#       title = "Cpne6") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
# # AXON AGE COPINE 6 ----
#
# p_axon_age_cpne6 <- plotCounts(
#  dds,
#  get_ensembl_gene_id("Cpne6"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(fraction == "IP" &
#           cohort != "C3") %>%
#  ggplot(aes(x = age, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = age)) +
#  facet_grid(cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Age and Region",
#       y = "Expression Count",
#       title = "Cpne6") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
# # AXON AGE NMNAT2 ----
#
# p_axon_age_nmnat2 <- plotCounts(
#  dds,
#  get_ensembl_gene_id("Nmnat2"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(fraction == "IP" &
#           cohort != "C3") %>%
#  ggplot(aes(x = age, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = age)) +
#  facet_grid(cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Age and Region",
#       y = "Expression Count",
#       title = "Nmnat2") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
#
# p_axon_fraction_nmnat2 <- plotCounts(
#  dds,
#  get_ensembl_gene_id("Nmnat2"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(count > 1) %>%
#  # filter(fraction == "IP" &
#  #          cohort != "C3") %>%
#  ggplot(aes(x = fraction, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = fraction)) +
#  facet_grid(cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Fraction and Region",
#       y = "Expression Count",
#       title = "Sarm1") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
#
# # AXON AGE SARM1 ----
#
# p_axon_age_sarm1 <- plotCounts(
#  dds,
#  get_ensembl_gene_id("Sarm1"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(count > 1) %>%
#  filter(fraction == "IP" &
#           cohort != "C3") %>%
#  ggplot(aes(x = age, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = age)) +
#  facet_grid(cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Age and Region",
#       y = "Expression Count",
#       title = "Sarm1") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
#
# p_axon_fraction_sarm1 <- plotCounts(
#  dds,
#  get_ensembl_gene_id("Sarm1"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(count > 1) %>%
#  # filter(fraction == "IP" &
#  #          cohort != "C3") %>%
#  ggplot(aes(x = fraction, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = fraction)) +
#  facet_grid(cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Fraction and Region",
#       y = "Expression Count",
#       title = "Sarm1") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
#
#
# # AXON AGE PTGES2 ----
#
# p_axon_age_ptges2 <- plotCounts(
#  dds,
#  get_ensembl_gene_id("Ptges2"),
#  c("cohort", "compartment", "age", "region", "fraction"),
#  returnData = T
# ) %>%
#  filter(fraction == "IP" &
#           cohort != "C3") %>%
#  ggplot(aes(x = age, y = count)) +
#  geom_violin() +
#  geom_quasirandom(aes(colour = age)) +
#  facet_grid(cols = vars(region), scales = "free_y") +
#  panel_border() +
#  labs(x = "Age and Region",
#       y = "Expression Count",
#       title = "Ptges2") +
#  scale_color_d3() +
#  scale_y_continuous(limits = c(0, NA)) +
#  theme(legend.position = "none")
#
#
#
#
# # ----
# # ----
# # ----
# # MB GENOTYPE Plan ----
# # Cohort 1 and Cohort 3 contain OVX samples
# # Cohort 1 Yound and Old samples were collected at different times, so there is a
# # confounding batch effect when comparing across these groups.
# # Cohort 3 was collected at the same time and sequenced short and long read
# # (7-10 replicates per group in short read, 3 replicates per group long read).
# # The cohort 3 short read data suffers because of the low input preparation.
# # So the main dataset used for genotype comparison is Cohort 3 long read.
# #
# # The first comparison to make is between Young OVX and Young WT animals.
# # No differences are found here, and the MA plot shows this (very little movement
# # from 0). Therefore, these samples are combined to a Young ALL group (n of 6)
# # for comparisons with Old samples.
# # Differences between Old WT and Young ALL samples are due to age. The list of
# # DEGs under this comparison is combined with the meta-analysis list of ageing
# # DEGs, to understand what component of Old OVX vs Young ALL DEGs share an
# # ageing component.
# # Note that ageing genes need to be restricted to genes that are ENRICHED
# # or SPLICED by TRAP.
# # Genotype genes between OLD OVX and YOUNG also need to be restricted, because
# # of library composition change.
# # Genotype OLD OVX vs OLD WT may not need restriction.
# # Finally, for the Old OVX vs Young ALL comparison (3 vs 6), DEGs are filtered,
# # requiring a pvalue of < 0.25(?) in a comparison of Old OVX vs Old WT,
# # so that genes that show no difference to Old WT samples are not mistakenly
# # categorised as OVX-related, when they are just ageing-related.
#
#
# # MB GENOTYPE - LIBRARY COMPOSITION CHECK ----
#
# p_mb_ip_genotype_library_composition <-
#   bind_rows(C1_MB_IP_MARKERS,
#                                              C3_MB_IP_MARKERS,
#                                              C3_MB_IP_ONT_MARKERS) %>%
#   mutate(marker = ifelse(gene %in% c("Th", "Slc6a3", "Ddc", "Slc18a2"),
#                          "Dopaminergic", "Other")) %>%
#   # group_by(age, gene) %>%
#   # summarise(count = mean(count)) %>%
#   ggplot(aes(x = genotype,
#              y = log2(count + 1),
#              colour = marker,
#              group = gene)) +
#   geom_smooth(method = "lm") +
#   facet_grid(rows = vars(cohort),
#              col = vars(age)) +
#   panel_border() +
#   scale_color_d3() +
#   labs(x = "Age",
#        y = "Log2 Expression Count",
#        colour = "Cell Type Marker") +
#   theme(legend.position = "top")
#
#
# # MB GENOTYPE - Make dds ONT objects --------------------------------------
# # make a young and old dds
# dds_ONT_MB_IP_YOUNG <- dds_ONT[,colData(dds_ONT)$age == "YOUNG"] %>% filter_zeros()
# dds_ONT_MB_IP_OLD <- dds_ONT[,colData(dds_ONT)$age == "OLD"] %>% filter_zeros()
# # MB GENOTYPE Young OVX vs Young WT subset common genes ----
# # filter genes for common detected
# dds_C1_MB_IP_YOUNG_FILTER <- filter_genes(dds_C1_MB_IP_YOUNG,
#                                         grouping = c("genotype", "gene_id"))
#
# # dds_C3_MB_IP_YOUNG_FILTER <- filter_genes(dds_C3_MB_IP_YOUNG,
# #                                           grouping = c("genotype", "gene_id"))
#
# dds_ONT_MB_IP_YOUNG_FILTER <- filter_genes(dds_ONT_MB_IP_YOUNG,
#                                       grouping = c("genotype", "gene_id"))
#
#
# # take subset: Going with ONT data!
#
# genes_MB_IP_YOUNG_GENOTYPE <- intersect(dds_C1_MB_IP_YOUNG_FILTER,
#                                       dds_ONT_MB_IP_YOUNG_FILTER)
#
# # MB genotype - YOUNG DESeq2 -----------------------------------------------------
# # This section confirms that the only gene DE between Young OVX and Young WT
# # Is Snca/SNCA
#
# dds_C1_MB_IP_YOUNG_GENOTYPE <- dds_C1_MB_IP_YOUNG[rownames(dds_C1_MB_IP_YOUNG) %in% genes_MB_IP_YOUNG_GENOTYPE,]
# dds_ONT_MB_IP_YOUNG_GENOTYPE <- dds_ONT_MB_IP_YOUNG[rownames(dds_ONT_MB_IP_YOUNG) %in% genes_MB_IP_YOUNG_GENOTYPE,]
#
# dds_C1_MB_IP_YOUNG_GENOTYPE@design <- ~ collection + genotype
# dds_ONT_MB_IP_YOUNG_GENOTYPE@design <- ~ collection + genotype
#
# colData(dds_C1_MB_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_C1_MB_IP_YOUNG_GENOTYPE))
# colData(dds_ONT_MB_IP_YOUNG_GENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_YOUNG_GENOTYPE))
#
# dds_C1_MB_IP_YOUNG_GENOTYPE <- DESeq(
#   dds_C1_MB_IP_YOUNG_GENOTYPE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# dds_ONT_MB_IP_YOUNG_GENOTYPE <- DESeq(
#   dds_ONT_MB_IP_YOUNG_GENOTYPE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# plotDispEsts(dds_C1_MB_IP_YOUNG_GENOTYPE)
# plotDispEsts(dds_ONT_MB_IP_YOUNG_GENOTYPE)
#
# res_C1_MB_YOUNG_IP_GENOTYPE <- DESeq2::results(
#   dds_C1_MB_IP_YOUNG_GENOTYPE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C1_MB_YOUNG_IP_GENOTYPE <- lfcShrink(
#   dds_C1_MB_IP_YOUNG_GENOTYPE,
#   res = res_C1_MB_YOUNG_IP_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   type = "ashr",
#   parallel = T
# )
#
# res_ONT_MB_YOUNG_IP_GENOTYPE <- DESeq2::results(
#   dds_ONT_MB_IP_YOUNG_GENOTYPE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_YOUNG_IP_GENOTYPE <- lfcShrink(
#   dds_ONT_MB_IP_YOUNG_GENOTYPE,
#   res = res_ONT_MB_YOUNG_IP_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_C1_MB_YOUNG_IP_GENOTYPE)
# summary(res_ONT_MB_YOUNG_IP_GENOTYPE)
# DESeq2::plotMA(res_ONT_MB_YOUNG_IP_GENOTYPE)
#
# MB_GENOTYPE_YOUNG_META <-
#   as_tibble(res_C1_MB_YOUNG_IP_GENOTYPE, rownames = "ensembl_gene_id") %>%
#   left_join(
#     as_tibble(res_ONT_MB_YOUNG_IP_GENOTYPE, rownames = "ensembl_gene_id"),
#     by = "ensembl_gene_id",
#     suffix = c("_C1", "_C3")
#   ) %>%
#   left_join(anno) %>%
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
#   # filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
#   drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
#   rowwise() %>%
#   mutate(
#     pvalue_C3 = ifelse(
#       sign(log2FoldChange_C3) == sign(log2FoldChange_C1),
#       pvalue_C3,
#       1 - pvalue_C3)
#   ) %>%
#   mutate(across(starts_with("pvalue"),
#                 ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
#   mutate(
#     # sumlog = sumlog(c_across(pvalue_C1:pvalue_C3))$p,
#     # calculate sumlog
#     sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
#     # calculate stouffers
#     log2FoldChange = mean(c(
#       log2FoldChange_C1,
#       log2FoldChange_C3
#     ))
#   ) %>%
#   ungroup() %>%
#   mutate(
#     # sumlog_adj = p.adjust(sumlog, method = "fdr"),
#     # correction for multiple comparisons
#     sumz_adj = p.adjust(sumz, method = "fdr"),
#     conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
#   # mutate(score = -log10(sumz_adj) * abs(log2FoldChange)) %>%
#   # arrange(desc(score), desc(log2FoldChange)) %>%
#   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_enriched = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)
#
# # # make a simple MB GENOTYPE OLD META table
# # MB_GENOTYPE_YOUNG_META_SIMPLE <- MB_GENOTYPE_YOUNG_META %>%
# #   select(external_gene_name,
# #          description,
# #          "baseMean" = baseMean_C1,
# #          log2FoldChange,
# #          starts_with("pval"),
# #          sumz_adj,
# #          conflict,
# #          score,
# #          mb_translated,
# #          axon_enriched)
#
# # # filter for genes that are at least pval < 0.1 in every cohort
# # # I am including cohort 1 here because there are genes that agree across cohorts
# # # and it helps to have that power of support
# MB_GENOTYPE_YOUNG_META_DE <- MB_GENOTYPE_YOUNG_META %>%
#   filter(sumz_adj < 0.01 &
#            !conflict &
#            pvalue_C1 < 0.1 &
#            pvalue_C3 < 0.1) %>%
#   arrange(sumz_adj)
#
# # MB GENOTYPE Old OVX vs Old WT subset common genes ----
# # filter genes for common detected
# dds_C1_MB_IP_OLD_FILTER <- filter_genes(dds_C1_MB_IP_OLD,
#                                           grouping = c("genotype", "gene_id"))
#
# # dds_C3_MB_IP_OLD_FILTER <- filter_genes(dds_C3_MB_IP_OLD,
# #                                           grouping = c("genotype", "gene_id"))
#
# dds_ONT_MB_IP_OLD_FILTER <- filter_genes(dds_ONT_MB_IP_OLD,
#                                            grouping = c("genotype", "gene_id"))
#
#
# # take subset: Going with ONT data!
#
# genes_MB_IP_OLD_GENOTYPE <- intersect(dds_C1_MB_IP_OLD_FILTER,
#                                         dds_ONT_MB_IP_OLD_FILTER)
#
# # MB genotype - OLD DESeq2 -----------------------------------------------------
# # This section is actually to produce an MA plot of Old OVX vs Old WT
# # To show that there is a broad expression effect of aged OVX compared to aged WT
# # whereas in young OVX vs WT, there is little effect
#
# dds_C1_MB_IP_OLD_GENOTYPE <- dds_C1_MB_IP_OLD[rownames(dds_C1_MB_IP_OLD) %in% genes_MB_IP_OLD_GENOTYPE,]
# dds_ONT_MB_IP_OLD_GENOTYPE <- dds_ONT_MB_IP_OLD[rownames(dds_ONT_MB_IP_OLD) %in% genes_MB_IP_OLD_GENOTYPE,]
#
# dds_C1_MB_IP_OLD_GENOTYPE@design <- ~ collection + genotype
# dds_ONT_MB_IP_OLD_GENOTYPE@design <- ~ collection + genotype
#
# colData(dds_C1_MB_IP_OLD_GENOTYPE) <- droplevels(colData(dds_C1_MB_IP_OLD_GENOTYPE))
# colData(dds_ONT_MB_IP_OLD_GENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_OLD_GENOTYPE))
#
# dds_C1_MB_IP_OLD_GENOTYPE <- DESeq(
#   dds_C1_MB_IP_OLD_GENOTYPE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# dds_ONT_MB_IP_OLD_GENOTYPE <- DESeq(
#   dds_ONT_MB_IP_OLD_GENOTYPE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# plotDispEsts(dds_C1_MB_IP_OLD_GENOTYPE)
# plotDispEsts(dds_ONT_MB_IP_OLD_GENOTYPE)
#
# res_C1_MB_IP_OLD_GENOTYPE <- DESeq2::results(
#   dds_C1_MB_IP_OLD_GENOTYPE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_C1_MB_IP_OLD_GENOTYPE <- lfcShrink(
#   dds_C1_MB_IP_OLD_GENOTYPE,
#   res = res_C1_MB_IP_OLD_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   type = "ashr",
#   parallel = T
# )
#
# res_ONT_MB_IP_OLD_GENOTYPE <- DESeq2::results(
#   dds_ONT_MB_IP_OLD_GENOTYPE,
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_OLD_GENOTYPE <- lfcShrink(
#   dds_ONT_MB_IP_OLD_GENOTYPE,
#   res = res_ONT_MB_IP_OLD_GENOTYPE,
#   contrast = c("genotype", "OVX", "WT"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_C1_MB_IP_OLD_GENOTYPE)
# summary(res_ONT_MB_IP_OLD_GENOTYPE)
# DESeq2::plotMA(res_ONT_MB_IP_OLD_GENOTYPE)
#
# MB_GENOTYPE_OLD_META <-
#   as_tibble(res_C1_MB_IP_OLD_GENOTYPE, rownames = "ensembl_gene_id") %>%
#   left_join(
#     as_tibble(res_ONT_MB_IP_OLD_GENOTYPE, rownames = "ensembl_gene_id"),
#     by = "ensembl_gene_id",
#     suffix = c("_C1", "_C3")
#   ) %>%
#   left_join(anno) %>%
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
#   # filter(external_gene_name != "1") %>% # remove genes with an external_gene_name of "1"
#   drop_na() %>% # remove genes where they are excluded in cohort 3 due to low counts
#   rowwise() %>%
#   mutate(
#     pvalue_C3 = ifelse(
#       sign(log2FoldChange_C3) == sign(log2FoldChange_C1),
#       pvalue_C3,
#       1 - pvalue_C3)
#   ) %>%
#   mutate(across(starts_with("pvalue"),
#                 ~ ifelse(.x == 1, .x - 1e-10, .x))) %>% # replace 1 pvalues with almost 1
#   mutate(
#     # sumlog = sumlog(c_across(pvalue_C1:pvalue_C3))$p,
#     # calculate sumlog
#     sumz = sumz(c_across(pvalue_C1:pvalue_C3))$p,
#     # calculate stouffers
#     log2FoldChange = mean(c(
#       log2FoldChange_C1,
#       log2FoldChange_C3
#     ))
#   ) %>%
#   ungroup() %>%
#   mutate(
#     # sumlog_adj = p.adjust(sumlog, method = "fdr"),
#     # correction for multiple comparisons
#     sumz_adj = p.adjust(sumz, method = "fdr"),
#     conflict = sign(log2FoldChange_C1) != sign(log2FoldChange_C3)) %>% # state whether there is a conflict in l2fc direction
#   # mutate(score = -log10(sumz_adj) * abs(log2FoldChange)) %>%
#   # arrange(desc(score), desc(log2FoldChange)) %>%
#   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#          axon_enriched = ensembl_gene_id %in% AXON_FRACTION_ENRICHED_AGNOSTIC_GENES)
#
# # # make a simple MB GENOTYPE OLD META table
# # MB_GENOTYPE_OLD_META_SIMPLE <- MB_GENOTYPE_OLD_META %>%
# #   select(external_gene_name,
# #          description,
# #          "baseMean" = baseMean_C1,
# #          log2FoldChange,
# #          starts_with("pval"),
# #          sumz_adj,
# #          conflict,
# #          score,
# #          mb_translated,
# #          axon_enriched)
#
# # filter for genes that are at least pval < 0.1 in every cohort
# MB_GENOTYPE_OLD_META_DE <- MB_GENOTYPE_OLD_META %>%
#   filter(sumz_adj < 0.01 &
#            !conflict &
#            pvalue_C1 < 0.1 &
#            pvalue_C3 < 0.1) %>%
#   arrange(sumz_adj)
#
# n_MB_GENOTYPE_OLD_NGENES <- MB_GENOTYPE_OLD_META_DE %>% nrow
#
#
#
# # FIGURE A: OVX YOUNG and OLD MA ----
# plot_data <- MB_GENOTYPE_OLD_META %>%
#   mutate(signif = external_gene_name %in% MB_GENOTYPE_OLD_META_DE$external_gene_name) %>%
#   arrange(signif) %>%
#   mutate(age = "OLD") %>%
#   bind_rows({
#     MB_GENOTYPE_YOUNG_META %>%
#       mutate(signif = external_gene_name %in% MB_GENOTYPE_YOUNG_META_DE$external_gene_name) %>%
#       arrange(signif)  %>%
#       mutate(age = "YOUNG")
#   })
#
# labels <- plot_data %>%
#   filter(ifelse(age == "OLD",
#                 external_gene_name %in% MB_GENOTYPE_OLD_META_DE$external_gene_name,
#                 external_gene_name %in% MB_GENOTYPE_YOUNG_META_DE$external_gene_name)) %>%
#   select(external_gene_name,
#          baseMean_C1,
#          log2FoldChange,
#          mb_translated,
#          age) %>%
#   mutate(fontface = ifelse(mb_translated, "bold", "plain"),
#          label_colour = mb_translated,
#          external_gene_name = ifelse(external_gene_name == "Snca", "SNCA", external_gene_name))
#
# p_mb_genotype_ma <- plot_data %>%
#   # mutate(age = factor(age, levels = c("old", "young"))) %>%
#   # arrange(age) %>%
#   ggplot(aes(x = baseMean_C1,
#              y = log2FoldChange)) +
#   geom_point(aes(fill = signif),
#              shape = 21,
#              colour = "black") +
#   geom_point(data = {plot_data %>% filter(signif)},
#              size = 2,
#              colour = pal_d3()(3)[3]) +
#   geom_label_repel(data = labels,
#                    aes(label = external_gene_name,
#                        fontface = fontface,
#                        colour = label_colour)) +
#   scale_x_log10() +
#   scale_y_continuous(limits = c(-1, 1),
#                      oob = scales::squish) +
#   scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
#   scale_color_manual(values = c("black", pal_d3()(3)[3])) +
#   theme(legend.position = "none") +
#   labs(x = "Expression Counts",
#        y = "Log2 Fold Change"
#        ) +
#   facet_wrap(vars(factor(age, levels = c("YOUNG", "OLD"))),
#              nrow = 1) +
#   panel_border()
# # ----
# # ----
# # ----
# # SNCA OVX: Transcripts ----
#
# # get transcript external gene names
# anno_SNCA <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "ensembl_transcript_id",
#                  "external_transcript_name"),
#   filters = "external_gene_name",
#   values = "SNCA",
#   mart = ensembl_human,
#   uniqueRows = TRUE
# )
#
# anno_Snca <- getBM(
#   attributes = c("ensembl_gene_id",
#                  "ensembl_transcript_id",
#                  "external_transcript_name"),
#   filters = "external_gene_name",
#   values = "Snca",
#   mart = ensembl,
#   uniqueRows = TRUE
# )
#
# anno_SNCA_joined <- bind_rows(anno_SNCA, anno_Snca)
#
# # TRAP transcripts
# files <-
#   list.files(
#     "/zfs/analysis/trap/active/testing_env/snca/input_lr",
#     pattern = "quant.sf",
#     recursive =  T,
#     full.names = T
#   )
#
# all(file.exists(files)) # everything exists
#
# names(files) <-
#   str_extract(files, "(?<=input_lr\\/)[:alnum:]+(?=\\/quant.sf)")
#
# snca_trap <- lapply(files, function(x){
#   read_delim(x,
#              delim = "\t") %>%
#     mutate(ensembl_transcript_id = str_extract(Name, "[:alnum:]+[\\.][:digit:]+(?=\\|)")) %>%
#     select(ensembl_transcript_id,
#            everything(),
#            -Name)
# }) %>%
#   bind_rows(.id = "barcode") %>%
#   left_join(colData(dds_ONT),
#             copy = T) %>%
#   mutate(ensembl_transcript_id = substr(ensembl_transcript_id, 1, 15)) %>%
#   filter(ensembl_transcript_id %in% anno_SNCA_joined$ensembl_transcript_id) %>%
#   mutate(species = ifelse(str_detect(ensembl_transcript_id, "ENST"),
#                           "SNCA", "Snca"),
#          genotype = factor(ifelse(genotype == "WT", "WT", "SNCA-OVX"),
#                            levels = c("WT", "SNCA-OVX")))
#
# # set the transcript order
# snca_order <- snca_trap %>%
#   filter(barcode == "barcode01") %>%
#   pull(ensembl_transcript_id)
#
# # plot the OVX effect gene-level
# p_ovx_gene_level <- snca_trap %>%
#   group_by(sample_name, genotype, species) %>%
#   summarise(reads = sum(NumReads)) %>%
#   ggplot(aes(x = genotype,
#              y = reads,
#              fill = species)) +
#   geom_quasirandom(shape = 21,
#                    size = 2) +
#   facet_wrap(vars(species)) +
#   panel_border() +
#   scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
#   labs(x = "Genotype",
#        y = "Counts") +
#   theme(legend.position = "none")
#
# # p_ovx_transcript_level <- snca_trap %>%
# #   group_by(species, sample_name) %>%
# #   mutate(ensembl_transcript_id = factor(paste0("T", row_number()),
# #                                         levels = paste0("T", row_number()))) %>%
# # ggplot(aes(x = ensembl_transcript_id,
# #              y = TPM,
# #            fill = species)) +
# #   geom_quasirandom(shape = 21,
# #                    size = 2) +
# #   facet_grid(rows = vars(genotype),
# #              cols = vars(species),
# #              scales = "free_x",
# #              space = "free_x") +
# #   panel_border() +
# #   labs(x = "Transcript",
# #        y = "Counts") +
# #   scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
# #   theme(legend.position = "none") +
# #   theme(axis.text.x = element_text(angle = 90,
# #                                    vjust = 0.5))
#
# # compare with human counts from ips neurons
#
# SNCA_snca_transcripts <- c(anno_SNCA$ensembl_transcript_id,
#                            anno_Snca$ensembl_transcript_id)
#
# files <-
#   list.files(
#     "/zfs/analysis/rwm_transcriptomics/data/combined_kallisto/final_counts_20200207",
#     pattern = "abundance.tsv",
#     recursive =  T,
#     full.names = T
#   )
#
# names(files) <-
#   str_extract(files, "(?<=final_counts_20200207\\/).+(?=\\/abundance.tsv)")
#
# snca_ips <- lapply(files, function(x){
#   read_delim(x,
#              delim = "\t") %>%
#     mutate(ensembl_transcript_id = substr(target_id, 1, 15)) %>%
#     filter(ensembl_transcript_id %in% SNCA_snca_transcripts) %>%
#     select(ensembl_transcript_id,
#            everything(),
#            -target_id)
# }) %>%
#   bind_rows(.id = "sample_name") %>%
#   mutate(species = "SNCA")
#
# # compare with human counts from ips neurons - nicheterwitz
#
# files <-
#   list.files(
#     "/zfs/analysis/bm_rna/2021/nichterwitz",
#     pattern = "quant.sf",
#     recursive =  T,
#     full.names = T
#   )
#
# names(files) <-
#   str_extract(files, "(?<=nichterwitz\\/).+(?=\\/quant.sf)")
#
# snca_lcm <- lapply(files, function(x){
#   read_delim(x,
#              delim = "\t") %>%
#     mutate(ensembl_transcript_id = str_extract(Name, "[:alnum:]+(?=[\\.][:digit:]+)")) %>%
#     select(ensembl_transcript_id,
#            everything(),
#            -Name) %>%
#     filter(ensembl_transcript_id %in% SNCA_snca_transcripts)
# })  %>%
#   bind_rows(.id = "sample_name") %>%
#   mutate(species = "SNCA")
#
# # compare with human counts from ips neurons - parkinnen
#
# files <-
#   list.files(
#     "/zfs/analysis/bm_rna/2021/parkinnen/",
#     pattern = "abundance.tsv",
#     recursive =  T,
#     full.names = T
#   )
#
# names(files) <-
#   str_extract(files, "(?<=parkinnen\\/).+(?=\\/abundance.tsv)")
#
# snca_lcm <- lapply(files, function(x){
#   read_delim(x,
#              delim = "\t") %>%
#     mutate(ensembl_transcript_id = str_extract(target_id, "[:alnum:]+(?=[\\.][:digit:]+)")) %>%
#     select(ensembl_transcript_id,
#            everything(),
#            -target_id) %>%
#     filter(ensembl_transcript_id %in% SNCA_snca_transcripts)
# })  %>%
#   bind_rows(.id = "sample_name") %>%
#   mutate(species = "SNCA")
#
#
#
# # plot tpm per transcript
# p_snca_wt_ovx_ips_lcm <- bind_rows(
#   {
#     snca_trap %>%
#       select(ensembl_transcript_id,
#              species,
#              sample_name,
#              genotype,
#              counts = "TPM")
#   }, {
#     snca_ips[order(factor(snca_ips$ensembl_transcript_id, levels = snca_order)),] %>%
#       mutate(species = "SNCA",
#              genotype = "Human iPSC\nDopamine Neurons") %>%
#       select(ensembl_transcript_id,
#              species,
#              sample_name,
#              genotype,
#              counts = "tpm")
#   }, {
#     snca_lcm[order(factor(snca_lcm$ensembl_transcript_id, levels = snca_order)),] %>%
#       mutate(species = "SNCA",
#              genotype = "Human LCM\nDopamine Neurons") %>%
#       select(ensembl_transcript_id,
#              species,
#              sample_name,
#              genotype,
#              counts = "tpm")
#   }
# ) %>%
#   arrange(sample_name) %>%
#   group_by(species, sample_name) %>%
#   left_join(anno_SNCA_joined) %>%
#   mutate(ensembl_transcript_id = factor(paste0("T", row_number()),
#                                         levels = paste0("T", row_number()))) %>%
#   mutate(genotype = factor(genotype,
#                            levels = c("WT",
#                                       "SNCA-OVX",
#                                       "Human iPSC\nDopamine Neurons",
#                                       "Human LCM\nDopamine Neurons"))) %>%
#   # filter(counts > 20) %>%
#   ggplot(aes(x = external_transcript_name,
#              # y = log2(counts+1),
#              y = counts,
#              fill = species)) +
#   geom_quasirandom(
#     # shape = 21,
#                    size = 1
#     ) +
#   facet_grid(rows = vars(genotype),
#              cols = vars(species),
#              scales = "free_x",
#              space = "free_x") +
#   panel_border() +
#   labs(x = "Transcript",
#        y = expression(Log[2] ~ Expression ~ Count)) +
#   scale_fill_manual(values = pal_d3()(3)[c(1, 3)]) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90,
#                                    vjust = 0.5)) +
#   scale_y_continuous(limits = c(0, 500),
#                      oob = scales::squish)
#
#
#
#
# # ----
# # ----
# # ----
# # MB GENOTYPE OLD OVX vs everything else - subset common genes ----
# # filter genes for common detected
# # Can't compare OLD OVX to YOUNG WT in C1 without encountering batch issue!
# # dds_C1_MB_IP_FILTER <- filter_genes(dds_C1_MB_IP,
# #                                         grouping = c("genotype", "age", "gene_id"))
#
# dds_C3_MB_IP_AGEGENOTYPE_FILTER <- filter_genes(dds_C3_MB_IP,
#                                     grouping = c("genotype", "age", "gene_id"))
#
# dds_ONT_MB_IP_AGEGENOTYPE_FILTER <- filter_genes(dds_ONT,
#                                      grouping = c("genotype", "age", "gene_id"))
#
# # take subset: Going with ONT data!
# dds_ONT_MB_IP_AGEGENOTYPE <- dds_ONT[rownames(dds_ONT) %in% dds_ONT_MB_IP_AGEGENOTYPE_FILTER,]
#
# # create new group
# colData(dds_ONT_MB_IP_AGEGENOTYPE)$age_genotype <- factor(paste(colData(dds_ONT_MB_IP_AGEGENOTYPE)$age,
#                                                                          colData(dds_ONT_MB_IP_AGEGENOTYPE)$genotype,
#                                                                          sep = "_"),
#                                                                    levels = c("YOUNG_WT", "YOUNG_OVX", "OLD_WT", "OLD_OVX"))
#
# dds_ONT_MB_IP_AGEGENOTYPE@design <- ~ collection + age_genotype
#
# colData(dds_ONT_MB_IP_AGEGENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_AGEGENOTYPE))
#
# dds_ONT_MB_IP_AGEGENOTYPE <- DESeq(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   # sfType = "poscounts",
#   # fitType = "local",
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# plotDispEsts(dds_ONT_MB_IP_AGEGENOTYPE)
#
# resultsNames(dds_ONT_MB_IP_AGEGENOTYPE)
#
# # YOUNG OVX vs YOUNG WT - Just Snca ----
# res_ONT_MB_IP_AGEGENOTYPE_YOUNGOVX_VS_YOUNGWT <- DESeq2::results(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   contrast = c("age_genotype", "YOUNG_OVX", "YOUNG_WT"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_AGEGENOTYPE_YOUNGOVX_VS_YOUNGWT <- lfcShrink(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   res = res_ONT_MB_IP_AGEGENOTYPE_YOUNGOVX_VS_YOUNGWT,
#   contrast = c("age_genotype", "YOUNG_OVX", "YOUNG_WT"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_ONT_MB_IP_AGEGENOTYPE_YOUNGOVX_VS_YOUNGWT)
#
# # res_ONT_MB_IP_AGEGENOTYPE_YOUNGOVX_VS_YOUNGWT %>%
# #   as_tibble(rownames = "ensembl_gene_id") %>%
# #   filter(padj < 0.01) %>%
# #   left_join(anno) %>%
# #   mutate(mb_translated = ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
# #          mb_age = ensembl_gene_id %in% MB_AGE_DE_GENES) %>% View
#
# # OLD WT/OLD OVX vs YOUNG ALL ----
#
# # There is no major difference between OVX and WT in YOUNG samples, so for OLD OVX/WT vs YOUNG comparisons,
# # merging YOUNG WT AND YOUNG OVX into just YOUNG category
#
# colData(dds_ONT_MB_IP_AGEGENOTYPE)$age_genotype2 <- factor(ifelse(colData(dds_ONT_MB_IP_AGEGENOTYPE)$age_genotype %in% c("YOUNG_WT", "YOUNG_OVX"),
#                                                                            "YOUNG",
#                                                                            as.character(colData(dds_ONT_MB_IP_AGEGENOTYPE)$age_genotype)),
#                                                                     levels = c("YOUNG", "OLD_WT", "OLD_OVX"))
#
# # recalculate dispersion
# dds_ONT_MB_IP_AGEGENOTYPE@design <- ~ collection + age_genotype2
#
# colData(dds_ONT_MB_IP_AGEGENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_AGEGENOTYPE))
#
# dds_ONT_MB_IP_AGEGENOTYPE <- DESeq(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   # sfType = "poscounts",
#   # fitType = "local",
#   minReplicatesForReplace = Inf,
#   parallel = T
# )
#
# # OLD WT vs YOUNG ALL
# res_ONT_MB_IP_AGEGENOTYPE_OLDWT_VS_YOUNGALL <- DESeq2::results(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   contrast = c("age_genotype2", "OLD_WT", "YOUNG"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_AGEGENOTYPE_OLDWT_VS_YOUNGALL <- lfcShrink(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   res = res_ONT_MB_IP_AGEGENOTYPE_OLDWT_VS_YOUNGALL,
#   contrast = c("age_genotype2", "OLD_WT", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_ONT_MB_IP_AGEGENOTYPE_OLDWT_VS_YOUNGALL)
# # DESeq2::plotMA(res_ONT_MB_IP_AGEGENOTYPE_OLDWT_VS_YOUNGALL)
#
# # Build list of ONT OLDWT vs YOUNGALL genes: Have to decide pvalue cutoff
# # going with 0.1 for now
# MB_ONT_OLDWT_VS_YOUNGALL_GENES <- res_ONT_MB_IP_AGEGENOTYPE_OLDWT_VS_YOUNGALL %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(padj < 0.1) %>%
#   pull(ensembl_gene_id)
#
# # OLD OVX vs YOUNG ALL ----
# res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL <- DESeq2::results(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   contrast = c("age_genotype2", "OLD_OVX", "YOUNG"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL <- lfcShrink(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   res = res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL,
#   contrast = c("age_genotype2", "OLD_OVX", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL)
# # DESeq2::plotMA(res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL)
#
# # Filter OLD OVX vs YOUNG ALL by ONT-specific ONT OLD OVX vs OLD WT ----
# # Using the subsetted dds_ONT object - important, because with the full dds_ONT object,
# # YOUNG OVX vs YOUNG WT dispersion gets in the way (e.g Snca is DE between YOUNG OVX and YOUNG WT)
#
# # There are a lot of DEGs at padj < 0.1, but only Snca at padj < 0.01
# # So the OLD OVX vs YOUNG ALL comparison is the most interesting
# # Need to find the component of these genes that are involved in ageing and that are not
# # And whether the ageing ones are worsened by OVX ageing.
#
# pval_ONT_MB_IP_OLD_GENOTYPE <- res_ONT_MB_IP_OLD_GENOTYPE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   select(ensembl_gene_id,
#          pvalue_ONT_OLDOVX_VS_OLDWT = pvalue)
#
# # AGE: ONT OLD VS YOUNG -----
#
# dds_ONT_MB_IP_AGEGENOTYPE@design <- ~ collection + age
#
# colData(dds_ONT_MB_IP_AGEGENOTYPE) <- droplevels(colData(dds_ONT_MB_IP_AGEGENOTYPE))
#
# dds_ONT_MB_IP_AGEGENOTYPE <- DESeq(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   minReplicatesForReplace = Inf,
#   parallel = TRUE
# )
#
# res_ONT_MB_IP_AGE <- DESeq2::results(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   contrast = c("age", "OLD", "YOUNG"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_AGE <- lfcShrink(
#   dds_ONT_MB_IP_AGEGENOTYPE,
#   res = res_ONT_MB_IP_AGE,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
#
# summary(res_ONT_MB_IP_AGE)
#
# # Get the list of ageing DEGs from all OLD ONT samples vs ALL YOUNG ONT samples
# MB_ONT_AGE_GENES <- res_ONT_MB_IP_AGE %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(padj < 0.1) %>%
#   pull(ensembl_gene_id)
#
#
# # OLD OVX VS YOUNG ALL CONTINUED ----
# # Rank genes and filter by OLD OVX vs OLD WT difference to leave
# # Just genes that in OLD OVX samples are different to ALL YOUNG and
# # also OLD WT samples
# MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL <- res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(padj < 0.01) %>%
#   left_join(anno) %>%
#   mutate(mb_translated = ifelse(ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#                               "Translated",
#                               ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                      "Depleted",
#                                      "Unchanged")),
#          mb_age = ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_META$ensembl_gene_id,
#          ont_wt_age = ensembl_gene_id %in% MB_ONT_OLDWT_VS_YOUNGALL_GENES,
#          ont_age = ensembl_gene_id %in% MB_ONT_AGE_GENES,
#          outcome = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
#   left_join(pval_ONT_MB_IP_OLD_GENOTYPE)
#   # filter(
#   #   # mb_translated == "Translated" &
#   #   #        !mb_age &
#   #   #        !ont_wt_age
#   #   pvalue_ONT_OLDOVX_VS_OLDWT < 0.2
#   # )
#
#
# # ----
# # ----
# # ----
# # Compare the lfc of ageing genes between OVX and WT ----
#
#
# # Aged WT vs Young ALL
# # subset samples
# dds_ONT_MB_IP_AGE_NO_OLDOVX <- dds_ONT[,colData(dds_ONT)$age == "YOUNG" |
#                                          colData(dds_ONT)$age == "OLD" &
#                                          colData(dds_ONT)$genotype == "WT"]
# # process results
# dds_ONT_MB_IP_AGE_NO_OLDOVX@design <- ~ collection + age
# colData(dds_ONT_MB_IP_AGE_NO_OLDOVX) <- droplevels(colData(dds_ONT_MB_IP_AGE_NO_OLDOVX))
# dds_ONT_MB_IP_AGE_NO_OLDOVX <- DESeq(
#   dds_ONT_MB_IP_AGE_NO_OLDOVX,
#   minReplicatesForReplace = Inf,
#   parallel = T
# )
# res_ONT_MB_IP_AGE_NO_OLDOVX <- DESeq2::results(
#   dds_ONT_MB_IP_AGE_NO_OLDOVX,
#   contrast = c("age", "OLD", "YOUNG"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_AGE_NO_OLDOVX <- lfcShrink(
#   dds_ONT_MB_IP_AGE_NO_OLDOVX,
#   res = res_ONT_MB_IP_AGE_NO_OLDOVX,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
# # plot enrichment of results
# res_ONT_MB_IP_AGE_NO_OLDOVX %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   filter(padj < 0.01) %>%
#   left_join(enrichment_status) %>%
#   mutate(outcome = ifelse(log2FoldChange < 0, "Down", "Up")) %>%
#   ggplot(aes(x = outcome,
#              fill = enrichment)) +
#   geom_bar(colour = "black") +
#   scale_fill_d3()
#
#
# # Aged OVX vs Young ALL
# # subset samples
# dds_ONT_MB_IP_AGE_OLDOVX <- dds_ONT[,colData(dds_ONT)$age == "YOUNG" |
#                                          colData(dds_ONT)$age == "OLD" &
#                                          colData(dds_ONT)$genotype == "OVX"]
# # process results
# dds_ONT_MB_IP_AGE_OLDOVX@design <- ~ collection + age
# colData(dds_ONT_MB_IP_AGE_OLDOVX) <- droplevels(colData(dds_ONT_MB_IP_AGE_OLDOVX))
# dds_ONT_MB_IP_AGE_OLDOVX <- DESeq(
#   dds_ONT_MB_IP_AGE_OLDOVX,
#   minReplicatesForReplace = Inf,
#   parallel = T
# )
# res_ONT_MB_IP_AGE_OLDOVX <- DESeq2::results(
#   dds_ONT_MB_IP_AGE_OLDOVX,
#   contrast = c("age", "OLD", "YOUNG"),
#   cooksCutoff = Inf,
#   filterFun = ihw,
#   parallel = T
# )
# res_ONT_MB_IP_AGE_OLDOVX <- lfcShrink(
#   dds_ONT_MB_IP_AGE_OLDOVX,
#   res = res_ONT_MB_IP_AGE_OLDOVX,
#   contrast = c("age", "OLD", "YOUNG"),
#   type = "ashr",
#   parallel = T
# )
# # plot enrichment of results
# res_ONT_MB_IP_AGE_OLDOVX %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   filter(padj < 0.01) %>%
#   left_join(enrichment_status) %>%
#   mutate(outcome = ifelse(log2FoldChange < 0, "Down", "Up")) %>%
#   ggplot(aes(x = outcome,
#              fill = enrichment)) +
#   geom_bar(colour = "black") +
#   scale_fill_d3()
#
# # combine the results of ageing from both genotypes
# res_ONT_AGE_GENOTYPE_SPLIT <- bind_rows({
#   res_ONT_MB_IP_AGE_NO_OLDOVX %>%
#     as_tibble(rownames = "ensembl_gene_id") %>%
#     left_join(anno) %>%
#     filter(padj < 0.1) %>%
#     left_join(enrichment_status) %>%
#     mutate(outcome = ifelse(log2FoldChange < 0, "Down", "Up")) %>%
#     mutate(genotype = "WT")
# }, {
#   res_ONT_MB_IP_AGE_OLDOVX %>%
#     as_tibble(rownames = "ensembl_gene_id") %>%
#     left_join(anno) %>%
#     filter(padj < 0.1) %>%
#     left_join(enrichment_status) %>%
#     mutate(outcome = ifelse(log2FoldChange < 0, "Down", "Up")) %>%
#     mutate(genotype = "OVX")
# })
#
# p_age_genotype_outcome <- res_ONT_AGE_GENOTYPE_SPLIT %>%
#   ggplot(aes(x = outcome)) +
#   geom_bar(aes(fill = enrichment),
#            colour = "black") +
#   scale_fill_d3() +
#   facet_wrap(vars(genotype)) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   theme(legend.position = c(0.8, 0.8)) +
#   geom_text(data = . %>%
#               group_by(outcome, genotype) %>%
#               tally(),
#             aes(label = paste("Total:", n),
#                 y = n),
#             vjust = -0.5) +
#   labs(x = "Change in Age",
#        y = "Number of Genes",
#        fill = "Enrichment")
#
#
# p_age_genotype_lfc <- res_ONT_AGE_GENOTYPE_SPLIT %>%
#   filter(ensembl_gene_id %in% {res_ONT_MB_IP_AGE_NO_OLDOVX %>%
#       as_tibble(rownames = "ensembl_gene_id") %>%
#       left_join(anno) %>%
#       filter(padj < 0.1) %>%
#       pull(ensembl_gene_id)}) %>%
#   group_by(ensembl_gene_id) %>%
#   filter(dplyr::n() == 2) %>%
#   mutate(score = -log10(padj) * log2FoldChange,
#          genotype = ifelse(genotype == "OVX", "SNCA-OVX", "WT")) %>%
#   ggplot(aes(x = genotype,
#              y = log2FoldChange)) +
#   # geom_point() +
#   geom_quasirandom() +
#   geom_hline(yintercept = 0,
#              linetype = "dotted") +
#   labs(x = "Genotype",
#        y = "Log2 Fold Change in Age")
#
#
# summary(res_ONT_MB_IP_OLD_GENOTYPE, alpha = 0.1)
#
# summary(res_C1_MB_IP_OLD_GENOTYPE, alpha = 0.1)
#
# summary(res_ONT_MB_IP_OLD_GENOTYPE) %>%
#   filter(enrichment == "Soma and Axon" &
#            padj < 0.1 &
#            log2FoldChange > 0) %>%
#   mutate(ageing = ensembl_gene_id %in% MB_AGE_DE_SPECIFIC_UP$ensembl_gene_id) %>% View
#
#
# res_ONT_MB_IP_AGEGENOTYPE_YOUNGOVX_VS_YOUNGWT %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   left_join(anno) %>%
#   left_join(enrichment_status) %>%
#   filter(padj < 0.1) %>%
#   filter(log2FoldChange > 0) %>%
#   # filter(enrichment == "Soma only" &
#   #          log2FoldChange > 0) %>%
#   View
#
#
# # ----
# # ----
# # ----
# # TEST ----
#
#   MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
#     ggplot(aes(x = mb_translated)) +
#     geom_bar()
#
#   MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
#     filter(mb_translated == "Translated") %>%
#     ggplot(aes(x = outcome)) +
#     geom_bar()
#
#   MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
#     filter(mb_translated == "Translated" &
#              outcome == "Down") %>% View
#
#   MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
#     filter(mb_translated == "Translated" &
#              outcome == "Up") %>%
#     left_join(publications_all_genes) %>%
#     left_join(publications_all_genes_axon) %>% View
#
# # FIGURE B: OLD OVX VS YOUNG ALL ONT: Enrichment breakdown ----
#
#
#
# plot_data <- res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
#   as_tibble(rownames = "ensembl_gene_id") %>%
#   filter(ensembl_gene_id %in% MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL$ensembl_gene_id) %>%
#   left_join(anno) %>%
#   mutate(mb_translated = ifelse(ensembl_gene_id %in% MB_FRACTION_TRANSLATED_GENES,
#                               "Translated",
#                               ifelse(ensembl_gene_id %in% MB_FRACTION_DEPLETED_GENES,
#                                      "Depleted",
#                                      "Unchanged")),
#          # mb_age = ensembl_gene_id %in% MB_AGE_DE_GENES,
#          ont_wt_age = ensembl_gene_id %in% MB_ONT_OLDWT_VS_YOUNGALL_GENES,
#          ont_age = ensembl_gene_id %in% MB_ONT_AGE_GENES)
#
# p_mb_genotype_oldovx_vs_youngall_enrichment <- plot_data %>%
#   filter(padj < 0.01) %>%
#   ggplot(aes(x = mb_translated,
#              fill = mb_translated)) +
#   geom_bar(colour = "black") +
#   scale_fill_d3() +
#   labs(x = "Somal Enrichment",
#        y = "Number of Genes") +
#   theme(legend.position = "none",
#         panel.grid.major.x = element_blank()) +
#   geom_text(aes( label = ..count..,
#                  y= ..count.. ), stat= "count", vjust = -.5)
#
# # FIGURE C: OLD OVX VS YOUNG ALL ONT: Ageing breakdown ----
#
# p_mb_genotype_oldovx_vs_youngall_ageing <- plot_data %>%
#   mutate(ageing = ifelse(ont_wt_age | ont_age,
#                          "Differentially Expressed", "Unchanged")) %>%
#   ggplot(aes(x = ageing,
#              fill = ageing)) +
#   geom_bar(colour = "black") +
#   scale_fill_d3() +
#   labs(x = "OVX Genes in Ageing",
#        y = "Number of Genes") +
#   theme(legend.position = "none",
#         panel.grid.major.x = element_blank()) +
#   geom_text(aes( label = ..count..,
#                  y= ..count.. ), stat= "count", vjust = -.5)
#
# # FIGURE D: OLD OVX vs YOUNG ALL: MB translated example ----
#
#
# p_mb_agegenotype_mb_translated_example <- view_counts(dds_ONT,
#             "Pfkp") %>%
#   left_join(colData(dds_ONT),
#             copy = T) %>%
#   ggplot(aes(x = genotype,
#              y = count,
#              fill = genotype)) +
#   geom_quasirandom(colour = "black",
#                    size = 4,
#                    shape = 21) +
#   scale_fill_d3() +
#   facet_wrap(vars(age)) +
#   labs(x = "Genotype",
#        y = "Expression Count") +
#   theme(legend.position = "none",
#         panel.grid.major.x = element_blank())
#
# # FIGURE E: OLD OVX vs YOUNG ALL: Non-enriched example ----
#
# p_mb_agegenotype_mb_unchanged_example <- view_counts(dds_ONT,
#                                                     "Apod") %>%
#   left_join(colData(dds_ONT),
#             copy = T) %>%
#   ggplot(aes(x = genotype,
#              y = count,
#              fill = genotype)) +
#   geom_quasirandom(colour = "black",
#                    size = 4,
#                    shape = 21) +
#   scale_fill_d3() +
#   facet_wrap(vars(age)) +
#   labs(x = "Genotype",
#        y = "Expression Count") +
#   theme(legend.position = "none",
#         panel.grid.major.x = element_blank())
#
#
#
#
# # MB OLD OVX vs YOUNG ALL: GO ----
# # temp <- MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
# #   filter(mb_translated != "Translated") %>%
# #   arrange(padj) %>%
# #   pull(ensembl_gene_id) %>%
# #   gost(organism = "mmusculus",
# #        ordered_query = T)
# # temp$result %>%
# #   arrange(p_value) %>%
# #   filter(source != "TF") %>% View
#
# # MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
# #   filter(mb_translated != "Translated") %>%
# #   arrange(padj) %>% View
#
# # view_counts(dds_ONT,
# #             "Tmem176b") %>%
# #   left_join(colData(dds_ONT),
# #             copy = T) %>%
# #   ggplot(aes(x = genotype,
# #              y = count,
# #              fill = genotype)) +
# #   geom_quasirandom(colour = "black",
# #                    size = 4,
# #                    shape = 21) +
# #   scale_fill_d3() +
# #   facet_wrap(vars(age)) +
# #   labs(x = "Genotype",
# #        y = "Expression Count") +
# #   theme(legend.position = "none",
# #         panel.grid.major.x = element_blank())
# # View Counts ----
#
#
# # view_counts(dds_ONT,
# #             "Apod") %>%
# #   left_join(colData(dds_ONT),
# #             copy = T) %>%
# #   ggplot(aes(x = genotype,
# #              y = count,
# #              fill = genotype)) +
# #   geom_quasirandom(colour = "black",
# #                    size = 4,
# #                    shape = 21) +
# #   scale_fill_d3() +
# #   facet_wrap(vars(age))
# #
# # view_counts(dds_C1_MB_IP,
# #             "Apod") %>%
# #   left_join(colData(dds_C1_MB_IP),
# #             copy = T) %>%
# #   ggplot(aes(x = genotype,
# #              y = count,
# #              fill = genotype)) +
# #   geom_quasirandom(colour = "black",
# #                    size = 4,
# #                    shape = 21) +
# #   scale_fill_d3() +
# #   facet_wrap(vars(age))
# #
# # view_counts(dds_C3_MB_IP,
# #             "Apod") %>%
# #   left_join(colData(dds_C3_MB_IP),
# #             copy = T) %>%
# #   ggplot(aes(x = genotype,
# #              y = count,
# #              fill = genotype)) +
# #   geom_quasirandom(colour = "black",
# #                    size = 4,
# #                    shape = 21) +
# #   scale_fill_d3() +
# #   facet_wrap(vars(age))
#
#
# # ----
# # ----
# # ----
# #
# # ----
# # ----
# # ----
# # MB GENOTYPE - fgsea ----
#
#
# # barplot(sort(res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL$log2FoldChange, decreasing = T))
# #
# # res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
# #   as_tibble() %>%
# #   arrange(log2FoldChange) %>%
# #   mutate(rank = row_number()) %>%
# #   ggplot(aes(x = rank,
# #              y = log2FoldChange)) +
# #   geom_point(colour = "black")
# #
# # # get entrez IDs
# # anno_fgsea <- getBM(
# #   attributes = c("ensembl_gene_id",
# #                  "entrezgene_id"),
# #   filters = "ensembl_gene_id",
# #   values = rownames(res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL),
# #   mart = ensembl,
# #   uniqueRows = TRUE
# # )
# #
# # # # load gmt
# # #
# # pathways <- fgsea::gmtPathways("R/objects/msigdb.v7.4.symbols.gmt.txt")
# # #
# # pathways <- fgsea::gmtPathways("R/objects/c2.cp.kegg.v7.4.symbols.gmt.txt")
# # # pathways <- fgsea::gmtPathways("R/objects/c5.go.v7.4.symbols.gmt.txt")
# # # pathways <- fgsea::gmtPathways("R/objects/c5.go.bp.v7.4.symbols.gmt.txt")
# # # pathways <- fgsea::gmtPathways("R/objects/c5.go.mf.v7.4.symbols.gmt.txt")
# #
# # ranks <- res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
# #   as_tibble(rownames = "ensembl_gene_id") %>%
# #   mutate(score = -log10(padj) * log2FoldChange)
# # ranks <- res_ONT_MB_IP_AGEGENOTYPE_OLDOVX_VS_YOUNGALL %>%
# #   as_tibble(rownames = "ensembl_gene_id") %>%
# #   mutate(rank = log2FoldChange) %>%
# #   # filter(!conflict) %>%
# #   inner_join(anno_human) %>%
# #   select(hsapiens_homolog_associated_gene_name,
# #          rank) %>%
# #   filter(hsapiens_homolog_associated_gene_name != "") %>%
# #   distinct(across(-rank), .keep_all = T) %>%
# #   deframe
# #
# # names(ranks)
# #
# # MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL_FGSEA <- fgsea(pathways,
# #                            ranks,
# #                            eps = 0,
# #                            minSize = 10,
# #                            maxSize = 500,
# #                            nPermSimple = 10000,
# #                            BPPARAM = MulticoreParam(),
# #                            nproc = 22)
# #
# # MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL_FGSEA_KEGG <- fgsea(pathways,
# #                                                  ranks,
# #                                                  eps = 0,
# #                                                  minSize = 10,
# #                                                  maxSize = 500,
# #                                                  nPermSimple = 10000,
# #                                                  BPPARAM = MulticoreParam(),
# #                                                  nproc = 22)
# #
# # MB_AGEGENOTYPE_OLDOVX_VS_YOUNGALL_FGSEA_KEGG %>% View
# #
# # plotEnrichment(pathways[["KEGG_LYSOSOME"]],
# #                ranks)
# #
# # leading_edges <- MB_GENOTYPE_OLD_FGSEA$leadingEdge
# # names(leading_edges) <- MB_GENOTYPE_OLD_FGSEA$pathway
# #
# # MB_GENOTYPE_OLD_META %>%
# #   left_join(anno_human) %>%
# #   filter(hsapiens_homolog_associated_gene_name %in% leading_edges[["KEGG_OXIDATIVE_PHOSPHORYLATION"]]) %>% View
# #
# # view_counts(dds_C1_MB_OLD_IP,
# #             "Rtn1") %>%
# #   left_join(colData(dds_C1_MB_OLD_IP), copy = T) %>%
# #   ggplot(aes(x = genotype,
# #              y = count)) +
# #   geom_point()
# #
# # view_counts(dds_ont_MB_OLD,
# #             "Rtn1") %>%
# #   left_join(colData(dds_ont), copy = T) %>%
# #   ggplot(aes(x = genotype,
# #              y = count)) +
# #   geom_point()

# ----
# ----
# ----

