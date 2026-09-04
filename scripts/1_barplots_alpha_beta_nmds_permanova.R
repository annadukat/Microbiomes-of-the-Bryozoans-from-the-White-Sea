library(qiime2R)
library(phyloseq)
library(tidyverse)
library(vegan)
library(decontam)
library(picante)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(UpSetR)
library(ggVennDiagram)
library(reshape2)  

options(stringsAsFactors = FALSE)

# Input files
asv_file      <- "table.qza"
taxonomy_file <- "taxonomy.qza"
tree_file     <- "rooted_tree.qza"
metadata_file <- "metadata.csv"

outdir <- "1_barplots_alpha_beta_nmds_permanova"
if (!dir.exists(outdir)) dir.create(outdir)

#  Read files
asv16 <- read_qza(asv_file)
asv_table <- asv16$data

taxonomy16 <- read_qza(taxonomy_file)
tree16 <- read_qza(tree_file)

metadata <- read.csv(metadata_file, sep = ";", header = TRUE, check.names = FALSE)
colnames(metadata)[1] <- "sample.id"
metadata$sample.id <- trimws(metadata$sample.id)
rownames(metadata) <- metadata$sample.id

# Identify PCR controls (k1-k4)
metadata$is.pcr_control <- grepl("^PCR_control(?:_|$)", metadata$group, ignore.case = TRUE)

# Build phyloseq object
ps_raw <- phyloseq(
  otu_table(asv_table, taxa_are_rows = TRUE),
  sample_data(metadata))

tax_table(ps_raw) <- tax_table(as.matrix(parse_taxonomy(taxonomy16$data)))
phy_tree(ps_raw) <- phy_tree(tree16$data)
print(ps_raw)

# Decontamination
otu_raw <- as(otu_table(ps_raw), "matrix")
if(taxa_are_rows(ps_raw)) otu_raw <- t(otu_raw)

contam_df <- isContaminant(otu_raw,
                           neg = sample_data(ps_raw)$is.pcr_control,
                           method = "prevalence",
                           threshold = 0.5)

contaminants <- rownames(contam_df)[contam_df$contaminant == TRUE]
cat("Removing", length(contaminants), "contaminant ASVs\n")

ps_clean <- prune_taxa(!taxa_names(ps_raw) %in% contaminants, ps_raw)

# Remove PCR control
ps_clean <- subset_samples(ps_clean, !is.pcr_control)
print(ps_clean)

# Filter singletons (ASVs present in <2 samples)
ps_clean <- filter_taxa(ps_clean, function(x) sum(x > 0) >= 2, TRUE)
print(ps_clean)

# Relative abundance (for barplots, NMDS, and PERMANOVA)
ps_rel <- transform_sample_counts(ps_clean, function(x) x / sum(x))

# Grouping variables and sample type for plot
sample_data(ps_rel)$group_orig <- sample_data(ps_rel)$group

sample_data(ps_rel)$SampleType <- ifelse(grepl("^Control_", sample_data(ps_rel)$group_orig),
                                         "Environmental_control", "Sample")

# Extract Species name (remove "Control_" prefix and trailing _number or _number.number)
tmp <- gsub("^Control_", "", sample_data(ps_rel)$group_orig)
tmp <- gsub("_\\d+(\\.\\d+)?$", "", tmp)
sample_data(ps_rel)$Species <- tmp

# Use the original sample ID (rowname) as x-axis label
sample_data(ps_rel)$SampleID <- rownames(sample_data(ps_rel))

# Verify species and sample IDs
head(sample_data(ps_rel)[, c("SampleID", "group_orig", "Species")])

# TAXONOMY BARPLOTS

prepare_bar_data <- function(ps_obj, rank, top_n = 20) {
  ps_glom <- tax_glom(ps_obj, taxrank = rank)
  top <- names(sort(taxa_sums(ps_glom), decreasing = TRUE))[1:min(top_n, ntaxa(ps_glom))]
  ps_top <- prune_taxa(top, ps_glom) %>%
    transform_sample_counts(function(x) 100 * x / sum(x))
  
  df <- psmelt(ps_top)
  colnames(df)[colnames(df) == rank] <- "Taxon"
  
  # Fix renamed sample variables
  if("sample_SampleID" %in% colnames(df)) df$SampleID <- df$sample_SampleID
  if("sample_Species" %in% colnames(df)) df$Species <- df$sample_Species
  if("sample_SampleType" %in% colnames(df)) df$SampleType <- df$sample_SampleType
  
  df <- df[, c("SampleID", "Abundance", "Taxon", "Species", "SampleType")]
  return(df)}

# Prepare genus data (discards ASVs that lacks a genus assignment, which is important since we have a lot of NA in the microbiome)
df_genus <- prepare_bar_data(ps_rel, "Genus", top_n = 20)

# Get unique species
species_list <- unique(df_genus$Species)

# Loop visualization
for(sp in species_list) {
  df_sp <- df_genus %>%
    filter(Species == sp) %>%
    mutate(SampleID = factor(SampleID, 
                             levels = unique(SampleID[order(SampleType, SampleID)])))
  
  # Calculate dynamic width: 0.8 inches per sample, min 10, max 30
  n_samples <- length(unique(df_sp$SampleID))
  plot_width <- max(10, min(30, n_samples * 0.8))
  
  p_genus <- ggplot(df_sp, aes(x = SampleID, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    ylab("Relative abundance (%)") + xlab("Sample ID") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(df_sp$Taxon)))) +
    theme_classic(base_size = 16) +  # bigger base font
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          plot.title = element_text(size = 18, face = "bold"),
          strip.text = element_text(size = 14, face = "bold")) +
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle(paste(sp, "- 20 most abundant genera"))
  
  print(p_genus)
  ggsave(file.path(outdir, paste0(sp, "_top20_genus.png")), p_genus,
         width = plot_width, height = 10, dpi = 600, limitsize = FALSE)}

# ALPHA DIVERSITY (rarefied) 

# Rarefaction to the min library size
rare_depth <- min(sample_sums(ps_clean))
ps_rare <- rarefy_even_depth(
  ps_clean,
  sample.size = rare_depth,
  rngseed = 123,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE)

# Metadata

meta <- data.frame(sample_data(ps_rare))
meta$SampleID <- rownames(meta)

meta$SampleType <- ifelse(
  grepl("^Control_", meta$group),
  "Environmental_control",
  "Sample")

meta$Species <- NA
is_sample <- meta$SampleType == "Sample"
meta$Species[is_sample] <- gsub("_\\d+(\\.\\d+)?$", "", meta$group[is_sample])

is_ctrl <- meta$SampleType == "Environmental_control"
ctrl_clean <- gsub("^Control_", "", meta$group[is_ctrl])
meta$Species[is_ctrl] <- gsub("_\\d+(\\.\\d+)?$", "", ctrl_clean)

verification <- meta[, c("SampleID", "group", "Species", "SampleType")]
print(verification[1:min(20, nrow(verification)), ])
print(table(meta$Species, meta$SampleType))

# Alpha diversity calculation
# Base metrics (Observed & Shannon)
alpha <- estimate_richness(ps_rare, measures = c("Observed", "Shannon"))

# Add Chao1
chao1_metrics <- estimate_richness(ps_rare, measures = "Chao1")
alpha$Chao1 <- chao1_metrics$Chao1

# Add Faith's PD
otu_mat <- t(otu_table(ps_rare))
tree <- phy_tree(ps_rare)
pd_out <- pd(otu_mat, tree, include.root = FALSE)
alpha$PD <- pd_out$PD

# SampleID
alpha$SampleID <- rownames(alpha)

# Align metadata
meta <- meta[match(alpha$SampleID, meta$SampleID), ]
stopifnot(all(meta$SampleID == alpha$SampleID))
alpha$Species <- meta$Species
alpha$SampleType <- meta$SampleType
alpha$CompareGroup <- alpha$Species

# Wilcoxon tests for all metrics
species_list <- unique(na.omit(alpha$Species))
alpha_stats <- data.frame()

for (sp in species_list) {
  df <- na.omit(alpha[alpha$Species == sp,
                      c("Species", "SampleType",
                        "Shannon", "Observed", "Chao1", "PD")])
  if (length(unique(df$SampleType)) < 2) next
  
  wt_shannon  <- wilcox.test(Shannon ~ SampleType, data = df)
  wt_observed <- wilcox.test(Observed ~ SampleType, data = df)
  wt_chao1    <- wilcox.test(Chao1 ~ SampleType, data = df)
  wt_pd       <- wilcox.test(PD ~ SampleType, data = df)
  
  alpha_stats <- rbind(alpha_stats, data.frame(
    Species = sp,
    pvalue_Shannon  = wt_shannon$p.value,
    pvalue_Observed = wt_observed$p.value,
    pvalue_Chao1    = wt_chao1$p.value,
    pvalue_PD       = wt_pd$p.value))}

# Summary statistics
summary_stats <- alpha %>%
  group_by(Species, SampleType) %>%
  summarise(
    n = n(),
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_sd   = sd(Observed, na.rm = TRUE),
    Observed_se   = Observed_sd / sqrt(n),
    Shannon_mean  = mean(Shannon, na.rm = TRUE),
    Shannon_sd    = sd(Shannon, na.rm = TRUE),
    Shannon_se    = Shannon_sd / sqrt(n),
    Chao1_mean    = mean(Chao1, na.rm = TRUE),
    Chao1_sd      = sd(Chao1, na.rm = TRUE),
    Chao1_se      = Chao1_sd / sqrt(n),
    PD_mean       = mean(PD, na.rm = TRUE),
    PD_sd         = sd(PD, na.rm = TRUE),
    PD_se         = PD_sd / sqrt(n),
    .groups = "drop")

write.csv(summary_stats,
          file.path(outdir, "alpha_diversity_summary_stats_all_metrics.csv"),
          row.names = FALSE)
print(summary_stats)

# Plotting function
plot_alpha <- function(data, metric, y_label, label_df = NULL) {
  p <- ggplot(data, aes(x = CompareGroup, y = .data[[metric]], fill = SampleType)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(0.7)) +
    geom_point(shape = 21, size = 2.5, color = "black",
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Species", y = y_label)
  
  if (!is.null(label_df)) {
    p <- p + geom_text(data = label_df,
                       aes(x = Species, y = y_pos, label = label),
                       inherit.aes = FALSE, size = 3.5, vjust = 0)}
  return(p)}

# Combined plots using alpha_long and label_data

alpha_long <- alpha %>%
  select(SampleID, Species, SampleType, all_of(c("Shannon", "Observed", "Chao1", "PD"))) %>%
  pivot_longer(cols = c("Shannon", "Observed", "Chao1", "PD"), names_to = "Metric", values_to = "Value")

# Plot all 4 metrics

p_all_four <- function(label_df = NULL) {
  p <- alpha_long %>%
    ggplot(aes(x = Species, y = Value, fill = SampleType)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(0.7)) +
    geom_point(shape = 21, size = 2.5, color = "black",
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold")) +
    labs(x = "Species", y = "Value")
  
  if (!is.null(label_df)) {
    p <- p + geom_text(data = label_df,
                       aes(x = Species, y = y_pos, label = label),
                       inherit.aes = FALSE, size = 3.5, vjust = 0)}
  return(p)}

p_all_four_wo <- p_all_four(NULL)

ggsave(file.path(outdir, "Alpha_diversity_All_four_metrics_no_labels.png"),
       p_all_four_wo, width = 10, height = 8, dpi = 600)

# NMDS and PERMANOVA (Bray‑Curtis on relative abundance)

# Ensure SampleType and Species are correctly set in ps_rel
sample_data(ps_rel)$group_orig <- sample_data(ps_rel)$group  # backup
sample_data(ps_rel)$SampleType <- ifelse(grepl("^Control_", sample_data(ps_rel)$group_orig),
                                         "Environmental_control", "Sample")

tmp <- gsub("^Control_", "", sample_data(ps_rel)$group_orig)
tmp <- gsub("_\\d+(\\.\\d+)?$", "", tmp)
sample_data(ps_rel)$Species <- tmp

# NMDS (on relative abundance three dimensions)

set.seed(123)

otu_nmds <- as(
  otu_table(ps_rel),
  "matrix")

if (taxa_are_rows(ps_rel)) {
  otu_nmds <- t(otu_nmds)}

# One 3-dimensional NMDS
nmds <- metaMDS(
  log1p(otu_nmds),
  distance = "bray",
  k = 3,
  trymax = 1000,
  autotransform = FALSE,
  trace = TRUE)

# Extract all three NMDS coordinates
ord_df <- as.data.frame(
  scores(
    nmds,
    display = "sites",
    choices = 1:3))

colnames(ord_df)[1:3] <- c(
  "NMDS1",
  "NMDS2",
  "NMDS3")

ord_df$SampleID <- rownames(ord_df)

# Add metadata without changing sample order
meta_df <- data.frame(
  sample_data(ps_rel),
  check.names = FALSE)

meta_df$SampleID <- rownames(meta_df)

metadata_index <- match(
  ord_df$SampleID,
  meta_df$SampleID)

ord_df$Species <-
  meta_df$Species[metadata_index]

ord_df$SampleType <-
  meta_df$SampleType[metadata_index]

ord_df$SampleType <-
  meta_df$SampleType[metadata_index]

print(table(ord_df$Species, ord_df$SampleType, useNA = "ifany"))

# Three pairwise projections
nmds_projection_df <- bind_rows(
  ord_df %>%
    transmute(
      SampleID,
      Species,
      SampleType,
      AxisX = NMDS1,
      AxisY = NMDS2,
      Projection = "NMDS1 vs NMDS2"),
  
  ord_df %>%
    transmute(
      SampleID,
      Species,
      SampleType,
      AxisX = NMDS1,
      AxisY = NMDS3,
      Projection = "NMDS1 vs NMDS3"),
  
  ord_df %>%
    transmute(
      SampleID,
      Species,
      SampleType,
      AxisX = NMDS2,
      AxisY = NMDS3,
      Projection = "NMDS2 vs NMDS3"))

nmds_projection_df$Projection <- factor(
  nmds_projection_df$Projection,
  levels = c(
    "NMDS1 vs NMDS2",
    "NMDS1 vs NMDS3",
    "NMDS2 vs NMDS3"))

# Use common limits so the geometric scale is comparable
axis_limits <- range(
  c(
    nmds_projection_df$AxisX,
    nmds_projection_df$AxisY),
  finite = TRUE)

axis_padding <- 0.05 * diff(axis_limits)

axis_limits <- c(
  axis_limits[1] - axis_padding,
  axis_limits[2] + axis_padding)

# Combined plot
p_nmds_3projections <- ggplot(
  nmds_projection_df,
  aes(
    x = AxisX,
    y = AxisY,
    color = Species,
    fill = Species)) +
  geom_point(
    aes(shape = SampleType),
    size = 4,
    alpha = 0.9) +
  facet_wrap(
    ~ Projection,
    nrow = 1) +
  scale_x_continuous(
    limits = axis_limits) +
  scale_y_continuous(
    limits = axis_limits) +
  scale_shape_manual(
    values = c(
      "Environmental_control" = 17,
      "Sample" = 16)) +
  scale_color_brewer(
    palette = "Set2") +
  scale_fill_brewer(
    palette = "Set2") +
  coord_equal() +
  theme_bw(base_size = 14) +
  labs(
    title = "Three-dimensional NMDS of microbial communities",
    subtitle = paste0(
      "Three projections of one NMDS solution; ",
      "Bray-Curtis; stress = ",
      round(nmds$stress, 3)),
    x = NULL,
    y = NULL,
    color = "Species",
    fill = "Species",
    shape = "Sample type") +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 13,
      face = "bold"),
    panel.spacing = grid::unit(
      1,
      "lines"))

print(p_nmds_3projections)

ggsave(
  file.path(
    outdir,
    "NMDS_3D_three_projections_combined.png"),
  p_nmds_3projections,
  width = 18,
  height = 6,
  dpi = 600)

# Prepare data for PERMANOVA

# OTU/ASV table: samples x ASVs
otu_perma <- as(otu_table(ps_rel), "matrix")

if (taxa_are_rows(ps_rel)) {
  otu_perma <- t(otu_perma)
}

# Metadata
meta_perma <- data.frame(
  sample_data(ps_rel),
  check.names = FALSE)

meta_perma <- meta_perma[rownames(otu_perma), , drop = FALSE]

stopifnot(
  identical(rownames(otu_perma), rownames(meta_perma)))

# Convert grouping variables to factors
meta_perma$Species <- factor(meta_perma$Species)
meta_perma$SampleType <- factor(meta_perma$SampleType)

# Bray-Curtis distance
dist_bc <- vegdist(
  otu_perma,
  method = "bray")

# Check
print(table(
  meta_perma$Species,
  meta_perma$SampleType,
  useNA = "ifany"))

# Global PERMANOVA (Species × SampleType)
set.seed(123)

global_res <- adonis2(
  dist_bc ~ Species * SampleType,
  data = meta_perma,
  permutations = 999)

print(global_res)

write.csv(
  as.data.frame(global_res),
  file.path(outdir, "PERMANOVA_global_Species_x_SampleType.csv"),
  row.names = TRUE)

# Pairwise PERMANOVA (SampleType within each Species)
pairwise_results <- data.frame()

for (sp in unique(meta_perma$Species)) {
  
  idx <- which(meta_perma$Species == sp)
  
  if (length(unique(meta_perma$SampleType[idx])) < 2) next
  
  dist_sub <- as.dist(as.matrix(dist_bc)[idx, idx])
  meta_sub <- meta_perma[idx, , drop = FALSE]
  
  set.seed(123)
  
  res <- adonis2(
    dist_sub ~ SampleType,
    data = meta_sub,
    permutations = 999)
  
  pairwise_results <- rbind(
    pairwise_results,
    data.frame(
      Species = as.character(sp),
      Df = res$Df[1],
      SumOfSqs = res$SumOfSqs[1],
      R2 = res$R2[1],
      F = res$F[1],
      pvalue_raw = res$`Pr(>F)`[1]))}

# BONFERRONI CORRECTION
pairwise_results$pvalue_adjusted <- p.adjust(
  pairwise_results$pvalue_raw,
  method = "bonferroni")

print(pairwise_results)

write.csv(
  pairwise_results,
  file.path(
    outdir,
    "PERMANOVA_pairwise_per_species_Bonferroni.csv"),
  row.names = FALSE)

# Pairwise distances, Bray-Curtis and Jaccard

# Use relative abundance data
ps_beta <- ps_rel

# OTU table (samples x ASVs) – relative abundances
otu_beta <- as(otu_table(ps_beta), "matrix")
if(taxa_are_rows(ps_beta)) otu_beta <- t(otu_beta)

# Metadata
meta_beta <- data.frame(sample_data(ps_beta))[, c("Species", "SampleType")]
meta_beta$SampleID <- rownames(meta_beta)

# Distance matrices
dist_bc <- vegdist(otu_beta, method = "bray")
dist_jaccard <- vegdist(otu_beta > 0, method = "jaccard", binary = TRUE)

dist_list <- list(
  BrayCurtis = as.matrix(dist_bc),
  Jaccard = as.matrix(dist_jaccard))

# Species and control sample identifiers
species_list <- unique(meta_beta$Species[meta_beta$SampleType == "Sample"])
species_list <- species_list[!is.na(species_list)]
control_samples <- rownames(meta_beta)[meta_beta$SampleType == "Environmental_control"]

get_pairwise_df_for_species <- function(sp, dist_mat, meta_beta, control_samples) {
  sp_samples <- rownames(meta_beta)[
    meta_beta$Species == sp & meta_beta$SampleType == "Sample"]
  
  if (length(sp_samples) < 2) return(NULL)
  
  other_samples <- rownames(meta_beta)[
    meta_beta$SampleType == "Sample" &
      meta_beta$Species != sp &
      !is.na(meta_beta$Species)]
  
  within_vals <- dist_mat[
    sp_samples, sp_samples, drop = FALSE][upper.tri(dist_mat[sp_samples, sp_samples, drop = FALSE])]
  
  between_vals <- if (length(other_samples) > 0) {
    as.vector(dist_mat[sp_samples, other_samples, drop = FALSE])} else numeric(0)
  
  control_vals <- if (length(control_samples) > 0) {
    as.vector(dist_mat[sp_samples, control_samples, drop = FALSE])} else numeric(0)
  
  within_control_vals <- if (length(control_samples) >= 2) {
    control_matrix <- dist_mat[control_samples, control_samples, drop = FALSE]
    control_matrix[upper.tri(control_matrix)]} else numeric(0)
  
  bind_rows(
    data.frame(Distance = within_vals, Type = "Within same species"),
    data.frame(
      Distance = between_vals,
      Type = rep("Between species", length(between_vals))),
    data.frame(
      Distance = control_vals,
      Type = rep("Sample vs Control", length(control_vals))),
    data.frame(
      Distance = within_control_vals,
      Type = rep("Within controls", length(within_control_vals))))}

# Combine all distances
plot_df <- bind_rows(lapply(names(dist_list), function(metric_name) {
  bind_rows(lapply(species_list, function(sp) {
    df <- get_pairwise_df_for_species(
      sp, dist_list[[metric_name]], meta_beta, control_samples
    )
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>%
      mutate(
        Species = sp,
        Metric = metric_name)}))}))

plot_df <- plot_df %>%
  mutate(
    Type = factor(
      Type,
      levels = c(
        "Within same species",
        "Between species",
        "Sample vs Control",
        "Within controls")),
    Metric = factor(
      Metric,
      levels = c("BrayCurtis", "Jaccard"),
      labels = c("Bray–Curtis", "Jaccard")),
    Species_label = gsub("_", " ", Species))

# Calculate and save medians
median_wide <- plot_df %>%
  group_by(Species, Metric, Type) %>%
  summarise(Median = median(Distance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(Species, Metric),
    names_from = Type,
    values_from = Median)

print(median_wide)

write.csv(
  median_wide,
  file.path(outdir, "Median_pairwise_distances.csv"),
  row.names = FALSE)

# Combined plot only
p_combined <- ggplot(
  plot_df,
  aes(x = Type, y = Distance, fill = Type)) +
  geom_boxplot(
    width = 0.7,
    outlier.shape = NA,
    alpha = 0.75) +
  geom_jitter(
    width = 0.15,
    alpha = 0.35,
    size = 1) +
  facet_grid(Metric ~ Species_label) +
  scale_fill_manual(values = c(
    "Within same species" = "lightblue",
    "Between species" = "salmon",
    "Sample vs Control" = "lightgreen",
    "Within controls" = "orange")) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Pairwise beta-diversity distances",
    x = NULL,
    y = "Dissimilarity",
    fill = "Comparison") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text.x = element_text(face = "italic"),
    strip.text.y = element_text(face = "bold"),
    legend.position = "top")

print(p_combined)

ggsave(file.path(outdir, "BrayCurtis_Jaccard.png"),
  p_combined,
  width = 14,
  height = 8,
  dpi = 600)
