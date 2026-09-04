library(qiime2R)
library(decontam)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load curated maAsLin2 results
top_enriched <- read.csv("results/3_maaslin2/Curated_Enriched_in_Sample.csv",
                         stringsAsFactors = FALSE)

top_depleted <- read.csv("results/3_maaslin2/Curated_Depleted_in_Sample.csv",
                         stringsAsFactors = FALSE)

dir.create("4_maaslin2_volanoplot_abundanceplot/abundance", recursive = TRUE, showWarnings = FALSE)

# Load ASV table
asv_qza <- read_qza("table.qza")
asv_raw <- asv_qza$data

# samples x ASVs
feature_table <- as.data.frame(t(asv_raw))

metadata <- read.csv(
  "metadata.csv",
  sep = ";",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE)

# Align samples
common <- intersect(
  rownames(metadata),
  rownames(feature_table))

metadata <- metadata[common, , drop = FALSE]
feature_table <- feature_table[common, , drop = FALSE]

# PCR controls identification
metadata$is.pcr_control <- grepl("^PCR_control(?:_|$)", metadata$group, ignore.case = TRUE)

# Decontamination
contam_df <- isContaminant(
  as.matrix(feature_table),
  neg = metadata$is.pcr_control,
  method = "prevalence",
  threshold = 0.5)

contaminants <- rownames(contam_df)[
  contam_df$contaminant == TRUE]

feature_table_clean <- feature_table[
  ,
  !colnames(feature_table) %in% contaminants,
  drop = FALSE]

# Remove PCR-control samples before the prevalence filter
metadata_clean <- metadata[!metadata$is.pcr_control, , drop = FALSE]
feature_table_clean <- feature_table_clean[rownames(metadata_clean), , drop = FALSE]

# Remove singletons (ASVs present in <2 samples)
keep_asvs <- colSums(feature_table_clean > 0) >= 2

feature_table_clean <- feature_table_clean[, keep_asvs, drop = FALSE]

# Create SampleType AND Species
metadata_clean$SampleType <- ifelse(
  grepl("^Control_", metadata_clean$group),
  "Environmental_control",
  "Sample")

species_temp <- gsub(
  "^Control_", "", metadata_clean$group)

species_temp <- gsub("_\\d+(\\.\\d+)?$", "", species_temp)

metadata_clean$Species <- species_temp

# Keep both host samples and environmental controls for abundance calculations
metadata_analysis <- metadata_clean
feature_analysis <- feature_table_clean[rownames(metadata_analysis), , drop = FALSE]
feature_analysis <- as.data.frame(lapply(feature_analysis, as.numeric),
                                  row.names = rownames(feature_analysis))

print(table(metadata_analysis$Species, metadata_analysis$SampleType))

# Convert counts to relative abundance (%)
sample_totals <- rowSums(feature_analysis)

# Remove samples with zero reads
keep_nonzero <- sample_totals > 0

feature_analysis <- feature_analysis[
  keep_nonzero,
  ,
  drop = FALSE]

metadata_analysis <- metadata_analysis[
  rownames(feature_analysis),
  ,
  drop = FALSE]

sample_totals <- rowSums(feature_analysis)
relative_abundance <- sweep(feature_analysis, 1, sample_totals, "/") * 100

# Host-only table retained for abundance ranks
metadata_host <- metadata_analysis[metadata_analysis$SampleType == "Sample", , drop = FALSE]
feature_host <- feature_analysis[rownames(metadata_host), , drop = FALSE]

# Rank ASVs by abundance within each species (1 = most abundant)

calculate_abundance_rank <- function(
    results_df,
    feature_host,
    metadata_host) {
  
  results_df$Mean_count_host <- NA_real_
  results_df$Abundance_rank <- NA_integer_
  
  for (sp in unique(results_df$Species)) {
    rows_sp <- results_df$Species == sp
    asvs <- results_df$ASV[rows_sp]
    asvs_clean <- sub("^X", "", asvs)
    feature_asvs_clean <- sub(
      "^X",
      "",
      colnames(feature_host))
    matches <- match(
      asvs_clean,
      feature_asvs_clean)
    valid <- !is.na(matches)
    
    if (!any(valid)) next
    
    sample_ids <- rownames(metadata_host)[
      metadata_host$Species == sp]
    
    if (length(sample_ids) == 0) next
    
    abundance_sp <- feature_host[
      sample_ids,
      matches[valid],
      drop = FALSE]
    
    # Mean read counts per ASV
    mean_counts <- colMeans(
      abundance_sp,
      na.rm = TRUE)
    
    # Rank: 1 = highest abundance
    ranks <- rank(
      -mean_counts,
      ties.method = "min")
    
    valid_indices <- which(rows_sp)[valid]
    
    results_df$Mean_count_host[
      valid_indices] <- mean_counts
    
    results_df$Abundance_rank[
      valid_indices] <- ranks
  }
  
  return(results_df)
}

# Mean abundance in hosts and corresponding environmental controls

calculate_mean_abundance <- function(
    results_df,
    relative_abundance,
    metadata_analysis) {
  
  results_df$Mean_abund_host <- NA_real_
  results_df$Mean_abund_control <- NA_real_
  
  for (sp in unique(results_df$Species)) {
    
    rows_sp <- results_df$Species == sp
    
    asvs <- results_df$ASV[rows_sp]
    
    # Remove possible MaAsLin2 "X" prefix
    asvs_clean <- sub("^X", "", asvs)
    
    feature_asvs_clean <- sub(
      "^X",
      "",
      colnames(relative_abundance))
    
    matches <- match(
      asvs_clean,
      feature_asvs_clean)
    
    valid <- !is.na(matches)
    
    if (!any(valid)) {
      warning(
        "No ASVs found for species: ",
        sp
      )
      next}
    
    valid_indices <- which(rows_sp)[valid]
    
    host_ids <- rownames(metadata_analysis)[
      metadata_analysis$Species == sp & metadata_analysis$SampleType == "Sample"]
    
    control_ids <- rownames(metadata_analysis)[
      metadata_analysis$Species == sp & metadata_analysis$SampleType == "Environmental_control"]
    
    if (length(host_ids) > 0) {
      results_df$Mean_abund_host[valid_indices] <- colMeans(
        relative_abundance[host_ids, matches[valid], drop = FALSE], na.rm = TRUE)}
    
    if (length(control_ids) > 0) {
      results_df$Mean_abund_control[valid_indices] <- colMeans(
        relative_abundance[control_ids, matches[valid], drop = FALSE], na.rm = TRUE)}
  }
  
  return(results_df)
}

# Enriched ASVs

top_enriched <- calculate_mean_abundance(
  top_enriched,
  relative_abundance,
  metadata_analysis)

enrich_abund <- top_enriched %>%
  select(
    Species,
    ASV,
    BestName,
    log2FC,
    FDR_q_value,
    Mean_abund_host,
    Mean_abund_control) %>%
  arrange(
    Species,
    desc(log2FC))

write.csv(
  enrich_abund,
  "4_maaslin2_volanoplot_abundanceplot/abundance/Enriched_ASVs_with_mean_abundance.csv",
  row.names = FALSE)

# Depleted ASVs

top_depleted <- calculate_mean_abundance(
  top_depleted,
  relative_abundance,
  metadata_analysis)

deplet_abund <- top_depleted %>%
  select(
    Species,
    ASV,
    BestName,
    log2FC,
    FDR_q_value,
    Mean_abund_host,
    Mean_abund_control) %>%
  arrange(
    Species,
    desc(abs(log2FC)))

write.csv(
  deplet_abund,
  "4_maaslin2_volanoplot_abundanceplot/abundance/Depleted_ASVs_with_mean_abundance.csv",
  row.names = FALSE)

print(head(enrich_abund, 10))

print(head(deplet_abund, 10))

# Calculate abundance rank among ALL ASVs within each host species

all_asv_rank <- data.frame()

for (sp in unique(metadata_host$Species)) {
  
  sample_ids <- rownames(metadata_host)[
    metadata_host$Species == sp]
  
  abundance_sp <- feature_host[
    sample_ids,, drop = FALSE]
  
  mean_counts <- colMeans(
    abundance_sp,
    na.rm = TRUE)
  
  rank_df <- data.frame(
    ASV = sub("^X", "", names(mean_counts)),
    Species = sp,
    Mean_count_host = mean_counts,
    Abundance_rank = rank(
      -mean_counts,
      ties.method = "min"))
  
  all_asv_rank <- rbind(
    all_asv_rank,
    rank_df)}

# Rank versus enrichment/depletion plots

enriched_top <- top_enriched %>%
  distinct() %>%
  mutate(Taxon_label = ifelse(!is.na(Genus), Genus, BestName)) %>%
  group_by(Species) %>%
  arrange(desc(abs(log2FC))) %>%
  slice_head(n = 10) %>%
  ungroup()

depleted_top <- top_depleted %>%
  distinct() %>%
  mutate(Taxon_label = ifelse(!is.na(Genus), Genus, BestName)) %>%
  group_by(Species) %>%
  arrange(desc(abs(log2FC))) %>%
  slice_head(n = 10) %>%
  ungroup()

print(range(enriched_top$Abundance_rank, na.rm = TRUE))
print(range(depleted_top$Abundance_rank, na.rm = TRUE))

p_enriched <- ggplot(
  enriched_top,
  aes(x = Abundance_rank, y = log2FC, colour = Taxon_label, shape = Class)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Taxon_label), size = 3, show.legend = FALSE, max.overlaps = 20) +
  scale_x_reverse() +
  facet_wrap(~Species, scales = "free_x") +
  theme_bw() +
  labs(x = "Abundance rank among host ASVs\n(1 = most abundant)",
       y = "log2FC enrichment", colour = "Taxon", shape = "Class",
       title = "Top 10 enriched ASVs (by |log2FC|)") +
  guides(colour = "none")

ggsave("4_maaslin2_volanoplot_abundanceplot/Enriched_rank_vs_log2FC_top10.png",
       p_enriched, width = 12, height = 8, dpi = 300)

p_depleted <- ggplot(
  depleted_top,
  aes(x = Abundance_rank, y = log2FC, colour = Taxon_label, shape = Class)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Taxon_label), size = 3, show.legend = FALSE, max.overlaps = 20) +
  scale_x_reverse() +
  facet_wrap(~Species, scales = "free_x") +
  theme_bw() +
  labs(x = "Abundance rank among host ASVs\n(1 = most abundant)",
       y = "log2FC depletion", colour = "Taxon", shape = "Class",
       title = "Top 10 depleted ASVs (by |log2FC|)") +
  guides(colour = "none")

ggsave("4_maaslin2_volanoplot_abundanceplot/Depleted_rank_vs_log2FC_top10.png",
       p_depleted, width = 12, height = 8, dpi = 300)

# VOLCANO PLOT 

# Load tables
enriched <- read.csv("4_maaslin2_volanoplot_abundanceplot/abundance/Enriched_ASVs_with_rank.csv",
                     stringsAsFactors = FALSE)
depleted <- read.csv(
  "4_maaslin2_volanoplot_abundanceplot/abundance/Depleted_ASVs_with_rank.csv",
  stringsAsFactors = FALSE)

volcano_data <- bind_rows(
  enriched %>%
    distinct() %>%
    mutate(Direction = "Enriched in Sample"),
  
  depleted %>%
    distinct() %>%
    mutate(Direction = "Enriched in Control")) %>%
  mutate(
    # Calculate -log10(q), protecting against q-values equal to zero
    neg_log10_q = -log10(
      pmax(FDR_q_value, .Machine$double.xmin)),
    
    # Prefer genus; otherwise use BestName and finally ASV
    Taxon = case_when(
      !is.na(Genus) & Genus != "" ~ Genus,
      !is.na(BestName) & BestName != "" ~ BestName,
      TRUE ~ ASV),
    
    # Give missing classes a visible category
    Class = ifelse(
      is.na(Class) | Class == "",
      "Unclassified",
      Class))

# Select the top 10 enriched and depleted ASVs per host species
volcano_top20 <- volcano_data %>%
  group_by(Species, Direction) %>%
  slice_max(
    order_by = abs(log2FC),
    n = 10,
    with_ties = FALSE) %>%
  ungroup()

# Volcano plot
p_volcano <- ggplot(
  volcano_top20,
  aes(
    x = log2FC,
    y = neg_log10_q)) +
  geom_point(
    aes(
      color = Direction,
      shape = Class),
    size = 3.5,
    alpha = 0.85) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40") +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40") +
  geom_text_repel(
    aes(label = Taxon),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE) +
  facet_wrap(
    ~Species,
    scales = "free",
    ncol = 2) +
  scale_color_manual(
    values = c(
      "Enriched in Sample" = "#79CDCD",
      "Enriched in Control" = "#FF7256"),
    labels = c(
      "Enriched in Sample" = "Top 10 enriched",
      "Enriched in Control" = "Top 10 depleted")) +
  labs(
    title = "Top enriched and depleted MaAsLin2 associations",
    subtitle = "Top 10 in each direction per species; q < 0.05",
    x = "log2FC",
    y = expression(-log[10]("FDR q-value")),
    color = "Association",
    shape = "Class") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank())

p_volcano

ggsave(
  "4_maaslin2_volanoplot_abundanceplot/Volcano_Top10_Enriched_Top10_Depleted.png",
  p_volcano,
  width = 15,
  height = 10,
  dpi = 300)
