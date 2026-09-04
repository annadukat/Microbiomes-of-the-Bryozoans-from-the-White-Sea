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

outdir <- "2_repeated_rarefaction"
if (!dir.exists(outdir)) dir.create(outdir)

# Read files
asv16 <- read_qza(asv_file)
asv_table <- asv16$data

taxonomy16 <- read_qza(taxonomy_file)
tree16 <- read_qza(tree_file)

metadata <- read.csv(
  metadata_file,
  sep = ";",
  header = TRUE,
  check.names = FALSE)

colnames(metadata)[1] <- "sample.id"
metadata$sample.id <- trimws(metadata$sample.id)
rownames(metadata) <- metadata$sample.id

# Identify PCR controls
metadata$is.pcr_control <- grepl(
  "^PCR_control(?:_|$)",
  metadata$group,
  ignore.case = TRUE)

# Build phyloseq object
ps_raw <- phyloseq(
  otu_table(asv_table, taxa_are_rows = TRUE),
  sample_data(metadata))

tax_table(ps_raw) <- tax_table(
  as.matrix(parse_taxonomy(taxonomy16$data)))

phy_tree(ps_raw) <- phy_tree(tree16$data)

print(ps_raw)

# Decontamination
otu_raw <- as(otu_table(ps_raw), "matrix")
if (taxa_are_rows(ps_raw)) {
  otu_raw <- t(otu_raw)}

contam_df <- isContaminant(
  otu_raw,
  neg = sample_data(ps_raw)$is.pcr_control,
  method = "prevalence",
  threshold = 0.5)

contaminants <- rownames(contam_df)[contam_df$contaminant == TRUE]

ps_clean <- prune_taxa(
  !taxa_names(ps_raw) %in% contaminants,
  ps_raw)

ps_analysis <- ps_clean

# Remove PCR controls
keep_samples <- !sample_data(ps_analysis)$is.pcr_control

ps_analysis <- prune_samples(keep_samples, ps_analysis)

# Remove singletons (occurring in fewer than 2 samples)
ps_analysis <- filter_taxa(
  ps_analysis,
  function(x) sum(x > 0) >= 2,
  TRUE)

# Remove empty samples and taxa
ps_analysis <- prune_samples(sample_sums(ps_analysis) > 0, ps_analysis)
ps_analysis <- prune_taxa(taxa_sums(ps_analysis) > 0, ps_analysis)

# Verify removals
remaining_pcr <- sample_data(ps_analysis)$is.pcr_control
stopifnot(!any(remaining_pcr))

# Alpha diversity at different depths
n_reps <- 100
n_points <- 20

ps_rarefaction <- ps_analysis
rare_depth <- min(sample_sums(ps_rarefaction))

# Extract OTU matrix (samples x taxa)
otu_mat <- as(otu_table(ps_rarefaction), "matrix")
if (taxa_are_rows(ps_rarefaction)) {
  otu_mat <- t(otu_mat)}
sample_ids <- rownames(otu_mat)

# Get tree and prune to common taxa for PD
tree_alpha <- phy_tree(ps_rarefaction)
common_taxa <- intersect(colnames(otu_mat), tree_alpha$tip.label)
otu_mat <- otu_mat[, common_taxa, drop = FALSE]
tree_alpha <- prune_taxa(common_taxa, tree_alpha)

alpha_results <- vector("list", length(sample_ids))
names(alpha_results) <- sample_ids

set.seed(123)

for (sample_id in sample_ids) {
  
  counts <- otu_mat[sample_id, , drop = TRUE]
  total_reads <- sum(counts)
  cat("Processing:", sample_id, "| library size =", total_reads, "\n")
  
  # Generate rarefaction depths
  depths <- round(
    seq(
      from = max(10, ceiling(total_reads * 0.10)),
      to = total_reads,
      length.out = n_points))
  depths <- sort(unique(c(depths, rare_depth)))
  depths <- depths[depths <= total_reads & depths >= 1]
  
  sample_results <- vector("list", length(depths))
  
  for (d in seq_along(depths)) {
    
    depth <- depths[d]
    results <- vector("list", n_reps)
    
    for (rep in seq_len(n_reps)) {
      
      # Rarefy
      rarefied_counts <- vegan::rrarefy(
        matrix(counts, nrow = 1),
        sample = depth)
      
      # Keep as matrix with column names
      colnames(rarefied_counts) <- colnames(otu_mat)
      
      # Alpha metrics (using the matrix)
      observed <- sum(rarefied_counts[1, ] > 0)
      shannon <- vegan::diversity(rarefied_counts, index = "shannon")[1]
      chao1 <- vegan::estimateR(rarefied_counts[1, ])["S.chao1"]
      
      # Faith's PD – use the matrix directly
      if (!is.null(tree_alpha) && Ntip(tree_alpha) > 0) {
        pd_result <- picante::pd(rarefied_counts, tree_alpha, include.root = FALSE)
        pd_value <- pd_result$PD[1]} else {
          pd_value <- NA}
      
      results[[rep]] <- data.frame(
        SampleID = sample_id,
        Depth = depth,
        Replicate = rep,
        Observed = observed,
        Shannon = shannon,
        Chao1 = as.numeric(chao1),
        PD = pd_value)}
    
    sample_results[[d]] <- bind_rows(results)}
  
  alpha_results[[sample_id]] <- bind_rows(sample_results)}

alpha_raw <- bind_rows(alpha_results)

write.csv(alpha_raw,
  file.path(outdir, "alpha_diversity_different_sequencing_depths_raw.csv"),
  row.names = FALSE)

# Summary
print(head(alpha_raw))

# Plot rarefaction curves

library(tidyverse)
library(viridis)

alpha_raw <- read.csv(
  file.path(outdir, "alpha_diversity_different_sequencing_depths_raw.csv"),
  header = TRUE)

# Read and prepare metadata
plot_metadata <- read.csv(
  metadata_file,
  sep = ";",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE)

colnames(plot_metadata)[1] <- "SampleID"

plot_metadata <- plot_metadata %>%
  transmute(
    SampleID = trimws(SampleID),
    Group = trimws(group),
    Species = str_remove(trimws(group), "_[0-9]+(?:\\.[0-9]+)?$"))

# Summarise rarefaction replicates
alpha_summary <- alpha_raw %>%
  group_by(SampleID, Depth) %>%
  summarise(
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_sd = sd(Observed, na.rm = TRUE),
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_sd = sd(Shannon, na.rm = TRUE),
    Chao1_mean = mean(Chao1, na.rm = TRUE),
    Chao1_sd = sd(Chao1, na.rm = TRUE),
    PD_mean = mean(PD, na.rm = TRUE),
    PD_sd = sd(PD, na.rm = TRUE),
    .groups = "drop")

# Reshape and add species information
alpha_long <- alpha_summary %>%
  pivot_longer(
    cols = c(Observed_mean, Shannon_mean, Chao1_mean, PD_mean),
    names_to = "Metric_code",
    values_to = "Mean") %>%
  mutate(
    SD = case_when(
      Metric_code == "Observed_mean" ~ Observed_sd,
      Metric_code == "Shannon_mean" ~ Shannon_sd,
      Metric_code == "Chao1_mean" ~ Chao1_sd,
      Metric_code == "PD_mean" ~ PD_sd),
    Metric = recode(
      Metric_code,
      Observed_mean = "Observed richness",
      Shannon_mean = "Shannon diversity",
      Chao1_mean = "Chao1 richness",
      PD_mean = "Faith's PD")) %>%
  left_join(
    plot_metadata %>% distinct(SampleID, .keep_all = TRUE),
    by = "SampleID")

# Stop if any samples could not be matched to metadata
unmatched_samples <- alpha_long %>%
  filter(is.na(Species)) %>%
  distinct(SampleID)

if (nrow(unmatched_samples) > 0) {
  print(unmatched_samples, n = Inf)
  stop("Some samples could not be matched to the metadata.")}

# Create one rarefaction plot
make_rarefaction_plot <- function(plot_data, species_name) {
  sample_ids <- sort(unique(plot_data$SampleID))
  sample_colours <- viridis(length(sample_ids))
  names(sample_colours) <- sample_ids
  display_name <- str_replace_all(species_name, "_", " ")
  
  ggplot(
    plot_data,
    aes(
      x = Depth,
      y = Mean,
      colour = SampleID,
      fill = SampleID,
      group = SampleID)) +
    geom_ribbon(
      aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD),
      alpha = 0.15,
      colour = NA,
      na.rm = TRUE) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_vline(
      xintercept = rare_depth,
      linetype = "dashed",
      colour = "red",
      alpha = 0.6,
      linewidth = 0.8) +
    annotate(
      "text",
      x = rare_depth + 200,
      y = Inf,
      label = paste("Rarefaction depth =", rare_depth),
      vjust = 1.5,
      hjust = 0,
      colour = "red",
      size = 3.5,
      fontface = "bold") +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = sample_colours) +
    scale_fill_manual(values = sample_colours) +
    labs(
      title = paste("Alpha diversity rarefaction curves:", display_name),
      subtitle = paste0("Mean ± SD across ", n_reps, " rarefaction replicates"),
      x = "Sequencing depth (number of reads)",
      y = "Diversity metric value",
      colour = "Sample",
      fill = "Sample") +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40"))
}

# Create and save one plot per species
species_plots <- alpha_long %>%
  split(.$Species) %>%
  imap(~make_rarefaction_plot(.x, .y))

iwalk(species_plots, function(species_plot, species_name) {
  safe_name <- str_replace_all(species_name, "[^A-Za-z0-9_-]", "_")
  
  ggsave(
    file.path(outdir, paste0("rarefaction_", safe_name, ".png")),
    species_plot,
    width = 10,
    height = 8,
    dpi = 300)
})

# Colors shared across all samples
sample_ids <- sort(unique(alpha_long$SampleID))
sample_colours <- viridis(length(sample_ids))
names(sample_colours) <- sample_ids

p_species <- ggplot(
  alpha_long,
  aes(
    x = Depth,
    y = Mean,
    colour = SampleID,
    fill = SampleID,
    group = SampleID)) +
  geom_ribbon(
    aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD),
    alpha = 0.15,
    colour = NA,
    na.rm = TRUE) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_vline(
    xintercept = rare_depth,
    linetype = "dashed",
    colour = "red",
    alpha = 0.6,
    linewidth = 0.8) +
  facet_grid(
    Metric ~ Species,
    scales = "free_y",
    labeller = labeller(
      Species = function(x) str_replace_all(x, "_", " "))) +
  scale_colour_manual(values = sample_colours) +
  scale_fill_manual(values = sample_colours) +
  labs(
    title = "Alpha diversity rarefaction curves",
    subtitle = paste0("Mean ± SD across ", n_reps, " rarefaction replicates"),
    x = "Sequencing depth (number of reads)",
    y = "Diversity metric value",
    colour = "Sample",
    fill = "Sample") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40"))

ggsave(
  file.path(outdir, "rarefaction_all_species_faceted.png"),
  p_species,
  width = 15,
  height = 12,
  dpi = 300)

print(p_species)
