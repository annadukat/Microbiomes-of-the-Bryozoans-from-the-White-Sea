library(qiime2R)
library(Maaslin2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(phyloseq)
library(decontam)

# Load data
asv_qza <- read_qza("table.qza")
tax_qza <- read_qza("taxonomy.qza")
metadata <- read.csv("metadata.csv", sep = ";", row.names = 1, check.names = FALSE)

asv_raw <- asv_qza$data
asv_ids <- rownames(asv_raw)

taxonomy <- tax_qza$data
rownames(taxonomy) <- asv_ids
taxonomy$Feature.ID <- asv_ids

# Feature table (samples x ASVs)
feature_table <- as.data.frame(t(asv_raw))

# Align samples
common <- intersect(rownames(metadata), rownames(feature_table))
metadata <- metadata[common, , drop = FALSE]
feature_table <- feature_table[common, , drop = FALSE]

# PCR controls identification
metadata$is.pcr_control <- grepl("^PCR_control(?:_|$)", metadata$group, ignore.case = TRUE)

# Build  phyloseq object
# OTU table: taxa are rows, samples are columns
otu_mat <- as.matrix(feature_table)          # samples x ASVs
otu_mat <- t(otu_mat)                        # now taxa x samples

# Taxonomy
tax_mat <- as.matrix(taxonomy[, "Taxon", drop = FALSE])
rownames(tax_mat) <- rownames(taxonomy)

ps_raw <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(tax_mat))

# Decontamination
otu_raw <- as(otu_table(ps_raw), "matrix")
if(taxa_are_rows(ps_raw)) otu_raw <- t(otu_raw)

contam_df <- isContaminant(
  otu_raw,
  neg = sample_data(ps_raw)$is.pcr_control,
  method = "prevalence",
  threshold = 0.5)

contaminants <- rownames(contam_df)[contam_df$contaminant == TRUE]

ps_clean <- prune_taxa(!taxa_names(ps_raw) %in% contaminants, ps_raw)

# Remove PCR control samples
ps_clean <- subset_samples(ps_clean, !is.pcr_control)

#  Remove singletons (ASVs present in < 2 samples)
ps_clean <- filter_taxa(ps_clean, function(x) sum(x > 0) >= 2, TRUE)

# Extract cleaned feature table and metadata for downstream analysis
feature_table <- as(otu_table(ps_clean), "matrix")

if(taxa_are_rows(ps_clean)) feature_table <- t(feature_table)
feature_table <- as.data.frame(feature_table)

metadata <- as(sample_data(ps_clean), "data.frame")

# Create SampleType and Species (now on the cleaned metadata)
metadata$SampleType <- ifelse(grepl("^Control_", metadata$group), "Environmental_control", "Sample")
species_temp <- gsub("^Control_", "", metadata$group)
species_temp <- gsub("_\\d+(\\.\\d+)?$", "", species_temp)
metadata$Species <- species_temp

# Verify that sample ordering is consistent
stopifnot(all(rownames(metadata) == rownames(feature_table)))

# Map ASV ID to label
get_best_name <- function(tax_string, asv_id) {
  if (is.na(tax_string) || tax_string == "") return(asv_id)
  if (grepl("chloroplast|mitochondria", tax_string, ignore.case = TRUE)) return("__EXCLUDED__")
  
  if (grepl("g__", tax_string)) {
    genus <- sub(".*g__", "", tax_string)
    genus <- strsplit(genus, ";")[[1]][1]
    genus <- trimws(genus)
    if (genus != "" && !grepl("uncultured", genus)) {
      return(paste(genus, asv_id, sep = " | "))}}
  if (grepl("f__", tax_string)) {
    family <- sub(".*f__", "", tax_string)
    family <- strsplit(family, ";")[[1]][1]
    family <- trimws(family)
    if (family != "" && !grepl("uncultured", family)) {
      return(paste("Family:", family, "|", asv_id))}}
  if (grepl("o__", tax_string)) {
    order <- sub(".*o__", "", tax_string)
    order <- strsplit(order, ";")[[1]][1]
    order <- trimws(order)
    if (order != "" && !grepl("uncultured", order)) {
      return(paste("Order:", order, "|", asv_id))}}
  if (grepl("c__", tax_string)) {
    class <- sub(".*c__", "", tax_string)
    class <- strsplit(class, ";")[[1]][1]
    class <- trimws(class)
    if (class != "" && !grepl("uncultured", class)) {
      return(paste("Class:", class, "|", asv_id))}}
  if (grepl("p__", tax_string)) {
    phylum <- sub(".*p__", "", tax_string)
    phylum <- strsplit(phylum, ";")[[1]][1]
    phylum <- trimws(phylum)
    if (phylum != "" && !grepl("uncultured", phylum)) {
      return(paste("Phylum:", phylum, "|", asv_id))}}
  if (grepl("d__", tax_string)) {
    domain <- sub(".*d__", "", tax_string)
    domain <- strsplit(domain, ";")[[1]][1]
    if (domain != "") return(paste("Domain:", domain, "|", asv_id))}
  return(asv_id)}

asv_to_label <- setNames(
  mapply(get_best_name, taxonomy$Taxon, rownames(taxonomy)),
  rownames(taxonomy))

# Loop over all species
q_threshold <- 0.05
species_list <- unique(metadata$Species)
all_results_combined <- list()

for (sp in species_list) {
  
  # Subset
  idx <- which(metadata$Species == sp)
  meta_sub <- metadata[idx, , drop = FALSE]
  feat_sub <- feature_table[idx, , drop = FALSE]
  
  # Filter low‑prevalence ASVs (present in ≥ 2 samples)
  feat_filt <- feat_sub[, colSums(feat_sub > 0) >= 2, drop = FALSE]
  
  # Prepare metadata for MaAsLin2
  maas_meta <- data.frame(SampleType = meta_sub$SampleType)
  rownames(maas_meta) <- rownames(meta_sub)
  maas_meta$SampleType <- factor(maas_meta$SampleType, levels = c("Environmental_control", "Sample"))
  
  outdir_sp <- file.path("3_maaslin2", gsub(" ", "_", sp))
  dir.create(outdir_sp, recursive = TRUE, showWarnings = FALSE)
  
  # Run MaAsLin2
  Maaslin2(
    input_data = feat_filt,
    input_metadata = maas_meta,
    output = outdir_sp,
    fixed_effects = "SampleType",
    reference = "SampleType,Environmental_control",
    standardize = FALSE,
    min_prevalence = 0,
    min_abundance = 0,
    plot_heatmap = FALSE,
    plot_scatter = FALSE)
  
  # Read results and apply labels
  res <- read.delim(file.path(outdir_sp, "all_results.tsv"), sep = "\t", stringsAsFactors = FALSE)
  res$ASV <- sub("^X", "", res$feature)
  res$Label <- asv_to_label[res$ASV]
  res$Label[is.na(res$Label)] <- res$ASV[is.na(res$Label)]
  res <- res[!grepl("^__EXCLUDED__", res$Label), ]
  res$Species <- sp
  all_results_combined[[sp]] <- res
  
  # Dot plot (only significant)
  res_sig <- res[res$qval < q_threshold, ]
  if (nrow(res_sig) > 0) {
    res_sig <- res_sig[order(res_sig$coef, decreasing = TRUE), ]
    res_sig$neglog10p <- -log10(res_sig$pval)
    p_sig <- ggplot(res_sig, aes(x = coef, y = reorder(Label, coef))) +
      geom_point(aes(size = neglog10p, color = coef > 0)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("TRUE" = "red3", "FALSE" = "steelblue"),
                         labels = c("More in Control", "More in Sample")) +
      labs(title = paste0(sp, " - Significant ASVs (q < ", q_threshold, ")"),
           x = "Coefficient", y = "Taxon", size = expression(-log[10](p)),
           color = "Direction") +
      theme_bw() + theme(axis.text.y = element_text(size = 8))
    
    height_sig <- max(6, nrow(res_sig) * 0.2)
    ggsave(file.path(outdir_sp, "dotplot_significant.png"), p_sig, width = 10, height = height_sig, dpi = 300)
  }
}

# Combine results from all species
combined_df <- do.call(rbind, all_results_combined)
combined_df <- combined_df[order(combined_df$Species, combined_df$qval), ]

# Function to extract ranks
parse_taxonomy_ranks <- function(tax_string) {
  ranks <- c("Kingdom" = "d__", "Phylum" = "p__", "Class" = "c__", 
             "Order" = "o__", "Family" = "f__", "Genus" = "g__", "Species" = "s__")
  result <- setNames(rep(NA, length(ranks)), names(ranks))
  if (is.na(tax_string) || tax_string == "") return(result)
  parts <- strsplit(tax_string, ";")[[1]]
  for (part in parts) {
    part <- trimws(part)
    for (rank_name in names(ranks)) {
      prefix <- ranks[rank_name]
      if (grepl(prefix, part)) {
        value <- sub(prefix, "", part)
        value <- trimws(value)
        if (value != "" && !grepl("^uncultured|^unclassified|^Incertae", value, ignore.case = TRUE)) {
          result[rank_name] <- value}
        break}}}
  return(result)}

# Parse all taxonomy strings
tax_parsed <- as.data.frame(t(sapply(taxonomy$Taxon, parse_taxonomy_ranks)))
tax_parsed$ASV <- rownames(taxonomy)  # assuming rownames are ASV IDs

# Function to get the best available name from taxonomy columns
get_best_name <- function(df) {
  ranks <- c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  best <- apply(df[, ranks, drop = FALSE], 1, function(row) {
    for (r in ranks) {
      val <- row[r]
      if (!is.na(val) && val != "" && !grepl("^uncultured|^unclassified|^Incertae", val, ignore.case = TRUE)) {
        return(paste0(r, ": ", val))}}
    return("Unclassified")})
  return(best)}

# Species vs Control
species_list <- unique(metadata$Species)
all_species_results <- list()

for (sp in species_list) {
  folder <- file.path("3_maaslin2", gsub(" ", "_", sp))
  res_file <- file.path(folder, "all_results.tsv")
  
  if (!file.exists(res_file)) {
    warning("File not found: ", res_file)
    next}
  
  res <- read.delim(res_file, sep = "\t", stringsAsFactors = FALSE)
  res$ASV <- sub("^X", "", res$feature)
  res$Species <- sp
  
  # Merge with taxonomy
  res <- merge(res, tax_parsed, by = "ASV", all.x = TRUE)
  
  # Rename conflict columns
  if ("Species.x" %in% colnames(res)) {
    res$Species <- res$Species.x
    res$Species.x <- NULL
  }
  if ("Species.y" %in% colnames(res)) {
    res$Tax_Species <- res$Species.y
    res$Species.y <- NULL
  }
  
  # Filter significant
  res_sig <- res[res$qval < q_threshold & !is.na(res$qval), ]
  
  if (nrow(res_sig) > 0) {
    res_sig$Direction <- ifelse(res_sig$coef > 0, "Enriched in Sample", "Enriched in Control")
    all_species_results[[sp]] <- res_sig}}

if (length(all_species_results) > 0) {
  combined_species <- do.call(rbind, all_species_results)
  
  cols_to_select <- c("Species", "ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                      "coef", "pval", "qval", "Direction")
  
  table_species <- combined_species %>%
    dplyr::select(dplyr::any_of(cols_to_select)) %>%
    dplyr::arrange(Species, qval) %>%
    dplyr::rename(log2FC = coef, p_value = pval, FDR_q_value = qval)
  
  table_species$log2FC <- round(table_species$log2FC, 3)
  table_species$p_value <- formatC(table_species$p_value, format = "e", digits = 3)
  table_species$FDR_q_value <- formatC(table_species$FDR_q_value, format = "e", digits = 3)
  
  # Add BestName
  table_species$BestName <- get_best_name(table_species)
  
  # Curation: Exclude chloroplasts, mitochondria, vague taxa
  
  # Convert FDR to numeric before filtering
  table_species$FDR_q_value <- as.numeric(table_species$FDR_q_value)
  
  # Taxa to remove
  exclude_pattern <- paste(
    "Chloroplast",
    "Mitochondria",
    "Kingdom:\\s*Bacteria",
    "^\\s*Unclassified\\s*$",
    "^\\s*Unassigned\\s*$",
    sep = "|")
  
  # Curate and separate results
  enriched <- table_species %>%
    filter(
      Direction == "Enriched in Sample",
      !is.na(BestName),
      !grepl(exclude_pattern, BestName, ignore.case = TRUE))
  
  negative <- table_species %>%
    filter(
      Direction == "Enriched in Control",
      !is.na(BestName),
      !grepl(exclude_pattern, BestName, ignore.case = TRUE))
  
  # Export complete curated tables
  write.csv(
    enriched,
    "3_maaslin2/Curated_Enriched_in_Sample.csv",
    row.names = FALSE)
  
  write.csv(
    negative,
    "3_maaslin2/Curated_Depleted_in_Sample.csv",
    row.names = FALSE)
  
  # Counts per species and direction
  counts <- bind_rows(enriched, negative) %>%
    count(Species, Direction, name = "Count") %>%
    tidyr::pivot_wider(
      names_from = Direction,
      values_from = Count,
      values_fill = 0) %>%
    mutate(
      Total = `Enriched in Sample` + `Enriched in Control`)
  
  write.csv(
    counts,
    "3_maaslin2/Counts_Enriched_Depleted.csv",
    row.names = FALSE)
  
  # Filter by q
  enriched_sig <- enriched %>% filter(FDR_q_value < q_threshold)
  negative_sig <- negative %>% filter(FDR_q_value < q_threshold)
  
  cat("Enriched (q <", q_threshold, "):", nrow(enriched_sig), "\n")
  cat("Depleted (q <", q_threshold, "):", nrow(negative_sig), "\n")
  
  if (nrow(enriched_sig) == 0 && nrow(negative_sig) == 0) {
    cat("\nNo ASVs pass the q-threshold. Here are the q-values after conversion:\n")
    print(summary(table_species$FDR_q_value))
    stop("Try increasing q_threshold or check curation filters.")}
  
  # Top 30 per species/direction by absolute log2FC
  top_enriched <- enriched_sig %>%
    group_by(Species) %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = 30) %>%
    ungroup()
  
  top_depleted <- negative_sig %>%
    group_by(Species) %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = 30) %>%
    ungroup()
  
  # Combine
  top_combined <- bind_rows(
    top_enriched %>% mutate(Direction = "Enriched in Sample"),
    top_depleted %>% mutate(Direction = "Enriched in Control"))
  
  # Prepare for plotting
  top_combined <- top_combined %>%
    mutate(Label = paste0(BestName, " (", Species, ")")) %>%
    group_by(Species, Direction) %>%
    mutate(order_val = rank(log2FC, ties.method = "first")) %>%
    ungroup()
  
  # Separate dotplots for Enriched and Depleted
  
  # Enriched plot
  plot_enriched <- top_enriched %>%
    group_by(Species) %>%
    mutate(Label = paste0(BestName, " (", Species, ")"),
           Label = forcats::fct_reorder(Label, log2FC, .desc = FALSE)) %>%   # < highest at top
    ungroup()
  
  p_enriched <- ggplot(plot_enriched, aes(x = log2FC, y = Label)) +
    geom_point(aes(size = -log10(FDR_q_value), colour = Species)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_grid(Species ~ ., scales = "free_y", space = "free_y") +
    scale_colour_brewer(palette = "Set2") +
    labs(title = paste0("Enriched in Sample – Top ASVs (q < ", q_threshold, ")"),
         x = "log₂ Fold Change", y = "Taxon", size = expression(-log[10](FDR~q))) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 7),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave("3_maaslin2/Dotplot_Enriched.png", p_enriched,
         width = 14, height = max(8, nrow(plot_enriched) * 0.12), dpi = 300)
  
  # Depleted plot
  plot_depleted <- top_depleted %>%
    group_by(Species) %>%
    mutate(Label = paste0(BestName, " (", Species, ")"),
           Label = forcats::fct_reorder(Label, log2FC, .desc = FALSE)) %>%   # highest at top
    ungroup()
  
  p_depleted <- ggplot(plot_depleted, aes(x = log2FC, y = Label)) +
    geom_point(aes(size = -log10(FDR_q_value), colour = Species)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_grid(Species ~ ., scales = "free_y", space = "free_y") +
    scale_colour_brewer(palette = "Set2") +
    labs(title = paste0("Depleted in Sample (Enriched in Control) – Top ASVs (q < ", q_threshold, ")"),
         x = "log₂ Fold Change", y = "Taxon", size = expression(-log[10](FDR~q))) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 7),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave("3_maaslin2/Dotplot_Depleted.png", p_depleted,
         width = 14, height = max(8, nrow(plot_depleted) * 0.12), dpi = 300)
}
