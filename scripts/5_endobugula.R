library(qiime2R)
library(phyloseq)
library(decontam)
library(Biostrings)
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)

asv_file      <- "table.qza"
taxonomy_file <- "taxonomy.qza"
tree_file     <- "rooted_tree.qza"
seq_file      <- "aligned_rep_seq.qza"
metadata_file <- "metadata.csv"

outdir <- "5_endobugula"
if (!dir.exists(outdir)) dir.create(outdir)

# Lead files

asv16 <- read_qza(asv_file)
asv_table <- asv16$data

taxonomy16 <- read_qza(taxonomy_file)
tree16 <- read_qza(tree_file)

# Load aligned sequences
aligned_qza <- read_qza("aligned_rep_seq.qza")

# EXTRACT the unmasked DNAStringSet
seqs_with_gaps <- aligned_qza$data@unmasked

# Remove the gaps (-) to get ASV sequences
raw_seqs <- DNAStringSet(gsub("-", "", as.character(seqs_with_gaps)))

# Preserve the ASV IDs from the alignment
names(raw_seqs) <- names(seqs_with_gaps)  # This gets the hash IDs

metadata <- read.csv(metadata_file, sep = ";", header = TRUE, check.names = FALSE)
colnames(metadata)[1] <- "sample.id"
metadata$sample.id <- trimws(metadata$sample.id)
rownames(metadata) <- metadata$sample.id

# PCR controls identification
metadata$is.pcr_control <- grepl("^PCR_control(?:_|$)", metadata$group, ignore.case = TRUE)

# Build phyloseq object
ps_raw <- phyloseq(
  otu_table(asv_table, taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(as.matrix(parse_taxonomy(taxonomy16$data))))

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

# Extract genus
genus_name <- "Candidatus_Endobugula"  

# Subset phyloseq to only this genus
ps_genus <- subset_taxa(ps_clean, Genus == genus_name)

# Check if any ASVs belong to this genus
if(ntaxa(ps_genus) == 0) {
  stop("Genus '", genus_name, "' not found! Check your taxonomy spelling.")} else {
  cat("Found", ntaxa(ps_genus), "ASVs for genus", genus_name, "\n")}

# Extract and save the matching sequences
genus_ids <- taxa_names(ps_genus)
genus_seqs <- raw_seqs[intersect(genus_ids, names(raw_seqs))]
if(length(genus_seqs) == 0) stop("Endobugula ASV IDs were not found in aligned_rep_seq.qza.")
writeXStringSet(genus_seqs, file.path(outdir, paste0(genus_name, "_ASVs.fasta")))

# Identification of samples where ASVs originated from
# ADD THE MISSING GROUPING VARIABLES TO ps_clean
sample_data(ps_clean)$group_orig <- sample_data(ps_clean)$group

# SampleType (Environment control vs real sample)
sample_data(ps_clean)$SampleType <- ifelse(
  grepl("^Control_", sample_data(ps_clean)$group_orig),
  "Environmental_control", 
  "Sample")

# Extract Species name (remove "Control_" prefix and trailing _number)
tmp <- gsub("^Control_", "", sample_data(ps_clean)$group_orig)
tmp <- gsub("_\\d+(\\.\\d+)?$", "", tmp)
sample_data(ps_clean)$Species <- tmp

# Verify it worked
print(colnames(sample_data(ps_clean)))

# Find Endobugula ASVs
ids_found <- rownames(tax_table(ps_clean))[grep("Endobugula", tax_table(ps_clean)[, "Genus"], ignore.case = TRUE)]
cat("Found", length(ids_found), "Endobugula ASVs\n")

if(length(ids_found) == 0) stop("No Endobugula found.")

# Extract abundance matrix
abundance_matrix <- as(otu_table(ps_clean), "matrix")
endo_abund <- abundance_matrix[ids_found, , drop = FALSE]

# Melt to long format (this is the safe way)
long_df <- melt(endo_abund, 
                varnames = c("ASV_ID", "Sample_ID"), 
                value.name = "Count")

# Remove zeros
long_df <- long_df[long_df$Count > 0, ]

# Add metadata (NOW these columns EXIST in ps_clean)
sam_data <- data.frame(sample_data(ps_clean))
long_df$Host_Species <- sam_data[long_df$Sample_ID, "Species"]
long_df$Original_Group <- sam_data[long_df$Sample_ID, "group_orig"]
long_df$Sample_Type <- sam_data[long_df$Sample_ID, "SampleType"]

# Calculate Relative Abundance (%)
sample_totals <- colSums(abundance_matrix)
long_df$Relative_Abundance_Percent <- (long_df$Count / sample_totals[long_df$Sample_ID]) * 100
long_df$Relative_Abundance_Percent <- round(long_df$Relative_Abundance_Percent, 4)

# Sort and print
long_df <- long_df[order(long_df$Host_Species, long_df$Sample_ID), ]

# Summary per host species
summary_df <- long_df %>%
  group_by(Host_Species) %>%
  summarise(
    Total_Reads = sum(Count),
    Num_ASVs = n_distinct(ASV_ID),
    Avg_Rel_Abund = mean(Relative_Abundance_Percent))
print(summary_df)

# Calculare procentage in each sample

# Reshape to wide format: rows = ASV_ID, columns = Sample_ID, values = Relative_Abundance_Percent
percent_wide <- dcast(long_df, ASV_ID ~ Sample_ID, value.var = "Relative_Abundance_Percent", fill = 0)

# Add the host species and taxonomy for each ASV (optional but useful)
# Get taxonomy for these ASVs
tax_subset <- tax_table(ps_clean)[ids_found, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
tax_df <- as.data.frame(tax_subset)
tax_df$ASV_ID <- rownames(tax_df)

# Merge with percent_wide
percent_wide <- merge(tax_df, percent_wide, by = "ASV_ID", all.y = TRUE)

# Reorder columns: put ASV_ID and taxonomy first, then samples
sample_cols <- colnames(percent_wide)[!colnames(percent_wide) %in% c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
percent_wide <- percent_wide[, c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", sample_cols)]

# Save to CSV
write.csv(percent_wide, "5_endobugula/Endobugula_ASV_Percent_per_Sample.csv", row.names = FALSE)

# BARPLOT OF Endobugula READS PER SAMPLE

# Aggregate total Endobugula reads per sample (summing across ASVs)
endo_per_sample <- long_df %>%
  group_by(Sample_ID, Host_Species, Sample_Type) %>%
  summarise(Total_Endobugula_Reads = sum(Count), .groups = "drop")

# Sort by Sample_Type and Total_Endobugula_Reads for a clean plot
endo_per_sample <- endo_per_sample %>%
  arrange(Sample_Type, Total_Endobugula_Reads)

# Keep the order for the plot (samples will appear from lowest to highest reads)
endo_per_sample$Sample_ID <- factor(endo_per_sample$Sample_ID, 
                                    levels = endo_per_sample$Sample_ID)
print(endo_per_sample)

write.csv(endo_per_sample, "5_endobugula/Endobugula_Reads_per_Sample_Exact.csv", row.names = FALSE)

# Create horizontal barplot with exact numbers on the bars

p <- ggplot(endo_per_sample, aes(x = Sample_ID, y = Total_Endobugula_Reads, fill = Sample_Type)) +
  geom_col() +
  geom_text(aes(label = Total_Endobugula_Reads), 
            hjust = -0.2,
            size = 3.5,
            color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Environmental_control" = "#F39C12", "Sample" = "#3498DB")) +
  labs(title = "Total Endobugula reads per sample",
       x = "Sample ID",
       y = "Endobugula reads (raw count)",
       fill = "Sample Type") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Print the plot
print(p)

# Save the plot
ggsave("5_endobugula/Endobugula_Reads_per_Sample_Barplot.png", p, width = 12, height = 10, dpi = 300)

#Plot reads per 2 hosts

# Define the two ASV IDs and their host species
asv_info <- list(
  "9be493e9aea5ea556ea35a3ff8fc61c8" = list(
    species = "Aquiloniella scabra",
    short   = "9be493"),
  "44cc17826fd7cbbd611d6ba70982a3e9" = list(
    species = "Rhamphostomella bilaminata",
    short   = "44cc17"))

# Store plots in a list
plot_list <- list()

for (asv in names(asv_info)) {
  
  # Subset to this ASV and samples with count > 0
  asv_data <- long_df[long_df$ASV_ID == asv & long_df$Count > 0, ]
  
  if (nrow(asv_data) == 0) {
    warning(paste("No data for ASV", asv))
    next
  }
  
  # Remove any samples whose ID starts with "P"
  asv_data <- asv_data[!grepl("^P", asv_data$Sample_ID, ignore.case = TRUE), ]
  
  # Split into Samples and Environmental Controls
  samples_df <- asv_data[asv_data$Sample_Type == "Sample", ]
  controls_df <- asv_data[asv_data$Sample_Type == "Environmental_control", ]
  
  # Sort Samples by ID ascending
  samples_df <- samples_df[order(as.character(samples_df$Sample_ID), decreasing = FALSE), ]
  
  # Sort Controls by ID ascending
  controls_df <- controls_df[order(as.character(controls_df$Sample_ID), decreasing = FALSE), ]
  
  # Combine data frames (order of rows doesn't matter for plotting)
  combined_df <- rbind(samples_df, controls_df)
  
  # Build the correct level order
  control_ids <- sort(unique(as.character(controls_df$Sample_ID)))
  sample_ids_desc <- rev(sort(unique(as.character(samples_df$Sample_ID))))
  desired_levels <- c(control_ids, sample_ids_desc)
  
  # Apply the levels
  combined_df$Sample_ID <- factor(combined_df$Sample_ID, levels = desired_levels)
  
  # Prepare labels
  species_name <- asv_info[[asv]]$species
  asv_short    <- asv_info[[asv]]$short
  
  # Create horizontal barplot
  p_asv <- ggplot(combined_df, aes(x = Sample_ID, y = Count, fill = Sample_Type)) +
    geom_col() +
    geom_text(aes(label = Count), 
              hjust = -0.2, 
              size = 3.5, 
              color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Environmental_control" = "#F39C12", "Sample" = "#3498DB")) +
    labs(title = paste(species_name, "– ASV", asv_short),
         x = NULL,
         y = "Read count",
         fill = "Sample Type") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          panel.grid.major.y = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Store the plot
  plot_list[[asv_short]] <- p_asv}

combined_plot <- (plot_list[["9be493"]] | plot_list[["44cc17"]]) +
  plot_annotation(
    title = "Ca. Endobugula abundance in host species",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))) &
  theme(legend.position = "bottom")

# Print and save
print(combined_plot)
ggsave("5_endobugula/Combined_ASV_Reads.png", combined_plot, width = 14, height = 10, dpi = 300)
