# Required libraries
# Install packages if not already installed
# install.packages(c("phyloseq", "vegan", "ape", "ggplot2", "reshape2", "dplyr", "ggpubr"))
library(phyloseq)
library(vegan)
library(ape)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

# Load OTU data for different taxonomy levels
load_taxonomy_data <- function(file_path) {
  read.csv(file_path, row.names = 1)
}

# Define file paths
file_base <- "C:/Users/Parth Doshi/Desktop/Project/B_GROUP_All/B_GROUP_All/"
file_paths <- list(
  Phylum = paste0(file_base, "Phylum-Level.csv"),
  Class = paste0(file_base, "Class-Level.csv"),
  Order = paste0(file_base, "Order-Level.csv"),
  Family = paste0(file_base, "Family-Level.csv"),
  Genus = paste0(file_base, "Genus-Level.csv"),
  Species = paste0(file_base, "Species-Level.csv")
)

# Load data into a list
taxonomy_data <- lapply(file_paths, load_taxonomy_data)

# Load the read statistics
read_stats <- read.csv(paste0(file_base, "ReadStatistics.csv"), row.names = 1)
total_reads <- read_stats$processed.read.count
names(total_reads) <- rownames(read_stats)

# Function to convert relative abundances to integer counts
convert_to_counts <- function(data, total_reads) {
  common_samples <- intersect(colnames(data), names(total_reads))
  data <- data[, common_samples]
  total_reads <- total_reads[common_samples]
  data <- sweep(data, 2, total_reads, "*")
  round(data)
}

# Apply the conversion to each taxonomy level
taxonomy_data <- lapply(taxonomy_data, convert_to_counts, total_reads)

# Calculate alpha diversity metrics
calculate_alpha_diversity <- function(data) {
  t_data <- t(data)
  diversity_data <- data.frame(
    Shannon_Diversity = diversity(t_data, index = "shannon"),
    Simpson_Evenness = diversity(t_data, index = "simpson"),
    Observed_Features = rowSums(t_data > 0)
  )
  return(diversity_data)
}

alpha_diversity <- lapply(taxonomy_data, calculate_alpha_diversity)

# Plotting function for alpha diversity
plot_alpha_diversity <- function(diversity_data, taxonomic_level, group, output_file) {
  diversity_data$SampleID <- rownames(diversity_data)
  diversity_melt <- melt(diversity_data, id.vars = "SampleID")
  
  plot <- ggplot(diversity_melt, aes(x = SampleID, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    labs(
      x = "Sample ID",
      y = "Alpha Diversity Metrics",
      title = paste("Alpha Diversity Metrics at", taxonomic_level, "Level for Group:", group)
    ) +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
  final_plot <- ggarrange(plot, ncol = 1, nrow = 1)
  ggsave(output_file, final_plot, width = 10, height = 8)
}

# Plot alpha diversity for each taxonomic level
output_files <- list(
  Phylum = "phylum_plot_B.png",
  Class = "class_plot_B.png",
  Order = "order_plot_B.png",
  Family = "family_plot_B.png",
  Genus = "genus_plot_B.png",
  Species = "species_plot_B.png"
)

mapply(plot_alpha_diversity, alpha_diversity, names(alpha_diversity), "B", output_files)

# Summarize diversity metrics
summarize_diversity <- function(alpha_diversity) {
  summary <- bind_rows(alpha_diversity, .id = "Taxonomic_Level") %>%
    group_by(Taxonomic_Level) %>%
    summarise(
      mean_shannon = mean(Shannon_Diversity, na.rm = TRUE),
      mean_simpson = mean(Simpson_Evenness, na.rm = TRUE),
      mean_observed = mean(Observed_Features, na.rm = TRUE)
    )
  return(summary)
}

diversity_summary <- summarize_diversity(alpha_diversity)
print(diversity_summary)

# Plot summarized diversity metrics
plot_diversity_summary <- function(diversity_summary) {
  diversity_melt <- melt(diversity_summary, id.vars = "Taxonomic_Level")
  
  ggplot(diversity_melt, aes(x = Taxonomic_Level, y = value, fill = Taxonomic_Level)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    labs(x = "Taxonomic Level", y = "Mean Diversity Metrics", title = "Mean Diversity Metrics by Taxonomic Level") +
    theme_minimal()
}

plot_diversity_summary(diversity_summary)

# Function to create a phyloseq object from OTU data
create_phyloseq <- function(otu_data) {
  otu_table <- otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)
  sample_data <- sample_data(data.frame(SampleID = colnames(otu_data), row.names = colnames(otu_data)))
  physeq <- phyloseq(otu_table, sample_data)
  return(physeq)
}

# Function to calculate beta diversity and perform PCoA for a given distance method
calculate_beta_diversity <- function(otu_data, method = "bray") {
  physeq <- create_phyloseq(otu_data)
  dist <- distance(physeq, method = method)
  pcoa <- ordinate(physeq, method = "PCoA", distance = dist)
  return(list(distance_matrix = dist, pcoa = pcoa, physeq = physeq))
}

# Function to plot PCoA with SampleID
plot_pcoa <- function(pcoa, physeq, title = "PCoA Plot") {
  sample_data_df <- data.frame(sample_data(physeq))
  pcoa_df <- data.frame(pcoa$vectors)
  pcoa_df$SampleID <- sample_data_df$SampleID
  
  eigenvalues <- pcoa$values$Eigenvalues
  variation_percent <- eigenvalues / sum(eigenvalues) * 100
  
  pcoa_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = SampleID)) + 
    geom_point(size = 3, alpha = 0.7) + 
    labs(
      title = title,
      x = paste0("PCoA Axis 1 (", round(variation_percent[1], 2), "%)"),
      y = paste0("PCoA Axis 2 (", round(variation_percent[2], 2), "%)")
    ) +
    theme_pubr() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 5)
    ) +
    scale_color_viridis_d() +
    guides(color = guide_legend(nrow = 8, byrow = TRUE))
  
  return(pcoa_plot)
}

# Calculate and plot beta diversity for species level using Bray-Curtis and Jaccard distances
bray_diversity_results <- calculate_beta_diversity(taxonomy_data$Species, method = "bray")
jaccard_diversity_results <- calculate_beta_diversity(taxonomy_data$Species, method = "jaccard")

# Extract components for Bray-Curtis
bray_dist <- bray_diversity_results$distance_matrix
bray_pcoa <- bray_diversity_results$pcoa
bray_physeq <- bray_diversity_results$physeq

# Extract components for Jaccard
jaccard_dist <- jaccard_diversity_results$distance_matrix
jaccard_pcoa <- jaccard_diversity_results$pcoa
jaccard_physeq <- jaccard_diversity_results$physeq

# Plot PCoA with SampleID for Bray-Curtis
bray_pcoa_plot <- plot_pcoa(bray_pcoa, bray_physeq, title = "PCoA Plot (Bray-Curtis) for Group B")
ggsave("Bray_curtic_pcoa_B.png", bray_pcoa_plot, width = 14, height = 12)

# Plot PCoA with SampleID for Jaccard
jaccard_pcoa_plot <- plot_pcoa(jaccard_pcoa, jaccard_physeq, title = "PCoA Plot (Jaccard) for Group B")
ggsave("jaccard_pcoa_B.png", jaccard_pcoa_plot, width = 14, height = 12)

# Print beta diversity results
print(bray_dist)
print(bray_pcoa)
print(bray_pcoa_plot)

print(jaccard_dist)
print(jaccard_pcoa)
print(jaccard_pcoa_plot)

# Display the plots
print(bray_pcoa_plot)
print(jaccard_pcoa_plot)
