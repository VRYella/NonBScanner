################################################################################
# NonBScanner R Example Script
# 
# This script demonstrates how to use NonBScanner from R
################################################################################

# Load required libraries
library(reticulate)

# Source the NonBScanner R interface
source("nonbscanner.R")

# Initialize NonBScanner
cat("Initializing NonBScanner...\n")
init_nonbscanner()

# ============================================================================
# Example 1: Analyze a single sequence
# ============================================================================
cat("\n=== Example 1: Single Sequence Analysis ===\n")

# Define a test sequence
test_seq <- "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTAAAAATTTTCGCGCGCGCGCGCACACACACACACACA"

# Analyze the sequence
motifs <- analyze_sequence(test_seq, "test_sequence")

# Display results
cat("\nFound", nrow(motifs), "motifs:\n")
print(head(motifs))

# ============================================================================
# Example 2: Analyze multiple sequences
# ============================================================================
cat("\n=== Example 2: Multiple Sequence Analysis ===\n")

# Define multiple sequences
sequences <- list(
  seq1 = "GGGTTAGGGTTAGGGTTAGGG",
  seq2 = "AAAAATTTTAAAAATTTT",
  seq3 = "CGCGCGCGCGCGCG",
  seq4 = "CACACACACACACACA"
)

# Analyze all sequences
all_motifs <- analyze_multiple_sequences(sequences)

cat("\nFound", nrow(all_motifs), "motifs across", length(sequences), "sequences\n")

# ============================================================================
# Example 3: Read and analyze FASTA file
# ============================================================================
cat("\n=== Example 3: FASTA File Analysis ===\n")

# Check if example FASTA exists
fasta_file <- "../example_all_motifs.fasta"

if (file.exists(fasta_file)) {
  # Read sequences
  fasta_sequences <- read_fasta(fasta_file)
  cat("Loaded", length(fasta_sequences), "sequence(s) from FASTA\n")
  
  # Analyze
  fasta_motifs <- analyze_multiple_sequences(fasta_sequences)
  cat("Detected", nrow(fasta_motifs), "motifs\n")
} else {
  cat("Example FASTA file not found. Skipping...\n")
}

# ============================================================================
# Example 4: Get motif classification information
# ============================================================================
cat("\n=== Example 4: Motif Classification Info ===\n")

info <- get_motif_info()
cat("NonBScanner Version:", info$version, "\n")
cat("Total Classes:", info$total_classes, "\n")
cat("Total Subclasses:", info$total_subclasses, "\n")

cat("\nAvailable Motif Classes:\n")
for (i in 1:length(info$classification)) {
  class_info <- info$classification[[i]]
  cat(sprintf("%2d. %s\n", i, class_info$name))
}

# ============================================================================
# Example 5: Generate summary statistics
# ============================================================================
cat("\n=== Example 5: Summary Statistics ===\n")

if (nrow(motifs) > 0) {
  summary_stats <- get_summary_stats(motifs)
  cat("\nSummary Statistics:\n")
  print(summary_stats)
}

# ============================================================================
# Example 6: Visualizations (requires ggplot2)
# ============================================================================
cat("\n=== Example 6: Visualizations ===\n")

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(motifs) > 0) {
  # Plot motif distribution
  cat("Creating motif distribution plot...\n")
  p1 <- plot_motif_distribution(motifs)
  print(p1)
  
  # Save plot
  ggsave("motif_distribution.png", p1, width = 10, height = 6)
  cat("Saved: motif_distribution.png\n")
  
  # Plot motif positions
  cat("Creating motif position plot...\n")
  p2 <- plot_motif_positions(motifs, sequence_length = nchar(test_seq))
  print(p2)
  
  # Save plot
  ggsave("motif_positions.png", p2, width = 10, height = 6)
  cat("Saved: motif_positions.png\n")
} else {
  cat("ggplot2 not available or no motifs detected. Skipping visualization...\n")
}

# ============================================================================
# Example 7: Export results to CSV
# ============================================================================
cat("\n=== Example 7: Export Results ===\n")

if (nrow(motifs) > 0) {
  # Export all motifs
  output_file <- "nonbscanner_results.csv"
  write.csv(motifs, output_file, row.names = FALSE)
  cat("Exported results to:", output_file, "\n")
  
  # Export summary
  summary_file <- "nonbscanner_summary.csv"
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("Exported summary to:", summary_file, "\n")
  
  # Export by class
  for (cls in unique(motifs$Class)) {
    if (cls != "NA") {
      class_motifs <- motifs[motifs$Class == cls, ]
      class_file <- paste0("nonbscanner_", gsub(" ", "_", cls), ".csv")
      write.csv(class_motifs, class_file, row.names = FALSE)
      cat("Exported", cls, "motifs to:", class_file, "\n")
    }
  }
}

cat("\n=== All Examples Completed ===\n")
