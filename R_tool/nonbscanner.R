################################################################################
# NonBScanner R Interface
# 
# Description: R wrapper for NonBScanner Python tool
#              Detects Non-B DNA motifs using the Python scanner
#
# Author: Dr. Venkata Rajesh Yella
# Version: 1.0
# License: MIT
################################################################################

#' Initialize NonBScanner
#'
#' Check Python availability and import required modules
#'
#' @return TRUE if successful, FALSE otherwise
#' @export
init_nonbscanner <- function() {
  # Check if reticulate is installed
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    message("Installing reticulate package...")
    install.packages("reticulate")
  }
  
  library(reticulate)
  
  # Check if Python is available
  if (!py_available()) {
    stop("Python is not available. Please install Python 3.8+")
  }
  
  # Try to import scanner module
  tryCatch({
    scanner <<- import("scanner")
    message("✓ NonBScanner initialized successfully")
    return(TRUE)
  }, error = function(e) {
    stop("Failed to import scanner module. Make sure scanner.py is in the Python path.\n", 
         "Error: ", e$message)
  })
}

#' Analyze DNA Sequence
#'
#' Detect Non-B DNA motifs in a single sequence
#'
#' @param sequence Character string containing DNA sequence
#' @param sequence_name Name identifier for the sequence (default: "sequence")
#' @return Data frame containing detected motifs
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' seq <- "GGGTTAGGGTTAGGGTTAGGG"
#' motifs <- analyze_sequence(seq, "test_seq")
#' print(motifs)
#' }
analyze_sequence <- function(sequence, sequence_name = "sequence") {
  if (!exists("scanner")) {
    stop("NonBScanner not initialized. Run init_nonbscanner() first.")
  }
  
  # Convert sequence to uppercase
  sequence <- toupper(sequence)
  
  # Call Python function
  motifs <- scanner$analyze_sequence(sequence, sequence_name)
  
  # Convert to R data frame
  df <- py_to_r(scanner$export_results_to_dataframe(motifs))
  
  return(df)
}

#' Analyze Multiple Sequences
#'
#' Detect Non-B DNA motifs in multiple sequences
#'
#' @param sequences Named list of sequences (names will be used as sequence IDs)
#' @param use_multiprocessing Logical, whether to use parallel processing (default: FALSE)
#' @return Data frame containing all detected motifs from all sequences
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' seqs <- list(
#'   seq1 = "GGGTTAGGGTTAGGGTTAGGG",
#'   seq2 = "AAAAATTTTAAAAATTTT"
#' )
#' motifs <- analyze_multiple_sequences(seqs)
#' print(motifs)
#' }
analyze_multiple_sequences <- function(sequences, use_multiprocessing = FALSE) {
  if (!exists("scanner")) {
    stop("NonBScanner not initialized. Run init_nonbscanner() first.")
  }
  
  # Convert sequences to uppercase
  sequences <- lapply(sequences, toupper)
  
  # Convert R list to Python dict
  py_sequences <- r_to_py(sequences)
  
  # Call Python function
  results <- scanner$analyze_multiple_sequences(py_sequences, 
                                               use_multiprocessing = use_multiprocessing)
  
  # Combine all results
  all_motifs <- list()
  for (seq_name in names(sequences)) {
    motifs <- results[[seq_name]]
    if (length(motifs) > 0) {
      all_motifs <- c(all_motifs, motifs)
    }
  }
  
  # Convert to data frame
  if (length(all_motifs) > 0) {
    df <- py_to_r(scanner$export_results_to_dataframe(all_motifs))
    return(df)
  } else {
    message("No motifs detected")
    return(data.frame())
  }
}

#' Read FASTA File
#'
#' Read sequences from a FASTA file
#'
#' @param filename Path to FASTA file
#' @return Named list of sequences
#' @export
#'
#' @examples
#' \dontrun{
#' sequences <- read_fasta("sequences.fasta")
#' }
read_fasta <- function(filename) {
  if (!file.exists(filename)) {
    stop("File not found: ", filename)
  }
  
  sequences <- list()
  current_name <- NULL
  current_seq <- c()
  
  lines <- readLines(filename)
  
  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, ">")) {
      # Save previous sequence
      if (!is.null(current_name)) {
        sequences[[current_name]] <- paste(current_seq, collapse = "")
      }
      # Start new sequence
      current_name <- strsplit(substring(line, 2), " ")[[1]][1]
      current_seq <- c()
    } else {
      current_seq <- c(current_seq, toupper(line))
    }
  }
  
  # Save last sequence
  if (!is.null(current_name)) {
    sequences[[current_name]] <- paste(current_seq, collapse = "")
  }
  
  return(sequences)
}

#' Analyze FASTA File
#'
#' Read and analyze sequences from a FASTA file
#'
#' @param filename Path to FASTA file
#' @param use_multiprocessing Logical, whether to use parallel processing (default: FALSE)
#' @return Data frame containing detected motifs
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' motifs <- analyze_fasta("sequences.fasta")
#' write.csv(motifs, "results.csv", row.names = FALSE)
#' }
analyze_fasta <- function(filename, use_multiprocessing = FALSE) {
  sequences <- read_fasta(filename)
  message("Loaded ", length(sequences), " sequence(s)")
  
  motifs <- analyze_multiple_sequences(sequences, use_multiprocessing)
  
  return(motifs)
}

#' Get Motif Classification Information
#'
#' Get information about supported motif classes and subclasses
#'
#' @return List containing classification information
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' info <- get_motif_info()
#' print(info)
#' }
get_motif_info <- function() {
  if (!exists("scanner")) {
    stop("NonBScanner not initialized. Run init_nonbscanner() first.")
  }
  
  info <- scanner$get_motif_classification_info()
  return(py_to_r(info))
}

#' Plot Motif Distribution
#'
#' Create a bar plot of motif class distribution
#'
#' @param motifs Data frame of detected motifs
#' @param title Plot title (default: "Non-B DNA Motif Distribution")
#' @return ggplot2 plot object
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' motifs <- analyze_sequence("GGGTTAGGGTTAGGGTTAGGG")
#' plot_motif_distribution(motifs)
#' }
plot_motif_distribution <- function(motifs, title = "Non-B DNA Motif Distribution") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Install with: install.packages('ggplot2')")
  }
  
  library(ggplot2)
  
  if (nrow(motifs) == 0) {
    message("No motifs to plot")
    return(NULL)
  }
  
  # Count by class
  class_counts <- as.data.frame(table(motifs$Class))
  names(class_counts) <- c("Class", "Count")
  
  # Create plot
  p <- ggplot(class_counts, aes(x = reorder(Class, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +
    labs(title = title, x = "Motif Class", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  return(p)
}

#' Plot Motif Positions
#'
#' Create a position plot showing motif locations
#'
#' @param motifs Data frame of detected motifs
#' @param sequence_length Length of the sequence (optional)
#' @param title Plot title (default: "Non-B DNA Motif Positions")
#' @return ggplot2 plot object
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' motifs <- analyze_sequence("GGGTTAGGGTTAGGGTTAGGG")
#' plot_motif_positions(motifs)
#' }
plot_motif_positions <- function(motifs, sequence_length = NULL, 
                                title = "Non-B DNA Motif Positions") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Install with: install.packages('ggplot2')")
  }
  
  library(ggplot2)
  
  if (nrow(motifs) == 0) {
    message("No motifs to plot")
    return(NULL)
  }
  
  # Convert Start and End to numeric
  motifs$Start <- as.numeric(as.character(motifs$Start))
  motifs$End <- as.numeric(as.character(motifs$End))
  
  # Add index for y-axis
  motifs$Index <- 1:nrow(motifs)
  
  # Create plot
  p <- ggplot(motifs, aes(y = Index)) +
    geom_segment(aes(x = Start, xend = End, yend = Index, color = Class),
                size = 3) +
    labs(title = title, x = "Position (bp)", y = "Motif Index", color = "Class") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right")
  
  if (!is.null(sequence_length)) {
    p <- p + xlim(0, sequence_length)
  }
  
  return(p)
}

#' Export Summary Statistics
#'
#' Generate summary statistics for detected motifs
#'
#' @param motifs Data frame of detected motifs
#' @return Data frame with summary statistics
#' @export
#'
#' @examples
#' \dontrun{
#' init_nonbscanner()
#' motifs <- analyze_sequence("GGGTTAGGGTTAGGGTTAGGG")
#' summary <- get_summary_stats(motifs)
#' print(summary)
#' }
get_summary_stats <- function(motifs) {
  if (nrow(motifs) == 0) {
    message("No motifs to summarize")
    return(data.frame())
  }
  
  summary <- data.frame(
    Total_Motifs = nrow(motifs),
    Unique_Classes = length(unique(motifs$Class)),
    Unique_Subclasses = length(unique(motifs$Subclass)),
    stringsAsFactors = FALSE
  )
  
  # Add class counts
  class_counts <- table(motifs$Class)
  for (cls in names(class_counts)) {
    summary[[paste0("Class_", cls)]] <- as.numeric(class_counts[cls])
  }
  
  return(summary)
}

# Package metadata
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("NonBScanner R Interface v1.0")
  packageStartupMessage("Non-B DNA Motif Detection Tool")
  packageStartupMessage("Author: Dr. Venkata Rajesh Yella")
  packageStartupMessage("")
  packageStartupMessage("Initialize with: init_nonbscanner()")
}
