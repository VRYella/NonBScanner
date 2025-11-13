################################################################################
# NonBScanner R Tool Installation Script
# 
# This script checks and installs all required dependencies for the R interface
################################################################################

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║       NonBScanner R Interface Installation                    ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Check R version
r_version <- getRversion()
cat("R version:", as.character(r_version), "\n")

if (r_version < "4.0.0") {
  warning("R version 4.0.0 or higher is recommended")
}

# Function to install package if not present
install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    cat("Installing", package_name, "...\n")
    install.packages(package_name, repos = "https://cran.r-project.org")
    cat("✓", package_name, "installed\n")
  } else {
    cat("✓", package_name, "already installed\n")
  }
}

# Required packages
cat("\n=== Checking Required R Packages ===\n")
required_packages <- c("reticulate", "ggplot2")

for (pkg in required_packages) {
  install_if_missing(pkg)
}

# Check Python availability
cat("\n=== Checking Python ===\n")
library(reticulate)

if (py_available()) {
  py_version <- py_config()$version
  cat("✓ Python available:", py_version, "\n")
  
  # Check Python version
  py_version_info <- py_run_string("import sys; print(sys.version_info[:2])")
  
  # Check required Python packages
  cat("\n=== Checking Python Packages ===\n")
  
  python_packages <- c("numpy", "pandas", "matplotlib", "seaborn", "scipy")
  
  for (pkg in python_packages) {
    tryCatch({
      py_run_string(paste0("import ", pkg))
      cat("✓", pkg, "installed\n")
    }, error = function(e) {
      cat("✗", pkg, "NOT installed\n")
      cat("  Install with: pip install", pkg, "\n")
    })
  }
  
} else {
  cat("✗ Python not available\n")
  cat("  Please install Python 3.8+ and ensure it's in your PATH\n")
}

# Test NonBScanner module
cat("\n=== Testing NonBScanner Module ===\n")

# Try to add parent directory to Python path
tryCatch({
  py_run_string("import sys; sys.path.insert(0, '..')")
  scanner <- import("scanner")
  cat("✓ NonBScanner module loaded successfully\n")
  
  # Get version info
  info <- scanner$get_motif_classification_info()
  cat("  Version:", py_to_r(info)$version, "\n")
  cat("  Classes:", py_to_r(info)$total_classes, "\n")
  
}, error = function(e) {
  cat("✗ Failed to load NonBScanner module\n")
  cat("  Error:", conditionMessage(e), "\n")
  cat("  Make sure scanner.py is in the parent directory\n")
})

# Installation summary
cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║       Installation Summary                                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("To use NonBScanner in R:\n\n")
cat("1. Source the interface:\n")
cat("   source('nonbscanner.R')\n\n")
cat("2. Initialize:\n")
cat("   init_nonbscanner()\n\n")
cat("3. Analyze sequences:\n")
cat("   motifs <- analyze_sequence('GGGTTAGGGTTAGGGTTAGGG', 'test')\n\n")
cat("4. See example_usage.R for more examples\n\n")

cat("Documentation: See README.md in this directory\n")
