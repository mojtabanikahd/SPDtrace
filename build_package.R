#!/usr/bin/env Rscript
# Build script for SPDtrace R package
# This script helps build and install the package

cat("Building SPDtrace R package...\n")

# Check if required packages are installed
required_packages <- c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cran.r-project.org")
}

# Load required packages
library(devtools)
library(roxygen2)

# Set working directory to package directory
package_dir <- getwd()
if (!file.exists("DESCRIPTION")) {
  stop("Please run this script from the package directory (where DESCRIPTION file is located)")
}

cat("Package directory:", package_dir, "\n")

# Clean previous builds
cat("Cleaning previous builds...\n")
if (dir.exists("man")) unlink("man", recursive = TRUE)
if (file.exists("NAMESPACE")) file.remove("NAMESPACE")

# Generate documentation
cat("Generating documentation...\n")
roxygen2::roxygenise()

# Check package
cat("Checking package...\n")
check_result <- devtools::check(quiet = TRUE)

if (length(check_result$errors) > 0) {
  cat("Package check failed with errors:\n")
  print(check_result$errors)
  stop("Package check failed")
}

if (length(check_result$warnings) > 0) {
  cat("Package check completed with warnings:\n")
  print(check_result$warnings)
} else {
  cat("Package check completed successfully!\n")
}

# Build package
cat("Building package...\n")
build_result <- devtools::build(quiet = TRUE)
cat("Package built successfully:", build_result, "\n")

# Install package
cat("Installing package...\n")
install_result <- devtools::install(quiet = TRUE)
cat("Package installed successfully!\n")

# Run tests
cat("Running tests...\n")
test_result <- devtools::test(quiet = TRUE)
cat("Tests completed!\n")

cat("\n=== Package Build Summary ===\n")
cat("Package built:", build_result, "\n")
cat("Package installed:", ifelse(install_result, "Yes", "No"), "\n")
cat("Tests passed:", test_result$n, "\n")
cat("Build completed successfully!\n")

# Instructions for use
cat("\n=== Next Steps ===\n")
cat("1. Load the package: library(SPDtrace)\n")
cat("2. View documentation: ?SPDtrace\n")
cat("3. Run vignette: vignette('SPDtrace-introduction')\n")
cat("4. Run examples: example(SPDtrace)\n")
