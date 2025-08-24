# Installation Guide for SPDtrace R Package

This guide provides step-by-step instructions for installing and building the SPDtrace R package.

## Prerequisites

### Required Software

1. **R** (version 4.0.0 or higher)
   - Download from: https://cran.r-project.org/
   - Install according to your operating system instructions

2. **C++ Compiler** (for Rcpp support)
   - **Windows**: Install Rtools from https://cran.r-project.org/bin/windows/Rtools/
   - **macOS**: Install Xcode Command Line Tools: `xcode-select --install`
   - **Linux**: Install build-essential: `sudo apt-get install build-essential`

### Required R Packages

Install the following packages in R:

```r
# Core packages for C++ compilation
install.packages("Rcpp")
install.packages("RcppArmadillo")

# Development tools
install.packages("devtools")
install.packages("roxygen2")
install.packages("testthat")

# Documentation and vignettes
install.packages("knitr")
install.packages("rmarkdown")

# Package dependencies
install.packages(c("dplyr", "stringr", "ggplot2", "ggpubr", "igraph", "Matrix"))
install.packages(c("rmatio", "readxl"))

# For Bioconductor packages (if needed)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GGMselect")
```

## Installation Methods

### Method 1: Install from Source (Recommended)

1. **Clone or download the package source**
   ```bash
   git clone https://github.com/mojtabanikahd/SPDtrace.git
   cd SPDtrace
   ```

2. **Build and install the package**
   ```r
   # In R, from the package directory
   library(devtools)
   install()
   ```

### Method 2: Use the Build Script

1. **Navigate to the package directory**
   ```bash
   cd SPDtrace_pkg
   ```

2. **Run the build script**
   ```r
   source("build_package.R")
   ```

### Method 3: Manual Build

1. **Generate documentation**
   ```r
   library(roxygen2)
   roxygenise()
   ```

2. **Check the package**
   ```r
   library(devtools)
   check()
   ```

3. **Build the package**
   ```r
   build()
   ```

4. **Install the package**
   ```r
   install()
   ```

## Verification

After installation, verify the package works correctly:

```r
# Load the package
library(SPDtrace)

# Check package information
packageVersion("SPDtrace")

# View available functions
ls("package:SPDtrace")

# Run basic example
?SPDtrace
```

## Troubleshooting

### Common Issues

1. **C++ compilation errors**
   - Ensure Rtools (Windows) or Xcode (macOS) is properly installed
   - Check that `Rcpp` package is installed and loaded
   - Verify C++ compiler is in your system PATH

2. **Package loading errors**
   - Update R to latest version
   - Install packages from CRAN: `install.packages("package_name")`
   - For Bioconductor packages: `BiocManager::install("package_name")`

3. **Memory issues**
   - Increase R memory limit: `memory.limit(size = 8000)`
   - Use garbage collection: `gc()`

4. **Permission errors**
   - Run R as administrator (Windows)
   - Check file permissions (Linux/macOS)

### Getting Help

- Check the package documentation: `?SPDtrace`
- View function examples: `example(SPDtrace)`
- Run vignette: `vignette("SPDtrace-introduction")`
- Report issues on the GitHub repository

## Next Steps

After successful installation:

1. **Read the vignette**: `vignette("SPDtrace-introduction")`
2. **Try the examples**: `example(SPDtrace)`
3. **Explore the functions**: `help(package = "SPDtrace")`
4. **Run simulations**: See `?run_first_scenario` and `?run_second_scenario`

## System Requirements

- **R**: >= 4.0.0
- **Memory**: >= 4GB RAM (8GB recommended for large networks)
- **Storage**: >= 100MB free space
- **Operating System**: Windows 10+, macOS 10.14+, or Linux (Ubuntu 18.04+)

## Support

For installation issues or questions:

- Check the troubleshooting section above
- Review the package documentation
- Open an issue on GitHub
- Contact the maintainer: nikahd@ce.sharif.ir
