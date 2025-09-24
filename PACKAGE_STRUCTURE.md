# SPDtrace R Package Structure

This document provides an overview of the complete SPDtrace R package structure and contents.

## Package Overview

**Package Name**: SPDtrace  
**Version**: 1.0.0  
**Title**: Efficient Learning of Differential Networks in Multi-Source Non-Paranormal Graphical Models  
**Authors**: Mojtaba Nikahd, Seyed Abolfazl Motahari  
**License**: GPL-3  

## Directory Structure

```
SPDtrace_pkg/
├── DESCRIPTION                 # Package metadata and dependencies
├── NAMESPACE                  # Package namespace and exports
├── README.md                  # Package overview and quick start
├── INSTALL.md                 # Installation instructions
├── PACKAGE_STRUCTURE.md       # This file
├── build_package.R            # Build and installation script
│
├── R/                         # R source code
│   ├── SPDtrace.R            # Main SPDtrace function and methods
│   ├── utilities.R           # Utility functions
│   └── data.R                # Data documentation
│
├── src/                       # C++ source code
│   ├── SolutionPath.cpp      # Core SPD-trace implementation
│   └── CrossFDTL.cpp         # CrossFDTL algorithm implementation
│
├── man/                       # Manual pages (generated)
│
├── inst/                      # Package data and resources
│   └── extdata/              # (removed in current refactor)
│
├── vignettes/                 # Package vignettes
│   └── SPDtrace-introduction.Rmd  # Introduction and usage guide
│
└── tests/                     # Package tests
    ├── testthat.R            # Test configuration
    └── testthat/             # Test files
        └── test-SPDtrace.R   # Unit tests for main functions
```

## Core Components

### 1. Main Functions (`R/SPDtrace.R`)

- `SPDtrace()`: Main function for differential network inference
- `print.SPDtrace_result()`: Print method for results
- `summary.SPDtrace_result()`: Summary method for results

### 2. Utility Functions (`R/utilities.R`)

- `kendall_tau_matrix_to_pearson_correlation_matrix()`: Convert correlation matrices
- `positive_semi_definite_maker()`: Ensure matrix positive semi-definiteness
- `weighted_rank_based_pearson_correlation_estimator()`: Estimate correlations
- `instability_evaluator_of_solution_paths()`: Evaluate solution stability
- `fisher_test()`: Perform Fisher's exact test for enrichment

### 3. [Removed] Simulation Functions

Simulation helpers and scenarios have been removed to streamline dependencies.

### 4. C++ Implementation (`src/`)

- `SolutionPath.cpp`: Core SPD-trace algorithm implementation
- `CrossFDTL.cpp`: CrossFDTL algorithm implementation

### 5. Package Data (`inst/extdata/`)

- Removed bundled datasets; users should supply their own data.

## Dependencies

### Required Packages
- `Rcpp` (>= 1.0.0): C++ integration
- `RcppArmadillo` (>= 0.10.0): Linear algebra operations
- `dplyr` (>= 1.0.0): Data manipulation
- `ggplot2` (>= 3.3.0): Plotting
- `igraph` (>= 1.2.0): Network analysis

### Suggested Packages
- `testthat` (>= 3.0.0): Testing framework
- `knitr` (>= 1.30): Vignette generation
- `rmarkdown` (>= 2.5): Documentation

## Key Features

1. Efficient Implementation: C++ backend for fast computation
2. Multiple Methods: Support for SPD-trace, CrossFDTL, and D-trace
3. Real Data Analysis: Includes ovarian cancer data for validation
4. Comprehensive Testing: Unit tests for all functions
5. Professional Documentation: Vignettes, examples, and help pages

## Usage Examples

### Basic Differential Network Inference
```r
library(SPDtrace)
result <- SPDtrace(CovA, CovB, sparsityLevel = 50)
print(result)
summary(result)
```

### Utility Functions
```r
# Convert correlation matrices
pearson_cor <- kendall_tau_matrix_to_pearson_correlation_matrix(tau_matrix)

# Ensure positive semi-definiteness
A_psd <- positive_semi_definite_maker(A)
```

## Building and Installation

### Quick Build
```r
source("build_package.R")
```

### Manual Build
```r
library(roxygen2)
roxygenise()
library(devtools)
install()
```

## Testing

Run the test suite:
```r
library(devtools)
test()
```

## Documentation

- Vignette: `vignette("SPDtrace-introduction")`
- Function Help: `?SPDtrace`
- Package Help: `help(package = "SPDtrace")`

## Citation

If you use this package in your research, please cite:

```bibtex
@article{Nikahd2025Differential,
  title     = {Efficient Learning of Differential Networks in Multi-Source Non-Paranormal Graphical Models},
  author    = {Nikahd, Mojtaba and Motahari, Seyed Abolfazl},
  year      = {2025},
  note      = {Manuscript under review},
  url       = {https://doi.org/10.48550/arXiv.2410.02496}
}
```

## Support and Contributing

- Issues: Report bugs on GitHub
- Contributions: Pull requests welcome
- Contact: nikahd@ce.sharif.ir
- Repository: https://github.com/mojtabanikahd/SPDtrace

## License

This project is licensed under the GPL-3 License - see the LICENSE file for details.
