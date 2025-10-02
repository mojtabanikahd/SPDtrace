# SPDtrace: Efficient Learning of Differential Networks

[![R-CMD-check](https://github.com/mojtabanikahd/SPDtrace/workflows/R-CMD-check/badge.svg)](https://github.com/mojtabanikahd/SPDtrace/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/SPDtrace)](https://CRAN.R-project.org/package=SPDtrace)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Overview

The **SPDtrace** package implements the SPD-trace method for efficient learning of differential networks in multi-source non-paranormal graphical models. This package provides tools for differential network inference, simulation studies, and real-world applications in genomics and other high-dimensional data analysis.

## Features

- **SPD-trace Method**: Core implementation of the proposed SPD-trace algorithm
- **Multiple Methods**: Support for SPD-trace, CrossFDTL, and D-trace methods
- **Simulation Tools**: Built-in functions for performance evaluation and robustness testing
- **Real Data Analysis**: Applicable to real datasets for biological validation
- **High Performance**: C++ implementation for efficient computation
- **Comprehensive Documentation**: Detailed vignettes and examples

## Installation

### From GitHub

```r
# Install devtools if you haven't already
if (!require(devtools)) install.packages("devtools")

# Install from GitHub
devtools::install_github("mojtabanikahd/SPDtrace")
```

### Dependencies

The package requires the following R packages:
- Rcpp (>= 1.0.0)
- RcppArmadillo (>= 0.10.0)
- dplyr (>= 1.0.0)
- ggplot2 (>= 3.3.0)
- igraph (>= 1.2.0)
- And several others (see DESCRIPTION file)

## Quick Start

```r
library(SPDtrace)

# Generate example data
set.seed(123)
n <- 100; p <- 20
Sigma_A <- diag(p) + 0.3 * (matrix(runif(p^2), p, p) > 0.8)
Sigma_A <- (Sigma_A + t(Sigma_A)) / 2; diag(Sigma_A) <- 1

Sigma_B <- Sigma_A + 0.2 * (matrix(runif(p^2), p, p) > 0.9)
Sigma_B <- (Sigma_B + t(Sigma_B)) / 2; diag(Sigma_B) <- 1

# Run SPD-trace
result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 50)

# View results
print(result)
summary(result)
```

## Core Functions

### Main Functions

- `SPDtrace()`: Main function for differential network inference
  

### Utility Functions

- `kendall_tau_matrix_to_pearson_correlation_matrix()`: Convert correlation matrices
- `positive_semi_definite_maker()`: Ensure matrix positive semi-definiteness
- `weighted_rank_based_pearson_correlation_estimator()`: Estimate correlations
- `instability_evaluator_of_solution_paths()`: Evaluate solution stability
 

## Documentation

- **Vignette**: `vignette("SPDtrace-introduction")`
- **Function Help**: `?SPDtrace`
- **Package Help**: `help(package = "SPDtrace")`

## Examples

### Basic Differential Network Inference

```r
# Construct two covariance matrices (example)
set.seed(1)
p <- 10
cov_matrix_A <- diag(p)
cov_matrix_B <- cov_matrix_A + 0.1 * (matrix(runif(p^2), p, p) > 0.9)
cov_matrix_B <- (cov_matrix_B + t(cov_matrix_B)) / 2; diag(cov_matrix_B) <- 1

# Run analysis
result <- SPDtrace(cov_matrix_A, cov_matrix_B, sparsityLevel = 100)

# Visualize results
plot(result$lambda_sequence, type = "l", main = "Lambda Sequence")
```

### Simulation Studies

```r
# Run performance comparison
results <- run_first_scenario(
  n_samples = 100, 
  n_variables = 20, 
  n_replicates = 10
)

# View performance summary
print(results$performance_summary)
```

 

## Real-World Applications

This package can be applied to public datasets (e.g., TCGA) for ovarian cancer analysis.

## Performance

- **Efficiency**: C++ implementation for fast computation
- **Scalability**: Handles networks with thousands of nodes
- **Memory**: Optimized memory usage for large datasets
- **Parallelization**: Support for parallel simulation studies

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

## Contributing

We welcome contributions! Please feel free to:

1. Report bugs and issues
2. Suggest new features
3. Submit pull requests
4. Improve documentation

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Mojtaba Nikahd
- **Email**: nikahd@ce.sharif.ir
- **GitHub**: [@mojtabanikahd](https://github.com/mojtabanikahd)
- **Paper**: [arXiv:2410.02496](https://doi.org/10.48550/arXiv.2410.02496)

## Acknowledgments

- The Cancer Genome Atlas (TCGA) for providing the ovarian cancer data
- Rcpp and RcppArmadillo developers for C++ integration
- The R community for package development tools and resources
