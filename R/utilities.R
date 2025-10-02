#' Convert Kendall's Tau Matrix to Pearson Correlation Matrix
#' 
#' @description
#' Converts a matrix of Kendall's tau correlation coefficients to Pearson correlation 
#' coefficients using the transformation: r = sin(π/2 * τ).
#' 
#' @param kendall_tau_matrix Matrix of Kendall's tau correlation coefficients
#' 
#' @return Matrix of Pearson correlation coefficients
#' 
#' @examples
#' # Create example Kendall's tau matrix
#' tau_matrix <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), 3, 3)
#' pearson_matrix <- kendall_tau_matrix_to_pearson_correlation_matrix(tau_matrix)
#' print(pearson_matrix)
#' 
#' @noRd
kendall_tau_matrix_to_pearson_correlation_matrix <- function(kendall_tau_matrix) {
  if (!is.matrix(kendall_tau_matrix)) {
    stop("Input must be a matrix")
  }
  sin((pi/2) * kendall_tau_matrix)
}

#' Make Matrix Positive Semi-Definite
#' 
#' @description
#' Ensures a symmetric matrix is positive semi-definite by adjusting diagonal elements
#' if necessary. This is useful for correlation matrices that may have negative eigenvalues.
#' 
#' @param A Symmetric matrix to be made positive semi-definite
#' 
#' @return Positive semi-definite matrix
#' 
#' @examples
#' # Create a matrix that might not be positive semi-definite
#' A <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
#' A_psd <- positive_semi_definite_maker(A)
#' eigen(A_psd)$values  # Should all be >= 0
#' 
#' @noRd
positive_semi_definite_maker <- function(A) {
  if (!is.matrix(A)) {
    stop("Input must be a matrix")
  }
  if (nrow(A) != ncol(A)) {
    stop("Input must be a square matrix")
  }
  
  # Check if matrix is symmetric
  if (!isSymmetric(A)) {
    warning("Matrix is not symmetric, making it symmetric")
    A <- (A + t(A)) / 2
  }
  
  # Get eigenvalues
  eigen_vals <- eigen(A)$values
  minimum_eigen_value <- min(eigen_vals)
  
  # Adjust if necessary
  if (minimum_eigen_value < 0) {
    A <- A - diag(minimum_eigen_value, nrow = nrow(A), ncol = ncol(A))
  }
  
  A
}

#' Weighted Rank-Based Pearson Correlation Estimator
#' 
#' @description
#' Estimates a Pearson correlation matrix from multiple datasets via a rank-based
#' procedure: (1) for each dataset (n × p), compute the Kendall's tau correlation
#' matrix; (2) integrate the tau matrices across datasets (weighted by their sample
#' sizes); (3) map the integrated tau to a Pearson correlation using
#' r = sin((π/2) * τ); and (4) if the resulting Pearson matrix is not positive
#' semi-definite, adjust it accordingly. This estimator is also practical for
#' validation and parameter-tuning workflows where a robust cross-dataset correlation
#' estimate is needed.
#' 
#' @param datasets A list of data matrices. Each matrix must have dimension n × p,
#'   where n is the number of samples and p is the number of features.
#' @param random_set_ratio Numeric in (0,1]. Proportion of samples used to compute the
#'   estimator. Defaults to 1 (all samples).
#' 
#' @return A positive semi-definite Pearson correlation matrix estimator.
#' 
#' @note Integrates information across datasets using Kendall’s tau for robustness,
#' then corrects the final Pearson correlation matrix for positive semi-definiteness
#' when necessary.
#' 
#' @examples
#' set.seed(123)
#' data1 <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
#' data2 <- matrix(rnorm(120 * 5), nrow = 120, ncol = 5)
#' data_list <- list(data1, data2)
#' 
#' # Estimate correlation matrix using all samples
#' weighted_rank_based_pearson_correlation_estimator(data_list)
#' 
#' # Estimate correlation matrix using 50% of samples
#' weighted_rank_based_pearson_correlation_estimator(data_list, random_set_ratio = 0.5)
#' 
#' @export
weighted_rank_based_pearson_correlation_estimator <- function(datasets, random_set_ratio = 1) {
  if (!is.list(datasets) || length(datasets) == 0) {
    stop("datasets must be a non-empty list")
  }
  
  if (random_set_ratio <= 0 || random_set_ratio > 1) {
    stop("random_set_ratio must be between 0 and 1")
  }
  
  number_of_all_cases <- sum(unlist(lapply(datasets, nrow)))
  number_of_nodes <- ncol(datasets[[1]])
  
  # Check all datasets have same number of columns
  for (i in seq_along(datasets)) {
    if (ncol(datasets[[i]]) != number_of_nodes) {
      stop("All datasets must have the same number of columns")
    }
  }
  
  empirical_kendall <- matrix(data = 0, nrow = number_of_nodes, ncol = number_of_nodes)
  
  for (p in seq_along(datasets)) {
    number_of_samples <- nrow(datasets[[p]])
    
    # Create random subset of samples
    indices <- sample(1:number_of_samples, 
                     size = round(random_set_ratio * number_of_samples), 
                     replace = FALSE)
    
    # Estimate empirical Kendall's tau using weighted mean
    empirical_kendall <- empirical_kendall +
      (number_of_samples / number_of_all_cases) * 
      cor(datasets[[p]][indices, ], method = "kendall")
  }
  
  # Convert to Pearson correlation
  empirical_pearson <- kendall_tau_matrix_to_pearson_correlation_matrix(empirical_kendall)
  
  # Make it positive semi-definite
  empirical_pearson <- positive_semi_definite_maker(empirical_pearson)
  
  empirical_pearson
}

