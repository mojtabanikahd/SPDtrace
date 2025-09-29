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
#' Estimates Pearson correlation using the mapping of empirical Kendall's tau from 
#' multiple datasets. This function is particularly useful for validation steps 
#' in parameter tuning.
#' 
#' @param datasets List of data matrices, where each matrix represents a dataset
#' @param random_set_ratio Ratio of samples to incorporate in estimation (default: 1)
#' 
#' @return Estimated Pearson correlation matrix
#' 
#' @examples
#' \dontrun{
#' # Create example datasets
#' set.seed(123)
#' n1 <- 100; n2 <- 80; p <- 10
#' dataset1 <- matrix(rnorm(n1 * p), n1, p)
#' dataset2 <- matrix(rnorm(n2 * p), n2, p)
#' 
#' datasets <- list(dataset1, dataset2)
#' 
#' # Estimate correlation
#' cor_matrix <- weighted_rank_based_pearson_correlation_estimator(datasets)
#' print(cor_matrix[1:5, 1:5])
#' }
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

#' Evaluate Instability of Solution Paths
#' 
#' @description
#' Evaluates the instability of solution paths across multiple random subsets.
#' This function is useful for assessing the robustness of differential network inference.
#' 
#' @param solution_paths List of solution paths from multiple random subsets
#' @param number_of_nodes Number of nodes in the network
#' 
#' @return Data frame with lambda values and corresponding instability measures
#' 
#' @examples
#' \dontrun{
#' # This function is typically used with results from SPDtrace
#' # Example usage would depend on the structure of solution_paths
#' }
#' 
#' @export
instability_evaluator_of_solution_paths <- function(solution_paths, number_of_nodes) {
  if (!is.list(solution_paths) || length(solution_paths) == 0) {
    stop("solution_paths must be a non-empty list")
  }
  
  if (number_of_nodes <= 0) {
    stop("number_of_nodes must be positive")
  }
  
  number_of_random_set <- length(solution_paths)
  
  instability_results <- data.frame()
  set_index <- rep(1, number_of_random_set)
  
  while (TRUE) {
    max_index <- NA
    max_lambda_value <- 0
    edge_frequency <- rep(0, number_of_nodes^2)
    number_of_results <- 0
    
    for (i in seq_along(solution_paths)) {
      edge_existence <- rep(0, number_of_nodes^2)
      
      if (set_index[i] <= length(solution_paths[[i]])) {
        if (set_index[i] > 1) {
          edge_existence[(solution_paths[[i]][[set_index[i] - 1]]$active_set + 1)] <- 1
        }
        edge_frequency <- edge_frequency + edge_existence
        number_of_results <- number_of_results + 1
        
        if (max_lambda_value < solution_paths[[i]][[set_index[i]]]$knots_lambdas) {
          max_index <- i
          max_lambda_value <- solution_paths[[i]][[set_index[i]]]$knots_lambdas
        }
      }
    }
    
    if (number_of_results < number_of_random_set) {
      break
    }
    
    if (max_lambda_value == 0) {
      break
    }
    
    # Calculate instability
    Adj <- matrix(edge_frequency / number_of_random_set, nrow = number_of_nodes)
    edge_existence_probability <- Adj[upper.tri(Adj)]
    instability <- sum(2 * (1 - edge_existence_probability) * edge_existence_probability) / 
                  choose(number_of_nodes, 2)
    
    instability_results <- rbind(instability_results,
                                data.frame(lambda = max_lambda_value,
                                          instability = instability))
    
    set_index[max_index] <- set_index[max_index] + 1
  }
  
  instability_results
}

#' Fisher's Exact Test for Gene Set Enrichment
#' 
#' @description
#' Performs Fisher's exact test to assess enrichment of hub genes in functional gene sets.
#' This is commonly used in biological network analysis.
#' 
#' @param all_genes Vector of all genes in the analysis
#' @param hub_genes Vector of hub genes (e.g., high-degree nodes)
#' @param functional_genes Vector of functional genes (e.g., known pathway genes)
#' 
#' @return List containing contingency table and test results
#' 
#' @examples
#' \dontrun{
#' # Example gene sets
#' all_genes <- paste0("Gene", 1:1000)
#' hub_genes <- paste0("Gene", sample(1:1000, 50))
#' functional_genes <- paste0("Gene", sample(1:1000, 100))
#' 
#' result <- fisher_test(all_genes, hub_genes, functional_genes)
#' print(result$p.value)
#' }
#' 
#' @export
fisher_test <- function(all_genes, hub_genes, functional_genes) {
  if (!is.vector(all_genes) || !is.vector(hub_genes) || !is.vector(functional_genes)) {
    stop("All inputs must be vectors")
  }
  
  functional_index <- all_genes %in% functional_genes
  hub_index <- all_genes %in% hub_genes
  
  # Create contingency table
  #                    In functional set | Not in functional set
  # In hub set         a                | b
  # Not in hub set     c                | d
  
  a <- sum(hub_index & functional_index)  # Hub genes in functional set
  b <- sum(hub_index & !functional_index) # Hub genes not in functional set
  c <- sum(!hub_index & functional_index) # Non-hub genes in functional set
  d <- sum(!hub_index & !functional_index) # Non-hub genes not in functional set
  
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  colnames(contingency_table) <- c("In functional set", "Not in functional set")
  rownames(contingency_table) <- c("In hub set", "Not in hub set")
  
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table, alternative = "greater")
  
  list(
    contingency_table = contingency_table,
    p.value = test_result$p.value,
    odds_ratio = test_result$estimate,
    conf.int = test_result$conf.int
  )
}
