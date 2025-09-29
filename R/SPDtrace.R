#' @noRd
# Private helper function to print edges from upper triangular matrix indices
convert_indices_to_edge <- function(upper_tri_indices, p, node_labels) {
  if (length(upper_tri_indices) == 0) {
    cat("No active edges found.\n")
    return(invisible(NULL))
  }
  
  cat("Active edges:\n")
  
  # Convert upper triangular indices to row/column pairs and print
  for (idx in upper_tri_indices) {
    # Convert upper triangular index to matrix coordinates
    # For upper triangular matrix, index i corresponds to row/col in upper triangle
    row <- ((idx-1) %/% p) + 1
    col <- ((idx-1) %% p) + 1
    
    # Ensure we're within bounds
    if (row < col && col <= p) {
      node1 <- if (!is.null(node_labels)) node_labels[row] else paste0("Node", row)
      node2 <- if (!is.null(node_labels)) node_labels[col] else paste0("Node", col)
      cat("  ", node1, "---", node2, "\n")
    }
  }
}
#' Sparse Differential Network Estimation via SPDtrace
#'
#' @description
#' SPDtrace is an efficient method for learning sparse structural changes (differential
#' network) between two classes of non-paranormal graphical models, by solving the
#' lasso-penalized D-trace loss function. It computes the exact solution path across
#' regularization parameters by exploiting the piecewise linearity of the solution path,
#' substantially reducing computational complexity in high-dimensional settings.
#'
#' @param CovA Estimated covariance matrix for class A.
#' @param CovB Estimated covariance matrix for class B.
#' @param sparsityLevel Integer, maximum number of edges to be identified in the
#'   differential network. Controls algorithm termination. Larger values increase runtime.
#' @param verbose Logical, whether to print progress logs.
#'
#' @return A list with the following components:
#'   \item{solution_path}{Object containing iteration logs, number of steps, intermediate
#'   results, and corresponding \eqn{\lambda} values.}
#'   \item{last_differential_network}{Final binary adjacency matrix (1 = differential edge present,
#'   0 = absent).}
#'   \item{lambda_sequence}{Vector of \eqn{\lambda} values (knots) where new edges are identified.}
#'   \item{call}{Function call, showing how the method was invoked.}
#'
#' @examples
#' \dontrun{
#' # Example with covariance matrices for two classes
#' set.seed(123)
#' p <- 20
#'
#' # Construct covariance for class A
#' Sigma_A <- diag(p) + 0.3 * (matrix(runif(p^2), p, p) > 0.8)
#' Sigma_A <- (Sigma_A + t(Sigma_A)) / 2; diag(Sigma_A) <- 1
#'
#' # Construct covariance for class B
#' Sigma_B <- Sigma_A + 0.2 * (matrix(runif(p^2), p, p) > 0.9)
#' Sigma_B <- (Sigma_B + t(Sigma_B)) / 2; diag(Sigma_B) <- 1
#'
#' # Run SPDtrace
#' result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 50, verbose = TRUE)
#'
#' # Inspect results
#' print(result$last_differential_network)
#' plot(result$lambda_sequence, type = "l", xlab = "Step", ylab = expression(lambda))
#'
#' # Alternative example using precision matrices for two classes
#' d <- 10
#' m <- 9
#'
#' # build the precision matrix of class A
#' Precision_A <- outer(1:d, 1:d, function(i, j) 0.5^(abs(i - j)))
#'
#' # build the precision matrix of class B
#' Precision_B <- Precision_A
#' Precision_A[abs(row(Precision_A) - col(Precision_B)) == m] <- 0.5
#'
#' # Compute the covariance matrix of both classes
#' CovA <- solve(Precision_A)
#' CovB <- solve(Precision_B)
#'
#' # Compute the differential network
#' SPDtrace(CovA = CovA, CovB = CovB, sparsityLevel = 1)
#' }
#'
#' @references
#' Nikahd, M., & Motahari, S. A. (2024). Efficient Learning of Differential Networks
#' in Multi-Source Non-Paranormal Graphical Models. Please cite this manuscript when
#' using the SPDtrace method.
#'
#' @export
SPDtrace <- function(CovA, CovB, sparsityLevel = NULL, verbose = TRUE) {
  
  # Input validation
  if (!is.matrix(CovA) || !is.matrix(CovB)) {
    stop("CovA and CovB must be matrices")
  }
  
  if (nrow(CovA) != ncol(CovA) || nrow(CovB) != ncol(CovB)) {
    stop("CovA and CovB must be square matrices")
  }
  
  if (nrow(CovA) != nrow(CovB)) {
    stop("CovA and CovB must have the same dimensions")
  }
  
  # Check for NA values
  if (any(is.na(CovA)) || any(is.na(CovB))) {
    stop("CovA and CovB cannot contain NA values")
  }
  
  if (is.null(sparsityLevel)) {
    sparsityLevel <- choose(nrow(CovA), 2) %/% 10  # Default to 10% of possible edges
  }
  
  if (sparsityLevel <= 0 || sparsityLevel > choose(nrow(CovA), 2)) {
    stop("sparsityLevel must be between 1 and the number of possible edges")
  }
  
  if (verbose) {
    cat("Matrix dimensions:", nrow(CovA), "x", ncol(CovA), "\n")
    cat("Sparsity level:", sparsityLevel, "\n")
  }
  
  sparsityLevel = 2*sparsityLevel
  result <- inference_Dtrace_solution_path(CovA, CovB, sparsityLevel)
  
  # Process results
  solution_path <- list()
  lambda_sequence <- numeric(0)
  sp = result$solution_path

  if (length(sp) > 0) {
    for (i in seq_along(sp)) {
      if (is.list(sp[[i]]) && "knots_lambdas" %in% names(sp[[i]])) {
        solution_path[[i]] <- sp[[i]]
        lambda_sequence[i] <- sp[[i]]$knots_lambdas
      }
    }
  }
  
  # Create differential network adjacency matrix
  p <- nrow(CovA)
  differential_network <- matrix(0, p, p)
  
  if (length(solution_path) > 0) {
    # Use the last solution in the path
    last_solution <- solution_path[[length(solution_path)]]
    if ("active_set" %in% names(last_solution)) {
      active_edges <- last_solution$active_set  # Convert to 1-based indexing
      for (idx in seq(active_edges)) {
        edge_idx <- active_edges[idx]
        if (edge_idx <= p^2) {
          row_idx <- as.integer((edge_idx) %% p) + 1
          col_idx <- as.integer((edge_idx) / p) + 1
          if (row_idx != col_idx) {
            differential_network[row_idx, col_idx] <- 1
            differential_network[col_idx, row_idx] <- 1
          }
        }
      }
    }
  }
  rownames(differential_network) <- rownames(CovA)
  colnames(differential_network) <- colnames(CovA)
  # Create output object
  output <- list(
    solution_path = solution_path,
    last_differential_network = differential_network,
    lambda_sequence = lambda_sequence,
    call = match.call()
  )
  
  class(output) <- "SPDtrace_result"
  
  if (verbose) {
    cat("Inference completed successfully!\n")
  }
  
  return(output)
}

#' Print method for SPDtrace results
#' 
#' @param x Object of class SPDtrace_result
#' @param ... Additional arguments passed to print
#' 
#' @export
print.SPDtrace_result <- function(x, ...) {
  cat("SPDtrace Result\n")
  cat("==============\n")
  cat("Method:", x$method, "\n")
  cat("Matrix dimensions:", nrow(x$last_differential_network), "x", ncol(x$last_differential_network), "\n")
  cat("Number of edges:", sum(x$last_differential_network[upper.tri(x$last_differential_network)]), "\n")
  cat("Lambda sequence length:", length(x$lambda_sequence), "\n")
  # Extract active edges from differential network matrix
  upper_tri_indices <- which(x$last_differential_network == 1)
  convert_indices_to_edge(upper_tri_indices, ncol(x$last_differential_network), colnames(x$last_differential_network))
  cat("\n")
  
  invisible(x)
}

#' Summary method for SPDtrace results
#' 
#' @param object Object of class SPDtrace_result
#' @param ... Additional arguments passed to summary
#' 
#' @export
summary.SPDtrace_result <- function(object, ...) {
  cat("SPDtrace Result Summary\n")
  cat("======================\n")
  cat("Method:", object$method, "\n")
  cat("Matrix dimensions:", nrow(object$last_differential_network), "x", ncol(object$last_differential_network), "\n")
  cat("Total possible edges:", choose(nrow(object$last_differential_network), 2), "\n")
  cat("Active edges:", sum(object$last_differential_network[upper.tri(object$last_differential_network)]), "\n")
  cat("Lambda sequence length:", length(object$lambda_sequence), "\n")
  # Extract active edges from differential network matrix
  upper_tri_indices <- which(object$last_differential_network[upper.tri(object$last_differential_network)] == 1)
  convert_indices_to_edge(upper_tri_indices, ncol(object$last_differential_network), colnames(object$last_differential_network))
  cat("\n")

  invisible(object)
}