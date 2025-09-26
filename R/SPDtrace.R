#' SPD-trace: Efficient Learning of Differential Networks
#' 
#' @description
#' Implements the SPD-trace method for efficient learning of differential networks 
#' in multi-source non-paranormal graphical models. This function provides the main 
#' interface for differential network inference.
#' 
#' @param CovA Covariance matrix for condition A (e.g., control group)
#' @param CovB Covariance matrix for condition B (e.g., treatment group)
#' @param sparsityLevel Maximum number of edges to include in the differential network
#' @param verbose Logical indicating whether to print progress information
#' 
#' @return A list containing:
#'   \item{solution_path}{List of solution path objects with knots and active sets}
#'   \item{differential_network}{Adjacency matrix of the inferred differential network}
#'   \item{lambda_sequence}{Sequence of regularization parameters}
#'   \item{method}{Method used for inference}
#'   \item{call}{Function call}
#' 
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' p <- 20
#' 
#' # Create two covariance matrices
#' Sigma_A <- diag(p) + 0.3 * (matrix(runif(p^2), p, p) > 0.8)
#' Sigma_A <- (Sigma_A + t(Sigma_A)) / 2
#' diag(Sigma_A) <- 1
#' 
#' Sigma_B <- Sigma_A + 0.2 * (matrix(runif(p^2), p, p) > 0.9)
#' Sigma_B <- (Sigma_B + t(Sigma_B)) / 2
#' diag(Sigma_B) <- 1
#' 
#' # Run SPD-trace
#' result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 50)
#' 
#' # View results
#' print(result$differential_network)
#' plot(result$lambda_sequence, type = "l")
#' }
#' 
#' @references
#' Nikahd, M., & Motahari, S. A. (2025). Efficient Learning of Differential Networks 
#' in Multi-Source Non-Paranormal Graphical Models. arXiv preprint arXiv:2410.02496.
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
  cat("Number of edges:", sum(x$last_differential_network) / 2, "\n")
  cat("Lambda sequence length:", length(x$lambda_sequence), "\n")
  cat("\n")
  
  if (length(x$lambda_sequence) > 0) {
    cat("Lambda range:", range(x$lambda_sequence), "\n")
  }
  
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
  cat("Active edges:", sum(object$last_differential_network) / 2, "\n")
  cat("Edge density:", round(sum(object$last_differential_network) / (nrow(object$last_differential_network)^2 - nrow(object$last_differential_network)), 4), "\n")
  
  if (length(object$lambda_sequence) > 0) {
    cat("\nLambda Statistics:\n")
    cat("  Min:", min(object$lambda_sequence), "\n")
    cat("  Max:", max(object$lambda_sequence), "\n")
    cat("  Mean:", mean(object$lambda_sequence), "\n")
    cat("  Median:", median(object$lambda_sequence), "\n")
  }
  
  invisible(object)
}
