#' Run First Simulation Scenario
#' 
#' @description
#' Runs the first simulation scenario comparing SPD-trace with APGD D-trace and CrossFDTL
#' on runtime and accuracy metrics.
#' 
#' @param n_samples Number of samples (default: 100)
#' @param n_variables Number of variables (default: 20)
#' @param n_replicates Number of simulation replicates (default: 10)
#' @param sparsity_levels Vector of sparsity levels to test (default: c(10, 20, 30))
#' @param verbose Logical indicating whether to print progress information
#' 
#' @return List containing simulation results and performance metrics
#' 
#' @examples
#' \dontrun{
#' # Run first simulation scenario
#' results <- run_first_scenario(n_samples = 100, n_variables = 20, n_replicates = 5)
#' 
#' # View results
#' print(results$performance_summary)
#' plot(results$precision_recall_plot)
#' }
#' 
#' @export
run_first_scenario <- function(n_samples = 100, n_variables = 20, n_replicates = 10, 
                              sparsity_levels = c(10, 20, 30), verbose = TRUE) {
  
  if (verbose) {
    cat("Running First Simulation Scenario\n")
    cat("================================\n")
    cat("Samples:", n_samples, "\n")
    cat("Variables:", n_variables, "\n")
    cat("Replicates:", n_replicates, "\n")
    cat("Sparsity levels:", paste(sparsity_levels, collapse = ", "), "\n\n")
  }
  
  # Initialize results storage
  results <- list()
  performance_data <- data.frame()
  
  for (sparsity in sparsity_levels) {
    if (verbose) cat("Testing sparsity level:", sparsity, "\n")
    
    for (rep in 1:n_replicates) {
      if (verbose && rep %% 5 == 0) cat("  Replicate:", rep, "/", n_replicates, "\n")
      
      # Generate synthetic data
      data_result <- generate_synthetic_data(n_samples, n_variables, sparsity)
      
      # Run SPD-trace
      start_time <- Sys.time()
      spdtrace_result <- tryCatch({
        SPDtrace(data_result$cov_A, data_result$cov_B, sparsity, method = "SPDtrace", verbose = FALSE)
      }, error = function(e) {
        list(differential_network = matrix(0, n_variables, n_variables), 
             lambda_sequence = numeric(0))
      })
      spdtrace_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Run CrossFDTL
      start_time <- Sys.time()
      crossfdtl_result <- tryCatch({
        SPDtrace(data_result$cov_A, data_result$cov_B, sparsity, method = "CrossFDTL", verbose = FALSE)
      }, error = function(e) {
        list(differential_network = matrix(0, n_variables, n_variables), 
             lambda_sequence = numeric(0))
      })
      crossfdtl_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Calculate performance metrics
      spdtrace_metrics <- calculate_performance_metrics(
        spdtrace_result$differential_network, 
        data_result$true_network
      )
      
      crossfdtl_metrics <- calculate_performance_metrics(
        crossfdtl_result$differential_network, 
        data_result$true_network
      )
      
      # Store results
      performance_data <- rbind(performance_data,
        data.frame(
          method = "SPDtrace",
          sparsity = sparsity,
          replicate = rep,
          precision = spdtrace_metrics$precision,
          recall = spdtrace_metrics$recall,
          f1_score = spdtrace_metrics$f1_score,
          runtime = spdtrace_time
        ),
        data.frame(
          method = "CrossFDTL",
          sparsity = sparsity,
          replicate = rep,
          precision = crossfdtl_metrics$precision,
          recall = crossfdtl_metrics$recall,
          f1_score = crossfdtl_metrics$f1_score,
          runtime = crossfdtl_time
        )
      )
    }
  }
  
  # Calculate summary statistics
  performance_summary <- performance_data %>%
    group_by(method, sparsity) %>%
    summarize(
      mean_precision = mean(precision, na.rm = TRUE),
      mean_recall = mean(recall, na.rm = TRUE),
      mean_f1 = mean(f1_score, na.rm = TRUE),
      mean_runtime = mean(runtime, na.rm = TRUE),
      sd_precision = sd(precision, na.rm = TRUE),
      sd_recall = sd(recall, na.rm = TRUE),
      sd_f1 = sd(f1_score, na.rm = TRUE),
      sd_runtime = sd(runtime, na.rm = TRUE),
      .groups = "drop"
    )
  
  results$performance_data <- performance_data
  results$performance_summary <- performance_summary
  results$parameters <- list(
    n_samples = n_samples,
    n_variables = n_variables,
    n_replicates = n_replicates,
    sparsity_levels = sparsity_levels
  )
  
  if (verbose) {
    cat("\nSimulation completed!\n")
    cat("Results summary:\n")
    print(performance_summary)
  }
  
  return(results)
}

#' Run Second Simulation Scenario
#' 
#' @description
#' Runs the second simulation scenario evaluating robustness of SPD-trace using 
#' heterogeneous vs. homogeneous datasets.
#' 
#' @param n_samples Number of samples (default: 100)
#' @param n_variables Number of variables (default: 20)
#' @param n_replicates Number of simulation replicates (default: 10)
#' @param heterogeneity_levels Vector of heterogeneity levels to test (default: c(0, 0.2, 0.4))
#' @param verbose Logical indicating whether to print progress information
#' 
#' @return List containing simulation results and robustness metrics
#' 
#' @examples
#' \dontrun{
#' # Run second simulation scenario
#' results <- run_second_scenario(n_samples = 100, n_variables = 20, n_replicates = 5)
#' 
#' # View results
#' print(results$robustness_summary)
#' }
#' 
#' @export
run_second_scenario <- function(n_samples = 100, n_variables = 20, n_replicates = 10,
                               heterogeneity_levels = c(0, 0.2, 0.4), verbose = TRUE) {
  
  if (verbose) {
    cat("Running Second Simulation Scenario\n")
    cat("=================================\n")
    cat("Samples:", n_samples, "\n")
    cat("Variables:", n_variables, "\n")
    cat("Replicates:", n_replicates, "\n")
    cat("Heterogeneity levels:", paste(heterogeneity_levels, collapse = ", "), "\n\n")
  }
  
  # Initialize results storage
  results <- list()
  robustness_data <- data.frame()
  
  for (heterogeneity in heterogeneity_levels) {
    if (verbose) cat("Testing heterogeneity level:", heterogeneity, "\n")
    
    for (rep in 1:n_replicates) {
      if (verbose && rep %% 5 == 0) cat("  Replicate:", rep, "/", n_replicates, "\n")
      
      # Generate heterogeneous data
      data_result <- generate_heterogeneous_data(n_samples, n_variables, heterogeneity)
      
      # Run SPD-trace on homogeneous subset
      start_time <- Sys.time()
      homogeneous_result <- tryCatch({
        SPDtrace(data_result$cov_homogeneous, data_result$cov_B, 
                sparsity = 20, method = "SPDtrace", verbose = FALSE)
      }, error = function(e) {
        list(differential_network = matrix(0, n_variables, n_variables), 
             lambda_sequence = numeric(0))
      })
      homogeneous_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Run SPD-trace on heterogeneous data
      start_time <- Sys.time()
      heterogeneous_result <- tryCatch({
        SPDtrace(data_result$cov_heterogeneous, data_result$cov_B, 
                sparsity = 20, method = "SPDtrace", verbose = FALSE)
      }, error = function(e) {
        list(differential_network = matrix(0, n_variables, n_variables), 
             lambda_sequence = numeric(0))
      })
      heterogeneous_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Calculate performance metrics
      homogeneous_metrics <- calculate_performance_metrics(
        homogeneous_result$differential_network, 
        data_result$true_network
      )
      
      heterogeneous_metrics <- calculate_performance_metrics(
        heterogeneous_result$differential_network, 
        data_result$true_network
      )
      
      # Store results
      robustness_data <- rbind(robustness_data,
        data.frame(
          data_type = "Homogeneous",
          heterogeneity = heterogeneity,
          replicate = rep,
          precision = homogeneous_metrics$precision,
          recall = homogeneous_metrics$recall,
          f1_score = homogeneous_metrics$f1_score,
          runtime = homogeneous_time
        ),
        data.frame(
          data_type = "Heterogeneous",
          heterogeneity = heterogeneity,
          replicate = rep,
          precision = heterogeneous_metrics$precision,
          recall = heterogeneous_metrics$recall,
          f1_score = heterogeneous_metrics$f1_score,
          runtime = heterogeneous_time
        )
      )
    }
  }
  
  # Calculate summary statistics
  robustness_summary <- robustness_data %>%
    group_by(data_type, heterogeneity) %>%
    summarize(
      mean_precision = mean(precision, na.rm = TRUE),
      mean_recall = mean(recall, na.rm = TRUE),
      mean_f1 = mean(f1_score, na.rm = TRUE),
      mean_runtime = mean(runtime, na.rm = TRUE),
      .groups = "drop"
    )
  
  results$robustness_data <- robustness_data
  results$robustness_summary <- robustness_summary
  results$parameters <- list(
    n_samples = n_samples,
    n_variables = n_variables,
    n_replicates = n_replicates,
    heterogeneity_levels = heterogeneity_levels
  )
  
  if (verbose) {
    cat("\nSimulation completed!\n")
    cat("Results summary:\n")
    print(robustness_summary)
  }
  
  return(results)
}

#' Generate Synthetic Data for Simulation
#' 
#' @description
#' Generates synthetic covariance matrices and true differential network for simulation studies.
#' 
#' @param n_samples Number of samples
#' @param n_variables Number of variables
#' @param sparsity Number of edges in the true differential network
#' 
#' @return List containing covariance matrices and true network
#' 
#' @keywords internal
generate_synthetic_data <- function(n_samples, n_variables, sparsity) {
  set.seed(123)
  
  # Generate base covariance matrix
  base_cov <- diag(n_variables) + 0.3 * (matrix(runif(n_variables^2), n_variables, n_variables) > 0.8)
  base_cov <- (base_cov + t(base_cov)) / 2
  diag(base_cov) <- 1
  
  # Generate differential network
  true_network <- matrix(0, n_variables, n_variables)
  edge_positions <- sample(which(upper.tri(true_network)), sparsity)
  true_network[edge_positions] <- 1
  true_network <- true_network + t(true_network)
  
  # Create second covariance matrix with differences
  cov_B <- base_cov + 0.2 * true_network
  cov_B <- (cov_B + t(cov_B)) / 2
  diag(cov_B) <- 1
  
  # Ensure positive definiteness
  cov_A <- positive_semi_definite_maker(base_cov)
  cov_B <- positive_semi_definite_maker(cov_B)
  
  list(
    cov_A = cov_A,
    cov_B = cov_B,
    true_network = true_network
  )
}

#' Generate Heterogeneous Data for Simulation
#' 
#' @description
#' Generates heterogeneous datasets for robustness testing.
#' 
#' @param n_samples Number of samples
#' @param n_variables Number of variables
#' @param heterogeneity_level Level of heterogeneity (0 = homogeneous, 1 = highly heterogeneous)
#' 
#' @return List containing homogeneous and heterogeneous covariance matrices
#' 
#' @keywords internal
generate_heterogeneous_data <- function(n_samples, n_variables, heterogeneity_level) {
  set.seed(123)
  
  # Base covariance matrix
  base_cov <- diag(n_variables) + 0.3 * (matrix(runif(n_variables^2), n_variables, n_variables) > 0.8)
  base_cov <- (base_cov + t(base_cov)) / 2
  diag(base_cov) <- 1
  
  # Create heterogeneous version
  noise_matrix <- matrix(rnorm(n_variables^2) * heterogeneity_level, n_variables, n_variables)
  noise_matrix <- (noise_matrix + t(noise_matrix)) / 2
  diag(noise_matrix) <- 0
  
  cov_heterogeneous <- base_cov + noise_matrix
  cov_heterogeneous <- positive_semi_definite_maker(cov_heterogeneous)
  
  # Create condition B covariance
  cov_B <- base_cov + 0.2 * (matrix(runif(n_variables^2), n_variables, n_variables) > 0.9)
  cov_B <- (cov_B + t(cov_B)) / 2
  diag(cov_B) <- 1
  cov_B <- positive_semi_definite_maker(cov_B)
  
  list(
    cov_homogeneous = base_cov,
    cov_heterogeneous = cov_heterogeneous,
    cov_B = cov_B,
    true_network = matrix(0, n_variables, n_variables)  # Placeholder
  )
}

#' Calculate Performance Metrics
#' 
#' @description
#' Calculates precision, recall, and F1-score for network inference results.
#' 
#' @param inferred_network Inferred differential network adjacency matrix
#' @param true_network True differential network adjacency matrix
#' 
#' @return List containing precision, recall, and F1-score
#' 
#' @keywords internal
calculate_performance_metrics <- function(inferred_network, true_network) {
  # Convert to binary (0/1)
  inferred_binary <- (inferred_network != 0) * 1
  true_binary <- (true_network != 0) * 1
  
  # Calculate metrics
  tp <- sum(inferred_binary & true_binary)
  fp <- sum(inferred_binary & !true_binary)
  fn <- sum(!inferred_binary & true_binary)
  
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1_score <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  
  list(
    precision = precision,
    recall = recall,
    f1_score = f1_score
  )
}
