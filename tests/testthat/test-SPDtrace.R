# Test file for SPDtrace package
# This file contains unit tests for the main functions

library(testthat)
library(SPDtrace)

# Test data setup
set.seed(123)
n <- 50
p <- 10

# Create test covariance matrices
Sigma_A <- diag(p) + 0.3 * (matrix(runif(p^2), p, p) > 0.8)
Sigma_A <- (Sigma_A + t(Sigma_A)) / 2
diag(Sigma_A) <- 1

Sigma_B <- Sigma_A + 0.2 * (matrix(runif(p^2), p, p) > 0.9)
Sigma_B <- (Sigma_B + t(Sigma_B)) / 2
diag(Sigma_B) <- 1

test_that("SPDtrace function works with valid inputs", {
  # Test basic functionality
  result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 20, verbose = FALSE)
  
  # Check output structure
  expect_type(result, "list")
  expect_true("SPDtrace_result" %in% class(result))
  expect_true(all(c("solution_path", "differential_network", "lambda_sequence", "call") %in% names(result)))
  
  # Check output types
  expect_true(is.list(result$solution_path))
  expect_true(is.matrix(result$differential_network))
  expect_true(is.numeric(result$lambda_sequence))
  
  # Check matrix dimensions
  expect_equal(nrow(result$differential_network), p)
  expect_equal(ncol(result$differential_network), p)
})

test_that("SPDtrace handles different methods", {
  # Test SPDtrace method
  result_spd <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 20, 
                         verbose = FALSE)
  
  # Test CrossFDTL method (if implemented)
  # result_cross <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 20, 
  #                          method = "CrossFDTL", verbose = FALSE)
  # expect_equal(result_cross$method, "CrossFDTL")
})

test_that("SPDtrace validates input parameters", {
  # Test invalid covariance matrices
  expect_error(SPDtrace(CovA = "invalid", CovB = Sigma_B, sparsityLevel = 20))
  expect_error(SPDtrace(CovA = Sigma_A, CovB = "invalid", sparsityLevel = 20))
  
  # Test non-square matrices
  non_square <- matrix(1:20, 4, 5)
  expect_error(SPDtrace(CovA = non_square, CovB = Sigma_B, sparsityLevel = 20))
  
  # Test mismatched dimensions
  Sigma_C <- diag(p + 1)
  expect_error(SPDtrace(CovA = Sigma_A, CovB = Sigma_C, sparsityLevel = 20))
  
  # Test invalid sparsity level
  expect_error(SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = -1))
  expect_error(SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = p^2 + 1))
  
  # Test invalid method
  expect_error(SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 20))
})

test_that("Utility functions work correctly", {
  # Test kendall_tau_matrix_to_pearson_correlation_matrix
  tau_matrix <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), 3, 3)
  pearson_matrix <- kendall_tau_matrix_to_pearson_correlation_matrix(tau_matrix)
  
  expect_true(is.matrix(pearson_matrix))
  expect_equal(dim(pearson_matrix), c(3, 3))
  expect_equal(pearson_matrix[1, 1], 1)  # Diagonal should be 1
  
  # Test positive_semi_definite_maker
  A <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
  A_psd <- positive_semi_definite_maker(A)
  
  expect_true(is.matrix(A_psd))
  expect_equal(dim(A_psd), c(2, 2))
  expect_true(all(eigen(A_psd)$values >= 0))  # Should be positive semi-definite
})


test_that("Print and summary methods work", {
  result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 20, verbose = FALSE)
  
  # Test print method
  expect_output(print(result), "SPDtrace Result")
  
  # Test summary method
  expect_output(summary(result), "SPDtrace Result Summary")
})

test_that("Fisher test function works correctly", {
  # Test data
  all_genes <- paste0("Gene", 1:1000)
  hub_genes <- paste0("Gene", sample(1:1000, 50))
  functional_genes <- paste0("Gene", sample(1:1000, 100))
  
  # Run test
  result <- fisher_test(all_genes, hub_genes, functional_genes)
  
  # Check output structure
  expect_true(is.list(result))
  expect_true(all(c("contingency_table", "p.value", "odds_ratio", "conf.int") %in% names(result)))
  
  # Check types
  expect_true(is.matrix(result$contingency_table))
  expect_true(is.numeric(result$p.value))
  expect_true(is.numeric(result$odds_ratio))
  expect_true(is.numeric(result$conf.int))
  
  # Check values
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(result$odds_ratio > 0)
  expect_true(length(result$conf.int) == 2)
})

# Test edge cases
test_that("Edge cases are handled correctly", {
  # Test with very small matrices
  tiny_A <- matrix(1, 2, 2)
  tiny_B <- matrix(1, 2, 2)
  
  result <- SPDtrace(CovA = tiny_A, CovB = tiny_B, sparsityLevel = 1, verbose = FALSE)
  expect_true(is.list(result))
  
  # Test with maximum sparsity
  max_sparsity <- choose(p, 2)
  result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = max_sparsity, verbose = FALSE)
  expect_true(is.list(result))
})

# Test error handling
test_that("Error handling works correctly", {
  # Test with NULL inputs
  expect_error(SPDtrace(CovA = NULL, CovB = Sigma_B, sparsityLevel = 20))
  expect_error(SPDtrace(CovA = Sigma_A, CovB = NULL, sparsityLevel = 20))
  
  # Test with NA values
  Sigma_NA <- Sigma_A
  Sigma_NA[1, 1] <- NA
  expect_error(SPDtrace(CovA = Sigma_NA, CovB = Sigma_B, sparsityLevel = 20))
})
