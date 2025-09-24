\VignetteIndexEntry{Introduction to SPDtrace}
\VignetteEngine{utils::Sweave}
\VignetteEncoding{UTF-8}

\title{Introduction to SPDtrace: Efficient Learning of Differential Networks}
\author{Mojtaba Nikahd, Seyed Abolfazl Motahari}
\date{\today}

\section{Introduction}

The \textbf{SPDtrace} package implements the SPD-trace method for efficient learning of differential networks in multi-source non-paranormal graphical models. This package provides tools for differential network inference, simulation studies, and real-world applications in genomics and other high-dimensional data analysis.

\section{Installation}

\begin{Schunk}
\begin{Sinput}
# Install from CRAN (when available)
install.packages("SPDtrace")

# Or install from GitHub
# devtools::install_github("yourusername/SPDtrace")
\end{Sinput}
\end{Schunk}

\section{Quick Start}

\subsection{Basic Usage}

The main function \texttt{SPDtrace()} provides a simple interface for differential network inference:

\begin{Schunk}
\begin{Sinput}
library(SPDtrace)

# Generate example data
set.seed(123)
n <- 100
p <- 20

# Create two covariance matrices
Sigma_A <- diag(p) + 0.3 * (matrix(runif(p^2), p, p) > 0.8)
Sigma_A <- (Sigma_A + t(Sigma_A)) / 2
diag(Sigma_A) <- 1

Sigma_B <- Sigma_A + 0.2 * (matrix(runif(p^2), p, p) > 0.9)
Sigma_B <- (Sigma_B + t(Sigma_B)) / 2
diag(Sigma_B) <- 1

# Run SPD-trace
result <- SPDtrace(CovA = Sigma_A, CovB = Sigma_B, sparsityLevel = 50)

# View results
print(result)
summary(result)
\end{Sinput}
\end{Schunk}

\subsection{Understanding the Output}

The \texttt{SPDtrace()} function returns an object of class \texttt{SPDtrace\_result} containing:

\begin{itemize}
\item \texttt{solution\_path}: List of solution path objects with knots and active sets
\item \texttt{differential\_network}: Adjacency matrix of the inferred differential network
\item \texttt{lambda\_sequence}: Sequence of regularization parameters
\item \texttt{method}: Method used for inference
\item \texttt{call}: Function call
\end{itemize}

\section{Core Functions}

\subsection{Main Inference Function}

\begin{Schunk}
\begin{Sinput}
# Basic usage
result <- SPDtrace(CovA, CovB, sparsityLevel = 30)

# With different methods
result_spd <- SPDtrace(CovA, CovB, sparsityLevel = 30, method = "SPDtrace")
result_cross <- SPDtrace(CovA, CovB, sparsityLevel = 30, method = "CrossFDTL")
result_dtrace <- SPDtrace(CovA, CovB, sparsityLevel = 30, method = "DTrace")
\end{Sinput}
\end{Schunk}

\subsection{Utility Functions}

The package provides several utility functions for data preprocessing and analysis:

\begin{Schunk}
\begin{Sinput}
# Convert Kendall's tau to Pearson correlation
tau_matrix <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), 3, 3)
pearson_matrix <- kendall_tau_matrix_to_pearson_correlation_matrix(tau_matrix)

# Ensure positive semi-definiteness
A_psd <- positive_semi_definite_maker(A)
\end{Sinput}
\end{Schunk}

\section{Advanced Usage}

\subsection{CrossFDTL Method}

For more advanced users, the package provides direct access to the CrossFDTL algorithm:

\begin{Schunk}
\begin{Sinput}
# Direct CrossFDTL usage
result <- CrossFDTL(CovA, CovB, lambda = 0.1, rho = 0.5, maxiter = 100)
\end{Sinput}
\end{Schunk}

\subsection{Solution Path Inference}

You can also use the DTrace solution path method:

\begin{Schunk}
\begin{Sinput}
# Solution path inference
result <- inference_Dtrace_solution_path(CovA, CovB, sparsityLevel = 50)
\end{Sinput}
\end{Schunk}

\section{Examples}

\subsection{Simulation Study}

Simulation helpers have been removed in this refactor to reduce dependencies.

\subsection{Performance Metrics}

Calculate performance metrics for your results:

\begin{Schunk}
\begin{Sinput}
# Calculate performance metrics
metrics <- calculate_performance_metrics(predicted_network, true_network)
print(metrics)
\end{Sinput}
\end{Schunk}

\section{Data Import}

Provide your own datasets and convert them to covariance matrices before calling \texttt{SPDtrace()}. For example, you can compute sample covariance from an expression matrix:

\begin{Schunk}
\begin{Sinput}
# Given two expression matrices X_A and X_B (samples x genes)
cov_A <- cov(X_A)
cov_B <- cov(X_B)
result <- SPDtrace(CovA = cov_A, CovB = cov_B, sparsityLevel = 50)
\end{Sinput}
\end{Schunk}

\section{Visualization}

Create publication-ready plots:

\begin{Schunk}
\begin{Sinput}
# Create network visualization
library(igraph)
g <- graph_from_adjacency_matrix(result$differential_network)
plot(g, vertex.size = 5, vertex.label = NA)

# Create performance comparison plots
library(ggplot2)
ggplot(performance_data, aes(x = method, y = accuracy)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Method Performance Comparison")
\end{Sinput}
\end{Schunk}

\section{Getting Help}

For more information and examples:

\begin{itemize}
\item Type \texttt{?SPDtrace} for the main function help
\item Type \texttt{?CrossFDTL} for the C++ function help
\item Check the package documentation: \texttt{help(package = "SPDtrace")}
\end{itemize}

