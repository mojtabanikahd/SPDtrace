#' Ovarian Cancer Gene Expression Data
#' 
#' @description
#' Gene expression data from The Cancer Genome Atlas (TCGA) for ovarian cancer samples.
#' This dataset contains expression profiles for platinum-sensitive and platinum-resistant tumors.
#' 
#' @format A data frame with 489 samples and 20,531 genes
#' 
#' @source The Cancer Genome Atlas (TCGA) - Ovarian Cancer Dataset
#' 
#' @examples
#' \dontrun{
#' # Load the data
#' data("ov_gene_exp")
#' 
#' # View dimensions
#' dim(ov_gene_exp)
#' 
#' # View first few rows and columns
#' ov_gene_exp[1:5, 1:5]
#' }
#' 
#' @references
#' The Cancer Genome Atlas Research Network. (2011). Integrated genomic analyses of ovarian carcinoma. 
#' Nature, 474(7353), 609-615.
#' 
#' @keywords datasets
"ov_gene_exp"

#' PI3K-Akt-mTOR Pathway Data
#' 
#' @description
#' Data related to the PI3K-Akt-mTOR signaling pathway in ovarian cancer.
#' This dataset contains pathway information and gene sets for biological validation.
#' 
#' @format A list containing pathway information and gene sets
#' 
#' @source Pathway databases and literature curation
#' 
#' @examples
#' \dontrun{
#' # Load the data
#' data("ov_pi3k_akt_mtor")
#' 
#' # View structure
#' str(ov_pi3k_akt_mtor)
#' }
#' 
#' @keywords datasets
"ov_pi3k_akt_mtor"

#' Weighted D-trace Matrix
#' 
#' @description
#' Weighted D-trace matrix for ovarian cancer differential network analysis.
#' This matrix is used in the StARS model selection procedure.
#' 
#' @format A matrix with dimensions matching the gene expression data
#' 
#' @source Computed from gene expression data using D-trace method
#' 
#' @examples
#' \dontrun{
#' # Load the data
#' data("ov_wdtrace_mat")
#' 
#' # View dimensions
#' dim(ov_wdtrace_mat)
#' 
#' # View summary statistics
#' summary(as.vector(ov_wdtrace_mat))
#' }
#' 
#' @keywords datasets
"ov_wdtrace_mat"

#' Platinum Resistant Genes
#' 
#' @description
#' List of genes known to be associated with platinum resistance in ovarian cancer.
#' This dataset is used for biological validation of inferred differential networks.
#' 
#' @format A data frame with gene information and resistance associations
#' 
#' @source Literature curation and pathway databases
#' 
#' @examples
#' \dontrun{
#' # Load the data
#' data("platinum_resistant_genes")
#' 
#' # View structure
#' str(platinum_resistant_genes)
#' 
#' # Count genes
#' nrow(platinum_resistant_genes)
#' }
#' 
#' @keywords datasets
"platinum_resistant_genes"
