# Load necessary libraries
library(MatrixEQTL)
library(dplyr)
library(ggplot2) # Optional, for plotting

# --- 1. Define Parameters for the Analysis ---

# Set a random seed for reproducibility of synthetic data
set.seed(123)

# Number of individuals (samples)
num_individuals <- 100

# Number of genes (expression features)
num_genes <- 1000

# Number of SNPs (genetic variants)
num_snps <- 5000

# Number of covariates (e.g., age, sex, principal components)
num_covariates <- 3

# Output file paths for results (will be created in your working directory)
output_file_cis <- "MatrixEQTL_cis_results.txt"
output_file_trans <- "MatrixEQTL_trans_results.txt"

# Only analyze associations for cis-acting SNPs (SNPs near the gene)
# Set to TRUE for cis-eQTLs, FALSE for trans-eQTLs (or both if combined)
# For this example, we'll run both separately.
use_only_cis <- TRUE # For cis-eQTL analysis
use_only_trans <- FALSE # For trans-eQTL analysis

# Distance for cis-eQTL analysis (SNPs within this distance from gene start)
# This is usually in base pairs.
cis_distance <- 1e6 # 1 MB (1,000,000 base pairs)

# Significance threshold for reported associations
pv_threshold_cis <- 1e-6 # P-value threshold for cis-eQTLs
pv_threshold_trans <- 1e-8 # P-value threshold for trans-eQTLs (stricter due to multiple testing)

# Error covariance matrix (optional, for paired data or complex designs)
# Set to numeric() for standard independent samples (most common)
error_cov_matrix <- numeric()

# --- 2. Generate Synthetic Data (Mimicking Real Data Structure) ---

cat("Generating synthetic data...\n")

# A. Gene Expression Data
# Rows are genes, columns are individuals.
# Values typically normalized counts (e.g., RPKM, FPKM, TPM, or VST/rlog from DESeq2/edgeR).
gene_expression <- matrix(
  rnorm(num_genes * num_individuals, mean = 5, sd = 2),
  nrow = num_genes,
  ncol = num_individuals
)
rownames(gene_expression) <- paste0("Gene_", 1:num_genes)
colnames(gene_expression) <- paste0("Individual_", 1:num_individuals)

# B. Genotype Data
# Rows are SNPs, columns are individuals.
# Values are typically 0, 1, or 2 (number of minor alleles).
# For simplicity, we'll generate random 0, 1, 2 genotypes.
genotype_data <- matrix(
  sample(0:2, num_snps * num_individuals, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
  nrow = num_snps,
  ncol = num_individuals
)
rownames(genotype_data) <- paste0("SNP_", 1:num_snps)
colnames(genotype_data) <- paste0("Individual_", 1:num_individuals)

# C. Covariate Data
# Rows are covariates, columns are individuals.
# These could be age, sex, batch effects, population stratification (e.g., first few PCs).
covariate_data <- matrix(
  rnorm(num_covariates * num_individuals),
  nrow = num_covariates,
  ncol = num_individuals
)
rownames(covariate_data) <- paste0("Covariate_", 1:num_covariates)
colnames(covariate_data) <- paste0("Individual_", 1:num_individuals)

# D. Gene Location Data (Required for cis-eQTLs)
# Data frame with gene ID, chromosome, and start position.
gene_locations <- data.frame(
  geneid = rownames(gene_expression),
  chr = sample(1:22, num_genes, replace = TRUE), # Assign random chromosomes
  s1 = sample(1:1e8, num_genes, replace = TRUE)  # Random start positions
)

# E. SNP Location Data (Required for cis-eQTLs)
# Data frame with SNP ID, chromosome, and position.
snp_locations <- data.frame(
  snpid = rownames(genotype_data),
  chr = sample(1:22, num_snps, replace = TRUE), # Assign random chromosomes
  pos = sample(1:1e8, num_snps, replace = TRUE) # Random positions
)

cat("Synthetic data generation complete.\n\n")

# --- 3. Prepare Data for MatrixEQTL ---
# MatrixEQTL uses a specific data object called 'SlicedData'.
# It can handle data stored in files or directly in R matrices.
# For this example, we'll use R matrices.

# Create SlicedData objects
snps <- SlicedData$new(genotype_data)
genes <- SlicedData$new(gene_expression)
cvrt <- SlicedData$new(covariate_data) # Can be NULL if no covariates

# Add gene and SNP locations for cis-eQTL analysis
# These data frames must have 'snpid', 'chr', 'pos' for SNPs and
# 'geneid', 'chr', 's1' (start position) for genes.
snps$CreateReplicationMatrix(snp_locations)
genes$CreateReplicationMatrix(gene_locations)

cat("SlicedData objects created.\n\n")

# --- 4. Run the eQTL Analysis with MatrixEQTL ---

# Model type: linear regression. Other options are 'linear' (default), 'ANOVA'.
# For standard eQTLs, 'linear' is appropriate.

# Run cis-eQTL analysis
cat("Running cis-eQTL analysis...\n")
me_cis <- Matrix_eQTL_main(
  snps = snps,
  gene = genes,
  cvrt = cvrt,
  output_file_name = output_file_cis,
  pvOutputThreshold = pv_threshold_cis,
  useModel = modelLINEAR, # Linear regression model
  errorCovariance = error_cov_matrix,
  verbose = TRUE,
  pvalue.hist = TRUE, # Generate p-value histogram data
  min.pv.by.genesnp = TRUE, # Store only the best SNP per gene, or best gene per SNP
  # Cis-eQTL specific parameters
  cisDist = cis_distance,
  cisOnly = TRUE # Crucial for cis-eQTL
)
cat("Cis-eQTL analysis complete. Results saved to:", output_file_cis, "\n\n")

# Run trans-eQTL analysis (reset cisOnly to FALSE)
# Note: For large datasets, trans-eQTLs can be computationally intensive.
# A more stringent p-value threshold is usually required for trans-eQTLs
# due to the massive number of tests.
cat("Running trans-eQTL analysis...\n")
me_trans <- Matrix_eQTL_main(
  snps = snps,
  gene = genes,
  cvrt = cvrt,
  output_file_name = output_file_trans,
  pvOutputThreshold = pv_threshold_trans,
  useModel = modelLINEAR,
  errorCovariance = error_cov_matrix,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  # Trans-eQTL specific parameter
  cisOnly = FALSE # Crucial for trans-eQTL
)
cat("Trans-eQTL analysis complete. Results saved to:", output_file_trans, "\n\n")


# --- 5. Inspect and Interpret Results ---

cat("Summary of cis-eQTL results:\n")
print(me_cis$cis$eqtls) # Access the results dataframe
cat("\n")

cat("Summary of trans-eQTL results:\n")
print(me_trans$trans$eqtls) # Access the results dataframe
cat("\n")

# Count significant associations
num_cis_eQTLs <- nrow(me_cis$cis$eqtls)
num_trans_eQTLs <- nrow(me_trans$trans$eqtls)

cat(paste("Number of significant cis-eQTLs found (p <", pv_threshold_cis, "):", num_cis_eQTLs, "\n"))
cat(paste("Number of significant trans-eQTLs found (p <", pv_threshold_trans, "):", num_trans_eQTLs, "\n"))


# Plot p-value histogram for cis-eQTLs (optional)
# This helps check for uniform p-value distribution under the null hypothesis
# and identify batch effects or other issues.
if (num_cis_eQTLs > 0) {
  hist(
    me_cis$cis$pvalue.hist$all,
    xlab = "P-value",
    main = "Cis-eQTL P-value Distribution",
    col = "darkblue",
    border = "white",
    breaks = 50
  )
} else {
  message("No cis-eQTLs found to plot histogram.")
}

# Plot p-value histogram for trans-eQTLs (optional)
if (num_trans_eQTLs > 0) {
  hist(
    me_trans$trans$pvalue.hist$all,
    xlab = "P-value",
    main = "Trans-eQTL P-value Distribution",
    col = "darkgreen",
    border = "white",
    breaks = 50
  )
} else {
  message("No trans-eQTLs found to plot histogram.")
}