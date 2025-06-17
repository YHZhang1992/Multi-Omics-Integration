# Load necessary libraries
library(WGCNA)
library(dplyr)
library(ggplot2)
library(corrplot)

# --- 0. Setup and Parallel Processing ---

# Allow WGCNA to use multiple CPU cores for faster computation
# Adjust the number of threads based on your system's capabilities.
# This should be run once per R session.
enableWGCNAThreads()
# Or specify a fixed number: enableWGCNAThreads(nThreads = 8)

# Set a random seed for reproducibility
set.seed(123)

cat("WGCNA threads enabled and seed set.\n\n")

# --- 1. Data Preparation: Loading and Preprocessing Multiple Datasets ---
# Consensus WGCNA requires multiple gene expression datasets.
# For demonstration, we'll generate synthetic data. In real scenarios,
# you would load your expression matrices (e.g., from CSV, tab-separated files).

# Define parameters for synthetic data
num_datasets <- 3      # Number of datasets (e.g., different cohorts/conditions)
num_genes <- 2000      # Number of genes (rows)
num_samples_per_dataset <- c(100, 120, 90) # Number of samples for each dataset

# Create lists to store expression data and clinical traits
multiExpr <- list()
multiTrait <- list()

cat("Generating synthetic gene expression and clinical trait data...\n")
for (i in 1:num_datasets) {
  # Simulate gene expression data: genes in rows, samples in columns
  # Values are often normalized counts (e.g., VST, rlog, log2(TPM+1))
  expr_matrix <- matrix(rnorm(num_genes * num_samples_per_dataset[i], mean = 10, sd = 2),
                        nrow = num_genes,
                        ncol = num_samples_per_dataset[i])
  rownames(expr_matrix) <- paste0("Gene_", 1:num_genes)
  colnames(expr_matrix) <- paste0("Sample_", i, "_", 1:num_samples_per_dataset[i])

  # Introduce some 'co-expression' patterns (simple simulation)
  # For example, a small block of genes co-expressing
  if (i == 1) { # Stronger module in dataset 1
    expr_matrix[1:50, ] <- expr_matrix[1:50, ] + rnorm(50 * num_samples_per_dataset[i], mean = 3, sd = 0.5)
  } else if (i == 2) { # Weaker module in dataset 2
    expr_matrix[1:50, ] <- expr_matrix[1:50, ] + rnorm(50 * num_samples_per_dataset[i], mean = 1, sd = 0.2)
  } else { # Another distinct module in dataset 3
    expr_matrix[100:150, ] <- expr_matrix[100:150, ] + rnorm(50 * num_samples_per_dataset[i], mean = 2, sd = 0.3)
  }

  multiExpr[[i]] <- list(data = as.data.frame(expr_matrix))

  # Simulate clinical trait data (e.g., disease status, age)
  # Rows are samples, columns are traits
  trait_data <- data.frame(
    Sample = colnames(expr_matrix),
    DiseaseStatus = sample(c(0, 1), num_samples_per_dataset[i], replace = TRUE, prob = c(0.6, 0.4)),
    Age = sample(20:70, num_samples_per_dataset[i], replace = TRUE)
  )
  rownames(trait_data) <- trait_data$Sample
  multiTrait[[i]] <- list(data = trait_data)
}

# Assign names to the list elements for clarity
names(multiExpr) <- paste0("Dataset", 1:num_datasets)
names(multiTrait) <- paste0("Dataset", 1:num_datasets)

cat("Synthetic data generation complete. Verifying dimensions...\n")
for (i in 1:num_datasets) {
  cat(paste("Dataset", i, " - Expression dimensions:", dim(multiExpr[[i]]$data)[1], "genes,", dim(multiExpr[[i]]$data)[2], "samples\n"))
  cat(paste("Dataset", i, " - Trait dimensions:", dim(multiTrait[[i]]$data)[1], "samples,", dim(multiTrait[[i]]$data)[2], "traits\n"))
}
cat("\n")

# --- 1.1. Data Check and Cleaning (Crucial Step for Real Data) ---
# In a real analysis, you MUST perform these steps:
# - Ensure all datasets have the SAME set of genes (rows). Filter/intersect if necessary.
# - Ensure sample IDs (columns) are unique across datasets if combined later.
# - Remove genes with too many missing values or low variance across samples.
# - Handle outliers (e.g., by removing samples or using robust normalization).
# - Match sample names between expression and trait data.

# For synthetic data, we assume genes are already matched and no missing values.
# Let's transpose expression data to have samples in rows and genes in columns,
# as required by WGCNA's 'multiExpr' structure for its internal processing.
# Also, remove the 'Sample' column from trait data if it's not a trait.
for (i in 1:num_datasets) {
  multiExpr[[i]]$data <- t(multiExpr[[i]]$data)
  multiTrait[[i]]$data <- multiTrait[[i]]$data %>% select(-Sample)
}

# Check dimensions again after transposition
cat("Dimensions after transposition (samples in rows, genes in columns):\n")
for (i in 1:num_datasets) {
  cat(paste("Dataset", i, " - Expression dimensions:", dim(multiExpr[[i]]$data)[1], "samples,", dim(multiExpr[[i]]$data)[2], "genes\n"))
}
cat("\n")

# --- 2. Network Construction Parameters: Choosing the Soft-Thresholding Power ($\beta$) ---
# The soft-thresholding power is crucial for constructing scale-free networks.
# It's usually determined by analyzing the scale-free topology fit index.
# For consensus WGCNA, we need a common power across all datasets or choose from a list.

# Test different powers for each dataset (optional, but good practice for real data)
# powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# sft_list <- list()
# for (i in 1:num_datasets) {
#   cat(paste("Picking soft threshold for Dataset", i, "...\n"))
#   sft_list[[i]] <- pickSoftThreshold(multiExpr[[i]]$data, powerVector = powers, verbose = 5)
#   # Plot the results for each dataset to visually inspect
#   # plot(sft_list[[i]]$fitIndices[, 1], -sign(sft_list[[i]]$fitIndices[, 3]) * sft_list[[i]]$fitIndices[, 2],
#   #      xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
#   #      type = "n", main = paste("Scale independence for Dataset", i));
#   # text(sft_list[[i]]$fitIndices[, 1], -sign(sft_list[[i]]$fitIndices[, 3]) * sft_list[[i]]$fitIndices[, 2],
#   #      labels = powers, cex = 0.9, col = "red");
#   # abline(h = 0.85, col = "red") # Good R^2 fit generally > 0.85
# }

# For consensus, a single power is usually chosen that works reasonably well across all datasets.
# Or, you can use pickSoftThreshold.fromList for a global fit.
# For this example, let's just pick a common power (e.g., 6)
# In real analysis, you would carefully examine the plots from pickSoftThreshold for each dataset.
softPower <- 6 # A common choice for unsigned networks. Adjust based on your data.

cat(paste("Chosen soft-thresholding power (beta):", softPower, "\n\n"))

# --- 3. Consensus Module Detection ---
# This is the core WGCNA step to find modules preserved across datasets.
# The `blockwiseConsensusModules` function performs all necessary steps:
# - Adjacency calculation for each dataset
# - Topological Overlap Matrix (TOM) calculation for each dataset
# - Consensus TOM calculation
# - Hierarchical clustering of consensus TOM
# - Dynamic Tree Cut to identify modules
# - Module merging
# - Consensus Module Eigengene (CME) calculation

# Define parameters for module detection
minModuleSize <- 30 # Minimum number of genes in a module
mergeCutHeight <- 0.25 # Module merging threshold (0.25 means modules with correlation > 0.75 are merged)
networkType <- "unsigned" # Use "unsigned", "signed", or "signed hybrid"
corType <- "pearson" # Correlation type ("pearson" or "bicor")
TOMType <- "signed" # Use "unsigned" or "signed" TOM. Should be consistent with networkType.
# For unsigned networks, it's common to use unsigned TOM or signed TOM (more specific).
# Here, using "signed" TOM means we differentiate between positive and negative correlations.

cat("Detecting consensus modules...\n")
bwnet <- blockwiseConsensusModules(
  multiExpr,
  power = softPower,
  minModuleSize = minModuleSize,
  mergeCutHeight = mergeCutHeight,
  networkType = networkType,
  corType = corType,
  TOMType = TOMType,
  pamRespectsDendro = FALSE, # Usually FALSE for blockwise
  saveTOMs = TRUE,           # Save TOMs for each dataset
  saveConsensusTOMs = TRUE,  # Save the consensus TOM
  saveTOMFileBase = "ConsensusTOM", # Base name for saving TOM files
  verbose = 3,               # Level of verbosity
  maxBlockSize = num_genes,  # Set to num_genes if you want to analyze all genes in one block (for smaller datasets)
                             # For larger datasets (>5000-10000 genes), blocks are necessary.
  checkTopology = TRUE
)
cat("Consensus module detection complete.\n\n")

# --- 4. Accessing and Summarizing Consensus Module Results ---

# Consensus module colors and labels
consensusColors <- bwnet$colors # Module assignments for each gene
consensusMEs <- bwnet$multiMEs  # Consensus Module Eigengenes

cat("Consensus module assignments:\n")
table(consensusColors)
cat("\n")

# Plotting the dendrogram and module colors for the consensus network
cat("Plotting consensus dendrogram and module colors...\n")
sizeGrWindow(12, 9) # Adjust window size for plotting
plotDendroAndColors(
  bwnet$dendrograms[[1]], # Use the dendrogram from the first block/dataset (often representative)
  consensusColors[bwnet$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Consensus Gene Dendrogram and Module Colors"
)
graphics.off() # Close the plot window

# --- 5. Relating Consensus Modules to Clinical Traits ---
# This step correlates the Consensus Module Eigengenes (CMEs) with your clinical traits.

cat("Relating consensus modules to clinical traits...\n")

# Prepare multiTrait for correlation (ensure numerical traits where appropriate)
# WGCNA expects traits to be columns in a data frame, samples in rows, matching multiExpr.
# If your traits are factors, convert them to numeric (e.g., 0/1 for disease status)

# Convert DiseaseStatus to numeric (0/1) for correlation
for (i in 1:num_datasets) {
  multiTrait[[i]]$data$DiseaseStatus <- as.numeric(multiTrait[[i]]$data$DiseaseStatus) - 1 # Ensure 0 and 1
}

# Calculate consensus ME-trait correlation
# First, ensure that the order of samples in MEs matches the order of samples in traits
# And that traits and MEs are within the same list structure.
multiCME_trait_cor <- list()
multiCME_trait_p <- list()

for (set in 1:num_datasets) {
  # Match samples
  current_MEs <- orderMEs(bwnet$multiMEs[[set]]$data)
  current_traits <- multiTrait[[set]]$data[rownames(current_MEs), ]

  # Calculate correlation
  modTraitCor <- cor(current_MEs, current_traits, use = "pairwise.complete.obs")
  modTraitP <- corPvalueStudent(modTraitCor, num_samples_per_dataset[set])

  multiCME_trait_cor[[set]] <- modTraitCor
  multiCME_trait_p[[set]] <- modTraitP
}

names(multiCME_trait_cor) <- paste0("Dataset", 1:num_datasets)
names(multiCME_trait_p) <- paste0("Dataset", 1:num_datasets)


# Display a summary of correlations and p-values for a chosen dataset (e.g., Dataset 1)
cat("\nConsensus Module-Trait Correlations (Dataset 1):\n")
print(multiCME_trait_cor[[1]])
cat("\nConsensus Module-Trait P-values (Dataset 1):\n")
print(multiCME_trait_p[[1]])


# Visualize module-trait relationships as a heatmap
cat("\nPlotting module-trait relationships heatmap...\n")
sizeGrWindow(10, 6)
# Combine correlations for plotting (e.g., concatenate them or pick one)
# For simplicity, let's plot for Dataset 1
textMatrix = paste(signif(multiCME_trait_cor[[1]], 2), "\n(",
                   signif(multiCME_trait_p[[1]], 1), ")", sep = "");
dim(textMatrix) = dim(multiCME_trait_cor[[1]])

labeledHeatmap(Matrix = multiCME_trait_cor[[1]],
               xLabels = names(multiTrait[[1]]$data),
               yLabels = names(multiCME_trait_cor[[1]]),
               ySymbols = names(multiCME_trait_cor[[1]]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = paste("Consensus Module-Trait Relationships (Dataset 1)"))
graphics.off()


# --- 6. Intramodular Connectivity (Module Membership) & Hub Gene Identification ---
# Gene significance (GS) measures correlation with a trait.
# Module membership (MM) measures correlation with the module eigengene.
# Hub genes are highly connected genes within their module (high MM).

cat("\nCalculating Module Membership and Hub Genes...\n")
# For a specific dataset (e.g., Dataset 1)
datExpr_set1 <- multiExpr[[1]]$data
datTrait_set1 <- multiTrait[[1]]$data

# Recalculate MEs for Dataset 1 with correct sample order
MEs_set1 <- moduleEigengenes(datExpr_set1, consensusColors)$eigengenes
geneModuleMembership <- as.data.frame(cor(datExpr_set1, MEs_set1, use = "pairwise.complete.obs"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr_set1)));

# Add gene names
names(geneModuleMembership) = paste("MM", names(MEs_set1), sep="");
names(MMPvalue) = paste("p.MM", names(MEs_set1), sep="");

# Gene significance for a trait (e.g., DiseaseStatus in Dataset 1)
traitNames <- names(datTrait_set1) # Get trait names
geneTraitSignificance <- as.data.frame(cor(datExpr_set1, datTrait_set1$DiseaseStatus, use = "pairwise.complete.obs"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr_set1)));

names(geneTraitSignificance) = paste("GS.", traitNames[1], sep="");
names(GSPvalue) = paste("p.GS.", traitNames[1], sep="");

# Combine MM, GS, and gene names
gene_info <- as.data.frame(rownames(datExpr_set1))
names(gene_info) <- "Gene"
gene_info$Module <- consensusColors
gene_info <- cbind(gene_info, geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue)

# Example: Find top hub genes for a specific module (e.g., "turquoise")
module <- "turquoise"
genes_in_module <- gene_info %>% filter(Module == module)
hub_genes_in_module <- genes_in_module %>%
  arrange(desc(MMturquoise)) %>% # Sort by module membership in descending order
  head(10) # Get top 10

cat(paste("\nTop 10 hub genes in the '", module, "' module (Dataset 1):\n", sep=""))
print(hub_genes_in_module)

# Plot correlation between Gene Significance and Module Membership
cat("\nPlotting GS vs MM for a representative module (e.g., turquoise)...\n")
# Choose a trait index and module color
trait <- "DiseaseStatus"
module <- "turquoise"
column <- match(module, names(MEs_set1));
traitColumn <- match(trait,names(datTrait_set1));
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[, column]),
                   abs(geneTraitSignificance[, traitColumn]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene Significance for", trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
graphics.off()

# --- 7. Exporting Network for External Visualization (Optional) ---
# You can export the network to formats compatible with Cytoscape.

# Select a module to export (e.g., "turquoise")
module_to_export <- "turquoise"
probes <- names(multiExpr[[1]]$data)
inModule <- (consensusColors == module_to_export);
modProbes <- probes[inModule];

# Select the corresponding topological overlap matrix (TOM)
# You need to load the saved TOM file from step 3 if you didn't keep it in memory
# For simplicity, we'll use the TOM from the first dataset's block for export.
# In a real scenario, you might construct a consensus TOM for export or use a specific dataset's TOM.
# If maxBlockSize was num_genes (as in this example), bwnet$TOMs[[1]] would contain the full TOM.
# If you run into memory issues with TOMs, you might need to reload specific ones.
# bwnet$TOMs contains a list of TOMs for each block. Here we assume one block.
# Let's get the TOM from the first block of the first dataset
TOM_matrix <- bwnet$consTOMs[[1]]$data # Access the consensus TOM for the first block

# Filter TOM for selected module
modTOM = TOM_matrix[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network as a tab-separated file for Cytoscape
cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", module_to_export, ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", module_to_export, ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02, # Only export edges above this TOM threshold
  nodeAttr = gene_info[match(modProbes, gene_info$Gene), c("Gene", "Module", "MMturquoise", "GS.DiseaseStatus")]
)

cat(paste("\nNetwork for module '", module_to_export, "' exported to CytoscapeInput-edges/nodes files.\n", sep=""))