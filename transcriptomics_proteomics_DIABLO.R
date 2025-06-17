# Install mixOmics if not already installed
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  install.packages("mixOmics")
}

library(mixOmics)

# Example data loading
transcriptomics <- read.csv("transcriptomics.csv", row.names = 1)
proteomics <- read.csv("proteomics.csv", row.names = 1)
# Ensure samples are rows, features are columns in both

# Load phenotype vector (factor)
Y <- factor(read.csv("phenotype.csv")$Class)

# Make sure sample orders match in both datasets and phenotype
all(rownames(transcriptomics) == rownames(proteomics))  # Should be TRUE
all(rownames(transcriptomics) == names(Y))              # Should be TRUE

# Optional: log-transform or normalize here if needed

data_list <- list(
  transcriptomics = transcriptomics,
  proteomics = proteomics
)

# design matrix defines strength of relationships between blocks (0 to 1).
design <- matrix(c(0, 1,
                   1, 0), 
                 ncol = 2, nrow = 2, byrow = TRUE)
colnames(design) <- rownames(design) <- names(data_list)

set.seed(123)  # for reproducibility

tune <- tune.block.splsda(
  X = data_list, 
  Y = Y, 
  ncomp = 2,
  test.keepX = list(transcriptomics = c(5, 10, 15), proteomics = c(5, 10, 15)),
  design = design,
  validation = 'Mfold',
  folds = 5,
  nrepeat = 10,
  dist = 'centroids.dist'
)

# Optimal number of features per block
list_keepX <- tune$choice.keepX

diablo_model <- block.splsda(
  X = data_list,
  Y = Y,
  ncomp = 2,
  keepX = list_keepX,
  design = design
)

# Sample plot showing discrimination on components
plotIndiv(diablo_model, legend = TRUE, title = 'DIABLO: Sample Plot')

# Circos plot to show relationships between features across datasets
circosPlot(diablo_model, cutoff = 0.7)

# Heatmap of selected features
plotHeatmap(diablo_model, ncomp = 1)

selected_features <- selectVar(diablo_model, comp = 1)

# Transcriptomics features selected
selected_features$transcriptomics$name

# Proteomics features selected
selected_features$proteomics$name

perf <- perf(diablo_model, validation = "Mfold", folds = 5, nrepeat = 10, progressBar = TRUE)

plot(perf)

# additionally for tuning

set.seed(123)
tune <- tune.block.splsda(
  X = data_list,
  Y = Y,
  ncomp = 2,
  test.keepX = list(transcriptomics = c(5,10,15), proteomics = c(5,10,15)),
  design = design,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  dist = "centroids.dist"
)

list_keepX <- tune$choice.keepX
print(tune$error.rate)       # Classification error rates per tested combination
print(tune$choice.keepX)     # Selected optimal number of features per dataset
print(tune$choice.ncomp)     # Selected optimal number of components
