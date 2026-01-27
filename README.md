# PFStatistics
`PFStatistics` is a personalized variant testing and selection framework that addresses how a variant impact on a phenotype varies continuously along some heterogeneity features, such as genetic ancestry. 

## Installation

To install from GitHub:
```r
# If devtools is not installed, first install it
install.packages("devtools")

# Then install PFStatistics from GitHub
devtools::install_github("juliefrwang/PFStatistics")
```

## Usage

### Loading the Package

After installation, load the package with:

```r
library(PFStatistics)
```

### Main Functions Overview

#### `get_importance_matrices`

This main function performs the core analysis by calculating the importance matrices for selected genetic variants.

**Usage:**

```r
results <- get_importance_matrices(
  genetic_variants = genetic_variants,
  genetic_variants_knockoff = genetic_variants_knockoff,
  additional_covariates = pcs,
  Z = eur, 
  y = y,
  n_folds = 10,  
  FDR_rate = 0.1  
)
```

**Arguments:**
- `genetic_variants`: Matrix of original genetic variants (SNPs), with dimensions n x p where n represents the number of samples (individuals), and p represents the number of SNPs.
- `genetic_variants_knockoff`: Matrix of knockoff genetic variants,  structured identically to genetic_variants with dimensions n x p, where each column is a knockoff version of the corresponding SNP in genetic_variants.
- `additional_covariates`: Matrix of additional covariates, with dimensions n x c, where c represents the number of covariates, it can be NULL.
- `Z`: Heterogeneity variable (e.g., estimated genetic ancestry or PCs),  dimensions n x m, where m is the number of heterogeneity variables
- `y`: Outcome variable, a vector of length n, representing the response variable (e.g., disease status or continuous phenotype) for each individual.
- `n_folds`: Number of folds for cross-validation.
- `FDR_rate`: Target false discovery rate for feature selection.

The function returns a list containing:
- `coefs`: Extracted model coefficients.
- `scaled_selection_matrix`: A matrix indicating scaled selection of SNPs.
- `selection_matrix`: A binary matrix indicating SNP selection.
- `W_statistic_matrix`: W-statistic matrix for SNPs.

#### `generate_knockoff_data`

This function generates knockoff variables for the given genetic variant matrix. Knockoff variables are used in the Lasso model to control false discovery rates.

**Usage:**

```r
knockoff_matrix <- generate_knockoff_data(genetic_variants_matrix)
```

**Arguments:**
- `genetic_variants_matrix`: Matrix containing SNP data where rows represent individuals (samples) ane columns represent genetic variants.


## Example

Below is an example of using `PFStatistics` to analyze SNP data:

```r
# Load SNP
snp_filepath <- "snp_data.csv"
snp_data <- read.csv(snp_filepath)

# Extract genetic_variants matrix and generate knockoff data
pcs <- snp_data[, c("PC1", "PC2", "PC3", "PC4")]
eur <- snp_data$EUR
y <- snp_data$AD
genetic_variants <- snp_data[, grepl("chr", colnames(snp_data))]  # Extract columns with "chr" in their names
genetic_variants_knockoff <- generate_knockoff_data(genetic_variants)

# Perform importance calculation
results <- get_importance_matrices(
  genetic_variants = genetic_variants,
  genetic_variants_knockoff = genetic_variants_knockoff,
  additional_covariates = pcs,
  Z = eur, # or other variables like 'pcs'
  y = y,
  n_folds = 10,
  FDR_rate = 0.1
)
```
## Contact