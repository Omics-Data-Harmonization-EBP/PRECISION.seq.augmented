# PRECISION.seq.augmented: An R Package for Performance Assessment of Depth Normalization Methods in MicroRNA Sequencing using Augmented Data


The PRECISION.seq.augmented R package offers a comprehensive framework for evaluating depth normalization methods in microRNA sequencing data analysis. 
It allows investigating normalization performance in both clustering and classification contexts and provides researchers with tools to assess how different normalization approaches affect analytical outcomes. 
Using AI-augmented miRNA-seq data derived from paired benchmark and test datasets, the package enables systematic comparison across controlled conditions with varying biological signal strengths and technical artifact magnitudes. 
PRECISION.seq.augmented implements multiple normalization techniques, clustering approaches, and classification algorithms, allowing researchers to identify optimal strategies for their specific analytical needs and reproduce findings from our publications. 
This package represents an essential resource for researchers seeking to maximize the reliability and reproducibility of insights derived from miRNA sequencing data.


## Installation

You can install the released version of *PRECISION.seq.augmented* directly from GitHub using `devtools`:

```R
devtools::install_github("Omics-Data-Harmonization-EBP/PRECISION.seq.augmented")
```

The R package `PoissonSeq` for PoissonSeq normalization was removed from CRAN, but you can install the archived version from GitHub:

```R
devtools::install_github("cran/PoissonSeq")
```

For successful installation, ensure all dependencies are properly installed. This package is based on R 4.2, and the following helper functions will install all required dependencies:

```R
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("BiocManager", "caret", "e1071", "glmnet", "pamr", "mclust", "cluster", "factoextra", "som", "curl", "digest"))

## from Bioconductor
Bioconductor.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE)
}
Bioconductor.packages(c("Biobase", "BiocGenerics", "edgeR", "EDASeq", "RUVSeq", "DESeq2", "preprocessCore", "sva"))

## from GitHub
devtools::install_github("cran/PoissonSeq")
```

## Main Functions

<!-- TODO Add link to documentation -->
The full *package documentation* with detailed function parameters and examples can be found on the [TODO: package documentation website](TODO).


### Data Access Functions
- `load.augmented.data()` - Load pre-generated augmented miRNA-seq datasets. Data will be downloaded from GitHub or loaded from storage if it was previously downloaded
- `cleanup.augmented.data()` - Remove cached augmented data
<!-- - `generate.scenario()` - Create customized datasets with user-defined parameters for biological signal and technical artifacts -->

### Core Object and Data Modulation
- `create.precision.cluster()` - Creates the main analysis object for clustering evaluation handling data and results
- `create.precision.classification()` - Creates the main analysis object for classification evaluation handling data and results
- `biological.effects()` - Apply biological effects with amplification factors
- `handling.effects()` - Modulate handling artifacts in test datasets

### Harmonization Methods
The package implements multiple data harmonization techniques applicable to both clustering and classification:

<!-- - `harmon.all()` - Apply all of the following harmonization methods sequentially -->
- `harmon.TC()` - Total Count normalization (scaling by library size)
- `harmon.UQ()` - Upper Quartile normalization (scaling by 75th percentile)
- `harmon.med()` - Median normalization (scaling by median count)
- `harmon.TMM()` - Trimmed Mean of M-values normalization (edgeR) 
- `harmon.DESeq()` - DESeq2 normalization (geometric mean approach)
- `harmon.PoissonSeq()` - PoissonSeq normalization (robust over-dispersed Poisson model)
- `harmon.QN()` - Quantile normalization
- `harmon.SVA()` - Surrogate Variable Analysis batch correction
- `harmon.RUVr()` - Remove Unwanted Variation (residuals)
- `harmon.RUVs()` - Remove Unwanted Variation (control samples)
- `harmon.RUVg()` - Remove Unwanted Variation (control genes)
- `harmon.ComBat.Seq()` - ComBat-seq batch effect adjustment

### Clustering Algorithms
The package implements multiple clustering approaches with various distance metrics:

<!-- - `cluster.all()` - Apply all of the following clustering methods sequentially -->
- `cluster.hc()` - Hierarchical clustering with multiple distance metrics (Euclidean distance, Pearson correlation, and Spearman correlation)
- `cluster.kmeans()` - K-means clustering with configurable starting points and iteration parameters
- `cluster.pam()` - Partitioning Around Medoids with the same distance metric options as hierarchical clustering (Euclidean, Pearson, Spearman).
- `cluster.som()` - Self-Organizing Maps for non-linear dimensionality reduction and clustering, particularly suited for high-dimensional data
- `cluster.mnm()` - Gaussian Mixture Model clustering with automated model selection using BIC.

### Classification Algorithms
The package implements multiple learning methods for sample classification to evaluate how normalization affects predictive performance across training and validation datasets:

<!-- - `classification.all()` - Apply all of the following clustering methods sequentially -->
- `classification.knn()` - k-Nearest Neighbor classification
- `classification.svm()` - Support Vector Machine classification
- `classification.pam()` - Prediction Analysis for Microarrays using nearest shrunken centroids, with built-in feature selection
- `classification.lasso()` -  Logistic Regression with LASSO regularization for automatic feature selection and reduced overfitting
- `classification.ranfor()` - Random Forest classification





<!-- 
## Basic Usage

```r
library(PRECISION.seq.augmented)

# Create analysis object
analysis <- create.precision.cluster(data = your_data, label = your_labels)

# Apply harmonization methods
analysis <- harmon.all(analysis)

# Perform clustering analysis
analysis <- cluster.all(analysis)

# Access results
ari_scores <- analysis@cluster.result
``` -->


<!-- ## Citation
TODO [Citation information to be added] -->


<!-- ### Key Features

TODO Update these

- **Comprehensive Analysis Pipeline**: Seamlessly integrate multiple harmonization, clustering, and classification approaches
- **Reproducible Workflows**: Standardized functions for consistent analysis
- **Performance Metrics**: Built-in calculation of Adjusted Rand Index (ARI) and Silhouette coefficients
- **Parallel Processing**: Support for multi-core processing to handle large datasets efficiently
- **Flexible Distance Metrics**: Multiple distance measures for clustering algorithms -->

<!-- ## Aims

The primary objectives of this package are to:

1. **Standardize RNA-seq data preprocessing** through multiple harmonization methods
2. **Provide comprehensive clustering analysis** with various algorithms and distance metrics
3. **Enable systematic performance evaluation** using established metrics (ARI, Silhouette)
4. **Facilitate comparative studies** of normalization and clustering approaches
5. **Support reproducible research** in genomics and bio-informatics -->
