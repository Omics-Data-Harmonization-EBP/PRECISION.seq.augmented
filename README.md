# PRECISION.seq.augmented: An R Package for Performance Assessment of Depth Normalization Methods in MicroRNA Sequencing using Augmented Data


PRECISION.seq.augmented is an R package designed for assessing the performance of depth normalization methods in microRNA sequencing using AI-augmented data.
The package provides a unified framework for applying multiple normalization techniques and clustering and classification algorithms, enabling researchers to systematically evaluate and compare different normalization methods.


## Installation

You can install the released version of *PRECISION.seq.augmented* directly from GitHub using `devtools` by:

```R
devtools::install_github(Omics-Data-Harmonization-EBP/PRECISION.seq.augmented")
```

The R package `PoissonSeq` for PoissonSeq normalization was removed from CRAN, but you can install the archived version from github using:

```R
devtools::install_github("cran/PoissonSeq")
```

If the package cannot be installed successfully, please ensure that the dependency packages are installed. This package is based on R 4.2, and the R codes for installing the dependent packages are:

```R
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("BiocManager", "magrittr", "tidyverse", "mclust", "aricode", "RSKC", "cluster", "factoextra", "pamr", "som"))

## from Bioconductor
Bioconductor.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE)
}
Bioconductor.packages(c("Biobase", "BiocGenerics","preprocessCore", "edgeR", "DESeq2", "affy", "sva", "RUVSeq", "EDASeq", "limma", "vsn"))
```

## Main Functions

<!-- TODO Add link to documentation -->
The full *package documentation* can be found [here](TODO).

### Core Object Creation
- `create.precision.cluster()` - Creates the main analysis object with data and labels

### Harmonization Methods
The package implements multiple data harmonization techniques:

<!-- - `harmon.all()` - Apply all of the following harmonization methods sequentially -->
- `harmon.TC()` - Total Count normalization
- `harmon.UQ()` - Upper Quartile normalization  
- `harmon.med()` - Median normalization
- `harmon.TMM()` - Trimmed Mean of M-values (edgeR) normalization
- `harmon.DESeq()` - DESeq2 normalization
- `harmon.PoissonSeq()` - PoissonSeq normalization
- `harmon.QN()` - Quantile normalization
- `harmon.QN.frozen()` - Quantile normalization with frozen parameters
- `harmon.sva()` - Surrogate Variable Analysis harmonization
- `harmon.RUVr()` - Remove Unwanted Variation (residuals)
- `harmon.RUVs()` - Remove Unwanted Variation (control samples)
- `harmon.RUVg()` - Remove Unwanted Variation (control genes)
- `harmon.ComBat.Seq()` - ComBat-seq batch effect adjustment

### Clustering Algorithms
Multiple clustering approaches with various distance metrics:

<!-- - `cluster.all()` - Apply all of the following clustering methods sequentially -->
- `cluster.hc()` - Hierarchical clustering (euclidean, pearson, spearman distances)
- `cluster.kmeans()` - K-means clustering
- `cluster.pam()` - Partitioning Around Medoids (euclidean, pearson, spearman distances)
- `cluster.som()` - Self-Organizing Maps
- `cluster.mnm()` - Gaussian Mixture model clustering

### Classification Algorithms
Multiple classification approaches:

<!-- - `classification.all()` - Apply all of the following clustering methods sequentially -->
- `classification.knn()` - k-Nearest Neighbor classification
- `classification.svm()` - Support Vector Machine classification
- `classification.lasso()` - Logistic Regression classification using LASSO
- `classification.ranfor()` - Random Forest classification

### Utility Functions
- `biological.effects()` - Apply biological effects with amplification factors



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
