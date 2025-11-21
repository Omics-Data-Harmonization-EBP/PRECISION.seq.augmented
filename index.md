# PRECISION.seq.augmented: An R Package for Performance Assessment of Depth Normalization Methods in MicroRNA Sequencing using Augmented Data

The PRECISION.seq.augmented R package offers a comprehensive framework
for evaluating depth normalization methods in microRNA sequencing data
analysis. It allows investigating normalization performance in both
clustering and classification contexts and provides researchers with
tools to assess how different normalization approaches affect analytical
outcomes. Using AI-augmented miRNA-seq data derived from paired
benchmark and test datasets, the package enables systematic comparison
across controlled conditions with varying biological signal strengths
and technical artifact magnitudes. PRECISION.seq.augmented implements
multiple normalization techniques, clustering approaches, and
classification algorithms, allowing researchers to identify optimal
strategies for their specific analytical needs and reproduce findings
from our publications. This package represents an essential resource for
researchers seeking to maximize the reliability and reproducibility of
insights derived from miRNA sequencing data.

## Installation

You can install the released version of *PRECISION.seq.augmented*
directly from GitHub using `devtools`:

``` r
devtools::install_github("Omics-Data-Harmonization-EBP/PRECISION.seq.augmented")
```

The R package `PoissonSeq` for PoissonSeq normalization was removed from
CRAN, but you can install the archived version from GitHub:

``` r
devtools::install_github("cran/PoissonSeq")
```

For successful installation, ensure all dependencies are properly
installed. This package is based on R 4.2, and the following helper
functions will install all required dependencies:

``` r
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("BiocManager", "caret", "e1071", "glmnet", "pamr", "mclust", "cluster", "factoextra", "som", "digest"))

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

The full *package documentation* with detailed function parameters and
examples can be found on the [package documentation
website](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/).

### Data Access Functions

- [`load.augmented.data()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/load.augmented.data.md) -
  Load pre-generated augmented miRNA-seq datasets. Data will be
  downloaded from GitHub or loaded from storage if it was previously
  downloaded
- [`cleanup.augmented.data()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cleanup.augmented.data.md) -
  Remove cached augmented data

### Core Object and Data Modulation

- [`create.precision.cluster()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/create.precision.cluster.md) -
  Creates the main analysis object for clustering evaluation handling
  data and results
- [`create.precision.classification()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/create.precision.classification.md) -
  Creates the main analysis object for classification evaluation
  handling data and results
- [`biological.effects()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/biological.effects.md) -
  Apply biological effects with amplification factors
- [`handling.effects()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/handling.effects.md) -
  Modulate handling artifacts in test datasets

### Harmonization Methods

The package implements multiple data harmonization techniques applicable
to both clustering and classification:

- [`harmon.all()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.all.md) -
  Apply all of the following harmonization methods sequentially
- [`harmon.TC()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.TC.md) -
  Total Count normalization (scaling by library size)
- [`harmon.UQ()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.UQ.md) -
  Upper Quartile normalization (scaling by 75th percentile)
- [`harmon.med()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.med.md) -
  Median normalization (scaling by median count)
- [`harmon.TMM()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.TMM.md) -
  Trimmed Mean of M-values normalization (edgeR)
- [`harmon.DESeq()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.DESeq.md) -
  DESeq2 normalization (geometric mean approach)
- [`harmon.PoissonSeq()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.PoissonSeq.md) -
  PoissonSeq normalization (robust over-dispersed Poisson model)
- [`harmon.QN()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.QN.md) -
  Quantile normalization
- [`harmon.SVA()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.SVA.md) -
  Surrogate Variable Analysis batch correction
- [`harmon.RUVr()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.RUVr.md) -
  Remove Unwanted Variation (residuals)
- [`harmon.RUVs()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.RUVs.md) -
  Remove Unwanted Variation (control samples)
- [`harmon.RUVg()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.RUVg.md) -
  Remove Unwanted Variation (control genes)
- [`harmon.ComBat.Seq()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.ComBat.Seq.md) -
  ComBat-seq batch effect adjustment

### Clustering Algorithms

The package implements multiple clustering approaches with various
distance metrics:

- [`cluster.all()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.all.md) -
  Apply all of the following clustering methods sequentially
- [`cluster.hc()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.hc.md) -
  Hierarchical clustering with multiple distance metrics (Euclidean
  distance, Pearson correlation, and Spearman correlation)
- [`cluster.kmeans()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.kmeans.md) -
  K-means clustering with configurable starting points and iteration
  parameters
- [`cluster.pam()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.pam.md) -
  Partitioning Around Medoids with the same distance metric options as
  hierarchical clustering (Euclidean, Pearson, Spearman).
- [`cluster.som()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.som.md) -
  Self-Organizing Maps for non-linear dimensionality reduction and
  clustering, particularly suited for high-dimensional data
- [`cluster.mnm()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.mnm.md) -
  Gaussian Mixture Model clustering with automated model selection using
  BIC.

### Classification Algorithms

The package implements multiple learning methods for sample
classification to evaluate how normalization affects predictive
performance across training and validation datasets:

- [`classification.all()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/classification.all.md) -
  Apply all of the following clustering methods sequentially
- [`classification.knn()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/classification.knn.md) -
  k-Nearest Neighbor classification
- [`classification.svm()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/classification.svm.md) -
  Support Vector Machine classification
- [`classification.pam()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/classification.pam.md) -
  Prediction Analysis for Microarrays using nearest shrunken centroids,
  with built-in feature selection
- [`classification.lasso()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/classification.lasso.md) -
  Logistic Regression with LASSO regularization for automatic feature
  selection and reduced overfitting
- [`classification.ranfor()`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/classification.ranfor.md) -
  Random Forest classification
