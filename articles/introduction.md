# Introduction to PRECISION.seq.augmented package

## Overview

The PRECISION.seq.augmented package provides a comprehensive framework
for clustering analysis of sequencing data through a two-stage pipeline:

1.  **Data Harmonization**: The first crucial step involves harmonizing
    the sequencing data to ensure consistency and compatibility across
    different samples. In this package, we have 9 different data
    harmonization methods, including:

- Trimmed mean of M-values (TMM)
- Total count harmonization (TC)
- Upper quantile harmonization (UQ)
- Median harmonization (med)
- Quantile harmonization (QU)
- DESeq harmonization (DESeq)
- PossionSeq harmonization (PossionSeq)
- Quantile normalization (QN)
- Surrogate Variable Analysis (SVA)
- RUV normalization (RUVg, RUVr, RUVs)
- ComBat-Seq harmonization (ComBat-Seq)

2.  **Clustering Analysis**: Following data harmonization, one can
    perform clustering analysis using the following clustering methods:

- K-means clustering (kmeans)
- Hierarchical clustering (hc)
- Self-organizing map (SOM)
- Gaussian mixture model (MNM)
- PAM clustering (pam.eucledian, pam.pearson, pam.spearman)

## Basic Usage

Here’s a basic example of how to perform clustering analysis:

``` r
# Step 0: Prepare the 'precision' object
example_cluster <- create.precision.cluster(data = data.test, label = data.group)

# Step 1: Harmonize via all the methods
example_cluster <- harmon.all(example_cluster)

# Step 2: Clustering via all the methods
example_cluster <- cluster.all(example_cluster)

# Step 4: Measure the consistency between the clustering labels and the true labels
cluster <- c("kmeans", "hc", "som", "mnm", "pam.euclidean", "pam.pearson", "pam.spearman")
harmon <- c("Raw", "TC", "UQ", "med", "TMM", "DESeq", "PoissonSeq", "QN", "RUVg", "RUVs", "RUVr")
ari_indexes <- data.frame(matrix(nrow = length(cluster), ncol = length(harmon)))
true_label <- as.factor(c(rep("MXF", 27), rep("PMFH", 27)))
for (i in seq_along(cluster)) {
  for (j in seq_along(harmon)) {
    eval(parse(text = paste0("est_cluster <-", "example_cluster@cluster.result", "$", cluster[i], "$", harmon[j])))
    ari_indexes[i, j] <- mclust::adjustedRandIndex(true_label, est_cluster)
  }
}
rownames(ari_indexes) <- cluster
colnames(ari_indexes) <- harmon
```
