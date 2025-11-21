# Apply all clustering methods to a precision object

This function applies all clustering methods to the training data
contained in a precision object. The clustering results are added to the
corresponding `cluster.result` slot in the object for each training
dataset.

## Usage

``` r
cluster.all(object, k = NULL)
```

## Arguments

- object:

  A
  [`precision`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing harmonized training data in the slot
  `harmon.train.data`.

- k:

  *Integer.* Number of clusters. If NULL or not given, uses the number
  of unique labels in training data.

## Value

Updated precision object with all clustering results added to the
`cluster.result` slot.

## Details

The following clustering methods are applied to the data:

- Hierarchical Clustering using multiple distance metrics (Euclidean,
  Pearson, and Spearman)
  ([`cluster.hc`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.hc.md))

- Gaussian Mixture Models (MNM)
  ([`cluster.mnm`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.mnm.md))

- K-Means Clustering
  ([`cluster.kmeans`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.kmeans.md))

- Self-Organizing Maps
  ([`cluster.som`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.som.md))

- Partitioning Around Medoids (Euclidean, Pearson, and Spearman)
  ([`cluster.pam`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cluster.pam.md))
