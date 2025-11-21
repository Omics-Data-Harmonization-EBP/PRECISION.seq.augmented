# Partitioning Around Medoids (PAM) clustering for harmonized data

This function performs Partitioning Around Medoids (PAM) clustering on
the harmonized training data using specified distance measures
(euclidean, pearson, or spearman). It adds the clustering results to the
`cluster.result` slot of the precision object.

## Usage

``` r
cluster.pam(object, k = NULL, distance = "euclidean")
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

- distance:

  Distance measure to use for clustering. Options are "euclidean",
  "pearson", or "spearman". Default is "euclidean".

## Value

Updated precision object with PAM clustering results added to the
`cluster.result` slot.
