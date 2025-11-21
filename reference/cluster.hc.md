# Hierarchical clustering for harmonized data

This function performs hierarchical clustering on the harmonized
training data using euclidean, pearson, or spearman distance measures.
It adds the clustering results to the `cluster.result` slot.

## Usage

``` r
cluster.hc(object, k = NULL, distance = "euclidean")
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

Updated precision object with hierarchical clustering results added to
the `cluster.result` slot.
