# K-means clustering for harmonized data

This function performs K-means clustering on the harmonized training
data and returns the clustering results.

## Usage

``` r
cluster.kmeans(object, k = NULL)
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

Updated precision object with k-means clustering results added to the
`cluster.result` slot.
