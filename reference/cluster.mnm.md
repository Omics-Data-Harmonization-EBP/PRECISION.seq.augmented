# Gaussian Mixture Model clustering for harmonized data

This function performs Gaussian Mixture Model (GMM) clustering on the
harmonized training data using the \`mclust\` package.

## Usage

``` r
cluster.mnm(object, k = NULL)
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

Updated precision object with Gaussian Mixture Model clustering results
added to the `cluster.result` slot.
