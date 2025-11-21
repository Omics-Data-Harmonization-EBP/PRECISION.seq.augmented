# Self-Organizing Map (SOM) clustering for harmonized data

This function performs Self-Organizing Map (SOM) clustering on the
harmonized training data. It uses the \`som\` package to create a SOM
model and returns the clustering results.

## Usage

``` r
cluster.som(object, k = NULL)
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

Updated precision object with SOM clustering results added to the
`cluster.result` slot.
