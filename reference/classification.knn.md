# K-Nearest Neighbors (KNN) Classification

This function applies KNN classification to the harmonized data in the
input
[precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
object containing. It supports two thresholding methods:
cross-validation to optimize the number of neighbors (k) or a fixed k
value of 1.

## Usage

``` r
classification.knn(object, threshold_method = "cv", kfold = 5, folds = NULL)
```

## Arguments

- object:

  A
  [precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing harmonized data. Must contain the slots
  `harmon.train.data` with harmonized training data and
  `harmon.test1.data` and `harmon.test2.data` with harmonized test data.

- threshold_method:

  A character string specifying the thresholding method. Use `"cv"` for
  cross-validation to determine the optimal k, or `"none"` to use a
  fixed k = 1.

- kfold:

  An integer specifying the number of folds for cross-validation. This
  parameter is only used if `threshold_method` is set to `"cv"`.

- folds:

  An optional list of pre-specified folds for cross-validation. If
  provided, these folds will be used instead of generating new ones.

## Value

The input object updated with KNN classification results added to the
`classification.result` slot, including predicted classes and associated
metrics.
