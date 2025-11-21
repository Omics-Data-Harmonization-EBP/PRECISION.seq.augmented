# Classification with LASSO Logistic Regression

This function applies LASSO classification to the harmonized data in the
input
[precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
object containing. It supports two thresholding methods:
cross-validation to optimize the lambda parameter or using all genes
without thresholding.

## Usage

``` r
classification.lasso(object, threshold_method = "cv", kfold = 5)
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
  cross-validation to determine the optimal k, or `"none"` to use all
  genes without thresholding. If `"none"` is used, a very small lambda
  value is applied to effectively include all genes in the model, which
  may be useful for datasets where all genes are relevant.

- kfold:

  An integer specifying the number of folds for cross-validation. This
  parameter is only used if `threshold_method` is set to `"cv"`.

## Value

The input object updated with LASSO classification results added to the
`classification.result` slot, including predicted classes and associated
metrics.
