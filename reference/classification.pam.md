# PAM Classification

This function applies classification using Prediction Analysis for
Microarrays (PAM) on the provided input object. It supports threshold
optimization using cross-validation or the use of all genes without
thresholding.

## Usage

``` r
classification.pam(object, threshold_method = "cv", vt.k = NULL, kfold = 5)
```

## Arguments

- object:

  A
  [precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing harmonized data. Must contain the slots
  `harmon.train.data` with harmonized training data and
  `harmon.test1.data` and `harmon.test2.data` with harmonized test data.

- threshold_method:

  A character string specifying the thresholding method. Options are: -
  `"cv"`: Use cross-validation to optimize the threshold. - `"none"`:
  Use all genes without applying a threshold.

- vt.k:

  A numeric vector of threshold values to evaluate during
  cross-validation. Only used if `threshold_method = "cv"`.

- kfold:

  An integer specifying the number of folds for cross-validation.

## Value

The input object updated with PAM classification results added to the
`classification.result` slot, including predicted classes and associated
metrics.
