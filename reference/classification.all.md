# Apply all classification methods to a precision object

This function applies all classification methods to the training data
contained in a precision object. The classification results are added to
the corresponding `classification.result` slot in the object for each
harmonized training and testing dataset.

## Usage

``` r
classification.all(object, threshold_method = "cv")
```

## Arguments

- object:

  A
  [`precision`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing harmonized training and testing data in the slots
  `harmon.train.data`, and `harmon.test1.data` and `harmon.test2.data`,
  respectively.

- threshold_method:

  A character string specifying the thresholding method. Options are: -
  `"cv"`: Use cross-validation to optimize the threshold. - `"none"`:
  Use all genes without applying a threshold.

## Value

Updated precision object with all classification results added to the
`classification.result` slot.
