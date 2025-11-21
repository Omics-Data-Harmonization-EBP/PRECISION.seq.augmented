# Apply Standard Quantile normalization to Precision Object

Quantile normalization for miRNA-Seq data using the preprocessCore
package. For more information, visit documentation on the
[preprocessCore Bioconductor
page](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html).

## Usage

``` r
harmon.QN(object)
```

## Arguments

- object:

  A
  [precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing raw data. Harmonization is applied to the slots
  `raw.train.data`, `raw.test1.data`, and `raw.test2.data` (if
  applicable), each containing a list with `data` and `label` elements.

## Value

A precision object with QN harmonization results added to the slots
`harmon.train.data`, `harmon.test1.data`, and `harmon.test2.data` (if
applicable).

## Details

This function applies quantile normalization to the raw data contained
in a precision object. The harmonization results are added to the
corresponding slots in the object for each dataset ("train", "test1",
"test2") if they exist.
