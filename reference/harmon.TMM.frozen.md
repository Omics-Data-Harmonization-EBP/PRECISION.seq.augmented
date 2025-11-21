# Apply TMM Harmonization with frozen parameters to Precision Object

Trimmed Mean of M-values (TMM) normalization for miRNA-Seq data using
the edgeR package. For more information, visit the documentation on the
[edgeR package
website](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

## Usage

``` r
harmon.TMM.frozen(object)
```

## Arguments

- object:

  A
  [precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing raw data. Harmonization is applied to the slots
  `raw.train.data`, `raw.test1.data`, and `raw.test2.data` (if
  applicable), each containing a list with `data` and `label` elements.

## Value

A precision object with TMM harmonization results added to the slots
`harmon.train.data`, `harmon.test1.data`, and `harmon.test2.data` (if
applicable).

## Details

This function applies TMM harmonization using frozen parameters to the
raw data contained in a precision object. The harmonization results are
added to the corresponding slots in the object for each dataset
("train", "test1", "test2") if they exist.
