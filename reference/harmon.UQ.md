# Perform Upper Quartile (UQ) Normalization

This function applies upper quartile normalization to miRNA-Seq data
contained within a precision object. The normalization scales the raw
read counts of each sample by the upper quartile value, ensuring
consistency across datasets.

## Usage

``` r
harmon.UQ(object)
```

## Arguments

- object:

  A
  [precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing raw data. Harmonization is applied to the slots
  `raw.train.data`, `raw.test1.data`, and `raw.test2.data` (if
  applicable), each containing a list with `data` and `label` elements.

## Value

A precision object with Upper Quartile Normalization results added to
the slots `harmon.train.data`, `harmon.test1.data`, and
`harmon.test2.data` (if applicable).

## Details

The function processes multiple datasets (e.g., "train", "test1",
"test2") if they are present in the precision object. The normalized
data is stored in the corresponding harmonization slots of the object.
