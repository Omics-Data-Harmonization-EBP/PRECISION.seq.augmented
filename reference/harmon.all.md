# Apply all harmonization methods to a precision object

This function applies all harmonization methods to the raw data
contained in a precision object. The harmonization results are added to
the corresponding slots in the object for each dataset ("train",
"test1", "test2") if they exist.

## Usage

``` r
harmon.all(object, add.raw = TRUE)
```

## Arguments

- object:

  A
  [precision](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
  object containing raw data. Harmonization is applied to the slots
  `raw.train.data`, `raw.test1.data`, and `raw.test2.data` (if
  applicable), each containing a list with `data` and `label` elements.

- add.raw:

  Logical. If TRUE, add the raw data to the harmonized data

## Value

A precision object with harmonization results added to the slots
`harmon.train.data`, `harmon.test1.data`, and `harmon.test2.data` (if
applicable).

## Details

The following harmonization methods are applied to the data:

- No Harmonization (if `add.raw = TRUE`)

- Total Count (TC)
  ([`harmon.TC`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.TC.md))

- Upper Quartile (UQ)
  ([`harmon.UQ`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.UQ.md))

- Median (median)
  ([`harmon.med`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.med.md))

- Trimmed Median of Means (TMM)
  ([`harmon.TMM`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.TMM.md))

- DESeq
  ([`harmon.DESeq`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.DESeq.md))

- PoissonSeq
  ([`harmon.PoissonSeq`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.PoissonSeq.md))

- Quantile Normalization (QN)
  ([`harmon.QN`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.QN.md))

- Remove Unwanted Variation (RUVg, RUVr, and RUVs)
  ([`harmon.RUVg`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.RUVg.md)),
  ([`harmon.RUVr`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.RUVr.md)),
  ([`harmon.RUVs`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/harmon.RUVs.md))
