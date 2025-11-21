# Create precision object for classification

Use this function to create a
[`precision`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/precision-class.md)
object for classification tasks. This object will contain raw training
and two sets of test data, along with their labels.

## Usage

``` r
create.precision.classification(
  traindata,
  testdata1,
  testdata2,
  trainlabel,
  testlabel1,
  testlabel2
)
```

## Arguments

- traindata:

  Training data matrix, (`p` by `n`), where `p` is the number of
  features and `n` is the number of samples.

- testdata1:

  First test data matrix, (`p` by `n`), where `p` is the number of
  features and `n` is the number of samples.

- testdata2:

  Second test data matrix, (`p` by `n`), where `p` is the number of
  features and `n` is the number of samples.

- trainlabel:

  Training data labels, a vector of length `n` (number of samples).

- testlabel1:

  First test data labels, a vector of length `n` (number of samples).

- testlabel2:

  Second test data labels, a vector of length `n` (number of samples).

## Value

A precision object for classification
