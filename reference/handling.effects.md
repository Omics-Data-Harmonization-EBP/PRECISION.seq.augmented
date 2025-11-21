# Add handling effects to benchmark data

This function simulates handling effects by modifying the clean input
data based on the differences between the benchmark and test data.

## Usage

``` r
handling.effects(clean.input, benchmark, test, group, d)
```

## Arguments

- clean.input:

  Numeric matrix of clean input data, (`p` by `n`), where `p` is the
  number of features and `n` is the number of samples.

- benchmark:

  Numeric matrix of benchmark data, (`p` by `n`), where `p` is the
  number of features and `n` is the number of samples.

- test:

  Numeric matrix of test data, (`p` by `n`), where `p` is the number of
  features and `n` is the number of samples.

- group:

  Factor or vector of character strings defining sample groups for each
  sample.

- d:

  Numeric effect strength factor

## Value

List containing modified data, group information, and effect size
