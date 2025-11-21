# Add biological effects to benchmark data

This function simulates biological effects on benchmark data by
amplifying the signal between *two* groups based on a specified
amplification factor. It modifies the input data and returns the
amplified data along with group information.

## Usage

``` r
biological.effects(benchmark, group, c)
```

## Arguments

- benchmark:

  Numeric matrix of expression values, (`p` by `n`), where `p` is the
  number of features and `n` is the number of samples.

- group:

  Factor or vector of character strings defining sample groups for each
  sample. Must have exactly two levels.

- c:

  Numeric amplification factor

## Value

List containing modified data, group information, and amplification
factor.
