# Clean Up Cached Augmented Data

This function deletes the cached augmented data file from both the
persistent storage and temporary directory.

## Usage

``` r
cleanup.augmented.data()
```

## Value

NULL. Used for side effects.

## See also

[`load.augmented.data`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/load.augmented.data.md)
for loading the data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download and cache the augmented data
augmented.data <- load.augmented.data()
cleanup.augmented.data()
} # }
```
