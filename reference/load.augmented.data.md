# Load augmented data from GitHub (with integrity check)

This function loads the augmented dataset. If the dataset was not
previously downloaded, it will download it from the package GitHub
release. It checks the integrity of the downloaded file using SHA256
hash to ensure the file is not corrupted or tampered with. It supports
both persistent storage in the user's R data directory and temporary
(R-session) storage. The temporary storage is useful for testing
purposes or when you do not want to keep the data after the R session
ends. When using temporary storage (`temp=TRUE`), the file will be
automatically removed when the R session ends (using
[`quit()`](https://rdrr.io/r/base/quit.html)). To delete the cached
data, use the
[`cleanup.augmented.data`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cleanup.augmented.data.md)
function. This function also allows for re-downloading the file even if
it exists, based on the \`force.redownload\` parameter.
`load.augmented.data()` returns the loaded R object, which is a list
containing the augmented benchmark and test datasets.

## Usage

``` r
load.augmented.data(temp = FALSE, force.redownload = FALSE)
```

## Arguments

- temp:

  Logical. If TRUE, store file in temporary directory and remove when R
  session ends. Default is FALSE (persistent storage).

- force.redownload:

  Logical. If TRUE, force re-download even if file exists and passes
  integrity check.

## Value

The loaded R object, which is a list containing the augmented benchmark
and test datasets.

## See also

[`cleanup.augmented.data`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/cleanup.augmented.data.md)
for deleting cached data.

## Examples

``` r
if (FALSE) { # \dontrun{

# First-time load: downloads & caches persistently
augmented.data <- load.augmented.data()

# Uses cached file if valid
augmented.data <- load.augmented.data()

# Force re-download
augmented.data <- load.augmented.data(force.redownload = TRUE)

# Store in tempdir() instead of persistent (will be cleaned up when session ends)
augmented.data <- load.augmented.data(temp = TRUE)

# Access benchmark and test datasets
benchmark <- augmented.data$benchmark
test <- augmented.data$test

# Manually cleanup cached data
cleanup.augmented.data()
} # }
```
