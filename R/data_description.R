#' MiRNA Sequencing Benchmark Data
#'
#' Myxofibrosarcoma (MXF) and pleomorphic malignant fibrous histiocytoma (PMFH)
#' are the two most common and aggressive subtypes of genetically complex soft
#' tissue sarcoma.
#' This dataset includes three libraries used for sequencing 54 individual tumor
#' samples. Library preparation and read capture were each processed by a single
#' experienced technician in one run.
#'
#'
#' @format A data frame with 1033 rows and 54 columns. Here are the examples of the column and row naming rule:
#' \describe{
#'   \item{MXF2516_D13}{One sample belonging to MXF and library D with sample ID MXF2516 and library ID D13.}
#'   \item{hsa-let-7a-2*}{Gene ID.}
#' }
"data.benchmark"

#' MiRNA Sequencing Test Data
#'
#' Myxofibrosarcoma (MXF) and pleomorphic malignant fibrous histiocytoma (PMFH)
#' are the two most common and aggressive subtypes of genetically complex soft
#' tissue sarcoma.
#' MiRNAs for the same 54 tumors used for the benchmark data were re-sequenced
#' using neither uniform handling nor balanced library assignment.
#' In this study these samples were sequenced in the order of sample collection and processed in multiple runs.
#'
#'
#' @format A dataframe with 1033 rows and 54 columns. Here are the examples of the column and row naming rule:
#' \describe{
#'   \item{MXF2516}{One sample belonging to MXF sample ID MXF2516.}
#'   \item{hsa-let-7a-2*}{Gene ID.}
#' }
"data.test"

#' MiRNA Information
#'
#' @format Information of the microRNAs in test and benchmark dataset, including
#' their names, exact sequence, length and GC content.
"data.miR.info"

#' Sample Labels of the test and benchmark data
#'
#' @format A set of labels for the test and benchmark samples, either "MXF"
#' (myxofibrosarcoma) or "PMFH" (pleomorphic malignant fibrous histiocytoma).
"data.group"


#' Load augmented data from GitHub (with integrity check)
#'
#' This function loads the augmented dataset. If the dataset was not previously downloaded,
#' it will download it from the package GitHub release.
#' It checks the integrity of the downloaded file using SHA256 hash to ensure
#' the file is not corrupted or tampered with.
#' It supports both persistent storage in the user's R data directory and
#' temporary (R-session) storage.
#' The temporary storage is useful for testing purposes or when you do not want
#' to keep the data after the R session ends. When using temporary storage (\code{temp=TRUE}),
#' the file will be automatically removed when the R session ends (using \code{quit()}).
#' To delete the cached data, use the \code{\link{cleanup.augmented.data}} function.
#' This function also allows for re-downloading the file even if it exists,
#' based on the `force.redownload` parameter.
#' \code{load.augmented.data()} returns the loaded R object, which is a list
#' containing the augmented benchmark and test datasets.
#'
#' @param temp Logical. If TRUE, store file in temporary directory and remove when R session ends.
#' Default is FALSE (persistent storage).
#' @param force.redownload Logical. If TRUE, force re-download even if file
#' exists and passes integrity check.
#' @return The loaded R object, which is a list containing the augmented
#' benchmark and test datasets.
#' @importFrom digest digest
#' @importFrom tools R_user_dir
#' @seealso \code{\link{cleanup.augmented.data}} for deleting cached data.
#' @examples \dontrun{
#'
#' # First-time load: downloads & caches persistently
#' augmented.data <- load.augmented.data()
#'
#' # Uses cached file if valid
#' augmented.data <- load.augmented.data()
#'
#' # Force re-download
#' augmented.data <- load.augmented.data(force.redownload = TRUE)
#'
#' # Store in tempdir() instead of persistent (will be cleaned up when session ends)
#' augmented.data <- load.augmented.data(temp = TRUE)
#'
#' # Access benchmark and test datasets
#' benchmark <- augmented.data$benchmark
#' test <- augmented.data$test
#'
#' # Manually cleanup cached data
#' cleanup.augmented.data()
#' }
#' @export
load.augmented.data <- function(temp = FALSE, force.redownload = FALSE) {
  ## Configuration

  # Filename and URL of the augmented datasets
  filename <- "MSKpair_augmented.rds"
  url <- "https://github.com/Omics-Data-Harmonization-EBP/PRECISION.seq.augmented/releases/download/Data/MSKpair_augmented.rds"
  expected.hash <- "9105e02ae4765b3b3561cd654777998c1244443067ed0981fab33130c0c214b4"

  persistent.dir <- tools::R_user_dir("PRECISION.seq.augmented", which = "data")
  temp.dir <- tempdir()

  primary.dir <- if (temp) temp.dir else persistent.dir
  dir.create(primary.dir, recursive = TRUE, showWarnings = FALSE)

  primary.file <- file.path(primary.dir, filename)
  persistent.file <- file.path(persistent.dir, filename)
  temp.file <- file.path(temp.dir, filename)

  # Helper function for integrity check
  is.valid.file <- function(path) {
    if (!file.exists(path)) {
      return(FALSE)
    }
    actual.hash <- digest::digest(file = path, algo = "sha256")
    identical(actual.hash, expected.hash)
  }

  ## Check if the file already exists and is valid
  if (!force.redownload) {
    # Check primary storage
    if (is.valid.file(primary.file)) {
      message("Loading data from ", primary.file)
      return(readRDS(primary.file))
    }
    # Check persistent storage, if valid copy to temp
    if (is.valid.file(persistent.file)) {
      message("Copying from persistent storage to ", primary.file)
      file.copy(persistent.file, primary.file, overwrite = TRUE)
      return(readRDS(primary.file))
    }
    # Check temporary storage, if valid copy to persistent
    if (is.valid.file(temp.file)) {
      message("Copying from temporary storage to ", primary.file)
      file.copy(temp.file, primary.file, overwrite = TRUE)
      return(readRDS(primary.file))
    }
  }

  ## Download the file if it doesn't exist or is invalid
  message("Downloading augmented data from GitHub Release...")
  download.file(url, destfile = primary.file, mode = "wb", method = "libcurl", quiet = FALSE)
  message("Download complete. Saved to: ", primary.file)
  message("To clean up cached data, use: cleanup.augmented.data()")

  # Verify hash after download
  if (!is.valid.file(primary.file)) {
    stop("Downloaded file failed integrity check. The file may be corrupted or tampered with.")
  }

  return(readRDS(primary.file))
}


#' Clean Up Cached Augmented Data
#'
#' This function deletes the cached augmented data file from both the
#' persistent storage and temporary directory.
#'
#' @return NULL. Used for side effects.
#' @importFrom tools R_user_dir
#' @seealso \code{\link{load.augmented.data}} for loading the data.
#' @examples
#' \dontrun{
#' # Download and cache the augmented data
#' augmented.data <- load.augmented.data()
#' cleanup.augmented.data()
#' }
#' @export
cleanup.augmented.data <- function() {
  filename <- "MSKpair_augmented.rds"
  dirs <- c(tools::R_user_dir("PRECISION.seq.augmented", which = "data"), tempdir())

  for (d in dirs) {
    f <- file.path(d, filename)
    if (file.exists(f)) {
      unlink(f, force = TRUE)
      message("Deleted: ", f)
    }
  }
}
