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
