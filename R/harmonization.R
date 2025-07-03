#' @import methods
#' @import sva
#' @import EDASeq
#' @import edgeR
#' @import RUVSeq
#' @import Biobase
#' @import BiocGenerics
#' @import tidyverse
#' @import DESeq2
#' @import preprocessCore
#' @import PoissonSeq
#' @importFrom magrittr %>%
NULL

##########################################################
# Precision Class
##########################################################

#' Precision Class
#'
#' @slot raw.train.data Training data and labels
#' @slot raw.test1.data First test dataset and labels
#' @slot raw.test2.data Second test dataset and labels
#' @slot harmon.train.data List of harmonized training data
#' @slot harmon.test1.data List of harmonized test1 data
#' @slot harmon.test2.data List of harmonized test2 data
#' @slot classification.result Classification results
#' @slot cluster.result Clustering results
#' @export
precision <- methods::setClass(
  "precision",
  slots = c(
    raw.train.data = "ANY", # For clustering, only consider train data
    raw.test1.data = "ANY",
    raw.test2.data = "ANY",
    harmon.train.data = "list",
    harmon.test1.data = "list",
    harmon.test2.data = "list",
    classification.result = "list",
    cluster.result = "list"
  )
)

#' Create precision object for classification
#'
#' @param traindata Training data matrix
#' @param testdata1 First test data matrix
#' @param testdata2 Second test data matrix
#' @param trainlabel Training data labels
#' @param testlabel1 First test data labels
#' @param testlabel2 Second test data labels
#' @return A precision object
#' @export
create.precision.classification <- function(traindata, testdata1, testdata2,
                                            trainlabel, testlabel1, testlabel2) {
  object <- methods::new(
    Class = "precision",
    raw.train.data = list(data = traindata, label = trainlabel),
    raw.test1.data = list(data = testdata1, label = testlabel1),
    raw.test2.data = list(data = testdata2, label = testlabel2)
  )
  return(object)
}

#' Create precision object for clustering
#'
#' @param data Data matrix
#' @param label Labels
#' @return A precision object
#' @export
create.precision.cluster <- function(data, label) {
  object <- methods::new(
    Class = "precision",
    raw.train.data = list(data = data, label = label)
  )
  return(object)
}

##########################################################
# Harmonization functions
##########################################################

#' Internal helper functions for TMM normalization
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom magrittr %>%
#' @noRd
.harmon.method.TMM <- function(raw, groups) {
  d <- edgeR::DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  ) %>%
    edgeR::calcNormFactors(method = "TMM", refColumn = 1)

  list(
    dat.harmonized = edgeR::cpm(d),
    scaling.factor = d$samples$norm.factors * d$samples$lib.size / 1.0e6
  )
}

#' Internal helper functions for TMM normalization (frozen parameters)
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom magrittr %>%
#' @noRd
.harmon.method.TMM.frozen <- function(train.data, train.group, test.data, test.group) {
  d <- edgeR::DGEList(
    counts = cbind(train.data[, 1], test.data),
    group = factor(c(train.group[1], test.group)),
    genes = rownames(test.data)
  ) %>%
    edgeR::calcNormFactors(method = "TMM", refColumn = 1)

  list(
    dat.harmonized = edgeR::cpm(d)[, -1],
    scaling.factor = (d$samples$norm.factors * d$samples$lib.size / 1.0e6)[-1]
  )
}

#' Apply TMM Harmonization to Precision Object
#'
#' Trimmed Mean of M-values (TMM) normalization for miRNA-Seq data using the edgeR package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR package website}.
#'
#' This function applies TMM harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with TMM harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.TMM <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$TMM <-
        .harmon.method.TMM(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Apply TMM Harmonization with frozen parameters to Precision Object
#'
#' Trimmed Mean of M-values (TMM) normalization for miRNA-Seq data using the edgeR package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR package website}.
#'
#' This function applies TMM harmonization using frozen parameters
#' to the raw data contained in a precision object.
#' The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with TMM harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.TMM.frozen <- function(object) {
  # Process training data
  object@harmon.train.data$TMM.frozen <- .harmon.method.TMM(
    object@raw.train.data$data,
    object@raw.train.data$label
  )

  # Process test datasets if they exist
  for (test_set in c("test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", test_set, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", test_set, ".data"))$TMM.frozen <-
        .harmon.method.TMM.frozen(
          object@raw.train.data$data,
          object@raw.train.data$label,
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper function for total count harmonization
#' @importFrom edgeR DGEList cpm
#' @noRd
.harmon.method.TC <- function(raw, groups) {
  d <- edgeR::DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  )

  list(
    dat.harmonized = edgeR::cpm(d, normalized.lib.sizes = FALSE),
    scaling.factor = d$samples$lib.size / 1e6
  )
}

#' Apply Total Count (TC) Normalization to Precision Object
#'
#' This function applies total count normalization to miRNA-Seq data
#' contained within a precision object. The normalization scales the raw
#' read counts of each sample by the total number of reads in the sample,
#' ensuring consistency across datasets.
#'
#' The function processes multiple datasets (e.g., "train", "test1", "test2")
#' if they are present in the precision object. The normalized data is stored
#' in the corresponding harmonization slots of the object.
#'
#' @param object A \link{precision} object. This object must contain raw data in
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and/or \code{raw.test2.data}.
#'   Each of these slots should include a \code{data} component (the dataset)
#'   and a \code{label} component (associated labels).
#' @return The input precision object with total count (TC) harmonization
#'   results added to the slots \code{harmon.train.data}, \code{harmon.test1.data},
#'   and/or \code{harmon.test2.data} under the \code{TC} component.
#' @export
harmon.TC <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$TC <-
        .harmon.method.TC(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper function for upper quartile harmonization
#' @importFrom edgeR DGEList
#' @noRd
.harmon.method.UQ <- function(raw, groups) {
  d <- edgeR::DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  )

  q.factor <- apply(d$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor / 1.0e6

  list(
    dat.harmonized = t(t(raw) / scaling.factor),
    scaling.factor = scaling.factor
  )
}

#' Perform Upper Quartile (UQ) Normalization
#'
#' This function applies upper quartile normalization to miRNA-Seq data
#' contained within a precision object. The normalization scales the raw
#' read counts of each sample by the upper quartile value, ensuring
#' consistency across datasets.
#'
#' The function processes multiple datasets (e.g., "train", "test1", "test2")
#' if they are present in the precision object. The normalized data is stored
#' in the corresponding harmonization slots of the object.
#'
#' @param object A \link{precision} object containing raw miRNA-Seq data.
#' @return A precision object with upper quartile (UQ) normalization results
#'   added to the harmonization slots (e.g., \code{harmon.train.data},
#'   \code{harmon.test1.data}, \code{harmon.test2.data}).
#' @export
harmon.UQ <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$UQ <-
        .harmon.method.UQ(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper function for median harmonization
#' @importFrom edgeR DGEList
#' @noRd
.harmon.method.med <- function(raw, groups) {
  d <- edgeR::DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  )

  m.factor <- apply(d$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1.0e6

  list(
    dat.harmonized = t(t(raw) / scaling.factor),
    scaling.factor = scaling.factor
  )
}

#' Apply Median Harmonization to Precision Object
#'
#' This function applies upper quartile normalization to miRNA-Seq data
#' contained within a precision object. The normalization scales the raw
#' read counts of each sample by the upper quartile value, ensuring
#' consistency across datasets.
#'
#' The function processes multiple datasets (e.g., "train", "test1", "test2")
#' if they are present in the precision object. The normalized data is stored
#' in the corresponding harmonization slots of the object.
#'
#' @param object A \link{precision} object. This object must contain raw data in
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and/or \code{raw.test2.data}.
#'   Each of these slots should include a \code{data} component (the dataset)
#'   and a \code{label} component (associated labels).
#' @return The input precision object with median harmonization results
#'   added to the slots \code{harmon.train.data}, \code{harmon.test1.data}, and/or
#'   \code{harmon.test2.data} under the \code{med} component.
#' @export
harmon.med <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$med <-
        .harmon.method.med(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper function for DESeq harmonization
#' @import DESeq2
#' @importFrom magrittr %>%
#' @noRd
.harmon.method.DESeq <- function(raw, groups) {
  # Prepare data for DESeq2
  colnames(raw) <- paste0("V", 1:(dim(raw)[2]))
  condition <- data.frame(
    SampleName = colnames(raw),
    Condition = factor(groups),
    row.names = colnames(raw)
  )

  # Create and process DESeq dataset
  d <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw,
    colData = condition,
    design = ~Condition
  ) %>%
    DESeq2::estimateSizeFactors()

  list(
    dat.harmonized = DESeq2::counts(d, normalized = TRUE),
    scaling.factor = DESeq2::sizeFactors(d)
  )
}

#' Apply DESeq Harmonization to Precision Object
#'
#' DESeq normalization for miRNA-Seq data using the DESeq2 package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2 package website}.
#'
#' This function applies DESeq harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with DESeq harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.DESeq <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$DESeq <-
        .harmon.method.DESeq(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper function for PoissonSeq harmonization
#' @import PoissonSeq
#' @noRd
.harmon.method.PoissonSeq <- function(raw) {
  scaling.factor <- PoissonSeq::PS.Est.Depth(raw)

  list(
    dat.harmonized = t(t(raw) / scaling.factor),
    scaling.factor = scaling.factor
  )
}

#' Apply PoissonSeq to Precision Object
#'
#' PoissonSeq normalization for miRNA-Seq data using the PoissonSeq package.
#' For more information, visit the
#' \href{https://github.com/cran/PoissonSeq}{PoissonSeq GitHub repository}.
#'
#' This function applies PoissonSeq harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with PoissonSeq harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.PoissonSeq <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$PoissonSeq <-
        .harmon.method.PoissonSeq(
          raw_data$data
        )
    }
  }
  object
}

#' Internal helper functions for quantile normalization
#' @importFrom preprocessCore normalize.quantiles
#' @noRd
.harmon.method.QN <- function(raw) {
  dat.harmonized <- preprocessCore::normalize.quantiles(as.matrix(raw))
  colnames(dat.harmonized) <- colnames(raw)
  rownames(dat.harmonized) <- rownames(raw)

  list(
    dat.harmonized = dat.harmonized
  )
}

#' Internal helper functions for quantile normalization (frozen parameters)
#' @importFrom preprocessCore normalize.quantiles
#' @noRd
.harmon.method.QN.frozen <- function(rawtrain, rawtest) {
  # Harmonize training data
  dat.harmonized <- preprocessCore::normalize.quantiles(as.matrix(rawtrain))
  colnames(dat.harmonized) <- colnames(rawtrain)
  rownames(dat.harmonized) <- rownames(rawtrain)

  # Apply frozen harmonization to test data
  ref.dis <- as.numeric(sort(dat.harmonized[, 1]))
  dat.harmonized.test <- apply(rawtest, 2, function(x) {
    ord <- rank(x)
    ref.dis[ord]
  })
  rownames(dat.harmonized.test) <- rownames(dat.harmonized)

  list(
    dat.harmonized = dat.harmonized.test
  )
}

#' Apply Standard Quantile normalization to Precision Object
#'
#' Quantile normalization for miRNA-Seq data using the preprocessCore package.
#' For more information, visit documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/preprocessCore.html}{preprocessCore Bioconductor page}.
#'
#' This function applies quantile normalization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with QN harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.QN <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$QN <-
        .harmon.method.QN(
          raw_data$data
        )
    }
  }
  object
}


#' Apply Standard Quantile normalization with frozen parameters
#'
#' Quantile normalization for miRNA-Seq data using the preprocessCore package.
#' For more information, visit documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/preprocessCore.html}{preprocessCore Bioconductor page}.
#'
#' This function applies quantile normalization with frozen parameters to the
#' raw data contained in a precision object. The harmonization results are added
#' to the corresponding slots in the object for each dataset
#' ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with QN harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.QN.frozen <- function(object) {
  # Process training data
  object@harmon.train.data$QN.frozen <- .harmon.method.QN(
    object@raw.train.data$data
  )

  # Process test datasets if they exist
  for (test_set in c("test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", test_set, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", test_set, ".data"))$QN.frozen <-
        .harmon.method.QN.frozen(
          object@raw.train.data$data,
          raw_data$data
        )
    }
  }
  object
}

#' Internal helper function for SVA
#' @import sva
#' @noRd
.harmon.method.SVA <- function(raw, groups) {
  # Filter out rows with zero sums
  dat.sva <- raw[rowSums(raw) > 0, ]

  # Create model matrices
  mod1 <- model.matrix(~groups)
  mod0 <- model.matrix(~1, data.frame(mod1))
  dat0 <- as.matrix(dat.sva)

  # Calculate surrogate variables
  n.sv <- sva::num.sv(dat0, mod1)
  invisible(capture.output(
    svseq <- sva::svaseq(dat0, mod1, mod0, n.sv = n.sv)$sv
  ))

  # Adjust data using surrogate variables
  adjust <- cbind(mod1, svseq)
  hat <- solve(t(adjust) %*% adjust) %*% t(adjust)
  beta <- (hat %*% t(dat.sva))
  P <- ncol(mod1)

  list(
    dat.harmonized = dat.sva - t(as.matrix(adjust[, -c(1:P)]) %*% beta[-c(1:P), ]),
    adjust.factor = svseq
  )
}

#' Apply SVA Harmonization to Precision Object
#'
#' SVA normalization for miRNA-Seq data using the sva package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/sva.html}{sva package website}.
#'
#' This function applies SVA harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with SVA harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.SVA <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$SVA <-
        .harmon.method.SVA(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper functions for RUVg harmonization
#' @import RUVSeq
#' @import EDASeq
#' @import edgeR
#' @importFrom BiocGenerics counts
#' @importFrom magrittr %>%
#' @noRd
.harmon.method.RUVg <- function(raw, groups) {
  # Ensure required packages are available
  if (!suppressMessages(require("Biobase"))) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  condition <- factor(groups)

  # Create expression set and design matrix
  set <- EDASeq::newSeqExpressionSet(
    as.matrix(raw),
    phenoData = data.frame(condition, row.names = colnames(raw))
  )
  design <- model.matrix(~condition,
    data = data.frame(condition, row.names = colnames(raw))
  )

  # Prepare DGE analysis
  y <- edgeR::DGEList(counts = counts(set), group = condition) %>%
    edgeR::calcNormFactors(method = "upperquartile")

  if (any(is.infinite(y$samples$norm.factors))) {
    y <- edgeR::DGEList(counts = counts(set), group = condition) %>%
      edgeR::calcNormFactors(method = "none")
  }

  # Estimate dispersions and fit model
  y <- y %>%
    edgeR::estimateGLMCommonDisp(design) %>%
    edgeR::estimateGLMTagwiseDisp(design)

  # Identify control genes
  fit <- edgeR::glmFit(y, design)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  top <- edgeR::topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15 * nrow(raw))]))]

  # Apply RUVg harmonization
  t <- RUVSeq::RUVg(set, spikes, k = 1)

  list(
    dat.harmonized = EDASeq::normCounts(t),
    adjust.factor = t$W
  )
}

#' Internal helper functions for RUVs harmonization
#' @import RUVSeq
#' @import EDASeq
#' @import edgeR
#' @importFrom BiocGenerics counts
#' @importFrom magrittr %>%
#' @noRd
.harmon.method.RUVs <- function(raw, groups) {
  # Ensure required packages are available
  if (!suppressMessages(require("Biobase"))) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  condition <- factor(groups)

  # Create expression set and design matrix
  set <- EDASeq::newSeqExpressionSet(
    as.matrix(raw),
    phenoData = data.frame(condition, row.names = colnames(raw))
  )
  design <- model.matrix(~condition,
    data = data.frame(condition, row.names = colnames(raw))
  )

  # Prepare DGE analysis
  y <- edgeR::DGEList(counts = counts(set), group = condition) %>%
    edgeR::calcNormFactors(method = "upperquartile")

  if (any(is.infinite(y$samples$norm.factors))) {
    y <- edgeR::DGEList(counts = counts(set), group = condition) %>%
      edgeR::calcNormFactors(method = "none")
  }

  # Estimate dispersions and fit model
  y <- y %>%
    edgeR::estimateGLMCommonDisp(design) %>%
    edgeR::estimateGLMTagwiseDisp(design)

  # Identify control genes and differences
  fit <- edgeR::glmFit(y, design)
  lrt <- edgeR::glmLRT(fit, coef = 2)
  top <- edgeR::topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15 * nrow(raw))]))]
  differences <- RUVSeq::makeGroups(condition)

  # Apply RUVs harmonization
  t <- RUVSeq::RUVs(set, rownames(raw), k = 1, differences)

  list(
    dat.harmonized = EDASeq::normCounts(t),
    adjust.factor = t$W
  )
}

#' Internal helper functions for RUVr harmonization
#' @import RUVSeq
#' @import EDASeq
#' @import edgeR
#' @importFrom BiocGenerics counts
#' @importFrom magrittr %>%
#' @noRd
.harmon.method.RUVr <- function(raw, groups) {
  # Ensure required packages are available
  if (!suppressMessages(require("Biobase"))) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  condition <- factor(groups)

  # Create expression set and design matrix
  set <- EDASeq::newSeqExpressionSet(
    as.matrix(raw),
    phenoData = data.frame(condition, row.names = colnames(raw))
  )
  design <- model.matrix(~condition, data = Biobase::pData(set))

  # Prepare DGE analysis
  y <- edgeR::DGEList(counts = counts(set), group = condition) %>%
    edgeR::calcNormFactors(method = "upperquartile")

  if (any(is.infinite(y$samples$norm.factors))) {
    y <- edgeR::DGEList(counts = counts(set), group = condition) %>%
      edgeR::calcNormFactors(method = "none")
  }

  # Estimate dispersions and get residuals
  y <- y %>%
    edgeR::estimateGLMCommonDisp(design) %>%
    edgeR::estimateGLMTagwiseDisp(design)

  fit <- edgeR::glmFit(y, design)
  residuals <- residuals(fit, type = "deviance")

  # Apply RUVr harmonization
  set <- EDASeq::betweenLaneNormalization(set, which = "upper")
  t <- RUVSeq::RUVr(set, rownames(raw), k = 1, residuals = residuals)

  list(
    dat.harmonized = EDASeq::normCounts(t),
    adjust.factor = t$W
  )
}

#' Apply RUVg Harmonization to Precision Object
#'
#' RUV normalization for miRNA-Seq data using the RUV package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/RUVSeq.html}{RUVSeq package website}.
#'
#' This function applies RUVg (gene-wise) harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with RUVg harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.RUVg <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$RUVg <-
        .harmon.method.RUVg(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Apply RUVs Harmonization to Precision Object
#'
#' RUV normalization for miRNA-Seq data using the RUV package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/RUVSeq.html}{RUVSeq package website}.
#'
#' This function applies RUVs (centered (technical) replicate/negative control
#' samples) harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with RUVs harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.RUVs <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$RUVs <-
        .harmon.method.RUVs(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Apply RUVr Harmonization to Precision Object
#'
#' RUV normalization for miRNA-Seq data using the RUV package.
#' For more information, visit the documentation on the
#' \href{https://bioconductor.org/packages/release/bioc/html/RUVSeq.html}{RUVSeq package website}.
#'
#' This function applies RUVr (using residuals) harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with RUVr harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.RUVr <- function(object) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$RUVr <-
        .harmon.method.RUVr(
          raw_data$data,
          raw_data$label
        )
    }
  }
  object
}

#' Internal helper function for ComBat-Seq normalization
#' @import sva
#' @noRd
.harmon.method.ComBat.Seq <- function(raw, batches) {
  # Apply ComBat-Seq normalization
  dat.harmonized <- sva::ComBat_seq(raw, batch = batches)

  list(
    dat.harmonized = dat.harmonized
  )
}

#' Apply ComBat-Seq to Precision Object
#'
#' ComBat-Seq batch effect adjustment for miRNA-Seq data using the ComBat-Seq method.
#' For more information, visit the documentation on the
#' \href{https://github.com/zhangyuqing/ComBat-seq}{ComBat-Seq documentation}.
#'
#' This function applies ComBat-Seq harmonization to the raw data contained in a
#' precision object. The harmonization results are added to the corresponding
#' slots in the object for each dataset ("train", "test1", "test2") if they exist.
#'
#' @param object A \link{precision} object containing raw data. The object must have
#'   slots named \code{raw.train.data}, \code{raw.test1.data}, and \code{raw.test2.data}
#'   (if applicable), each containing a list with \code{data} and \code{label} elements.
#' @return A precision object with ComBat-Seq harmonization results added to the
#'   slots \code{harmon.train.data}, \code{harmon.test1.data}, and \code{harmon.test2.data}
#'   (if applicable).
#' @export
harmon.ComBat.Seq <- function(object, batches) {
  # Process all datasets if they exist
  for (dataset in c("train", "test1", "test2")) {
    raw_data <- slot(object, paste0("raw.", dataset, ".data"))
    if (!is.null(raw_data)) {
      slot(object, paste0("harmon.", dataset, ".data"))$ComBat.Seq <-
        .harmon.method.ComBat.Seq(
          raw_data$data,
          batches
        )
    }
  }
  object
}
