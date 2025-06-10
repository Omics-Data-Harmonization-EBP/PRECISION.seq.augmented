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

# Internal helper functions for TMM normalization
.harmon.method.TMM <- function(raw, groups) {
  d <- DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  ) %>%
    calcNormFactors(method = "TMM", refColumn = 1)

  list(
    dat.harmonized = cpm(d),
    scaling.factor = d$samples$norm.factors * d$samples$lib.size / 1e6
  )
}

.harmon.method.TMM.frozen <- function(train.data, train.group, test.data, test.group) {
  d <- DGEList(
    counts = cbind(train.data[, 1], test.data),
    group = factor(c(train.group[1], test.group)),
    genes = rownames(test.data)
  ) %>%
    calcNormFactors(method = "TMM", refColumn = 1)

  list(
    dat.harmonized = cpm(d)[, -1],
    scaling.factor = (d$samples$norm.factors * d$samples$lib.size / 1e6)[-1]
  )
}

#' Apply standard TMM normalization
#'
#' @param object A precision object containing raw data
#' @return A precision object with TMM harmonization results added
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

#' Apply TMM normalization with frozen parameters
#'
#' @param object A precision object containing raw data
#' @return A precision object with TMM frozen harmonization results added
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

# Internal helper function for total count harmonization
.harmon.method.TC <- function(raw, groups) {
  d <- DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  )

  list(
    dat.harmonized = cpm(d, normalized.lib.sizes = FALSE),
    scaling.factor = d$samples$lib.size / 1e6
  )
}

#' Apply total count harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with total count (TC) harmonization results added
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

# Internal helper function for upper quartile harmonization
.harmon.method.UQ <- function(raw, groups) {
  d <- DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  )

  q.factor <- apply(d$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor / 1e6

  list(
    dat.harmonized = t(t(raw) / scaling.factor),
    scaling.factor = scaling.factor
  )
}

#' Apply upper quartile harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with upper quartile (UQ) harmonization results added
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

# Internal helper function for median harmonization
.harmon.method.med <- function(raw, groups) {
  d <- DGEList(
    counts = raw,
    group = factor(groups),
    genes = rownames(raw)
  )

  m.factor <- apply(d$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1e6

  list(
    dat.harmonized = t(t(raw) / scaling.factor),
    scaling.factor = scaling.factor
  )
}

#' Apply median harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with median harmonization results added
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

# Internal helper function for DESeq harmonization
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

#' Apply DESeq harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with DESeq harmonization results added
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

# Internal helper function for PoissonSeq harmonization
.harmon.method.PoissonSeq <- function(raw) {
  scaling.factor <- PoissonSeq::PS.Est.Depth(raw)

  list(
    dat.harmonized = t(t(raw) / scaling.factor),
    scaling.factor = scaling.factor
  )
}

#' Apply PoissonSeq harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with PoissonSeq harmonization results added
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

# Internal helper functions for quantile normalization
.harmon.method.QN <- function(raw) {
  dat.harmonized <- preprocessCore::normalize.quantiles(as.matrix(raw))
  colnames(dat.harmonized) <- colnames(raw)
  rownames(dat.harmonized) <- rownames(raw)

  list(
    dat.harmonized = dat.harmonized
  )
}

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

#' Apply standard quantile normalization
#'
#' @param object A precision object containing raw data
#' @return A precision object with quantile normalization (QN) results added
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

#' Apply quantile normalization with frozen parameters
#'
#' @param object A precision object containing raw data
#' @return A precision object with frozen quantile normalization (QN) results added
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

# Internal helper function for SVA
.harmon.method.SVA <- function(raw, groups) {
  # Filter out rows with zero sums
  dat.sva <- raw[rowSums(raw) > 0, ]

  # Create model matrices
  mod1 <- model.matrix(~groups)
  mod0 <- model.matrix(~1, data.frame(mod1))
  dat0 <- as.matrix(dat.sva)

  # Calculate surrogate variables
  n.sv <- num.sv(dat0, mod1)
  invisible(capture.output(
    svseq <- svaseq(dat0, mod1, mod0, n.sv = n.sv)$sv
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

#' Apply SVA harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with SVA harmonization results added
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

# Internal helper functions for RUV harmonization
.harmon.method.RUVg <- function(raw, groups) {
  condition <- factor(groups)

  # Create expression set and design matrix
  set <- newSeqExpressionSet(
    as.matrix(raw),
    phenoData = data.frame(condition, row.names = colnames(raw))
  )
  design <- model.matrix(~condition,
    data = data.frame(condition, row.names = colnames(raw))
  )

  # Prepare DGE analysis
  y <- DGEList(counts = DESeq2::counts(set), group = condition) %>%
    calcNormFactors(method = "upperquartile")

  if (any(is.infinite(y$samples$norm.factors))) {
    y <- DGEList(counts = DESeq2::counts(set), group = condition) %>%
      calcNormFactors(method = "none")
  }

  # Estimate dispersions and fit model
  y <- y %>%
    estimateGLMCommonDisp(design) %>%
    estimateGLMTagwiseDisp(design)

  # Identify control genes
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15 * nrow(raw))]))]

  # Apply RUVg harmonization
  t <- RUVg(set, spikes, k = 1)

  list(
    dat.harmonized = normCounts(t),
    adjust.factor = t$W
  )
}

.harmon.method.RUVs <- function(raw, groups) {
  condition <- factor(groups)

  # Create expression set and design matrix
  set <- newSeqExpressionSet(
    as.matrix(raw),
    phenoData = data.frame(condition, row.names = colnames(raw))
  )
  design <- model.matrix(~condition,
    data = data.frame(condition, row.names = colnames(raw))
  )

  # Prepare DGE analysis
  y <- DGEList(counts = DESeq2::counts(set), group = condition) %>%
    calcNormFactors(method = "upperquartile")

  if (any(is.infinite(y$samples$norm.factors))) {
    y <- DGEList(counts = DESeq2::counts(set), group = condition) %>%
      calcNormFactors(method = "none")
  }

  # Estimate dispersions and fit model
  y <- y %>%
    estimateGLMCommonDisp(design) %>%
    estimateGLMTagwiseDisp(design)

  # Identify control genes and differences
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15 * nrow(raw))]))]
  differences <- makeGroups(condition)

  # Apply RUVs harmonization
  t <- RUVs(set, rownames(raw), k = 1, differences)

  list(
    dat.harmonized = normCounts(t),
    adjust.factor = t$W
  )
}

.harmon.method.RUVr <- function(raw, groups) {
  condition <- factor(groups)

  # Create expression set and design matrix
  set <- newSeqExpressionSet(
    as.matrix(raw),
    phenoData = data.frame(condition, row.names = colnames(raw))
  )
  design <- model.matrix(~condition, data = pData(set))

  # Prepare DGE analysis
  y <- DGEList(counts = DESeq2::counts(set), group = condition) %>%
    calcNormFactors(method = "upperquartile")

  if (any(is.infinite(y$samples$norm.factors))) {
    y <- DGEList(counts = DESeq2::counts(set), group = condition) %>%
      calcNormFactors(method = "none")
  }

  # Estimate dispersions and get residuals
  y <- y %>%
    estimateGLMCommonDisp(design) %>%
    estimateGLMTagwiseDisp(design)

  fit <- glmFit(y, design)
  residuals <- residuals(fit, type = "deviance")

  # Apply RUVr harmonization
  set <- betweenLaneNormalization(set, which = "upper")
  t <- RUVr(set, rownames(raw), k = 1, residuals = residuals)

  list(
    dat.harmonized = normCounts(t),
    adjust.factor = t$W
  )
}

#' Apply RUVg harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with RUVg harmonization results added
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

#' Apply RUVs harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with RUVs harmonization results added
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

#' Apply RUVr harmonization
#'
#' @param object A precision object containing raw data
#' @return A precision object with RUVr harmonization results added
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

# Internal helper function for ComBat-Seq normalization
.harmon.method.ComBat.Seq <- function(raw, batches) {
  # Apply ComBat-Seq normalization
  dat.harmonized <- sva::ComBat_seq(raw, batch = batches)

  list(
    dat.harmonized = dat.harmonized
  )
}

#' Apply ComBat-Seq normalization
#'
#' @param object A precision object containing raw data
#' @param batches A vector of batch labels for the samples
#' @return A precision object with ComBat-Seq normalization results added
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
