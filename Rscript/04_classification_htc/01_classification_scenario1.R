library(parallel)
library(sva)
library(EDASeq)
library(edgeR)
library(RUVSeq)
library(Biobase)
library(BiocGenerics)
library(tidyverse)
library(mclust)
library(aricode)
library(RSKC)
library(cluster)
library(factoextra)
library(caret)
library(e1071)
library(PRECISION.seq.augmented)
load("MSKpair_300_classification.RData")

harmon.all <- function(object) {
  object@harmon.train.data$Raw$dat.harmonized <- object@raw.train.data$data
  object@harmon.test1.data$Raw$dat.harmonized <- object@raw.test1.data$data
  object@harmon.test2.data$Raw$dat.harmonized <- object@raw.test2.data$data
  object <- harmon.TC(object)
  object <- harmon.UQ(object)
  object <- harmon.med(object)
  object <- harmon.TMM(object)
  object <- harmon.DESeq(object)
  object <- harmon.PoissonSeq(object)
  object <- harmon.QN(object)
  object <- harmon.RUVr(object)
  object <- harmon.RUVs(object)
  object <- harmon.RUVg(object)
  return(object)
}


dirty.scenario <- function(c, d) {
  clean.datasets <- lapply(1:300, function(x) {
    biological.effects(
      benchmark_subsets[[x]]$data,
      benchmark_subsets[[x]]$label, c
    )
  })
  dirty.datasets <- lapply(1:300, function(x) {
    handling.effects(
      clean.datasets[[x]]$data, benchmark_subsets[[x]]$data, test_subsets[[x]]$data,
      clean.datasets[[x]]$group, d
    )
  })


  train.dataset <- lapply(dirty.datasets, function(x) {
    list(
      data = x$data[, c(1:100, 201:300)],
      group = x$group[c(1:100, 201:300)]
    )
  })
  test1.dataset <- lapply(dirty.datasets, function(x) {
    list(
      data = x$data[, c(101:200, 301:400)],
      group = x$group[c(101:200, 301:400)]
    )
  })
  test2.dataset <- lapply(clean.datasets, function(x) {
    list(
      data = x$data[, c(101:200, 301:400)],
      group = x$group[c(101:200, 301:400)]
    )
  })


  return(list(
    train.dataset = train.dataset,
    test1.dataset = test1.dataset,
    test2.dataset = test2.dataset
  ))
}

measure <- function(k, datalist, threshold_method) {
  analysis <- create.precision.classification(
    traindata = datalist$train.dataset[[k]]$data,
    testdata1 = datalist$test1.dataset[[k]]$data,
    testdata2 = datalist$test2.dataset[[k]]$data,
    trainlabel = datalist$train.dataset[[k]]$group,
    testlabel1 = datalist$test1.dataset[[k]]$group,
    testlabel2 = datalist$test2.dataset[[k]]$group
  )
  analysis <- harmon.all(analysis)

  # Run only the specified threshold method version
  suppressMessages({
    suppressWarnings({
      analysis <- classification.pam(analysis, threshold_method = threshold_method)
      analysis <- classification.knn(analysis, threshold_method = threshold_method)
      analysis <- classification.lasso(analysis, threshold_method = threshold_method)
      analysis <- classification.svm(analysis, threshold_method = threshold_method)
      analysis <- classification.ranfor(analysis, threshold_method = threshold_method)
    })
  })

  # Define base classification methods and harmonization methods
  base_methods <- c("pam", "knn", "lasso", "svm", "ranfor")
  harmon <- c(
    "Raw", "TC", "UQ", "med", "TMM", "DESeq", "PoissonSeq", "QN",
    "RUVr", "RUVs", "RUVg"
  )

  # Function to create and initialize a data frame with proper names
  init_df <- function() {
    df <- data.frame(matrix(nrow = length(base_methods), ncol = length(harmon)))
    colnames(df) <- harmon
    rownames(df) <- base_methods
    return(df)
  }

  # Initialize metrics for the specific threshold method
  metrics <- list(
    train = list(
      accuracy = init_df(),
      sensitivity = init_df(),
      specificity = init_df()
    ),
    test1 = list(
      accuracy = init_df(),
      sensitivity = init_df(),
      specificity = init_df()
    ),
    test2 = list(
      accuracy = init_df(),
      sensitivity = init_df(),
      specificity = init_df()
    )
  )

  # Calculate metrics
  for (i in 1:length(base_methods)) {
    for (j in 1:length(harmon)) {
      group.levels <- levels(factor(datalist$train.dataset[[k]]$group))
      # Get predicted labels
      pred_train <- factor(
        analysis@classification.result[[base_methods[i]]][[threshold_method]][[j]]$train_class,
        levels = group.levels
      )
      pred_test1 <- factor(
        analysis@classification.result[[base_methods[i]]][[threshold_method]][[j]]$test1_class,
        levels = group.levels
      )
      pred_test2 <- factor(
        analysis@classification.result[[base_methods[i]]][[threshold_method]][[j]]$test2_class,
        levels = group.levels
      )

      # Calculate confusion matrices
      cm_train <- table(
        Predicted = pred_train,
        Actual = datalist$train.dataset[[k]]$group
      )
      cm_test1 <- table(
        Predicted = pred_test1,
        Actual = datalist$test1.dataset[[k]]$group
      )
      cm_test2 <- table(
        Predicted = pred_test2,
        Actual = datalist$test2.dataset[[k]]$group
      )

      # Calculate metrics
      metrics$train$accuracy[i, j] <- sum(diag(cm_train)) / sum(cm_train)
      metrics$train$sensitivity[i, j] <- cm_train[2, 2] / sum(cm_train[, 2])
      metrics$train$specificity[i, j] <- cm_train[1, 1] / sum(cm_train[, 1])

      metrics$test1$accuracy[i, j] <- sum(diag(cm_test1)) / sum(cm_test1)
      metrics$test1$sensitivity[i, j] <- cm_test1[2, 2] / sum(cm_test1[, 2])
      metrics$test1$specificity[i, j] <- cm_test1[1, 1] / sum(cm_test1[, 1])

      metrics$test2$accuracy[i, j] <- sum(diag(cm_test2)) / sum(cm_test2)
      metrics$test2$sensitivity[i, j] <- cm_test2[2, 2] / sum(cm_test2[, 2])
      metrics$test2$specificity[i, j] <- cm_test2[1, 1] / sum(cm_test2[, 1])
    }
  }

  return(metrics)
}

dirty.data.analysis <- function(c, d) {
  simulated_data <- dirty.scenario(c, d)

  # First run "none" threshold method
  result_none <- mclapply(1:300, function(x) measure(k = x, simulated_data, "none"), mc.cores = 50)
  gc()

  # Then run "cv" threshold method
  result_cv <- mclapply(1:300, function(x) measure(k = x, simulated_data, "cv"), mc.cores = 50)
  gc()

  # Combine results
  result <- list(
    none = result_none,
    cv = result_cv
  )

  return(list(
    result = result,
    simulated_data = simulated_data
  ))
}

params <- list(
  c(-0.2, 1.2),
  c(-0.2, 1.5),
  c(-0.2, 2),
  c(0.2, 1.2),
  c(0.2, 1.5),
  c(0.2, 2),
  c(0.4, 1.2),
  c(0.4, 1.5),
  c(0.4, 2)

  # c(0.6, 1.2),
  # c(0.6, 1.5),
  # c(0.6, 2),
  #
  # c(1, 1.2),
  # c(1, 1.5),
  # c(1, 2),
  #
  # c(1.2, 1.2),
  # c(1.2, 1.5),
  # c(1.2, 2)
)


dirty_results_list <- list()
for (param in params) {
  key <- paste0("dirty_res_c", param[1], "d", param[2])
  dirty_results_list[[key]] <- dirty.data.analysis(param[1], param[2])
}

save(dirty_results_list, file = "classfication_sceanrio1_part1.RData")
