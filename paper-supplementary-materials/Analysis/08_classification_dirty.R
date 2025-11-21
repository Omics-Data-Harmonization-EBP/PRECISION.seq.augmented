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
library(penalizedSVM)
load("MSKpair_300_classification.RData")

# Define harmonization methods
harmon_methods <- c("Raw", "TC", "UQ", "med", "TMM", "TMM.frozen", "DESeq", "PoissonSeq", 
                    "QN", "QN.frozen", "RUVg", "RUVr", "RUVs")

# Define base classification methods and threshold methods
base_methods <- c("lasso", "pam", "knn",  "svm", "svm_fs",  "ranfor", "ranfor_fs")
threshold_methods <- c("none", "cv")

harmon.all <- function(object){
  suppressMessages({
    suppressWarnings({
      object@harmon.train.data$Raw$dat.harmonized <- object@raw.train.data$data
      object@harmon.test1.data$Raw$dat.harmonized <- object@raw.test1.data$data
      object@harmon.test2.data$Raw$dat.harmonized <- object@raw.test2.data$data
      object <- harmon.TC(object)
      object <- harmon.UQ(object)
      object <- harmon.med(object)
      object <- harmon.TMM(object)
      object <- harmon.TMM.frozen(object)
      object <- harmon.DESeq(object)
      object <- harmon.PoissonSeq(object)
      object <- harmon.QN(object)
      object <- harmon.QN.frozen(object)
      object <- harmon.RUVg(object)
      object <- harmon.RUVr(object)
      object <- harmon.RUVs(object)
    })
  })
  return(object)
}

measure <- function(k, datalist, threshold_method){
  analysis <- create.precision.classification(traindata = datalist$train.dataset[[k]]$data,
                                              testdata1 = datalist$test1.dataset[[k]]$data,
                                              testdata2 = datalist$test2.dataset[[k]]$data,
                                              trainlabel = datalist$train.dataset[[k]]$group,
                                              testlabel1 = datalist$test1.dataset[[k]]$group,
                                              testlabel2 = datalist$test2.dataset[[k]]$group)
  analysis <- harmon.all(analysis)
  
  # Run only the specified threshold method version
  suppressMessages({
    suppressWarnings({
      analysis <- classification.lasso(analysis, threshold_method = threshold_method)
      analysis <- classification.pam(analysis, threshold_method = threshold_method)
      analysis <- classification.knn(analysis, threshold_method = threshold_method)
      analysis <- classification.svm(analysis, threshold_method = threshold_method)
      analysis <- classification.svm(analysis, threshold_method = threshold_method, feature_selection = TRUE)
      analysis <- classification.ranfor(analysis, threshold_method = threshold_method)
      analysis <- classification.ranfor(analysis, threshold_method = threshold_method, feature_selection = TRUE)
    })
  })
  
  # Define base classification methods and harmonization methods
  base_methods <-  c("lasso", "pam", "knn",  "svm", "svm_with_fs",  "ranfor", "ranfor_with_fs")
  harmon <- 1:13
  
  # Store only predicted classes
  predictions <- list()
  
  for (method in base_methods) {
    predictions[[method]] <- list()
    for (harm in harmon) {
      predictions[[method]][[harm]] <- list(
        train = analysis@classification.result[[method]][[threshold_method]][[harm]]$train_class,
        test1 = analysis@classification.result[[method]][[threshold_method]][[harm]]$test1_class,
        test2 = analysis@classification.result[[method]][[threshold_method]][[harm]]$test2_class
      )
    }
  }
  
  return(predictions)
}

dirty.scenario <- function(c, d){
  clean.datasets <- lapply(1:300, function(x) biological.effects(benchmark_subsets[[x]]$data,
                                                                 benchmark_subsets[[x]]$label, c))
  dirty.datasets <- lapply(1:300, function (x)  handling.effects(clean.datasets[[x]]$data, benchmark_subsets[[x]]$data, test_subsets[[x]]$data,
                                                                 clean.datasets[[x]]$group, d))
  
  
  train.dataset <- lapply(dirty.datasets, function(x) list(data = x$data[, c(1:100, 201:300)],
                                                           group = x$group[c(1:100, 201:300)]))
  test1.dataset <- lapply(dirty.datasets, function(x) list(data = x$data[, c(101:200, 301:400)],
                                                           group = x$group[c(101:200, 301:400)]))
  test2.dataset <- lapply(clean.datasets, function(x) list(data = x$data[, c(101:200, 301:400)],
                                                           group = x$group[c(101:200, 301:400)]))
  
  
  return(list(train.dataset = train.dataset,
              test1.dataset = test1.dataset,
              test2.dataset = test2.dataset))
}


dirty.data.analysis <- function(c, d){
  simulated_data <- dirty.scenario(c, d)
  
  # First run "none" threshold method
  # result_none <- mclapply(1:300, function(x) measure(k = x, simulated_data, "none"), mc.cores = 50)
  # gc()
  
  # Then run "cv" threshold method
  result_cv <- mclapply(1:300, function(x) measure(k = x, simulated_data, "cv"), mc.cores = 70)
  gc()
  
  # Combine results
  result <- list(
  #  none = result_none,
    cv = result_cv
  )
  
  return(list(result = result, 
              simulated_data = simulated_data))
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
  c(0.4, 2),
  
  c(0.6, 1.2),
  c(0.6, 1.5),
  c(0.6, 2),
  
  c(1, 1.2),
  c(1, 1.5),
  c(1, 2),
  
  c(1.2, 1.2),
  c(1.2, 1.5),
  c(1.2, 2)
)

# Process one parameter set at a time
# for (i in seq_along(params)) {
#   param <- params[[i]]
# 
#   # Create filename based on parameters
#   filename <- sprintf("classification_scenario2_c%.1f_d%.1f.RData", param[1], param[2])
# 
#   # Process single parameter set
#   message(sprintf("Processing c=%.1f, d=%.1f", param[1], param[2]))
#   result <- dirty.data.analysis(param[1], param[2])
# 
#   # Save result immediately and clear memory
#   save(result, file = filename)
#   rm(result)
#   gc()
# 
#   message(sprintf("Completed and saved results for c=%.1f, d=%.1f", param[1], param[2]))
# }

for (i in seq_along(params)) {
  param <- params[[i]]
  
  # Create filename based on parameters
  filename <- sprintf("classification_scenario2_c%.1f_d%.1f.RData", param[1], param[2])
  
  # Process single parameter set
  message(sprintf("Loading c=%.1f, d=%.1f", param[1], param[2]))
  load(filename)
  null_detect <- which(sapply(result$result$cv, is.null)|sapply(result$result$cv, length) !=  7)
  result$result$cv[null_detect] <- mclapply(null_detect,
                                            function(x) measure(k = x, result$simulated_data, "cv"), mc.cores = 70)
  
  # Save result immediately and clear memory
  save(result, file = filename)
  rm(result)
  gc()
  
  message(sprintf("Completed and saved results for c=%.1f, d=%.1f", param[1], param[2]))
}
