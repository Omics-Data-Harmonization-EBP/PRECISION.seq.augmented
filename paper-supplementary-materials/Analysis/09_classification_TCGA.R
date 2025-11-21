suppressPackageStartupMessages({
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
})

source("RUV_without_effects.R") 
load("brca_LumA_vs_LumB_trainEarlyAB_testLate_with_reps.RData")  
stopifnot(is.list(datasets) && length(datasets) > 0)
harmon.all <- function(object){
  suppressMessages({
    suppressWarnings({
      object@harmon.train.data$Raw$dat.harmonized  <- object@raw.train.data$data
      object@harmon.test1.data$Raw$dat.harmonized  <- object@raw.test1.data$data
      object@harmon.test2.data$Raw$dat.harmonized  <- object@raw.test2.data$data
      
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
.preprocess.classification.data <- function(data) {
  lapply(data, function(x) {
    if (!any(x$dat.harmonized < 0)) log2(x$dat.harmonized + 1) else x$dat.harmonized
  })
}

.safe.exec <- function(expr) {
  tryCatch(expr, error = function(e) { message("Error: ", conditionMessage(e)); NULL })
}

.calc.error <- function(predicted, observed) {
  predicted <- as.character(predicted)
  observed  <- as.character(observed)
  mean(predicted != observed)
}

classification.lasso <- function(object, threshold_method = "cv", kfold = 5) {
  threshold_method <- match.arg(threshold_method, c("cv", "none"))
  
  processed_train_data <- .preprocess.classification.data(object@harmon.train.data)
  processed_test1_data <- .preprocess.classification.data(object@harmon.test1.data)
  processed_test2_data <- .preprocess.classification.data(object@harmon.test2.data)
  
  Xtr_list <- lapply(processed_train_data, t)  # samples x features
  Xte1_list <- lapply(processed_test1_data, t)
  Xte2_list <- lapply(processed_test2_data, t)
  
  ytr <- factor(object@raw.train.data$label)
  yte1 <- factor(object@raw.test1.data$label, levels = levels(ytr))
  yte2 <- factor(object@raw.test2.data$label, levels = levels(ytr))
  
  object@classification.result$lasso[[threshold_method]] <-
    lapply(seq_along(Xtr_list), function(i) {
      mdl <- .train.lasso(Xtr_list[[i]], ytr, threshold_method, kfold = kfold)
      p1  <- .predict.lasso(mdl, Xte1_list[[i]], yte1)
      p2  <- .predict.lasso(mdl, Xte2_list[[i]], yte2)
      list(train_class = mdl$pred, test1_class = p1$pred, test2_class = p2$pred)
    })
  
  object
}

.train.lasso <- function(x, y, threshold_method = "cv", kfold = 5, alpha = 1) {
  set.seed(42)
  if (threshold_method == "cv") {
    cv.fit <- glmnet::cv.glmnet(
      x = x, y = y, family = "multinomial",
      type.measure = "class", alpha = alpha, nfolds = kfold
    )
    lambda_value <- cv.fit$lambda.min
  } else {
    lambda_value <- 1e-6
    cv.fit <- glmnet::glmnet(x = x, y = y, family = "multinomial", alpha = alpha, lambda = lambda_value)
  }
  pred <- predict(cv.fit, newx = x, s = lambda_value, type = "class")
  if (is.matrix(pred)) pred <- drop(pred)
  list(pred = factor(pred, levels = levels(y)),
       mc = mean(pred != y),
       model = cv.fit,
       lambda = lambda_value,
       coefs = .safe.exec(coef(cv.fit, s = lambda_value)))
}

.predict.lasso <- function(model, x, y) {
  pred <- predict(model$model, newx = x, s = model$lambda, type = "class")
  if (is.matrix(pred)) pred <- drop(pred)
  prob <- predict(model$model, newx = x, s = model$lambda)
  list(pred = factor(pred, levels = levels(y)), mc = .calc.error(pred, y), prob = prob)
}

.has_two_classes <- function(y) length(unique(y)) >= 2

.stratified_random_split <- function(y, train_frac = 0.75, seed = 42, min_per_set = 1) {
  set.seed(seed)
  y <- factor(y)
  tr_idx <- integer(0); te_idx <- integer(0)
  for (lev in levels(y)) {
    ids <- which(y == lev)
    n   <- length(ids)
    ntr <- max(min_per_set, floor(train_frac * n))
    ntr <- min(n - min_per_set, ntr)  # leave ≥1 for test
    if (ntr < min_per_set) ntr <- max(1, n - min_per_set)
    tr  <- sample(ids, ntr)
    te  <- setdiff(ids, tr)
    tr_idx <- c(tr_idx, tr)
    te_idx <- c(te_idx, te)
  }
  list(train = sort(tr_idx), test = sort(te_idx), strategy = "random_stratified")
}

.pick_split <- function(meta, prefer_era = TRUE) {
  stopifnot(all(c("era","Labels") %in% colnames(meta)))
  if (prefer_era) {
    eras <- as.character(meta$era)
    train_idx <- which(eras %in% c("Early_partA","Early_partB"))
    test_idx  <- which(eras == "Late")
    
    ok_counts  <- (length(train_idx) > 0) && (length(test_idx) > 0)
    ok_classes <- .has_two_classes(meta$Labels[train_idx])
    
    if (ok_counts && ok_classes) {
      return(list(train = train_idx, test = test_idx, strategy = "era_EarlyAB_to_Late"))
    } else {
      message("Era‑based split not feasible (counts/classes). Falling back to stratified random.")
    }
  }
  .stratified_random_split(meta$Labels)
}

.run_one_dataset <- function(ds_name, ds, threshold_method = "cv") {
  expr <- as.matrix(ds$expr)
  meta <- ds$meta
  stopifnot(is.matrix(expr), ncol(expr) == nrow(meta))
  rownames(meta) <- colnames(expr)  
  keep <- rowSums(expr > 0) >= ceiling(0.01 * ncol(expr))
  expr <- expr[keep, , drop = FALSE]
  
  split <- .pick_split(meta, prefer_era = TRUE)
  tr <- split$train; te <- split$test
  ytr <- factor(meta$Labels[tr])
  yte <- factor(meta$Labels[te], levels = levels(ytr))
  
  analysis <- create.precision.classification(
    traindata  = expr[, tr, drop = FALSE],
    testdata1  = expr[, te, drop = FALSE],
    testdata2  = expr[, te, drop = FALSE],
    trainlabel = as.character(ytr),
    testlabel1 = as.character(yte),
    testlabel2 = as.character(yte)
  )
  
  analysis <- harmon.all(analysis)
  suppressMessages(suppressWarnings({
    analysis <- classification.lasso(analysis, threshold_method = threshold_method)
    analysis <- classification.pam(analysis,   threshold_method = threshold_method)
    analysis <- classification.knn(analysis,   threshold_method = threshold_method)
    analysis <- classification.svm(analysis,   threshold_method = threshold_method)
    analysis <- classification.svm(analysis, threshold_method = threshold_method, feature_selection = TRUE)
    analysis <- classification.ranfor(analysis, threshold_method = threshold_method)
    analysis <- classification.ranfor(analysis, threshold_method = threshold_method, feature_selection = TRUE)
  }))
  H <- length(analysis@classification.result$lasso[[threshold_method]])
  base_methods <- c("lasso","pam","knn","svm","svm_with_fs", "ranfor","ranfor_with_fs")
  predictions <- list()
  for (m in base_methods) {
    predictions[[m]] <- vector("list", H)
    for (h in seq_len(H)) {
      slot <- analysis@classification.result[[m]][[threshold_method]][[h]]
      if (is.null(slot)) next
      predictions[[m]][[h]] <- list(
        train = slot$train_class,
        test1 = slot$test1_class,
        test2 = slot$test2_class
      )
    }
  }
  
  vec_len <- function(x) if (is.null(x)) NA_integer_ else length(x)
  metrics <- list()
  parts <- strsplit(ds_name, "__", fixed = TRUE)[[1]]
  scen  <- parts[1]; pair <- parts[2]; repid <- parts[3]
  
  for (m in base_methods) {
    for (h in seq_len(H)) {
      pr <- predictions[[m]][[h]]
      if (is.null(pr)) next
      tr_pred <- factor(as.character(pr$train), levels = levels(ytr))
      te_pred <- factor(as.character(pr$test1), levels = levels(ytr))
      
      metrics[[length(metrics)+1L]] <- data.frame(
        dataset      = ds_name,
        scenario     = scen,
        pair         = pair,
        replicate    = repid,
        split        = split$strategy,
        method       = m,
        harmon_id    = h,
        N_train      = length(tr),
        N_test       = length(te),
        err_train    = .calc.error(tr_pred, ytr),
        err_test1    = .calc.error(te_pred, yte),
        stringsAsFactors = FALSE
      )
    }
  }
  metrics <- dplyr::bind_rows(metrics)
  
  list(
    name        = ds_name,
    split       = split,
    predictions = predictions,
    metrics     = metrics
  )
}

threshold_method <- "cv"   # or "none"

ds_names <- names(datasets)

N_CORES <- as.integer(Sys.getenv("MC_CORES", parallel::detectCores() - 1L))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

set.seed(42)  

runner <- function(nm) {
  message("Running: ", nm)
  tryCatch(
    .run_one_dataset(nm, datasets[[nm]], threshold_method = threshold_method),
    error = function(e) {
      message("Error in ", nm, ": ", conditionMessage(e))
      NULL
    }
  )
}

if (.Platform$OS.type == "windows") {
  warning("mclapply not supported on Windows; running sequentially.")
  results <- lapply(ds_names, runner)
} else {
  results <- parallel::mclapply(
    ds_names,
    runner,
    mc.cores      = N_CORES,
    mc.set.seed   = TRUE,   
    mc.preschedule = FALSE  
  )
}

names(results) <- ds_names
keep <- !vapply(results, is.null, logical(1))
results <- results[keep]

all_metrics <- dplyr::bind_rows(lapply(results, `[[`, "metrics"))
predictions_by_dataset <- lapply(results, `[[`, "predictions")

save(results, all_metrics, predictions_by_dataset,
     file = "classification_results_simulated_all_v5.RData")
