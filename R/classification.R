#' @import e1071
#' @import caret
#' @import ranger
#' @import pamr
#' @import glmnet
#' @importFrom stats predict coef
#' @importFrom utils data
NULL

### Internal utility functions


#' Try-catch wrapper for error handling
#' @keywords internal
#' @noRd
trycatch.func <- function(expr, msg = "") {
  out <- tryCatch({
    expr
  }, warning = function(cond) {
    message("warning: ")
    message(cond)
    # Choose a return value in case of warning
    return(NULL)
  }, error = function(cond) {
    message("error: ")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  }, finally = {
    message(msg)
  })
  return(out)
}

#' Calculate external error rate
#' @param pred.obj Predicted values
#' @param obs.grp Observed groups
#' @return Error rate
#' @noRd
tabulate.ext.err.func <- function(pred.obj, obs.grp) {
  return(1.0 - sum(diag(table(pred.obj, obs.grp))) / length(obs.grp))
}

#' Modified log2 transformation for harmonized data
#' @param x Data object containing harmonized data
#' @return Transformed data object
#' @noRd
.transform_data <- function(x) {
  if (!any(x$dat.harmonized < 0)) {
    log2(x$dat.harmonized + 1)
  } else {
    x$dat.harmonized
  }
}

##########################################################
# PAM
##########################################################

#' PAM Classification
#'
#' This function applies classification using Prediction Analysis for Microarrays (PAM)
#' on the provided input object. It supports threshold optimization using
#' cross-validation or the use of all genes without thresholding.
#'
#' @param object A \link{precision} object containing harmonized data.
#'  Must contain the slots \code{harmon.train.data} with harmonized training data
#'  and \code{harmon.test1.data} and \code{harmon.test2.data} with harmonized test data.
#' @param threshold_method A character string specifying the thresholding method.
#'   Options are:
#'   - \code{"cv"}: Use cross-validation to optimize the threshold.
#'   - \code{"none"}: Use all genes without applying a threshold.
#' @param vt.k A numeric vector of threshold values to evaluate during cross-validation.
#'   Only used if \code{threshold_method = "cv"}.
#' @param kfold An integer specifying the number of folds for cross-validation.
#' @return The input object updated with PAM classification results added to the
#'  \code{classification.result} slot, including predicted classes and
#'  associated metrics.
#' @import pamr
#' @export
classification.pam <- function(object, threshold_method = "cv", vt.k = NULL, kfold = 5) {
  # Validate threshold method
  threshold_method <- match.arg(threshold_method, c("cv", "none"))

  # Prepare datasets and labels
  datasets <- list(
    data = list(
      train = lapply(object@harmon.train.data, .transform_data),
      test1 = lapply(object@harmon.test1.data, .transform_data),
      test2 = lapply(object@harmon.test2.data, .transform_data)
    ),
    labels = list(
      train = object@raw.train.data$label,
      test1 = object@raw.test1.data$label,
      test2 = object@raw.test2.data$label
    )
  )

  # Run PAM classification for each harmonization method
  object@classification.result$pam[[threshold_method]] <- lapply(seq_along(datasets$data$train), function(i) {
    # Train model
    model <- pam.intcv(
      X = datasets$data$train[[i]],
      y = datasets$labels$train,
      threshold_method = threshold_method,
      vt.k = vt.k,
      kfold = kfold
    )

    # Get predictions for both test sets
    pred1 <- pam.predict(model, datasets$data$test1[[i]], datasets$labels$test1)
    pred2 <- pam.predict(model, datasets$data$test2[[i]], datasets$labels$test2)

    # Return metrics
    list(train_class = model$pred, test1_class = pred1$pred, test2_class = pred2$pred)
  })

  return(object)
}

#' Internal function for PAM cross-validation
#' @param X Training data matrix
#' @param y Labels for training data
#' @param threshold_method Character. Either "cv" for cross-validation optimized threshold or "none" for using all genes
#' @param vt.k Threshold values (only used if threshold_method = "cv")
#' @param kfold Number of folds for cross-validation
#' @import pamr
#' @noRd
pam.intcv <- function(X, y, threshold_method = "cv", vt.k = NULL, kfold = 5) {
  set.seed(42)
  start_time <- proc.time()

  data.pam <- list(
    x = as.matrix(X), y = factor(y),
    geneids = rownames(X), genenames = rownames(X)
  )

  if (threshold_method == "cv") {
    # Train initial model with a sequence of thresholds
    raw_model <- pamr::pamr.train(data.pam, threshold = vt.k)

    # Perform cross-validation to find best threshold
    cv_model <- pamr::pamr.cv(
      fit = raw_model,
      data = data.pam,
      nfold = kfold
    )

    # Find best threshold (one with minimum CV error)
    best.threshold <- cv_model$threshold[which.min(cv_model$error)]

    # Train final model with best threshold
    final_model <- pamr::pamr.train(data.pam, threshold = best.threshold)

    # Calculate error on training set
    pred <- pamr::pamr.predict(final_model, X, threshold = best.threshold, type = "class")
    mc <- mean(pred != y)
  } else {
    # Use all genes (no threshold)
    final_model <- pamr::pamr.train(data.pam, threshold = 0)

    # Calculate error on training set
    pred <- pamr::pamr.predict(final_model, X, threshold = 0, type = "class")
    mc <- mean(pred != y)
  }

  list(
    pred = pred,
    mc = mc,
    time = proc.time() - start_time,
    model = final_model,
    cfs = if (threshold_method == "cv") {
      trycatch.func(pamr::pamr.listgenes(final_model, data.pam, threshold = best.threshold))
    } else {
      NULL
    }
  )
}

#' Prediction function for PAM model
#' @param pam.intcv.model Trained PAM model
#' @param pred.obj New data matrix for prediction
#' @param pred.obj.group.id Group IDs for the new data
#' @return List containing predicted classes and misclassification counts
#' @noRd
#' @importFrom pamr pamr.predict
pam.predict <- function(pam.intcv.model, pred.obj, pred.obj.group.id) {
  # Use threshold = 0 if model was trained without threshold selection
  threshold <- if (length(pam.intcv.model$model$threshold) == 1) 0 else pam.intcv.model$model$threshold

  pred <- pamr::pamr.predict(
    pam.intcv.model$model,
    newx = pred.obj,
    threshold = threshold,
    type = "class"
  )

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)

  prob <- pamr::pamr.predict(
    pam.intcv.model$model,
    newx = pred.obj,
    threshold = threshold,
    type = "posterior"
  )

  return(list(pred = pred, mc = mc, prob = prob))
}


##########################################################
# KNN
##########################################################

#' K-Nearest Neighbors (KNN) Classification
#'
#' This function applies KNN classification to the harmonized data in the
#' input \link{precision} object containing.
#' It supports two thresholding methods: cross-validation to optimize the number of neighbors (k)
#' or a fixed k value of 1.
#'
#' @param object A \link{precision} object containing harmonized data.
#'  Must contain the slots \code{harmon.train.data} with harmonized training data
#'  and \code{harmon.test1.data} and \code{harmon.test2.data} with harmonized test data.
#' @param threshold_method A character string specifying the thresholding method.
#'  Use \code{"cv"} for cross-validation to determine the optimal k,
#'  or \code{"none"} to use a fixed k = 1.
#' @param kfold An integer specifying the number of folds for cross-validation.
#'   This parameter is only used if \code{threshold_method} is set to \code{"cv"}.
#' @param folds An optional list of pre-specified folds for cross-validation.
#'   If provided, these folds will be used instead of generating new ones.
#' @return The input object updated with KNN classification results added to the
#'  \code{classification.result} slot, including predicted classes and
#'  associated metrics.
#' @export
classification.knn <- function(object, threshold_method = "cv", kfold = 5, folds = NULL) {
  # Validate threshold method
  threshold_method <- match.arg(threshold_method, c("cv", "none"))

  # Prepare datasets and labels
  datasets <- list(
    data = list(
      train = lapply(object@harmon.train.data, .transform_data),
      test1 = lapply(object@harmon.test1.data, .transform_data),
      test2 = lapply(object@harmon.test2.data, .transform_data)
    ),
    labels = list(
      train = object@raw.train.data$label,
      test1 = object@raw.test1.data$label,
      test2 = object@raw.test2.data$label
    )
  )

  # Run KNN classification for each harmonization method
  object@classification.result$knn[[threshold_method]] <- lapply(seq_along(datasets$data$train), function(i) {
    # Train model with transposed data
    model <- knn.intcv(
      kfold = kfold,
      X = t(datasets$data$train[[i]]),
      y = datasets$labels$train,
      threshold_method = threshold_method
    )

    # Get predictions for both test sets
    pred1 <- knn.predict(model, t(datasets$data$test1[[i]]), datasets$labels$test1)
    pred2 <- knn.predict(model, t(datasets$data$test2[[i]]), datasets$labels$test2)

    # Return metrics
    list(train_class = model$pred, test1_class = pred1$pred, test2_class = pred2$pred)
  })

  return(object)
}

#' Internal function for KNN classification
#' @param kfold Number of folds for cross-validation
#' @param X Training data matrix
#' @param y Labels for training data
#' @param threshold_method Character. Either "cv" for cross-validation optimized
#' threshold or "none" for using all genes
#' @import caret
#' @noRd
knn.intcv <- function(kfold = 5, X, y, threshold_method = "cv") {
  start_time <- proc.time()
  set.seed(42)

  if (threshold_method == "cv") {
    # Use cross-validation to find optimal k
    ctrl <- caret::trainControl(method = "repeatedcv", repeats = 2, number = kfold)
    knn_model <- caret::train(
      x = data.matrix(X), y = factor(y),
      method = "knn",
      tuneLength = 5, # Try 5 different k values
      trControl = ctrl
    )
    pred <- predict(knn_model, newdata = data.matrix(X))
    mc <- 1 - max(knn_model$results$Accuracy)
  } else {
    # Use fixed k=1 without cross-validation
    ctrl <- caret::trainControl(method = "none")
    knn_model <- caret::train(
      x = data.matrix(X), y = factor(y),
      method = "knn",
      tuneGrid = data.frame(k = 1),
      trControl = ctrl
    )

    # Calculate error on training set
    pred <- predict(knn_model, newdata = data.matrix(X))
    mc <- mean(pred != y)
  }

  list(
    pred = pred,
    mc = mc,
    time = proc.time() - start_time,
    model = knn_model,
    cfs = NULL
  )
}

#' Prediction function for KNN model
#' @param knn.intcv.model Trained KNN model
#' @param pred.obj New data matrix for prediction
#' @param pred.obj.group.id Group IDs for the new data
#' @return List containing predicted classes and misclassification counts
#' @noRd
knn.predict <- function(knn.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(knn.intcv.model$model, newdata = data.matrix(pred.obj))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  list(pred = pred, mc = mc)
}

##########################################################
# SVM
##########################################################

#' Support Vector Machine (SVM) Classification
#'
#' This function applies SVM classification to the harmonized data in the
#' input \link{precision} object containing.
#' It supports two thresholding methods: cross-validation to optimize the cost parameter
#' or a fixed cost value of 1.
#'
#' @param object A \link{precision} object containing harmonized data.
#'  Must contain the slots \code{harmon.train.data} with harmonized training data
#'  and \code{harmon.test1.data} and \code{harmon.test2.data} with harmonized test data.
#' @param threshold_method A character string specifying the thresholding method.
#'  Use \code{"cv"} for cross-validation to determine the optimal k,
#'  or \code{"none"} to use a fixed k = 1.
#' @param kfold An integer specifying the number of folds for cross-validation.
#'   This parameter is only used if \code{threshold_method} is set to \code{"cv"}.
#' @return The input object updated with SVM classification results added to the
#'  \code{classification.result} slot, including predicted classes and
#'  associated metrics.
#' @export
classification.svm <- function(object, threshold_method = "cv", kfold = 5) {
  # Validate threshold method
  threshold_method <- match.arg(threshold_method, c("cv", "none"))

  # Transform all datasets
  dat.lists <- list(
    train = lapply(object@harmon.train.data, .transform_data),
    test1 = lapply(object@harmon.test1.data, .transform_data),
    test2 = lapply(object@harmon.test2.data, .transform_data)
  )

  # Get labels
  labels <- list(
    train = object@raw.train.data$label,
    test1 = object@raw.test1.data$label,
    test2 = object@raw.test2.data$label
  )

  # Run SVM classification for each harmonization method
  object@classification.result$svm[[threshold_method]] <- lapply(seq_along(dat.lists$train), function(i) {
    # Train model
    svm.intcv.model <- svm.intcv(
      kfold = kfold,
      X = t(dat.lists$train[[i]]),
      y = labels$train,
      threshold_method = threshold_method
    )

    # Get predictions for both test sets
    pred1 <- svm.predict(svm.intcv.model, t(dat.lists$test1[[i]]), labels$test1)
    pred2 <- svm.predict(svm.intcv.model, t(dat.lists$test2[[i]]), labels$test2)

    # Return metrics
    list(train_class = svm.intcv.model$pred, test1_class = pred1$pred, test2_class = pred2$pred)
  })

  return(object)
}

#' Internal function for SVM classification
#' @param kfold Number of folds for cross-validation
#' @param X Training data matrix
#' @param y Labels for training data
#' @param threshold_method Character. Either "cv" for cross-validation optimized
#' threshold or "none" for using all genes
#' @import e1071 tune.svm svm tune.control
#' @noRd
svm.intcv <- function(kfold = 5, X, y, threshold_method = "cv") {
  ptm <- proc.time()
  set.seed(42)

  if (threshold_method == "cv") {
    # Use cross-validation to tune parameters
    svm_tune <- e1071::tune.svm(
      x = data.matrix(X),
      y = factor(y),
      tunecontrol = e1071::tune.control(cross = kfold)
    )
    model <- svm_tune$best.model
  } else {
    # Use default parameters without tuning
    model <- e1071::svm(
      x = data.matrix(X),
      y = factor(y),
      kernel = "radial", # default kernel
      cost = 1, # default cost parameter
      scale = TRUE
    )
  }

  # Calculate error on training set
  mc <- mean(predict(model, data.matrix(X)) != factor(y))
  time <- proc.time() - ptm

  return(list(
    pred = predict(model, data.matrix(X)),
    mc = mc, time = time, model = model
  ))
}

#' Prediction function for SVM model
#' @param svm.intcv.model Trained SVM model
#' @param pred.obj New data matrix for prediction
#' @param pred.obj.group.id Group IDs for the new data
#' @return List containing predicted classes and misclassification counts
#' @noRd
svm.predict <- function(svm.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(svm.intcv.model$model, newdata = pred.obj)
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)

  return(list(pred = pred, mc = mc))
}

##########################################################
# LASSO
##########################################################


#' Classification with LASSO Logistic Regression
#'
#' This function applies LASSO classification to the harmonized data in the
#' input \link{precision} object containing.
#' It supports two thresholding methods: cross-validation to optimize the lambda
#' parameter or using all genes without thresholding.
#'
#' @param object A \link{precision} object containing harmonized data.
#'  Must contain the slots \code{harmon.train.data} with harmonized training data
#'  and \code{harmon.test1.data} and \code{harmon.test2.data} with harmonized test data.
#' @param threshold_method A character string specifying the thresholding method.
#'  Use \code{"cv"} for cross-validation to determine the optimal k,
#'  or \code{"none"} to use all genes without thresholding.
#'  If \code{"none"} is used, a very small lambda value is applied
#'  to effectively include all genes in the model, which may be useful for
#'  datasets where all genes are relevant.
#' @param kfold An integer specifying the number of folds for cross-validation.
#'   This parameter is only used if \code{threshold_method} is set to \code{"cv"}.
#' @return The input object updated with LASSO classification results added to the
#'  \code{classification.result} slot, including predicted classes and
#'  associated metrics.
#' @export
classification.lasso <- function(object, threshold_method = "cv", kfold = 5) {
  # Validate threshold method
  threshold_method <- match.arg(threshold_method, c("cv", "none"))

  # Transform all datasets
  dat.lists <- list(
    train = lapply(object@harmon.train.data, .transform_data),
    test1 = lapply(object@harmon.test1.data, .transform_data),
    test2 = lapply(object@harmon.test2.data, .transform_data)
  )

  # Get labels
  labels <- list(
    train = object@raw.train.data$label,
    test1 = object@raw.test1.data$label,
    test2 = object@raw.test2.data$label
  )

  # Run LASSO classification for each harmonization method
  object@classification.result$lasso[[threshold_method]] <- lapply(seq_along(dat.lists$train), function(i) {
    # Train model
    lasso.model <- lasso.intcv(
      kfold = kfold,
      X = t(dat.lists$train[[i]]),
      y = labels$train,
      threshold_method = threshold_method
    )

    # Get predictions for both test sets
    pred1 <- lasso.predict(lasso.model, t(dat.lists$test1[[i]]), labels$test1)
    pred2 <- lasso.predict(lasso.model, t(dat.lists$test2[[i]]), labels$test2)

    # Return metrics
    list(
      train_class = c(lasso.model$pred),
      test1_class = c(pred1$pred),
      test2_class = c(pred2$pred)
    )
  })

  return(object)
}

#' Internal function for lasso classification
#' @param kfold Number of folds for cross-validation
#' @param X Training data matrix
#' @param y Labels for training data
#' @param threshold_method Character. Either "cv" for cross-validation optimized
#' threshold or "none" for using all genes
#' @param seed Random seed for reproducibility
#' @param alp Numeric. The alpha parameter for glmnet (1 for LASSO, 0 for Ridge)
#' @import glmnet
#' @noRd
lasso.intcv <- function(kfold = 5, X, y, threshold_method = "cv", seed = 1, alp = 1) {
  ptm <- proc.time()
  set.seed(seed)

  if (threshold_method == "cv") {
    # Use cross-validation to find optimal lambda
    cv.fit <- glmnet::cv.glmnet(
      x = data.matrix(X),
      y = factor(y),
      family = "binomial",
      type.measure = "class",
      alpha = alp,
      nfold = kfold
    )
    lambda_value <- cv.fit$lambda.min
    pred <- predict(cv.fit, newx = data.matrix(X), s = lambda_value, type = "class")
    mc <- mean(pred != y)
  } else {
    # Use very small lambda to effectively include all genes
    tiny_lambda <- 1e-6
    cv.fit <- glmnet::glmnet(
      x = data.matrix(X),
      y = factor(y),
      family = "binomial",
      alpha = alp,
      lambda = tiny_lambda
    )
    lambda_value <- tiny_lambda

    # Calculate error on training set
    pred <- predict(cv.fit, newx = data.matrix(X), s = lambda_value, type = "class")
    mc <- mean(pred != y)
  }

  coefs <- if (threshold_method == "cv") {
    trycatch.func(coef(cv.fit, s = lambda_value))
  } else {
    trycatch.func(coef(cv.fit, s = lambda_value))
  }

  return(list(
    pred = pred,
    mc = mc,
    time = proc.time() - ptm,
    model = cv.fit,
    lambda = lambda_value,
    cfs = coefs
  ))
}

#' Prediction function for LASSO model
#' @param lasso.intcv.model Trained LASSO model
#' @param pred.obj New data matrix for prediction
#' @param pred.obj.group.id Group IDs for the new data
#' @return List containing predicted classes and misclassification counts
#' @noRd
lasso.predict <- function(lasso.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(
    lasso.intcv.model$model,
    newx = pred.obj,
    s = lasso.intcv.model$lambda,
    type = "class"
  )

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  prob <- predict(
    lasso.intcv.model$model,
    newx = pred.obj,
    s = lasso.intcv.model$lambda
  )

  return(list(pred = pred, mc = mc, prob = prob))
}

##########################################################
# Random Forest
##########################################################


#' Random Forest Classification
#'
#' This function applies Random Forest classification to the harmonized data in the
#' input \link{precision} object containing.
#' It supports two thresholding methods: cross-validation to optimize the tuning parameter
#' or using default parameters without tuning.
#'
#' @param object A \link{precision} object containing harmonized data.
#'  Must contain the slots \code{harmon.train.data} with harmonized training data
#'  and \code{harmon.test1.data} and \code{harmon.test2.data} with harmonized test data.
#' @param threshold_method A character string specifying the thresholding method.
#'  Use \code{"cv"} for cross-validation to determine the optimal k,
#'  or \code{"none"} to use default parameters without tuning.
#' @param kfold An integer specifying the number of folds for cross-validation.
#'   This parameter is only used if \code{threshold_method} is set to \code{"cv"}.
#' @return The input object updated with Random Forest classification results added to the
#'  \code{classification.result} slot, including predicted classes and
#'  associated metrics.
#' @export
classification.ranfor <- function(object, threshold_method = "cv", kfold = 5) {
  # Validate threshold method
  threshold_method <- match.arg(threshold_method, c("cv", "none"))

  # Transform all datasets
  dat.lists <- list(
    train = lapply(object@harmon.train.data, .transform_data),
    test1 = lapply(object@harmon.test1.data, .transform_data),
    test2 = lapply(object@harmon.test2.data, .transform_data)
  )

  # Get labels
  labels <- list(
    train = object@raw.train.data$label,
    test1 = object@raw.test1.data$label,
    test2 = object@raw.test2.data$label
  )

  # Run Random Forest classification for each harmonization method
  object@classification.result$ranfor[[threshold_method]] <- lapply(seq_along(dat.lists$train), function(i) {
    # Train model
    ranfor.model <- ranfor.intcv(
      kfold = kfold,
      X = t(dat.lists$train[[i]]),
      y = labels$train,
      threshold_method = threshold_method
    )

    # Get predictions for both test sets
    pred1 <- ranfor.predict(ranfor.model, t(dat.lists$test1[[i]]), labels$test1)
    pred2 <- ranfor.predict(ranfor.model, t(dat.lists$test2[[i]]), labels$test2)

    # Return metrics
    list(
      train_class = c(ranfor.model$pred),
      test1_class = c(pred1$pred),
      test2_class = c(pred2$pred)
    )
  })

  return(object)
}

#' Internal function for Random Forest classification
#' @param kfold Number of folds for cross-validation
#' @param X Training data matrix
#' @param y Labels for training data
#' @param threshold_method Character. Either "cv" for cross-validation optimized
#' threshold or "none" for using default parameters without tuning
#' @param seed Random seed for reproducibility
#' @import caret
#' @noRd
ranfor.intcv <- function(kfold = 5, X, y, threshold_method = "cv", seed = 1) {
  ptm <- proc.time()
  set.seed(seed)

  if (threshold_method == "cv") {
    # Use cross-validation to tune parameters
    control <- caret::trainControl(
      method = "cv",
      number = kfold,
      search = "random"
    )

    rf <- caret::train(
      x = data.matrix(X),
      y = factor(y),
      method = "ranger",
      metric = "Accuracy",
      tuneLength = 3,
      preProcess = c("center", "scale"),
      trControl = control
    )

    pred <- predict(rf, newdata = data.matrix(X))
    mc <- 1 - max(rf$results$Accuracy)
  } else {
    # Use default parameters without tuning
    control <- caret::trainControl(method = "none")

    rf <- caret::train(
      x = data.matrix(X),
      y = factor(y),
      method = "ranger",
      metric = "Accuracy",
      tuneGrid = data.frame(
        mtry = floor(sqrt(ncol(X))), # default mtry for classification
        splitrule = "gini",
        min.node.size = 1
      ),
      preProcess = c("center", "scale"),
      trControl = control
    )

    # Calculate error on training set
    pred <- predict(rf, newdata = data.matrix(X))
    mc <- mean(pred != y)
  }

  return(list(
    pred = pred,
    mc = mc,
    time = proc.time() - ptm,
    model = rf,
    cfs = NULL
  ))
}

#' Prediction function for Random Forest model
#' @param ranfor.intcv.model Trained Random Forest model
#' @param pred.obj New data matrix for prediction
#' @param pred.obj.group.id Group IDs for the new data
#' @return List containing predicted classes and misclassification counts
#' @noRd
ranfor.predict <- function(ranfor.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(ranfor.intcv.model$model, newdata = pred.obj)
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)

  return(list(pred = pred, mc = mc))
}
