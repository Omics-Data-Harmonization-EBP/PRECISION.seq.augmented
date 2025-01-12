#' @import e1071
#' @import caret
#' @import ranger
#' @import pamr
#' @import glmnet
#' @importFrom stats predict coef
#' @importFrom utils data
NULL

# Internal utility functions
#' @keywords internal
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

  }, finally={
    message(msg)

  })
  return(out)
}

#' Calculate external error rate
#' @param pred.obj Predicted values
#' @param obs.grp Observed groups
#' @return Error rate
tabulate.ext.err.func <- function(pred.obj, obs.grp) {
    return(1 - sum(diag(table(pred.obj, obs.grp)))/length(obs.grp))
}

##########################################################
# PAM
##########################################################

#' Main PAM classification function
#' @param object Input object containing harmonized data
#' @param vt.k Threshold values
#' @param n.k Number of threshold values
#' @param kfold Number of folds
#' @return Updated object with classification results
classification.pam <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL) {
    # Transform data with log2 if non-negative
    transform_data <- function(x) {
        if (!any(x$dat.harmonized < 0)) log2(x$dat.harmonized + 1) else x$dat.harmonized
    }
    
    # Prepare datasets and labels
    datasets <- list(
        data = list(
            train = lapply(object@harmon.train.data, transform_data),
            test1 = lapply(object@harmon.test1.data, transform_data),
            test2 = lapply(object@harmon.test2.data, transform_data)
        ),
        labels = list(
            train = object@raw.train.data$label,
            test1 = object@raw.test1.data$label,
            test2 = object@raw.test2.data$label
        )
    )
    
    # Run PAM classification for each harmonization method
    object@classification.result$pam <- lapply(seq_along(datasets$data$train), function(i) {
        # Train model
        model <- pam.intcv(
            X = datasets$data$train[[i]],
            y = datasets$labels$train,
            vt.k = vt.k,
            n.k = n.k,
            kfold = kfold,
            folds = folds
        )
        
        # Get predictions for both test sets
        pred1 <- pam.predict(model, datasets$data$test1[[i]], datasets$labels$test1)
        pred2 <- pam.predict(model, datasets$data$test2[[i]], datasets$labels$test2)
        
        # Return metrics
        c(mccv = model$mc, mcexternal1 = pred1$mc, mcexternal2 = pred2$mc)
    })
    
    return(object)
}

# Internal helper functions
new.pamr.cv <- function(fit, data, nfold = 5, ...) {
    x <- data$x[fit$gene.subset, fit$sample.subset]
    y <- factor(if (is.null(fit$newy)) data$y[fit$sample.subset] else fit$newy[fit$sample.subset])
    
    cv_results <- get("nsccv", envir = asNamespace("pamr"))(
        x = x, y = y, object = fit,
        folds = get("balanced.folds", envir = asNamespace("pamr"))(y, nfolds = nfold),
        survival.time = data$survival.time,
        censoring.status = data$censoring.status,
        ngroup.survival = fit$ngroup.survival,
        problem.type = fit$problem.type, ...
    )
    
    list2env(list(call = match.call(), newy = fit$newy, 
                  sample.subset = fit$sample.subset), cv_results)
    return(cv_results)
}

pam.intcv <- function(X, y, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL, seed = 1) {
    set.seed(seed)
    start_time <- proc.time()
    
    data.pam <- list(x = X, y = factor(y), 
                     geneids = rownames(X), genenames = rownames(X))
    
    fit.cv <- new.pamr.cv(
        fit = pamr::pamr.train(data.pam, threshold = vt.k, n.threshold = n.k),
        data = data.pam, nfold = kfold
    )
    
    best.threshold <- fit.cv$threshold[max(which(fit.cv$error == min(fit.cv$error)))]
    final_model <- pamr::pamr.train(data.pam, threshold = best.threshold, n.threshold = n.k)
    
    list(
        mc = min(fit.cv$error),
        time = proc.time() - start_time,
        model = final_model,
        cfs = trycatch.func(pamr::pamr.listgenes(final_model, data.pam, threshold = best.threshold))
    )
}

pam.predict <- function(pam.intcv.model, pred.obj, pred.obj.group.id) {
    threshold <- pam.intcv.model$model$threshold
    
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
    
    list(pred = pred, mc = mc, prob = prob)
}


##########################################################
# KNN
##########################################################

#' KNN Classification Function
#' @param object Input object containing harmonized data
#' @param vt.k Threshold values (not used in KNN but kept for consistency)
#' @param n.k Number of thresholds (not used in KNN but kept for consistency)
#' @param kfold Number of folds for cross-validation
#' @param folds Pre-specified folds (optional)
#' @return Updated object with KNN classification results
classification.knn <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL) {
    # Transform data with log2 if non-negative
    transform_data <- function(x) {
        if (!any(x$dat.harmonized < 0)) log2(x$dat.harmonized + 1) else x$dat.harmonized
    }
    
    # Prepare datasets and labels
    datasets <- list(
        data = list(
            train = lapply(object@harmon.train.data, transform_data),
            test1 = lapply(object@harmon.test1.data, transform_data),
            test2 = lapply(object@harmon.test2.data, transform_data)
        ),
        labels = list(
            train = object@raw.train.data$label,
            test1 = object@raw.test1.data$label,
            test2 = object@raw.test2.data$label
        )
    )
    
    # Run KNN classification for each harmonization method
    object@classification.result$knn <- lapply(seq_along(datasets$data$train), function(i) {
        # Train model with transposed data
        model <- knn.intcv(
            kfold = kfold,
            X = t(datasets$data$train[[i]]),
            y = datasets$labels$train
        )
        
        # Get predictions for both test sets
        pred1 <- knn.predict(model, t(datasets$data$test1[[i]]), datasets$labels$test1)
        pred2 <- knn.predict(model, t(datasets$data$test2[[i]]), datasets$labels$test2)
        
        # Return metrics
        c(mccv = model$mc, mcexternal1 = pred1$mc, mcexternal2 = pred2$mc)
    })
    
    return(object)
}

# Internal function for KNN cross-validation
knn.intcv <- function(kfold = 5, X, y, seed = 1) {
    start_time <- proc.time()
    set.seed(seed)
    
    ctrl <- trainControl(method = "repeatedcv", repeats = 2, number = kfold)
    knn_model <- train(
        x = data.matrix(X), y = factor(y),
        method = "knn", tuneLength = 5, trControl = ctrl
    )
    
    list(
        mc = 1 - max(knn_model$results$Accuracy),
        time = proc.time() - start_time,
        model = knn_model,
        cfs = NULL
    )
}

# Internal function for KNN prediction
knn.predict <- function(knn.intcv.model, pred.obj, pred.obj.group.id) {
    pred <- predict(knn.intcv.model$model, newdata = data.matrix(pred.obj))
    mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
    list(pred = pred, mc = mc)
}

##########################################################
# SVM
##########################################################

#' SVM Classification Function
#' @param object Input object containing harmonized data
#' @param kfold Number of folds for cross-validation
#' @return Updated object with SVM classification results
classification.svm <- function(object, kfold = 5) {
    # Transform data with log2 if non-negative
    transform_data <- function(x) {
        if (!any(x$dat.harmonized < 0)) {
            log2(x$dat.harmonized + 1)
        } else {
            x$dat.harmonized
        }
    }
    
    # Transform all datasets
    dat.lists <- list(
        train = lapply(object@harmon.train.data, transform_data),
        test1 = lapply(object@harmon.test1.data, transform_data),
        test2 = lapply(object@harmon.test2.data, transform_data)
    )
    
    # Get labels
    labels <- list(
        train = object@raw.train.data$label,
        test1 = object@raw.test1.data$label,
        test2 = object@raw.test2.data$label
    )
    
    # Run SVM classification for each harmonization method
    object@classification.result$svm <- lapply(seq_along(dat.lists$train), function(i) {
        # Train model
        svm.intcv.model <- svm.intcv(
            kfold = kfold,
            X = t(dat.lists$train[[i]]),
            y = labels$train
        )
        
        # Get predictions for both test sets
        pred1 <- svm.predict(svm.intcv.model, t(dat.lists$test1[[i]]), labels$test1)
        pred2 <- svm.predict(svm.intcv.model, t(dat.lists$test2[[i]]), labels$test2)
        
        # Return metrics
        c(mccv = svm.intcv.model$mc, mcexternal1 = pred1$mc, mcexternal2 = pred2$mc)
    })
    
    return(object)
}

# Internal SVM cross-validation function
svm.intcv <- function(kfold = 5, X, y, seed = 1) {
    ptm <- proc.time()
    set.seed(seed)

    svm_tune <- tune.svm(
        x = data.matrix(X),
        y = factor(y),
        tunecontrol = tune.control(cross = kfold)
    )
    
    mc <- mean(predict(svm_tune$best.model, data.matrix(X)) != factor(y))
    time <- proc.time() - ptm
    
    return(list(mc = mc, time = time, model = svm_tune$best.model))
}

# Internal SVM prediction function
svm.predict <- function(svm.intcv.model, pred.obj, pred.obj.group.id) {
    pred <- predict(svm.intcv.model$model, newdata = pred.obj)
    mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
    
    return(list(pred = pred, mc = mc))
}

##########################################################
# LASSO
##########################################################

#' LASSO Classification Function
#' @param object Input object containing harmonized data
#' @param kfold Number of folds for cross-validation
#' @return Updated object with LASSO classification results
classification.lasso <- function(object, kfold = 5) {
    # Transform data with log2 if non-negative
    transform_data <- function(x) {
        if (!any(x$dat.harmonized < 0)) {
            log2(x$dat.harmonized + 1)
        } else {
            x$dat.harmonized
        }
    }
    
    # Transform all datasets
    dat.lists <- list(
        train = lapply(object@harmon.train.data, transform_data),
        test1 = lapply(object@harmon.test1.data, transform_data),
        test2 = lapply(object@harmon.test2.data, transform_data)
    )
    
    # Get labels
    labels <- list(
        train = object@raw.train.data$label,
        test1 = object@raw.test1.data$label,
        test2 = object@raw.test2.data$label
    )
    
    # Run LASSO classification for each harmonization method
    object@classification.result$lasso <- lapply(seq_along(dat.lists$train), function(i) {
        # Train model
        lasso.model <- lasso.intcv(
            kfold = kfold,
            X = t(dat.lists$train[[i]]),
            y = labels$train
        )
        
        # Get predictions for both test sets
        pred1 <- lasso.predict(lasso.model, t(dat.lists$test1[[i]]), labels$test1)
        pred2 <- lasso.predict(lasso.model, t(dat.lists$test2[[i]]), labels$test2)
        
        # Return metrics
        c(mccv = lasso.model$mc, mcexternal1 = pred1$mc, mcexternal2 = pred2$mc)
    })
    
    return(object)
}

# Internal LASSO cross-validation function
lasso.intcv <- function(kfold = 5, X, y, seed = 1, alp = 1) {
    ptm <- proc.time()
    set.seed(seed)
    
    cv.fit <- glmnet::cv.glmnet(
        x = data.matrix(X),
        y = factor(y),
        family = "binomial",
        type.measure = "class",
        alpha = alp,
        nfold = kfold
    )
    
    mc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.1se)]
    coefs <- trycatch.func(coef(cv.fit, s = "lambda.1se"))
    
    return(list(mc = mc, time = proc.time() - ptm, model = cv.fit, cfs = coefs))
}

# Internal LASSO prediction function
lasso.predict <- function(lasso.intcv.model, pred.obj, pred.obj.group.id) {
    pred <- predict(
        lasso.intcv.model$model,
        newx = pred.obj,
        s = lasso.intcv.model$model$lambda.1se,
        type = "class"
    )
    
    mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
    prob <- predict(
        lasso.intcv.model$model,
        newx = pred.obj,
        s = lasso.intcv.model$model$lambda.1se
    )
    
    return(list(pred = pred, mc = mc, prob = prob))
}

##########################################################
# Random Forest
##########################################################

#' Random Forest Classification Function
#' @param object Input object containing harmonized data
#' @param kfold Number of folds for cross-validation
#' @return Updated object with Random Forest classification results
classification.ranfor <- function(object, fold = 5) {
    # Transform data with log2 if non-negative
    transform_data <- function(x) {
        if (!any(x$dat.harmonized < 0)) {
            log2(x$dat.harmonized + 1)
        } else {
            x$dat.harmonized
        }
    }
    
    # Transform all datasets
    dat.lists <- list(
        train = lapply(object@harmon.train.data, transform_data),
        test1 = lapply(object@harmon.test1.data, transform_data),
        test2 = lapply(object@harmon.test2.data, transform_data)
    )
    
    # Get labels
    labels <- list(
        train = object@raw.train.data$label,
        test1 = object@raw.test1.data$label,
        test2 = object@raw.test2.data$label
    )
    
    # Run Random Forest classification for each harmonization method
    object@classification.result$ranfor <- lapply(seq_along(dat.lists$train), function(i) {
        # Train model
        ranfor.model <- ranfor.intcv(
            kfold = kfold,
            X = t(dat.lists$train[[i]]),
            y = labels$train
        )
        
        # Get predictions for both test sets
        pred1 <- ranfor.predict(ranfor.model, t(dat.lists$test1[[i]]), labels$test1)
        pred2 <- ranfor.predict(ranfor.model, t(dat.lists$test2[[i]]), labels$test2)
        
        # Return metrics
        c(mccv = ranfor.model$mc, mcexternal1 = pred1$mc, mcexternal2 = pred2$mc)
    })
    
    return(object)
}

# Internal Random Forest cross-validation function
ranfor.intcv <- function(kfold = 5, X, y, seed = 1) {
    ptm <- proc.time()
    set.seed(seed)
    
    control <- trainControl(
        method = 'cv',
        number = 5,
        search = 'random'
    )
    
    rf <- train(
        x = data.matrix(X),
        y = factor(y),
        method = "ranger",
        metric = 'Accuracy',
        tuneLength = 3,
        preProcess = c("center", "scale"),
        trControl = control
    )
    
    return(list(
        mc = 1 - max(rf$results$Accuracy),
        time = proc.time() - ptm,
        model = rf,
        cfs = NULL
    ))
}

# Internal Random Forest prediction function
ranfor.predict <- function(ranfor.intcv.model, pred.obj, pred.obj.group.id) {
    pred <- predict(ranfor.intcv.model$model, newdata = pred.obj)
    mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
    
    return(list(pred = pred, mc = mc))
}


## classification methods for all
classification.all <- function(object) {
  object <- classification.pam(object)
  object <- classification.knn(object)
  object <- classification.lasso(object)
  object <- classification.svm(object)
  object <- classification.ranfor(object)

  return(object)
}
