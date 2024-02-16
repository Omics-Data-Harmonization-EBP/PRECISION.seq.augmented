
library(e1071)
library(caret)

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

tabulate.ext.err.func <- function(pred.obj, obs.grp)
  return(1 - sum(diag(table(pred.obj, obs.grp)))/length(obs.grp))

#pam

new.pamr.cv <- function (fit, data, nfold = 5, ...){
  x <- data$x[fit$gene.subset, fit$sample.subset]
  if (is.null(fit$newy)) {
    y <- factor(data$y[fit$sample.subset])
  }
  else {
    y <- factor(data$newy[fit$sample.subset])
  }
  this.call <- match.call()
  nsccv2 <- get("nsccv", envir = asNamespace("pamr"))
  balanced.folds <- get("balanced.folds", envir = asNamespace("pamr"))
  folds = balanced.folds(y, nfolds = nfold)
  junk <- nsccv2(x, y, object = fit, folds = folds,
                 survival.time = data$survival.time,
                 censoring.status = data$censoring.status,
                 ngroup.survival = fit$ngroup.survival,
                 problem.type = fit$problem.type,
                 ...) # changed here
  junk$call <- this.call
  junk$newy <- fit$newy
  junk$sample.subset <- fit$sample.subset
  return(junk)
}


pam.intcv <- function(X, y, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL, seed=1){

  ptm <- proc.time()
  set.seed(seed)
  data.pam  <- list(x = X, y = factor(y), geneids = rownames(X), genenames = rownames(X))
  fit.pam	<- pamr::pamr.train(data.pam, threshold=vt.k, n.threshold = n.k)
  fit.cv <-  new.pamr.cv(fit = fit.pam, data = data.pam, nfold = kfold)
  best.threshold <- fit.cv$threshold[max(which(fit.cv$error == min(fit.cv$error)))]

  mc <- fit.cv$error[which.min(fit.cv$error)]

  model <- pamr::pamr.train(data.pam, threshold = best.threshold, n.threshold = n.k)

  ## if nonzero == 0 (no feature selected)
  coefs <- trycatch.func(pamr::pamr.listgenes(model, data.pam, threshold = best.threshold))

  time <- proc.time() - ptm
  return(list(mc = mc, time = time, model = model, cfs = coefs))
}


pam.predict <- function(pam.intcv.model, pred.obj, pred.obj.group.id){
  pred <- pamr::pamr.predict(pam.intcv.model$model, newx = pred.obj,
                             threshold = pam.intcv.model$model$threshold,
                             type = "class")

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  prob <- pamr::pamr.predict(pam.intcv.model$model, newx = pred.obj,
                             threshold = pam.intcv.model$model$threshold,
                             type = "posterior")
  return(list(pred=pred, mc=mc, prob=prob))
}


classification.pam <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){

  dat.listtrain <- lapply(object@harmon.train.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})

  dat.listtest <- lapply(object@harmon.test.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})


  pam <- function (datatrain,datatest,labeltrain, labeltest, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){

    #X = trainobject@norm.data
    #y = trainobject@raw.data$label
    X = datatrain
    y = labeltrain

    pam.intcv.model=pam.intcv(X=X,y=y)
    mccv=pam.intcv.model$mc
    pam.pred <- pam.predict(pam.intcv.model = pam.intcv.model, pred.obj = datatest, pred.obj.group.id = labeltest)

    mcexternal <- pam.pred$mc

    return(c(mccv,mcexternal))
  }

  c1<-lapply(1:length(dat.listtrain),function (i) pam(dat.listtrain[[i]],dat.listtest[[i]],object@raw.train.data$label,object@raw.test.data$label))
  names(c1) <- c('TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVg','RUVs','RUVr')
  object@classification.result$pam <- c1
  return(object)
}

#knn

knn.intcv <- function(kfold = 5, X, y, seed=1){
  ptm <- proc.time()
  set.seed(seed)


  ctrl <- trainControl(method = "repeatedcv",
                       repeats = 3,
                       number = kfold)
  knn <- train(x = data.matrix((X)), y = factor(y),
               method = "knn",
               tuneLength = 9,
               trControl = ctrl)

  time <- proc.time() - ptm
  return(list(mc = 1 - max(knn$results$Accuracy), time = time, model = knn, cfs = NULL))
}

knn.predict <- function(knn.intcv.model, pred.obj, pred.obj.group.id){


  pred <- predict(knn.intcv.model$model, newdata = data.matrix((pred.obj)))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)

  return(list(pred=pred, mc=mc))
}

classification.knn <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){

  dat.listtrain <- lapply(object@harmon.train.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})

  dat.listtest <- lapply(object@harmon.test.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})

  knn <- function(datatrain,datatest,labeltrain, labeltest, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){
    knn.intcv.model=knn.intcv(kfold = kfold, X=t(datatrain), y=labeltrain)
    mccv=knn.intcv.model$mc
    pred.obj=t(datatest)
    pred.obj.group.id=labeltest
    pred <- knn.predict(knn.intcv.model,pred.obj,pred.obj.group.id)
    mcexternal=pred$mc

    return(c(mccv,mcexternal))
  }


  c1<-lapply(1:length(dat.listtrain),function (i) knn(dat.listtrain[[i]],dat.listtest[[i]],object@raw.train.data$label,object@raw.test.data$label))
  names(c1) <- c('TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVg','RUVs','RUVr')
  object@classification.result$knn <- c1
  return(object)
}




#svm

svm.intcv <- function(kfold = 5, X, y, seed=1){
  ptm <- proc.time()
  set.seed(seed)

  svm_tune = tune.svm(x = data.matrix((X)), y = factor(y), tunecontrol = tune.control(cross = kfold))

  time <- proc.time() - ptm
  return(list(mc = svm_tune$best.performance, time = time, model = svm_tune$best.model))
}

svm.predict <- function(svm.intcv.model, pred.obj, pred.obj.group.id){

  pred <- predict(svm.intcv.model$model, newdata = (pred.obj))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)

  return(list(pred=pred, mc=mc))
}

classification.svm <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){

  dat.listtrain <- lapply(object@harmon.train.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})

  dat.listtest <- lapply(object@harmon.test.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  svm <- function(datatrain,datatest,labeltrain, labeltest, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){
    svm.intcv.model=svm.intcv(kfold = kfold, X=t(datatrain), y=labeltrain)
    mccv=svm.intcv.model$mc
    pred.obj=t(datatest)
    pred.obj.group.id=labeltest
    pred <- svm.predict(svm.intcv.model,pred.obj,pred.obj.group.id)
    mcexternal=pred$mc

    return(c(mccv,mcexternal))
  }


  c1<-lapply(1:length(dat.listtrain),function (i) svm(dat.listtrain[[i]],dat.listtest[[i]],object@raw.train.data$label,object@raw.test.data$label))
  names(c1) <- c('TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVg','RUVs','RUVr')
  object@classification.result$svm <- c1
  return(object)
}





#lasso

lasso.intcv<-function(kfold = 5, X, y, seed = 1, alp = 1){
  ptm <- proc.time()
  set.seed(seed)

  cv.fit <- glmnet::cv.glmnet(x = data.matrix((X)), y = factor(y),
                              family = "binomial", type.measure = "class", alpha = alp, nfold = kfold)
  mc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.1se)]
  #best.lambda <- cv.fit$lambda.1se # can be extracted from cv.fit
  coefs <- trycatch.func(coef(cv.fit, s = "lambda.1se"))
  time <- proc.time() - ptm
  return(list(mc=mc, time=time, model=cv.fit, cfs=coefs))
}

lasso.predict<-function(lasso.intcv.model, pred.obj, pred.obj.group.id){
  pred <- predict(lasso.intcv.model$model, newx = (pred.obj),
                  s = lasso.intcv.model$model$lambda.1se,
                  type = "class")

  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  prob <- predict(lasso.intcv.model$model, newx = (pred.obj),
                  s = lasso.intcv.model$model$lambda.1se)
  return(list(pred=pred, mc=mc, prob=prob))
}

classification.lasso <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){

  dat.listtrain <- lapply(object@harmon.train.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})

  dat.listtest <- lapply(object@harmon.test.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  lasso <- function(datatrain,datatest,labeltrain, labeltest, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){
    lasso.intcv.model=lasso.intcv(kfold = kfold, X=t(datatrain), y=labeltrain)
    mccv=lasso.intcv.model$mc
    pred.obj=t(datatest)
    pred.obj.group.id=labeltest
    pred <- lasso.predict(lasso.intcv.model,pred.obj,pred.obj.group.id)
    mcexternal=pred$mc

    return(c(mccv,mcexternal))
  }

  c1<-lapply(1:length(dat.listtrain),function (i) lasso(dat.listtrain[[i]],dat.listtest[[i]],object@raw.train.data$label,object@raw.test.data$label))
  names(c1) <- c('TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVg','RUVs','RUVr')
  object@classification.result$lasso <- c1
  return(object)
}



#random forest

ranfor.intcv <- function(kfold = 5, X, y, seed=1){
  ptm <- proc.time()
  set.seed(seed)


  control <- trainControl(method='cv',
                          number=5,
                          search = 'random')

  rf <- train(x = data.matrix((X)), y = factor(y),
              method = 'rf',
              metric = 'Accuracy',
              tuneLength = 5,
              preProcess = c("center", "scale"),
              trControl = control)



  time <- proc.time() - ptm
  return(list(mc = 1 - max(rf$results$Accuracy), time = time, model = rf$finalModel, cfs = NULL))
}

ranfor.predict <- function(ranfor.intcv.model, pred.obj, pred.obj.group.id){

  pred <- predict(ranfor.intcv.model$model, newdata = (pred.obj))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)

  return(list(pred=pred, mc=mc))
}

classification.ranfor <- function(object, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){

  dat.listtrain <- lapply(object@harmon.train.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  
  dat.listtest <- lapply(object@harmon.test.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  ranfor <- function(datatrain,datatest,labeltrain, labeltest, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL){
    ranfor.intcv.model=ranfor.intcv(kfold = kfold, X=t(datatrain), y=labeltrain)
    mccv=ranfor.intcv.model$mc
    pred.obj=t(datatest)
    pred.obj.group.id=labeltest
    pred <- ranfor.predict(ranfor.intcv.model,pred.obj,pred.obj.group.id)
    mcexternal=pred$mc

    return(c(mccv,mcexternal))
  }

  c1<-lapply(1:length(dat.listtrain),function (i) ranfor(dat.listtrain[[i]],dat.listtest[[i]],object@raw.train.data$label,object@raw.test.data$label))
  names(c1) <- c('TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVg','RUVs','RUVr')
  object@classification.result$ranfor <- c1
  return(object)
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
