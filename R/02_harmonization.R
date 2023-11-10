library(sva)
library(EDASeq)
library(edgeR)
library(RUVSeq)
library(Biobase)
library(BiocGenerics)
library(tidyverse)

## Scaling based methods
harmon.method.TMM <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  d <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- d$samples$norm.factors * d$samples$lib.size / 1e6
  dat.harmonized <- cpm(d)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}
harmon.TMM <- function(object){
  object@harmon.train.data$TMM <- harmon.method.TMM(raw = object@raw.train.data$data,
                                                    groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$TMM <- harmon.method.TMM(raw = object@raw.test.data$data,
                                                     groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.TC <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  scaling.factor <- dat.DGE$samples$lib.size/1e6
  dat.harmonized <- cpm(dat.DGE, normalized.lib.sizes = F)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}
harmon.TC <- function(object){
  object@harmon.train.data$TC <- harmon.method.TC(raw = object@raw.train.data$data,
                                                  groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$TC <- harmon.method.TC(raw = object@raw.test.data$data,
                                                   groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.UQ <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor/1e6
  dat.harmonized <- t(t(raw)/scaling.factor)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}
harmon.UQ <- function(object){
  object@harmon.train.data$UQ <- harmon.method.UQ(raw = object@raw.train.data$data,
                                                  groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$UQ <- harmon.method.UQ(raw = object@raw.test.data$data,
                                                   groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.med <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor/1e6
  dat.harmonized <- t(t(raw)/scaling.factor)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}
harmon.med <- function(object){
  object@harmon.train.data$med <- harmon.method.med(raw = object@raw.train.data$data,
                                                    groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$med <- harmon.method.med(raw = object@raw.test.data$data,
                                                     groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.DESeq <- function(raw, groups) {
  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) = colnames(raw)
  dat.DGE <- DESeq2::DESeqDataSetFromMatrix(countData = raw, colData = condition, design = ~ Condition)
  dat.DGE <- DESeq2::estimateSizeFactors(dat.DGE)
  scaling.factor <- DESeq2::sizeFactors(dat.DGE)
  dat.harmonized <- DESeq2::counts(dat.DGE, normalized = T)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}
harmon.DESeq <- function(object){
  object@harmon.train.data$DESeq <- harmon.method.DESeq(raw = object@raw.train.data$data,
                                                        groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$DESeq <- harmon.method.DESeq(raw = object@raw.test.data$data,
                                                         groups = object@raw.test.data$label)
  }
  return(object)
}


harmon.method.PoissonSeq <- function(raw) {
  scaling.factor <- PoissonSeq::PS.Est.Depth(raw)
  dat.harmonized <- t(t(raw)/scaling.factor)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}
harmon.PoissonSeq <- function(object){
  object@harmon.train.data$PoissonSeq <- harmon.method.PoissonSeq(raw = object@raw.train.data$data)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$PoissonSeq <- harmon.method.PoissonSeq(raw = object@raw.test.data$data)
  }
  return(object)
}




## model-based methods
harmon.method.QN <- function(object){
  raw <- object@raw.data$data

  dat.harmonized <- preprocessCore::normalize.quantiles(as.matrix(raw))
  colnames(dat.harmonized) <- colnames(raw)
  rownames(dat.harmonized) <- rownames(raw)
  object@harmon.data$QN <- list(dat.harmonized = dat.harmonized)
  return(object)
}
harmon.QN <- function(object){
  object@harmon.train.data$QN <- harmon.method.QN(raw = object@raw.train.data$data,
                                                  groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$QN <- harmon.method.QN(raw = object@raw.test.data$data,
                                                   groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.SVA <- function(raw, groups){ ### methods need biological label
  dat.sva <- raw
  dat.sva <- dat.sva[rowSums(dat.sva) > 0,]

  mod1 <- model.matrix(~ groups)
  mod0 <- model.matrix(~ 1, data.frame(mod1))
  dat0 <- as.matrix(dat.sva)
  # svseq <- svaseq(dat0, mod1, mod0, n.sv = 1)$sv
  n.sv <- num.sv(dat0, mod1)
  invisible(capture.output(svseq <- svaseq(dat0, mod1, mod0, n.sv = n.sv)$sv))
  adjust <- cbind(mod1, svseq)
  hat <- solve(t(adjust) %*% adjust) %*% t(adjust)
  beta <- (hat %*% t(dat.sva))
  P <- ncol(mod1)
  dat.harmonized <- dat.sva - t(as.matrix(adjust[,-c(1:P)]) %*% beta[-c(1:P),])
  res <- list(dat.harmonized = dat.harmonized, adjust.factor = svseq)
  return(res)
}
harmon.SVA <- function(object){
  object@harmon.train.data$SVA <- harmon.method.SVA(raw = object@raw.train.data$data,
                                                    groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$SVA <- harmon.method.SVA(raw = object@raw.test.data$data,
                                                     groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.RUVg <- function(raw, groups) { ### methods need biological label
  dat.ruv <- raw

  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15*nrow(dat.ruv))]))]
  t <- RUVg(set, spikes, k = 1)
  dat.harmonized <- normCounts(t)
  res <- list(dat.harmonized = dat.harmonized, adjust.factor = t$W)
  return(res)
}
harmon.RUVg <- function(object){
  object@harmon.train.data$RUVg <- harmon.method.RUVg(raw = object@raw.train.data$data,
                                                      groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$RUVg <- harmon.method.RUVg(raw = object@raw.test.data$data,
                                                       groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.RUVs <- function(raw, groups) { ### methods need biological label
  dat.ruv <- raw

  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15*nrow(dat.ruv))]))]
  differences <- makeGroups(condition)
  controls <- rownames(dat.ruv)
  t <- RUVs(set, controls, k = 1, differences)
  dat.harmonized <- normCounts(t)
  res <- list(dat.harmonized = dat.harmonized, adjust.factor = t$W)
  return(res)
}
harmon.RUVs <- function(object){
  object@harmon.train.data$RUVs <- harmon.method.RUVs(raw = object@raw.train.data$data,
                                                      groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$RUVs <- harmon.method.RUVs(raw = object@raw.test.data$data,
                                                       groups = object@raw.test.data$label)
  }
  return(object)
}



harmon.method.RUVr <- function(raw, groups) { ### methods need biological label
  dat.ruv <- raw

  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15*nrow(dat.ruv))]))]
  design <- model.matrix(~ condition, data = pData(set))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type = "deviance")
  setUQ <- betweenLaneNormalization(set, which = "upper")
  controls <- rownames(dat.ruv)
  t <- RUVr(setUQ, controls, k = 1, res)
  dat.harmonized <- normCounts(t)
  res <- list(dat.harmonized = dat.harmonized, adjust.factor = t$W)
  return(res)
}
harmon.RUVr <- function(object){
  object@harmon.train.data$RUVr <- harmon.method.RUVr(raw = object@raw.train.data$data,
                                                      groups = object@raw.train.data$label)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$RUVr <- harmon.method.RUVr(raw = object@raw.test.data$data,
                                                       groups = object@raw.test.data$label)
  }
  return(object)
}


harmon.method.ComBat.Seq <- function(raw, batches){ ### methods need batch information
  dat.combat <- raw
  dat.harmonized <- sva::ComBat_seq(dat.combat, batch = batches)
  object@harmon.data$ComBat.Seq <- list(dat.harmonized = dat.harmonized) #### returns the count matrix
  return(object)
}
harmon.RUVr <- function(object, batches){
  object@harmon.train.data$ComBat.Seq <- harmon.method.ComBat.Seq(raw = object@raw.train.data$data,
                                                                  batches = batches)
  if(!is.null(obejct@raw.test.data)){
    object@harmon.test.data$ComBat.Seq <- harmon.method.ComBat.Seq(raw = object@raw.test.data$data,
                                                                   batches = batches)
  }
  return(object)
}


# ## normalization for all
# harmon.all <- function(object){
#   object@harmon.data$Raw$dat.harmonized <- object@raw.data$data
#   object <- harmon.TC(object)
#   object <- harmon.UQ(object)
#   object <- harmon.med(object)
#   object <- harmon.TMM(object)
#   object <- harmon.DESeq(object)
#   object <- harmon.PoissonSeq(object)
#   object <- harmon.QN(object)
#   object <- harmon.RUVg(object)
#   object <- harmon.RUVs(object)
#   object <- harmon.RUVr(object)
#
#   return(object)
# }
