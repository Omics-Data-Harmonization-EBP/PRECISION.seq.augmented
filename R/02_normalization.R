library(sva)
library(EDASeq)
library(edgeR)
library(RUVSeq)
library(Biobase)
library(BiocGenerics)
library(tidyverse)

## Scaling based methods
norm.TMM <- function(object) {
  raw <- object@raw.data$data
  groups <- object@raw.data$label

  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  d <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- d$samples$norm.factors * d$samples$lib.size / 1e6
  dat.normed <- cpm(d)
  object@norm.data$TMM <- list(dat.normed = dat.normed, scaling.factor = scaling.factor)
  return(object)
}

norm.TC <- function(object) {
  raw <- object@raw.data$data
  groups <- object@raw.data$label

  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  scaling.factor <- dat.DGE$samples$lib.size/1e6
  dat.normed <- cpm(dat.DGE, normalized.lib.sizes = F)
  object@norm.data$TC <- list(dat.normed = dat.normed, scaling.factor = scaling.factor)
  return(object)
}

norm.UQ <- function(object) {
  raw <- object@raw.data$data
  groups <- object@raw.data$label

  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor/1e6
  dat.normed <- t(t(raw)/scaling.factor)
  object@norm.data$UQ <- list(dat.normed = dat.normed, scaling.factor = scaling.factor)
  return(object)
}

norm.med <- function(object) {
  raw <- object@raw.data$data
  groups <- object@raw.data$label

  dat.DGE <- DGEList(counts = raw,
                     group = factor(groups),
                     genes = rownames(raw))
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor/1e6
  dat.normed <- t(t(raw)/scaling.factor)
  object@norm.data$med <- list(dat.normed = dat.normed, scaling.factor = scaling.factor)
  return(object)
}

norm.DESeq <- function(object) {
  raw <- object@raw.data$data
  groups <- object@raw.data$label

  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) = colnames(raw)
  dat.DGE <- DESeq2::DESeqDataSetFromMatrix(countData = raw, colData = condition, design = ~ Condition)
  dat.DGE <- DESeq2::estimateSizeFactors(dat.DGE)
  scaling.factor <- DESeq2::sizeFactors(dat.DGE)
  dat.normed <- DESeq2::counts(dat.DGE, normalized = T)
  object@norm.data$DESeq <- list(dat.normed = dat.normed, scaling.factor = scaling.factor)
  return(object)
}

norm.PoissonSeq <- function(object) {
  raw <- object@raw.data$data

  scaling.factor <- PoissonSeq::PS.Est.Depth(raw)
  dat.normed <- t(t(raw)/scaling.factor)
  object@norm.data$PoissonSeq <- list(dat.normed = dat.normed)
  return(object)
}


## model-based methods
norm.QN <- function(object){
  raw <- object@raw.data$data

  dat.normed <- preprocessCore::normalize.quantiles(raw)
  colnames(dat.normed) <- colnames(raw)
  rownames(dat.normed) <- rownames(raw)
  object@norm.data$QN <- list(dat.normed = dat.normed)
  return(object)
}

norm.SVA <- function(object){ ### methods need biological label
  dat.sva <- object@raw.data$data
  dat.sva <- dat.sva[rowSums(dat.sva) > 0,]
  groups <- object@raw.data$label

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
  dat.normed <- dat.sva - t(as.matrix(adjust[,-c(1:P)]) %*% beta[-c(1:P),])
  object@norm.data$SVA <- list(dat.normed = dat.normed, adjust.factor = svseq)
  return(object)
}

norm.RUVg <- function(object) { ### methods need biological label
  dat.ruv <- object@raw.data$data
  groups <- object@raw.data$label

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
  dat.normed <- normCounts(t)
  object@norm.data$RUVg <- list(dat.normed = dat.normed, adjust.factor = t$W)
  return(object)
}

norm.RUVs <- function(object) { ### methods need biological label
  dat.ruv <- object@raw.data$data
  groups <- object@raw.data$label

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
  dat.normed <- normCounts(t)
  object@norm.data$RUVs <- list(dat.normed = dat.normed, adjust.factor = t$W)
  return(object)
}

norm.RUVr <- function(object) { ### methods need biological label
  dat.ruv <- object@raw.data$data
  groups <- object@raw.data$label

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
  dat.normed <- normCounts(t)
  object@norm.data$RUVr <- list(dat.normed = dat.normed, adjust.factor = t$W)
  return(object)
}

norm.ComBat.Seq <- function(object, batches){ ### methods need batch information
  dat.combat <- object@raw.data$data
  dat.normed <- sva::ComBat_seq(dat.combat, batch = batches)
  object@norm.data$ComBat.Seq <- list(dat.normed = dat.normed) #### returns the count matrix
  return(object)
}


# ## normalization for all
# norm.all <- function(object){
#   object@norm.data$Raw$dat.normed <- object@raw.data$data
#   object <- norm.TC(object)
#   object <- norm.UQ(object)
#   object <- norm.med(object)
#   object <- norm.TMM(object)
#   object <- norm.DESeq(object)
#   object <- norm.PoissonSeq(object)
#   object <- norm.QN(object)
#   object <- norm.RUVg(object)
#   object <- norm.RUVs(object)
#   object <- norm.RUVr(object)
#
#   return(object)
# }
