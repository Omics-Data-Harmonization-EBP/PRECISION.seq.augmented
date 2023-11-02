library(mclust)
library(aricode)
library(RSKC)
library(cluster)
library(factoextra)
memory.limit(size=1800000)

## Kmeans
cluster.kmeans <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  object@cluster.result$kmeans <- lapply(dat.list, function(x) kmeans(x, k, nstart = 100)$cluster)
  return(object)
}

## HC
cluster.hc <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  hc <- function(data, group){
    dist_mat <- dist(data, method = 'euclidean')
    hctree <- hclust(dist_mat, method = "ward.D2")
    res <- cutree(hctree, group)
    return(res)
  }
  object@cluster.result$hc <- lapply(dat.list, function(x) hc(x, k))
  return(object)
}

## SOM
cluster.som <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  object@cluster.result$som <- lapply(dat.list, function(x) som::som(x, xdim = as.numeric(k), ydim = 1)$visual$x+1)
  return(object)
}

## MNM
cluster.mnm <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  object@cluster.result$mnm <- lapply(dat.list, function(x) mclust::Mclust(x, G = k, verbose = FALSE)$classification)
  return(object)
}

# Robust Sparse kmeans
library(RSKC)
GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  tss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)
  bcss.perfeature <- tss.perfeature-wcss.perfeature
  r <- bcss.perfeature

  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature, r=r))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature, r=r))
}

RSKC.permute <- function (x, K = NULL, nperms = 25, wbounds = NULL, alpha = 0.1, silent = TRUE,
                          nvals = 10, centers = NULL) {
  if (is.null(wbounds))
    wbounds <- exp(seq(log(12), log(sqrt(ncol(x)) * 0.9),
                       len = nvals))
  if (min(wbounds) <= 1)
    stop("Wbounds should be greater than 1, since otherwise only one weight will be nonzero.")
  if (length(wbounds) < 2)
    stop("Wbounds should be a vector of at least two elements.")
  if (is.null(K) && is.null(centers))
    stop("Must provide either K or centers.")
  if (!is.null(K) && !is.null(centers)) {
    if (nrow(centers) != K)
      stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if (nrow(centers) == K)
      K <- NULL
  }
  if (!is.null(centers) && ncol(centers) != ncol(x))
    stop("If centers is provided, then ncol(centers) must equal ncol(x).")
  permx <- list()
  nnonzerows <- NULL
  for (i in 1:nperms) {
    permx[[i]] <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) permx[[i]][, j] <- sample(x[, j])
  }
  tots <- NULL
  out <- lapply(wbounds, function(i) RSKC(d = x, ncl = K, alpha = alpha,  L1 = i))

  for (i in 1:length(out)) {
    nnonzerows <- c(nnonzerows, sum(out[[i]]$weights != 0))
    bcss <- GetWCSS(x, out[[i]]$labels)$bcss.perfeature
    tots <- c(tots, sum(out[[i]]$weight * bcss))
  }
  permtots <- matrix(NA, nrow = length(wbounds), ncol = nperms)
  for (k in 1:nperms) {
    if (!silent)
      cat("Permutation ", k, "of ", nperms,
          fill = TRUE)
    perm.out <- lapply(wbounds, function(i) RSKC(d = permx[[k]], ncl = K, alpha = alpha,  L1 = i))

    for (i in 1:length(perm.out)) {
      perm.bcss <- GetWCSS(permx[[k]], perm.out[[i]]$labels)$bcss.perfeature
      permtots[i, k] <- sum(perm.out[[i]]$weights * perm.bcss)
    }
  }
  gaps <- (log(tots) - apply(log(permtots), 1, mean))
  out <- list(tots = tots, permtots = permtots, nnonzerows = nnonzerows,
              gaps = gaps, sdgaps = apply(log(permtots), 1, sd), wbounds = wbounds,
              bestw = wbounds[which.max(gaps)])
  if (!silent)
    cat(fill = TRUE)
  class(out) <- "RSKC.permute"
  return(out)
}

cluster.rskmeans <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  km.perm = lapply(dat.list, function(x)
    tryCatch({
      RSKC.permute(x, K = k, silent = TRUE, wbounds = exp(seq(log(10), log(sqrt(ncol(x))*.7), len=10)))
    }, error = function(e){
      list(bestw = 12)
    }))

  object@cluster.result$rskmeans <- lapply(1:length(dat.list), function(i) RSKC(dat.list[[i]], ncl = k, alpha = 0.1, L1 = km.perm[[i]]$bestw, silent = TRUE)$labels)
  names(object@cluster.result$rskmeans) <- names(dat.list)
  return(object)
}

## pam
cluster.pam.euclidean <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  object@cluster.result$pam.euclidean <- lapply(dat.list, function(x) cluster::pam(x, 2, nstart = 100)$clustering)
  return(object)
}

cluster.pam.pearson <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  frame.cor <- lapply(dat.list, function(x) get_dist(x, method = "pearson"))
  object@cluster.result$pam.pearson <- lapply(frame.cor, function(x) cluster::pam(x, 2, nstart = 100)$clustering)
  return(object)
}

cluster.pam.spearman <- function(object, k = NULL){
  if(is.null(k)){k = length(unique(object@raw.data$label))}

  dat.list <- lapply(object@harmon.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  dat.list <- lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)

  frame.cor <- lapply(dat.list, function(x) get_dist(x, method = "spearman"))
  object@cluster.result$pam.spearman <- lapply(frame.cor, function(x) cluster::pam(x, 2, nstart = 100)$clustering)
  return(object)
}
