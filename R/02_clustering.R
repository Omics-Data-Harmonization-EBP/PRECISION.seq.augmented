#' @import mclust
#' @import aricode
#' @import RSKC
#' @import cluster
#' @import factoextra
#' @importFrom stats kmeans dist hclust cutree scale
#' @importFrom som som

#' Preprocess harmonized training data
#' @param harmon.train.data List of harmonized training data
#' @return List of preprocessed data matrices
#' @noRd
.preprocess.data <- function(harmon.train.data) {
  dat.list <- lapply(harmon.train.data, function(x) {
    if(any(x$dat.harmonized < 0)){
      x$dat.harmonized
    } else{
      log2(x$dat.harmonized + 1)
    }})
  lapply(dat.list, function(x) scale(t(x)) %>% t %>% na.omit %>% t)
}

#' Perform K-means clustering on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with k-means clustering results
#' @export
cluster.kmeans <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Perform kmeans clustering
  object@cluster.result$kmeans <- lapply(dat.list, 
    function(x) kmeans(x, k, nstart = 100)$cluster)
  
  return(object)
}

#' Perform hierarchical clustering on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with hierarchical clustering results
#' @export
cluster.hc <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Define internal hierarchical clustering function
  hc <- function(data, group) {
    dist_mat <- dist(data, method = "euclidean")
    hctree <- hclust(dist_mat, method = "ward.D2")
    cutree(hctree, group)
  }
  
  # Perform hierarchical clustering
  object@cluster.result$hc <- lapply(dat.list, function(x) hc(x, k))
  
  return(object)
}

#' Perform Self-Organizing Map (SOM) clustering on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with SOM clustering results
#' @export
cluster.som <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Perform SOM clustering
  object@cluster.result$som <- lapply(dat.list, function(x) {
    som::som(x, xdim = as.numeric(k), ydim = 1)$visual$x + 1
  })
  
  return(object)
}

#' Perform Gaussian Mixture Model clustering on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with Gaussian Mixture Model clustering results
#' @export
cluster.mnm <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Perform Gaussian Mixture Model clustering
  object@cluster.result$mnm <- lapply(dat.list, function(x) {
    mclust::Mclust(x, G = k, verbose = FALSE)$classification
  })
  
  return(object)
}

#' Perform PAM clustering with Euclidean distance on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with PAM clustering results using Euclidean distance
#' @export
cluster.pam.euclidean <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Perform PAM clustering with Euclidean distance
  object@cluster.result$pam.euclidean <- lapply(dat.list, function(x) {
    cluster::pam(x, k, nstart = 100)$clustering
  })
  
  return(object)
}

#' Perform PAM clustering with Pearson correlation distance on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with PAM clustering results using Pearson correlation distance
#' @export
cluster.pam.pearson <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Calculate Pearson correlation distance and perform PAM clustering
  frame.cor <- lapply(dat.list, function(x) factoextra::get_dist(x, method = "pearson"))
  object@cluster.result$pam.pearson <- lapply(frame.cor, function(x) {
    cluster::pam(x, k, nstart = 100)$clustering
  })
  
  return(object)
}

#' Perform PAM clustering with Spearman correlation distance on harmonized data
#' 
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' 
#' @return Updated object with PAM clustering results using Spearman correlation distance
#' @export
cluster.pam.spearman <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }
  
  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)
  
  # Calculate Spearman correlation distance and perform PAM clustering
  frame.cor <- lapply(dat.list, function(x) factoextra::get_dist(x, method = "spearman"))
  object@cluster.result$pam.spearman <- lapply(frame.cor, function(x) {
    cluster::pam(x, k, nstart = 100)$clustering
  })
  
  return(object)
}
