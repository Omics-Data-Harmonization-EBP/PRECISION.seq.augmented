#' @import mclust
#' @import aricode
#' @import RSKC
#' @import cluster
#' @import factoextra
#' @importFrom stats kmeans dist hclust cutree
#' @importFrom som som
#' @importFrom magrittr %>%
NULL

#' Preprocess harmonized training data
#' @param harmon.train.data List of harmonized training data
#' @return List of preprocessed data matrices
#' @noRd
.preprocess.data <- function(harmon.train.data) {
  dat.list <- lapply(harmon.train.data, function(x) {
    if (any(x$dat.harmonized < 0)) {
      x$dat.harmonized
    } else {
      log2(x$dat.harmonized + 1)
    }
  })
  lapply(dat.list, function(x) {
    scale(t(x)) %>%
      t() %>%
      na.omit() %>%
      t()
  })
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
  object@cluster.result$kmeans <- lapply(
    dat.list,
    function(x) kmeans(x, k, nstart = 100)$cluster
  )

  return(object)
}

#' Perform hierarchical clustering on harmonized data
#'
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' @param distance Character. Distance measure to use: "euclidean", "pearson", or "spearman"
#'
#' @return Updated object with hierarchical clustering results
#' @export
cluster.hc <- function(object, k = NULL, distance = "euclidean") {
  # Validate distance parameter
  distance <- match.arg(distance, c("euclidean", "pearson", "spearman"))

  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }

  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)

  # Define internal hierarchical clustering function
  hc <- function(data, group, dist_method) {
    if (dist_method == "euclidean") {
      dist_mat <- dist(data, method = "euclidean")
    } else {
      dist_mat <- factoextra::get_dist(data, method = dist_method)
    }
    hctree <- hclust(dist_mat, method = "ward.D2")
    cutree(hctree, group)
  }

  # Perform hierarchical clustering
  object@cluster.result$hc[[distance]] <- lapply(dat.list, function(x) hc(x, k, distance))

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

#' Perform PAM clustering with specified distance measure on harmonized data
#'
#' @param object An object containing harmonized training data
#' @param k Integer. Number of clusters. If NULL, uses the number of unique labels in training data
#' @param distance Character. Distance measure to use: "euclidean", "pearson", or "spearman"
#'
#' @return Updated object with PAM clustering results using specified distance measure
#' @export
cluster.pam <- function(object, k = NULL, distance = "euclidean") {
  # Validate distance parameter
  distance <- match.arg(distance, c("euclidean", "pearson", "spearman"))

  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }

  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)

  # Perform PAM clustering with specified distance measure
  if (distance == "euclidean") {
    object@cluster.result$pam[[distance]] <- lapply(dat.list, function(x) {
      cluster::pam(x, k, nstart = 100)$clustering
    })
  } else {
    frame.cor <- lapply(dat.list, function(x) factoextra::get_dist(x, method = distance))
    object@cluster.result$pam[[distance]] <- lapply(frame.cor, function(x) {
      cluster::pam(x, k, nstart = 100)$clustering
    })
  }

  return(object)
}
