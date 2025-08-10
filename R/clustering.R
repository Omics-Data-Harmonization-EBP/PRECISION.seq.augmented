#' Apply all clustering methods to a precision object
#'
#' This function applies all clustering methods to the training data
#' contained in a precision object.
#' The clustering results are added to the corresponding \code{cluster.result}
#' slot in the object for each training dataset.
#'
#' The following clustering methods are applied to the data:
#'   \itemize{
#'     \item Hierarchical Clustering using multiple distance metrics
#'       (Euclidean, Pearson, and Spearman)
#'       (\code{\link{cluster.hc}})
#'     \item Gaussian Mixture Models (MNM) (\code{\link{cluster.mnm}})
#'     \item K-Means Clustering (\code{\link{cluster.kmeans}})
#'     \item Self-Organizing Maps (\code{\link{cluster.som}})
#'     \item Partitioning Around Medoids (Euclidean, Pearson, and Spearman)
#'       (\code{\link{cluster.pam}})
#'   }
#'
#' @param object A \code{\link{precision}} object containing harmonized training data
#'  in the slot \code{harmon.train.data}.
#' @param k \emph{Integer.} Number of clusters. If NULL or not given,
#' uses the number of unique labels in training data.
#' @return Updated precision object with all clustering results added to
#' the \code{cluster.result} slot.
#' @export cluster.all
cluster.all <- function(object, k = NULL) {
  object <- cluster.hc(object, k = k, distance = "euclidean")
  object <- cluster.hc(object, k = k, distance = "pearson")
  object <- cluster.hc(object, k = k, distance = "spearman")
  object <- cluster.mnm(object, k = k)
  object <- cluster.kmeans(object, k = k)
  object <- cluster.som(object, k = k)
  object <- cluster.pam(object, k = k, distance = "euclidean")
  object <- cluster.pam(object, k = k, distance = "pearson")
  object <- cluster.pam(object, k = k, distance = "spearman")
  return(object)
}


#' Preprocess harmonized training data
#'
#' This function preprocesses the harmonized training data by applying modified
#' log2 transformation and scaling to each given dataset.
#'
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
  lapply(dat.list, function(x) t(na.omit(t(scale(t(x))))))
}

#' K-means clustering for harmonized data
#'
#' This function performs K-means clustering on the harmonized training data
#' and returns the clustering results.
#'
#' @param object A \code{\link{precision}} object containing harmonized training data
#'  in the slot \code{harmon.train.data}.
#' @param k \emph{Integer.} Number of clusters. If NULL or not given,
#' uses the number of unique labels in training data.
#' @return Updated precision object with k-means clustering results added to
#' the \code{cluster.result} slot.
#' @importFrom stats kmeans
#' @export cluster.kmeans
cluster.kmeans <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }

  # Preprocess data
  dat.list <- .preprocess.data(object@harmon.train.data)

  # Perform kmeans clustering
  object@cluster.result$kmeans <- lapply(
    dat.list,
    function(x) stats::kmeans(x, k, nstart = 100)$cluster
  )

  return(object)
}

#' Hierarchical clustering for harmonized data
#'
#' This function performs hierarchical clustering on the harmonized training data
#' using euclidean, pearson, or spearman distance measures.
#' It adds the clustering results to the \code{cluster.result} slot.
#'
#' @param object A \code{\link{precision}} object containing harmonized training data
#'  in the slot \code{harmon.train.data}.
#' @param k \emph{Integer.} Number of clusters. If NULL or not given,
#' uses the number of unique labels in training data.
#' @param distance Distance measure to use for clustering.
#' Options are "euclidean", "pearson", or "spearman". Default is "euclidean".
#'
#' @return Updated precision object with hierarchical clustering results added to the
#' \code{cluster.result} slot.
#' @importFrom factoextra get_dist
#' @importFrom stats dist hclust cutree
#' @export cluster.hc
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
      dist_mat <- stats::dist(data, method = "euclidean")
    } else {
      dist_mat <- factoextra::get_dist(data, method = dist_method)
    }
    hctree <- stats::hclust(dist_mat, method = "ward.D2")
    stats::cutree(hctree, group)
  }

  # Perform hierarchical clustering
  object@cluster.result$hc[[distance]] <- lapply(dat.list, function(x) hc(x, k, distance))

  return(object)
}


#' Self-Organizing Map (SOM) clustering for harmonized data
#'
#' This function performs Self-Organizing Map (SOM) clustering on the harmonized
#' training data. It uses the `som` package to create a SOM model and returns
#' the clustering results.
#'
#' @param object A \code{\link{precision}} object containing harmonized training data
#'  in the slot \code{harmon.train.data}.
#' @param k \emph{Integer.} Number of clusters. If NULL or not given,
#' uses the number of unique labels in training data.
#'
#' @return Updated precision object with SOM clustering results added to the
#' \code{cluster.result} slot.
#' @importFrom som som
#' @export cluster.som
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

#' Gaussian Mixture Model clustering for harmonized data
#'
#' This function performs Gaussian Mixture Model (GMM) clustering on the harmonized
#' training data using the `mclust` package.
#'
#' @param object A \code{\link{precision}} object containing harmonized training data
#'  in the slot \code{harmon.train.data}.
#' @param k \emph{Integer.} Number of clusters. If NULL or not given,
#' uses the number of unique labels in training data.
#'
#' @return Updated precision object with Gaussian Mixture Model clustering
#' results added to the \code{cluster.result} slot.
#' @import mclust
#' @export cluster.mnm
cluster.mnm <- function(object, k = NULL) {
  # Get k from labels if not specified
  if (is.null(k)) {
    k <- length(unique(object@raw.train.data$label))
  }

  # # Fix for error: 'could not find function "mclustBIC"'
  # mclustBIC <- mclust::mclustBIC

  # Reuse preprocessing function
  dat.list <- .preprocess.data(object@harmon.train.data)

  # Perform Gaussian Mixture Model clustering
  object@cluster.result$mnm <- lapply(dat.list, function(x) {
    mclust::Mclust(x, G = k, verbose = FALSE)$classification
  })

  return(object)
}

#' Partitioning Around Medoids (PAM) clustering for harmonized data
#'
#' This function performs Partitioning Around Medoids (PAM) clustering on the harmonized
#' training data using specified distance measures (euclidean, pearson, or spearman).
#' It adds the clustering results to the \code{cluster.result} slot of the precision
#' object.
#'
#' @param object A \code{\link{precision}} object containing harmonized training data
#'  in the slot \code{harmon.train.data}.
#' @param k \emph{Integer.} Number of clusters. If NULL or not given,
#'  uses the number of unique labels in training data.
#' @param distance Distance measure to use for clustering.
#' Options are "euclidean", "pearson", or "spearman". Default is "euclidean".
#'
#' @return Updated precision object with PAM clustering
#' results added to the \code{cluster.result} slot.
#' @importFrom cluster pam
#' @importFrom factoextra get_dist
#' @export cluster.pam
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
