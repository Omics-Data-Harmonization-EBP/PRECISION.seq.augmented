library(parallel)
library(edgeR)
library(tidyverse)
library(EDASeq)
library(RUVSeq)
library(mclust)
library(PRECISION.seq.augmented)
load("MSKpair_300_cluster.RData")

harmon.all <- function(object) {
  object@harmon.train.data$Raw$dat.harmonized <- object@raw.train.data$data
  object <- harmon.TC(object)
  object <- harmon.UQ(object)
  object <- harmon.med(object)
  object <- harmon.TMM(object)
  object <- harmon.DESeq(object)
  object <- harmon.PoissonSeq(object)
  object <- harmon.QN(object)
  object <- harmon.RUVr(object)
  object <- harmon.RUVs(object)
  object <- harmon.RUVg(object)
  return(object)
}
cluster.all <- function(object, k) {
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

measure <- function(k, datalist, num_cluster) {
  analysis <- create.precision.cluster(data = datalist[[k]]$data, label = datalist[[k]]$group)
  analysis <- harmon.all(analysis)
  analysis <- cluster.all(analysis, k = num_cluster)
  cluster <- c(
    "hc_euclidean", "hc_pearson", "hc_spearman",
    "mnm",
    "kmeans",
    "som",
    "pam_euclidean", "pam_pearson", "pam_spearman"
  )
  harmon <- c(
    "Raw", "TC", "UQ", "med", "TMM", "DESeq", "PoissonSeq", "QN",
    "RUVr", "RUVs", "RUVg"
  )

  ari_indexes <- data.frame(matrix(nrow = length(cluster), ncol = length(harmon)))
  silhouette_indexes <- data.frame(matrix(nrow = length(cluster), ncol = length(harmon)))
  true_label <- as.factor(c(rep("MXF", 100), rep("PMFH", 100)))

  for (i in 1:length(cluster)) {
    for (j in 1:length(harmon)) {
      if (startsWith(cluster[i], "hc_")) {
        distance <- sub("hc_", "", cluster[i])
        est_cluster <- analysis@cluster.result$hc[[distance]][[harmon[j]]]
      } else if (startsWith(cluster[i], "pam_")) {
        distance <- sub("pam_", "", cluster[i])
        est_cluster <- analysis@cluster.result$pam[[distance]][[harmon[j]]]
      } else {
        est_cluster <- analysis@cluster.result[[cluster[i]]][[harmon[j]]]
      }
      # Calculate ARI
      ari_indexes[i, j] <- mclust::adjustedRandIndex(true_label, est_cluster)
      # Calculate silhouette value
      if (startsWith(cluster[i], "hc_") || startsWith(cluster[i], "pam_")) {
        dist_matrix <- switch(sub(".*_", "", cluster[i]),
          "euclidean" = dist(t(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized)),
          "pearson" = as.dist(1 - cor(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized)),
          "spearman" = as.dist(1 - cor(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized, method = "spearman"))
        )
      } else {
        dist_matrix <- dist(t(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized))
      }
      sil <- cluster::silhouette(as.numeric(est_cluster), dist_matrix)
      silhouette_indexes[i, j] <- mean(sil[, "sil_width"])
    }
  }
  rownames(ari_indexes) <- cluster
  colnames(ari_indexes) <- harmon
  rownames(silhouette_indexes) <- cluster
  colnames(silhouette_indexes) <- harmon

  indexes <- list(
    ari = ari_indexes,
    silhouette = silhouette_indexes
  )
  return(indexes)
}


# scenario 3: misspecify the number of clusters
dirty.scenario <- function(c, d) {
  clean.datasets <- lapply(1:300, function(x) {
    biological.effects(
      benchmark_sub[[x]]$data,
      benchmark_sub[[x]]$label, c
    )
  })
  dirty.datasets <- lapply(1:300, function(x) {
    handling.effects(
      clean.datasets[[x]]$data, benchmark_sub[[x]]$data, test_sub[[x]]$data,
      clean.datasets[[x]]$group, d
    )
  })
  return(dirty.datasets)
}

cluster.misspecify.results <- function(num_cluster) {
  dirtydata <- dirty.scenario(c = 1.5, d = 1.2)
  indexessummary <- mclapply(1:300, function(k) {
    measure(k, datalist = dirtydata, num_cluster = num_cluster)
  }, mc.cores = 50)

  return(list(
    indexessummary = indexessummary,
    dirtydata = dirtydata
  ))
}

# misspecify_results_list <- list()
# for (param in c(2, 4, 6, 8)) {
#   key <- paste0("misspecify_res_k", param)
#   misspecify_results_list[[key]] <- cluster.misspecify.results(param)
#   gc()
# }
# save(misspecify_results_list, file = "clustering_sceanrio3.RData")

## more for server
load("clustering_sceanrio3.RData")
for (i in 2:4) {
  null_indices <- which(sapply(misspecify_results_list[[i]]$indexessummary, is.null))
  rerun_results <- mclapply(null_indices, measure,
    datalist = misspecify_results_list[[i]]$dirtydata,
    num_cluster = c(2, 4, 6, 8)[i], mc.cores = 50
  )
  misspecify_results_list[[i]]$indexessummary[null_indices] <- rerun_results
}
save(misspecify_results_list, file = "clustering_sceanrio3.RData")
