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
cluster.all <- function(object) {
  object <- cluster.hc(object, distance = "euclidean")
  object <- cluster.hc(object, distance = "pearson")
  object <- cluster.hc(object, distance = "spearman")
  object <- cluster.mnm(object)
  object <- cluster.kmeans(object)
  object <- cluster.som(object)
  object <- cluster.pam(object, distance = "euclidean")
  object <- cluster.pam(object, distance = "pearson")
  object <- cluster.pam(object, distance = "spearman")
  return(object)
}

measure <- function(k, datalist) {
  analysis <- create.precision.cluster(data = datalist[[k]]$data, label = datalist[[k]]$group)
  analysis <- harmon.all(analysis)
  analysis <- cluster.all(analysis)
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
  true_label <- as.factor(datalist[[k]]$group)

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



# scenario 5: sample size
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

subsample_sample_size_element <- function(df, sample_size) {
  set.seed(42)
  labels <- df[[1]]$group
  labels_MXF <- which(labels == "MXF")
  labels_PMFH <- which(labels == "PMFH")

  # Subsample
  subsample_MXF <- sample(labels_MXF, sample_size)
  subsample_PMFH <- sample(labels_PMFH, sample_size)

  # Combine subsamples
  subsample_indices <- c(subsample_MXF, subsample_PMFH)
  subsample_list <- lapply(df, function(x) {
    list(
      data = x$data[, subsample_indices],
      group = x$group[subsample_indices]
    )
  })

  return(subsample_list)
}
sample_size_results <- function(size) {
  dirtydata <- dirty.scenario(c = 1.5, d = 1.2)
  dirtydata_sample_size <- subsample_sample_size_element(dirtydata, size)
  indexessummary <- mclapply(1:300, measure, datalist = dirtydata_sample_size, mc.cores = 50)
  gc()

  return(list(
    dirtydata = dirtydata_sample_size,
    indexessummary = indexessummary
  ))
}

# sample_size_results_list <- list()
# for (param in c(27, 50, 100)) {
#   key <- paste0("sample_size_res_s", param)
#   sample_size_results_list[[key]] <- sample_size_results(param)
#   gc()
# }
# save(sample_size_results_list, file = "clustering_sceanrio5.RData")

## more for server
load("clustering_sceanrio5.RData")
for (i in 2:3) {
  null_indices <- which(sapply(sample_size_results_list[[i]]$indexessummary, is.null))
  rerun_results <- mclapply(null_indices, measure,
    datalist = sample_size_results_list[[i]]$dirtydata, mc.cores = 50
  )
  sample_size_results_list[[i]]$indexessummary[null_indices] <- rerun_results
}
save(sample_size_results_list, file = "clustering_sceanrio5.RData")
