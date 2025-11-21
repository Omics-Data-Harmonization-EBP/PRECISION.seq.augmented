library(parallel)
library(edgeR)
library(tidyverse)
library(EDASeq)
library(RUVSeq)
library(mclust)
library(cluster)
library(PRECISION.seq.augmented)
load("MSKpair_300_cluster.RData")
source("RUV_no_bio_effects.R")

harmon.all <- function(object){
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
cluster.all <- function(object){
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

measure <- function(k, datalist){
  analysis <- create.precision.cluster(data = datalist[[k]]$data, label = datalist[[k]]$group)
  analysis <- harmon.all(analysis)
  analysis <- cluster.all(analysis)
  cluster <- c('hc_euclidean', 'hc_pearson', 'hc_spearman', 
               'mnm', 
               'kmeans', 
               'som',
               'pam_euclidean', 'pam_pearson', 'pam_spearman')
  harmon <- c('Raw', 'TC', 'UQ', 'med', 'TMM', 'DESeq', 'PoissonSeq', 'QN', 
              'RUVr', 'RUVs', 'RUVg')
  
  ari_indexes <- data.frame(matrix(nrow=length(cluster), ncol=length(harmon)))
  silhouette_indexes <- data.frame(matrix(nrow=length(cluster), ncol=length(harmon)))
  true_label <- as.factor(c(rep('MXF',100), rep('PMFH',100)))
  
  for (i in 1:length(cluster)){
    for (j in 1:length(harmon)){
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
      ari_indexes[i,j] <- mclust::adjustedRandIndex(true_label, est_cluster)
      # Calculate silhouette value
      if (startsWith(cluster[i], "hc_") || startsWith(cluster[i], "pam_")) {
        dist_matrix <- switch(sub(".*_", "", cluster[i]),
                            "euclidean" = dist(t(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized)),
                            "pearson" = as.dist(1 - cor(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized)),
                            "spearman" = as.dist(1 - cor(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized, method = "spearman")))
      } else {
        dist_matrix <- dist(t(analysis@harmon.train.data[[harmon[j]]]$dat.harmonized))
      }
      sil <- cluster::silhouette(as.numeric(est_cluster), dist_matrix)
      silhouette_indexes[i,j] <- mean(sil[, "sil_width"])
    }
  }
  rownames(ari_indexes) <- cluster
  colnames(ari_indexes) <- harmon
  rownames(silhouette_indexes) <- cluster
  colnames(silhouette_indexes) <- harmon
  
  indexes <- list(ari = ari_indexes,
                 silhouette = silhouette_indexes)
  return(indexes)
}


# scenario 1: clean data analysis
clean.scenario <- function(c){
  clean.datasets <- lapply(1:300, function(x) biological.effects(benchmark_sub[[x]]$data,
                                                                 benchmark_sub[[x]]$label, c))
  return(clean.datasets)
}

clean.data.results <- function(c){
  cleandata <- clean.scenario(c)
  indexessummary <- mclapply(1:300, measure, datalist = cleandata, mc.cores = 50)
  return(list(indexessummary = indexessummary,
              cleandata = cleandata))
}

# clean_results_list <- list()
# for (param in c(0.2, 0.6, 1)) {
#   key <- paste0("clean_res_c", param)
#   clean_results_list[[key]] <- clean.data.results(param)
#   gc()
# }
# save(clean_results_list, file = "clustering_sceanrio1.RData")

## more for server
load("clustering_sceanrio1.RData")

null_indices <- which(sapply(clean_results_list[[3]]$indexessummary, is.null))
rerun_results <- mclapply(null_indices, measure,
                          datalist = clean_results_list[[3]]$cleandata, mc.cores = 50)
clean_results_list[[3]]$indexessummary[null_indices] <- rerun_results

save(clean_results_list, file = "clustering_sceanrio1.RData")


